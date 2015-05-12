import pulp
from collections import defaultdict
from itertools import izip
from src.abstractIlpSolving import SequenceGraphLpProblem
from src.helperFunctions import labels_from_kmer


class Block(object):
    """
    Represents one block (connected component). Each block may have anywhere from 1 to n paralogs represented in it.
    Initialized by passing one connected component subgraph from a UnitigGraph.

    TODO: this new implementation removes the concept of instances. It is not possible for a unitig to appear more than
    once in the source sequence. If it does, we will have no idea where it came from. I re-wrote UnitigGraph to fail
    to build a graph if this happens.
    """
    def __init__(self, subgraph, min_ploidy=0, max_ploidy=4):
        self._size = len(subgraph.source_kmers)
        self._variables = {}
        self.paralogs = None
        for a, b in subgraph.edges_iter():
            if 'positions' in subgraph.edge[a][b]:
                self.paralogs = sorted(subgraph.edge[a][b]['positions'].iterkeys())
                break
        assert self.paralogs is not None
        self.adjusted_count = None
        # each block gets its own trash bin - a place for extra kmer counts to go
        self.trash = pulp.LpVariable(str(id(self)), lowBound=0)
        # adjust the kmer set to remove kmers flagged as bad
        self.kmers = set()
        for k in subgraph.kmers:
            l, r = labels_from_kmer(k)
            if 'bad' not in subgraph.edge[l][r]:
                self.kmers.add(k)
        for para in subgraph.paralogs.iterkeys():
            if len(self.kmers) > 0:
                self._variables[para] = pulp.LpVariable("{}_{}".format(para, str(id(self))), lowBound=min_ploidy,
                                                        upBound=max_ploidy, cat="Integer")
            else:
                self._variables[para] = None
        self.variable_map = defaultdict(dict)
        for k in subgraph.source_kmers:
            l, r = labels_from_kmer(k)
            for para, pos in subgraph.edge[l][r]['positions'].iteritems():
                self.variable_map[para][pos] = self._variables[para]

    def __len__(self):
        return self._size

    @property
    def variables(self):
        return [v for v in self._variables.itervalues() if v is not None]


class KmerIlpModel(SequenceGraphLpProblem):
    """
    Represents a kmer UnitigGraph model to infer copy number of highly
    similar paralogs using ILP.

    normalizing is the average kmer count in a region with known copy number of 2.

    This model can be fine tuned with the following parameters:
    breakpoint_penalty: determines how much we allow neighboring variables to differ.
        Larger values restricts breakpoints.
    data_penalty: how much do we trust the data? Larger values locks the variables more
        strongly to the data.
    trash_penalty: How much do we want to favor dumping extra kmer counts into the trash bin?
    expected_value_penalty: How much do we want to favor the expected copy number total per block?
    """
    def __init__(self, UnitigGraph, breakpoint_penalty=1, data_penalty=1, trash_penalty=1, expected_value_penalty=1,
                 infer_c=None, infer_d=None, default_ploidy=2):
        SequenceGraphLpProblem.__init__(self)
        self.graph = UnitigGraph
        self.breakpoint_penalty = breakpoint_penalty
        self.data_penalty = data_penalty
        self.trash_penalty = trash_penalty
        self.expected_value_penalty = expected_value_penalty
        self.expected_ploidy = {}
        for para in UnitigGraph.paralogs.iterkeys():
            if para == "Notch2NL-C" and infer_c is not None:
                self.expected_ploidy[para] = infer_c
            elif para == "Notch2NL-D" and infer_d is not None:
                self.expected_ploidy[para] = infer_d
            else:
                self.expected_ploidy[para] = default_ploidy
        # list of Block objects
        self.blocks = []
        self.block_map = defaultdict(list)
        self._build_blocks(UnitigGraph, infer_c, infer_d)

    def _build_blocks(self, UnitigGraph, infer_c, infer_d):
        """
        Requires a UnitigGraph, which is a networkx UnitigGraph built over the genome region of interest.
        See unitigGraph.py
        """
        # build the blocks
        for subgraph in UnitigGraph.connected_component_iter():
            b = Block(subgraph)
            self.blocks.append(b)

        # build the block map
        for b in self.blocks:
            for para in b.variable_map:
                for pos, var in b.variable_map[para].iteritems():
                    self.block_map[para].append([pos, var, b])

        # sort the block maps by start positions
        for para in self.block_map:
            self.block_map[para] = sorted(self.block_map[para], key=lambda x: x[0])

        # now we tie the blocks together
        for para in self.block_map:
            # filter out all blocks without variables (no kmers)
            variables = [var for pos, var, block in self.block_map[para] if var is not None]
            for i in xrange(1, len(variables)):
                var_a, var_b = variables[i - 1], variables[i]
                if var_a != var_b:
                    self.constrain_approximately_equal(var_a, var_b, penalty=self.breakpoint_penalty)

        # tie each block sum to be approximately equal to the expected value subject to expected_value_penalty
        for block in self.blocks:
            exp = sum(self.expected_ploidy[p] for p in block.paralogs)
            self.constrain_approximately_equal(exp, sum(block.variables), penalty=self.expected_value_penalty)

        # now we force all Notch2 variables to be equal to 2
        for pos, var, block in self.block_map["Notch2"]:
            self.add_constraint(var == 2)

        # if we have previously inferred C/D copy numbers, set those values
        if infer_c is not None:
            for pos, var, block in self.block_map["Notch2NL-C"]:
                self.add_constraint(var == infer_c)
        if infer_d is not None:
            for pos, var, block in self.block_map["Notch2NL-D"]:
                self.add_constraint(var == infer_d)

        # finally, constrain all trash bins to be as close to zero as possible
        for b in self.blocks:
            self.constrain_approximately_equal(b.trash, 0, penalty=self.trash_penalty)

    def introduce_data(self, kmer_counts, normalizing):
        """
        Introduces data to this ILP kmer model. For this, the input is assumed to be a dict
        representing the results of kmer counting a WGS dataset (format seq:count)
        """
        for block in self.blocks:
            if len(block) == 0:
                continue
            count = sum(kmer_counts.get(k, 0) for k in block.kmers)
            adjusted_count = (2.0 * count) / (len(block.kmers) * normalizing)
            block.adjusted_count = adjusted_count
            self.constrain_approximately_equal(adjusted_count, sum(block.variables + [block.trash]),
                                               penalty=self.data_penalty)

    def report_results(self):
        result = defaultdict(list)
        prev_var, prev_raw_val = None, None
        for para in self.block_map:
            for pos, var, block in self.block_map[para]:
                raw_val = block.adjusted_count / len(block.variables)
                if var is not None:
                    result[para].append([pos, pulp.value(var), raw_val])
                elif prev_var is not None:
                    result[para].append([pos, prev_var, prev_raw_val])
                prev_var = var
                prev_raw_val = raw_val
        return result