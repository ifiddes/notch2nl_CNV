import pulp
from collections import defaultdict
from src.abstractIlpSolving import SequenceGraphLpProblem
from src.helperFunctions import labels_from_kmer


class Block(object):
    """
    Represents one block (connected component). Each block may have anywhere from 1 to n paralogs represented in it.
    Initialized by passing one connected component subgraph from a UnitigGraph. Creates a ILP variable for each
    instance of a paralog in this unitig.
    """
    def __init__(self, subgraph, min_ploidy=0, max_ploidy=4):
        self._size = len(subgraph.source_kmers)
        self.variable_map = {}
        self.adjusted_count = None
        # each block gets its own trash bin - a place for extra kmer counts to go
        self.trash = pulp.LpVariable(str(id(self)), lowBound=0)
        # adjust the kmer set to remove kmers flagged as bad
        self.kmers = set()
        for k in subgraph.kmers:
            l, r = labels_from_kmer(k)
            if 'bad' not in subgraph.edge[l][r]:
                self.kmers.add(k)
            else:
                if k in subgraph.source_kmers:
                    self._size -= 1
        for para in subgraph.paralogs:
            for start, stop in subgraph.paralogs[para]:
                if len(self.kmers) > 0:
                    self.variable_map[(para, start, stop)] = pulp.LpVariable("{}_{}".format(para, start),
                                                                             lowBound=min_ploidy, upBound=max_ploidy,
                                                                             cat="Integer")
                else:
                    self.variable_map[(para, start, stop)] = None

    def __len__(self):
        return self._size

    def variable_iter(self):
        for (para, start, stop), variable in self.variable_map.iteritems():
            yield para, start, stop, variable

    @property
    def variables(self):
        return [v for v in self.variable_map.itervalues() if v is not None]


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
        # maps Block objects to the paralogs they contain
        self.block_map = {x: [] for x in UnitigGraph.paralogs.iterkeys()}
        self._build_blocks(UnitigGraph, infer_c, infer_d)

    def _build_blocks(self, UnitigGraph, infer_c, infer_d):
        """
        Requires a UnitigGraph, which is a networkx UnitigGraph built over the genome region of interest.
        See unitigGraph.py
        """
        # build the blocks, don't tie them together yet
        for subgraph in UnitigGraph.connected_component_iter():
            b = Block(subgraph)
            self.blocks.append(b)

        # build the block map, which relates a block to its position
        for block in self.blocks:
            for para, start, stop, variable in block.variable_iter():
                self.block_map[para].append([start, stop, variable, block])

        # sort the block maps by start positions
        for para in self.block_map:
            self.block_map[para] = sorted(self.block_map[para], key=lambda x: x[0])

        # now we tie the blocks together
        for para in self.block_map:
            # filter out all blocks without variables (no kmers)
            variables = [var for start, stop, var, block in self.block_map[para] if var is not None]
            for i in xrange(1, len(variables)):
                var_a, var_b = variables[i - 1], variables[i]
                self.constrain_approximately_equal(var_a, var_b, penalty=self.breakpoint_penalty)

        # tie each block sum to be approximately equal to the expected value subject to expected_value_penalty
        for block in self.blocks:
            exp = sum(self.expected_ploidy[p] for p, st, so, v in block.variable_iter() if v is not None)
            self.constrain_approximately_equal(exp, sum(block.variables), penalty=self.expected_value_penalty)

        # now we force all Notch2 variables to be equal to 2
        if "Notch2" in self.block_map:
            for start, stop, var, block in self.block_map["Notch2"]:
                self.add_constraint(var == 2)

        # if we have previously inferred C/D copy numbers, set those values
        if infer_c is not None:
            for start, stop, var, block in self.block_map["Notch2NL-C"]:
                self.add_constraint(var == infer_c)

        if infer_d is not None:
            for start, stop, var, block in self.block_map["Notch2NL-D"]:
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

    def report_normalized_raw_data_map(self):
        """
        Reports copy number for each ILP variable, once. If a block lacks a variable, reports the previous value.
        format [position, span, value]
        """
        copy_map = defaultdict(list)
        for para in self.block_map:
            for i in xrange(len(self.block_map[para])):
                start, stop, var, block = self.block_map[para][i]
                start += self.graph.paralogs[para][0]
                stop += self.graph.paralogs[para][0]
                if var is not None:
                    copy_map[para].append([start, stop, block.adjusted_count / len(block.variables)])
                    prev_var = block.adjusted_count / len(block.variables)
                else:
                    copy_map[para].append([start, stop, prev_var])
        return copy_map

    def report_copy_map(self):
        """
        Reports the raw counts seen at each variable. This is normalized by the number of variables in the block.
        """
        copy_map = defaultdict(list)
        for para in self.block_map:
            for i in xrange(len(self.block_map[para])):
                start, stop, var, block = self.block_map[para][i]
                start += self.graph.paralogs[para][0]
                stop += self.graph.paralogs[para][0]
                if var is not None:
                    copy_map[para].append([start, stop, pulp.value(var)])
                    prev_var = pulp.value(var)
                else:
                    copy_map[para].append([start, stop, prev_var])
        return copy_map