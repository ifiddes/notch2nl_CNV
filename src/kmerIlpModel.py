import pulp
from collections import defaultdict
from src.abstractIlpSolving import SequenceGraphLpProblem
from src.helperFunctions import remove_label


class Block(object):
    """
    Represents one block (connected component).
    Each block may have anywhere from 1 to n paralogs represented in
    it. Initialized by passing in all of the nodes from one UnitigGraph
    CC. Stores a mapping between each paralog this CC represents
    to the start and stop positions in the source sequence.
    Also stores all of the kmers in this block as a set.
    """
    def __init__(self, subgraph, min_ploidy=0, max_ploidy=4):
        self.variables = {}
        self.subgraph = subgraph
        self.adjusted_count = None
        # each block gets its own trash bin - a place for extra kmer counts to go
        self.trash = pulp.LpVariable(str(id(self)), lowBound=0)
        # find all sequence edges
        sequence_edges = [x for x in subgraph.edges() if remove_label(x[0]) == remove_label(x[1])]
        # find all valid kmers
        self.kmers = {remove_label(a) for a, b in sequence_edges if 'bad' not in subgraph[a][b]}
        begin_a, begin_b = sequence_edges[0]
        # now we want to find the first and last nodes
        # loop over paralogs
        for p in subgraph[begin_a][begin_b]['positions'].iterkeys():
            # pick a random sequence edge to start at (networkx stores these unordered)
            start_a, start_b = begin_a, begin_b
            # loop over the instances of this paralog
            for i in xrange(len(subgraph[start_a][start_b]['positions'][p])):
                # find the start and stop node for each instance of this paralog
                s = subgraph[start_a][start_b]['positions'][p][i]
                # loop over all remaining sequence edges and update the start as necessary
                for new_a, new_b in sequence_edges[1:]:
                    new_s = subgraph[new_a][new_b]['positions'][p][i]
                    if new_s < s:
                        s = new_s
                        start_a, start_b = new_a, new_b
                if len(self.kmers) > 0:
                    self.variables[(p, s)] = pulp.LpVariable("{}_{}".format(p, s), lowBound=min_ploidy,
                                                             upBound=max_ploidy, cat="Integer")
                else:
                    self.variables[(p, s)] = None

    def kmer_iter(self):
        """iterator over all kmers in this block"""
        for kmer in self.kmers:
            yield kmer

    def variable_iter(self):
        """iterator over variables and related values in this block"""
        for (para, start), variable in self.variables.iteritems():
            yield para, start, variable

    def get_variables(self):
        """returns all LP variables in this block"""
        return [x for x in self.variables.values() if x is not None]

    def get_kmers(self):
        """returns set of all kmers in block"""
        return self.kmers

    def get_trash(self):
        """returns the trash variable for this block"""
        return self.trash


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
    diploid_penalty: How much do we want to favor copy number of default_ploidy?
    trash_penalty: How much do we want to favor dumping extra kmer counts into the trash bin?
    expected_value_penalty: How much do we want to favor the expected copy number total per block?
    """
    def __init__(self, UnitigGraph, normalizing, breakpoint_penalty=15, data_penalty=1, diploid_penalty=1,
                 trash_penalty=1, expected_value_penalty=1, infer_c=None, infer_d=None, default_ploidy=2):
        SequenceGraphLpProblem.__init__(self)
        self.graph = UnitigGraph
        # list of Block objects
        self.blocks = []
        # maps Block objects to the paralogs they contain
        self.block_map = {x[0]: [] for x in UnitigGraph.paralogs}
        # maps the genome start position of each paralog for plotting later
        self.offset_map = {x[0]: int(x[1]) for x in UnitigGraph.paralogs}
        self.normalizing = normalizing
        self.breakpoint_penalty = breakpoint_penalty
        self.data_penalty = data_penalty
        self.trash_penalty = trash_penalty
        self.diploid_penalty = diploid_penalty
        self.expected_value_penalty = expected_value_penalty
        self.default_ploidy = default_ploidy
        self.infer_c = infer_c
        self.infer_d = infer_d
        # prepare expected_ploidy_dict, using previously inferred values for C/D if they were passed
        self.expected_ploidy_dict = {x[0]: self.default_ploidy for x in self.graph.paralogs}
        if infer_c is not None:
            self.expected_ploidy_dict["Notch2NL-C"] = infer_c
        if infer_d is not None:
            self.expected_ploidy_dict["Notch2NL-D"] = infer_d
        self.build_blocks(UnitigGraph)

    def build_blocks(self):
        """
        Requires a UnitigGraph, which is a networkx UnitigGraph built over the genome region of interest.
        See unitigGraph.py
        """
        # build the blocks, don't tie them together yet
        for subgraph in self.graph.connected_component_iter():
            b = Block(subgraph)
            self.blocks.append(b)

        # build the block map, which relates a block to its position
        for block in self.blocks:
            for para, start, variable in block.variable_iter():
                self.block_map[para].append([start, variable, block])

        # sort the block maps by start positions
        for para in self.block_map:
            self.block_map[para] = sorted(self.block_map[para], key=lambda x: x[0])

        # now we tie the blocks together
        for para in self.block_map:
            # filter out all blocks without variables (no kmers)
            variables = [var for start, var, block in self.block_map[para] if var is not None]
            for i in xrange(1, len(variables)):
                var_a, var_b = variables[i - 1], variables[i]
                self.constrain_approximately_equal(var_a, var_b, penalty=self.breakpoint_penalty)

        # tie each variable to be approximately equal to copy number 2 subject to the diploid_penalty constraint
        for block in self.blocks:
            for para, start, variable in block.variable_iter():
                self.constrain_approximately_equal(self.default_ploidy, variable, penalty=self.diploid_penalty)

        # tie each block sum to be approximately equal to the expected value subject to expected_value_penalty
        for block in self.blocks:
            exp = sum(self.expected_ploidy_dict[p] for p, s, v in block.variable_iter() if v is not None)
            self.constrain_approximately_equal(exp, sum(block.get_variables()), penalty=self.expected_value_penalty)

        # now we force all Notch2 variables to be equal to 2
        if "Notch2" in self.block_map:
            for start, var, block in self.block_map["Notch2"]:
                self.add_constraint(var == 2)

        # if we have previously inferred C/D copy numbers, set those values
        if self.infer_c is not None:
            for start, var, block in self.block_map["Notch2NL-C"]:
                self.add_constraint(var == self.infer_c)

        if self.infer_d is not None:
            for start, var, block in self.block_map["Notch2NL-D"]:
                self.add_constraint(var == self.infer_d)

        # finally, constrain all trash bins to be as close to zero as possible
        for b in self.blocks:
            self.constrain_approximately_equal(b.get_trash(), 0, penalty=self.trash_penalty)

    def introduce_data(self, kmer_counts):
        """
        Introduces data to this ILP kmer model. For this, the input is assumed to be a dict
        representing the results of kmer counting a WGS dataset (format seq:count)
        """
        for block in self.blocks:
            if len(block.kmers) == 0:
                continue
            count = sum(kmer_counts[k] * self.graph.weights[k] for k in block.get_kmers())
            adjusted_count = (1.0 * count) / (len(block.kmers) * self.normalizing)
            block.adjusted_count = adjusted_count
            self.constrain_approximately_equal(adjusted_count, sum(block.get_variables() + [block.get_trash()]),
                                               penalty=self.data_penalty)

    def report_normalized_raw_data_map(self):
        """
        Reports copy number for each ILP variable, once. If a block lacks a variable, reports the previous value.
        format [position, span, value]
        """
        copy_map = defaultdict(list)
        for para in self.block_map:
            offset = self.offset_map[para]
            for i in xrange(len(self.block_map[para]) - 1):
                start, var, block = self.block_map[para][i]
                span = self.block_map[para][i + 1][0] - start
                if var is not None:
                    copy_map[para].append([start + offset, span, block.adjusted_count / len(block.get_variables())])
                    prev_var = block.adjusted_count / len(block.get_variables())
                else:
                    copy_map[para].append([start + offset, span, prev_var])
            final_start, final_var, final_block = self.block_map[para][-1]
            final_span = self.graph.sizes[para] - final_start
            if final_var is not None:
                copy_map[para].append([final_start + offset, final_span, block.adjusted_count / len(block.get_variables())])
            else:
                copy_map[para].append([final_start + offset, final_span, prev_var])
        return copy_map

    def report_copy_map(self):
        """
        Reports the raw counts seen at each variable. This is normalized by the number of variables in the block.
        """
        copy_map = defaultdict(list)
        for para in self.block_map:
            offset = self.offset_map[para]
            for i in xrange(len(self.block_map[para]) - 1):
                start, var, block = self.block_map[para][i]
                span = self.block_map[para][i + 1][0] - start
                if var is not None:
                    copy_map[para].append([start + offset, span, pulp.value(var)])
                    prev_var = pulp.value(var)
                else:
                    copy_map[para].append([start + offset, span, prev_var])
            final_start, final_var, final_block = self.block_map[para][-1]
            final_span = self.graph.sizes[para] - final_start
            if final_var is not None:
                copy_map[para].append([final_start + offset, final_span, pulp.value(var)])
            else:
                copy_map[para].append([final_start + offset, final_span, prev_var])
        return copy_map