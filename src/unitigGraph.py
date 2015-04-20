import networkx as nx
from collections import defaultdict
from src.helperFunctions import remove_label, labels_from_kmer, labels_from_node, labeled_kmer_iter, strandless
from jobTree.src.bioio import reverseComplement as reverse_complement


class UnitigGraph(nx.Graph):
    """
    This UnitigGraph uses a sequence edge/adjacency edge paradigm to represent both strands without implementing
    a full bidirected DeBruijn graph. This paradigm works as follows:
    Each (strandless) kmer is assigned two nodes in the graph, a left and right node. A undirected edge connects these
    two nodes called a sequence edge. Adjacency edges are undirected edges that connect left or right nodes to their
    neighbor in the input sequence. In this way, you can enter or exit a kmer from either side, creating a bidirected
    graph. Any path through this graph must alternate between sequence edges and adjacency edges.

    Strandless kmers are defined as whichever comes first lexicographically, the kmer or its reverse complement. This
    is so that we can properly represent input sequencing data which we do not have strand information for.

    To finish representing a UnitigGraph, this resulting graph must be pruned. Any node which has more than one
    adjacency edge exiting or entering the node has all adjacency edges removed.
    """

    def __init__(self, kmer_size=49):
        nx.Graph.__init__(self)
        self.kmer_size = kmer_size
        self.has_sequences = False
        self.has_normalizing = False
        self.is_pruned = False
        self.is_built = False
        self.paralogs = []
        self.kmers = set()
        self.normalizing_kmers = set()
        self.sizes = {}
        self.weights = {}

    def add_normalizing(self, seq):
        """
        Adds normalizing kmers to the graph. These kmers are from a region of notch2
        that is after the duplication breakpoint, and so has exactly two copies in everyone.

        These kmers are stored as whichever strand comes first lexicographically.
        This is how jellyfish in the -C mode will report the kmer.
        """
        for i in xrange(len(seq) - self.kmer_size + 1):
            s = strandless(seq[i:i + self.kmer_size].upper())
            if "N" in s:
                continue
            self.normalizing_kmers.add(s)
        self.has_normalizing = True

    def construct_ref_nodes(self, name, offset, seq):
        """
        Constructs nodes for the reference portion of the graph. Name should be the paralog name,
        while offset is the position away from the start of the chromosome for that paralog.
        """
        self.paralogs.append([name, offset])
        self.sizes[name] = len(seq)
        for i in xrange(len(seq) - self.kmer_size + 1):
            kmer = strandless(seq[i:i + self.kmer_size].upper())
            if "N" in kmer:
                continue
            self.kmers.add(kmer)
            l, r = labels_from_kmer(kmer)
            if self.has_node(l) is not True and self.has_node(r) is not True:
                # should not be possible to have just left or just right
                assert not (self.has_node(l) or self.has_node(r))
                self.add_node(l)
                self.add_node(r)
                self.add_edge(l, r, count=1, positions=defaultdict(list))
                self.edge[l][r]['positions'][name].append(i)
            else:
                self.edge[l][r]['count'] += 1
                self.edge[l][r]['positions'][name].append(i)
            assert len(self) % 2 == 0

    def construct_individual_nodes(self, seq):
        """
        Adds additional nodes for individual to find SNP bubbles.
        """
        for i in xrange(len(seq) - self.kmer_size + 1):
            kmer = strandless(seq[i:i + self.kmer_size].upper())
            if "N" in kmer:
                continue
            self.kmers.add(kmer)
            l, r = labels_from_kmer(kmer)
            if self.has_node(l) is not True and self.has_node(r) is not True:
                # should not be possible to have just left or just right
                assert not (self.has_node(l) or self.has_node(r))
                self.add_node(l)
                self.add_node(r)
                self.add_edge(l, r)
            assert len(self) % 2 == 0

    def construct_adjacencies(self, seq, source_seq=True):
        """
        Constructs adjacency edges.
        """
        prev = seq[:self.kmer_size].upper()
        prev_strandless = strandless(prev)
        for i in xrange(1, len(seq) - self.kmer_size + 1):
            prev_size = len(self)
            kmer = seq[i:i + self.kmer_size].upper()
            if "N" in kmer or "N" in prev:
                continue
            kmer_strandless = strandless(kmer)
            if prev == prev_strandless:
                # exiting right side of previous kmer
                if kmer == kmer_strandless:
                    # entering left side of next kmer
                    l, r = prev + "_R", kmer + "_L"
                    self.add_edge(l, r)
                    if source_seq is True:
                        self.edge[l][r]['source'] = True
                else:
                    # entering right side of next kmer
                    l, r = prev + "_R", reverse_complement(kmer) + "_R"
                    self.add_edge(l, r)
                    if source_seq is True:
                        self.edge[l][r]['source'] = True
            else:
                # exiting left side of previous kmer
                if kmer == kmer_strandless:
                    # entering left side of next kmer
                    l, r = reverse_complement(prev) + "_L", kmer + "_L"
                    self.add_edge(l, r)
                    if source_seq is True:
                        self.edge[l][r]['source'] = True
                else:
                    # entering right side of next kmer
                    l, r = reverse_complement(prev) + "_L", reverse_complement(kmer) + "_R"
                    self.add_edge(l, r)
                    if source_seq is True:
                        self.edge[l][r]['source'] = True
            assert prev_size == len(self)
            prev = kmer
            prev_strandless = kmer_strandless
        self.has_sequences = True

    def prune_source_edges(self):
        """
        Creates unitig graphs out of the graph by removing all adjacency edges which fit the following rules:
        1) adjacent nodes have multiple non self-loop edges
        2) adjacent nodes have different counts
        3) these edges came from the source sequence and not from individual reads
        """
        assert self.has_sequences
        assert self.is_pruned is False
        for n in self.nodes_iter():
            this_source_degree = self.source_degree(n)
            if this_source_degree > 1:
                # remove all source adjacency edges from this node
                for a, b in self.edges(n):
                    if a == b:
                        # edge case: is this edge a self loop?
                        continue
                    elif remove_label(a) == remove_label(b):
                        # never remove sequence edges
                        continue
                    if 'source' in self.edge[a][b]:
                        # only remove source edges (for now)
                        self.remove_edge(a, b)
            elif this_source_degree == 1:
                # make sure we aren't leaving edges connecting sequences with different counts
                for a, b in self.edges(n):
                    if a == b:
                        # edge case: is this edge a self loop?
                        continue
                    elif remove_label(a) == remove_label(b):
                        # never remove sequence edges
                        continue
                    if 'source' in self.edge[a][b]:
                        # only remove source edges (for now)
                        # find the sequence edges associated with a and b and make sure they have the same count
                        dest_a, dest_b = labels_from_node(b)
                        source_a, source_b = labels_from_node(a)
                        if self.edge[source_a][source_b]['count'] != self.edge[dest_a][dest_b]['count']:
                            self.remove_edge(a, b)

    def prune_individual_edges(self):
        """
        Removes individual edges that join unitigs together. Bubble edges are defined as edges which have a degree of
        at least 3.
        """
        to_remove = []
        for subgraph in self.connected_component_iter():
            # determine if this subgraph has a cycle - we know that this should be a tree, so if any node has
            # degree greater than 2 there is a cycle coming off of it. We ignore self loops (as long as they are
            # source)
            self_loop_nodes = [a for a, b in self.selfloop_edges() if 'source' in self.edge[a][b]]
            bubble_nodes = [n for n in subgraph.nodes_iter() if subgraph.degree(n) > 2 and n not in self_loop_nodes]
            # we count the number of source sequences each bubble node is attached to
            source_sequences = {}
            for n in bubble_nodes:
                a, b = labels_from_node(n)
                source_sequences[n] = tuple(subgraph.edge[a][b]['positions'].keys())
            if len(source_sequences) == 0:
                # this subgraph has no anchors - can't know which paralog, if any, this came from. Remove this subgraph.
                for n in subgraph.nodes_iter():
                    self.remove_node(n)
                continue
            if len(source_sequences) == 1:
                # there are no bridges that need resolving here
                continue
            largest = max(len(x) for x in source_sequences)
            for n, seqs in source_sequences.iteritems():
                if seqs == largest:
                    for a, b in self.edges(n):
                        if a == b:
                            # edge case: is this edge a self loop?
                            continue
                        elif remove_label(a) == remove_label(b):
                            # never remove sequence edges
                            continue
                        elif 'source' in self.edge[a][b]:
                            continue
                        to_remove.append([a, b])
                        subgraph.remove_edge(a, b)
            # this subgraph should now contain at least 2 unitigs - verify that this is true
            for new_subgraph in subgraph.connected_component_iter(internal=True):
                source_sequences = {tuple(new_subgraph.edge[a][b]['positions'].keys()) for a, b in
                                    new_subgraph.edges_iter() if 'positions' in new_subgraph.edge[a][b]}
                assert len(source_sequences) == 1, source_sequences
        for a, b in to_remove:
            self.remove_edge(a, b)

    def source_degree(self, n):
        """
        finds the degree of a node EXCLUDING edges that are not from the source sequence
        """
        return len([x for x in self.adj[n] if 'source' in self.edge[n][x]])

    def finish_build(self, graphviz=False):
        """
        Finishes building the graph.

        If graphviz is true, adds a label tag to each sequence edge to improve understandability in graphviz
        """
        # assert self.is_pruned is True and self.has_sequences is True
        # assert self.is_built is False
        if graphviz is True:
            self_loops = set(self.selfloop_edges())
            for edge in self.edges():
                if edge in self_loops:
                    continue
                l, r = edge
                self.node[l]['fontsize'] = self.node[r]['fontsize'] = 10
                self.edge[l][r]['penwidth'] = 3
                if remove_label(l) == remove_label(r):
                    self.edge[l][r]['fontsize'] = 8
                    # make a fancy label for this sequence edge if its from a source sequence
                    if 'positions' in self.edge[l][r]:
                        self.edge[l][r]["label"] = " - ".join(
                            [": ".join([y, ", ".join([str(x) for x in self.edge[l][r]['positions'][y]])]) for y in
                             self.edge[l][r]['positions']]) + "\\ncount: " + str(self.edge[l][r]["count"])
                    self.edge[l][r]["color"] = "purple"
                elif 'source' in self.edge[l][r]:
                    self.edge[l][r]['color'] = "blue"
                else:
                    self.edge[l][r]['color'] = "green"
        self.paralogs = sorted(self.paralogs, key=lambda x: x[0])
        # default weight for each kmer is 2.0 (for diploid) unless derived weights are used (use weight_kmers())
        self.weights = {x: 2.0 for x in self.kmers}
        self.is_built = True

    def connected_component_iter(self, internal=False):
        """
        Yields connected components. Internal is used for pruning edges within this object and should not be set to
        True outside of this setup.
        """
        if internal is False:
            assert self.is_built is True
        for subgraph in nx.connected_component_subgraphs(self):
            yield subgraph

    def flag_nodes(self, kmer_iter):
        """
        Iterates over a kmer_iter and flags nodes as being bad.
        This is used to flag nodes whose kmer is represented elsewhere
        in the genome, so that we won't count it later.
        """
        assert self.is_built is True
        for k in kmer_iter:
            k = k.rstrip()
            assert k in self.kmers
            self.edge[k + "_L"][k + "_R"]['bad'] = True

    def weight_kmers(self, weight_dict):
        """
        Takes a python dictionary mapping k1mers to an empirically derived
        weight. Applies a weight tag to each k1mer in the graph.
        """
        assert self.is_built is True
        for k, w in weight_dict.iteritems():
            assert k in self.kmers
            self.edge[k + "_L"][k + "_R"]['weight'] = w
