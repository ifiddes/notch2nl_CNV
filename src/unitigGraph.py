import networkx as nx
from collections import defaultdict
from src.helperFunctions import remove_label, create_labels, strandless
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

    def __init__(self, kmer_size=50):
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
        self.construct_nodes(name, seq)

    def construct_nodes(self, name, seq):
        """
        constructs left and right nodes for each kmer, adding a sequence edge between nodes.
        """
        for i in xrange(len(seq) - self.kmer_size + 1):
            kmer = strandless(seq[i:i + self.kmer_size].upper())
            if "N" in kmer:
                continue
            self.kmers.add(kmer)
            l = kmer + "_L"
            r = kmer + "_R"
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

    def construct_adjacencies(self, seq):
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
                    self.add_edge(prev + "_R", kmer + "_L")
                else:
                    # entering right side of next kmer
                    self.add_edge(prev + "_R", reverse_complement(kmer) + "_R")
            else:
                # exiting left side of previous kmer
                if kmer == kmer_strandless:
                    # entering left side of next kmer
                    self.add_edge(reverse_complement(prev) + "_L", kmer + "_L")
                else:
                    # entering right side of next kmer
                    self.add_edge(reverse_complement(prev) + "_L", reverse_complement(kmer) + "_R")
            assert prev_size == len(self)
            prev = kmer
            prev_strandless = kmer_strandless
        self.has_sequences = True

    def prune_graph(self):
        """
        Creates unitig graphs out of the graph by removing all adjacency edges which fit the following rules:
        1) adjacent nodes have multiple non self-loop edges
        2) adjacent nodes have different counts
        """
        assert self.has_sequences
        assert self.is_pruned is False
        to_remove = []
        for n in self.nodes():
            if self.degree(n) > 2:
                for a, b in self.edges(n):
                    if a == b or remove_label(a) != remove_label(b):
                        to_remove.append([a, b])
            elif self.degree(n) == 2:
                # make sure we aren't leaving edges connecting sequences with different counts
                for a, b in self.edges(n):
                    if remove_label(a) != remove_label(b):
                        dest_a, dest_b = create_labels(b)
                        source_a, source_b = create_labels(a)
                        if self.edge[source_a][source_b]['count'] != self.edge[dest_a][dest_b]['count']:
                            to_remove.append([a, b])
        for a, b in to_remove:
            if self.has_edge(a, b):
                self.remove_edge(a, b)
        self.is_pruned = True

    def finish_build(self, graphviz=False):
        """
        Finishes building the graph.

        If graphviz is true, adds a label tag to each sequence edge to improve understandability in graphviz
        """
        assert self.is_pruned is True and self.has_sequences is True
        assert self.is_built is False
        if graphviz is True:
            self_loops = set(self.selfloop_edges())
            for edge in self.edges():
                if edge not in self_loops and remove_label(edge[0]) == remove_label(edge[1]):
                    # make a fancy label for this sequence edge
                    l = remove_label(edge[0]) + "\\n" + " - ".join(
                        [": ".join([y, ", ".join([str(x) for x in self.edge[edge[0]][edge[1]]['positions'][y]])]) for y
                         in self.edge[edge[0]][edge[1]]['positions']]) + "\\ncount: " + str(
                        self.edge[edge[0]][edge[1]]["count"])
                    self.edge[edge[0]][edge[1]]["label"] = l
                else:
                    self.edge[edge[0]][edge[1]]['penwidth'] = 2

        self.paralogs = sorted(self.paralogs, key=lambda x: x[0])
        # default weight for each kmer is 2.0 (for diploid) unless derived weights are used (use weight_kmers())
        self.weights = {x: 2.0 for x in self.kmers}
        assert len(self.kmers.intersection(self.normalizing_kmers)) == 0
        self.is_built = True

    def connected_component_iter(self):
        """
        Yields connected components.
        """
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
