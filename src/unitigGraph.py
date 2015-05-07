import networkx as nx
import cPickle as pickle
from itertools import groupby
from operator import itemgetter
from collections import defaultdict, OrderedDict
from src.helperFunctions import remove_label, labels_from_kmer, labels_from_node, canonical
from jobTree.src.bioio import reverseComplement as reverse_complement


class UnitigGraph(nx.Graph):
    """
    This UnitigGraph uses a sequence edge/adjacency edge paradigm to represent both strands without implementing
    a full bidirected DeBruijn graph. This paradigm works as follows:
    Each (canonical) kmer is assigned two nodes in the graph, a left and right node. A undirected edge connects these
    two nodes called a sequence edge. Adjacency edges are undirected edges that connect left or right nodes to their
    neighbor in the input sequence. In this way, you can enter or exit a kmer from either side, creating a bidirected
    graph. Any path through this graph must alternate between sequence edges and adjacency edges.

    canonical kmers are defined as whichever comes first lexicographically, the kmer or its reverse complement. This
    is so that we can properly represent input sequencing data which we do not have strand information for.

    To finish representing a UnitigGraph, this resulting graph must be pruned. Any node which has more than one
    adjacency edge exiting or entering the node has all adjacency edges removed.

    Derived should be set to False for only the parent graph. Derived graphs (created by the copying function of
    networkx generally when finding connected subgraphs) will not retain the masking, source edge or normalization
    information.
    """

    def __init__(self, kmer_size=49, derived=True):
        nx.Graph.__init__(self)
        self.kmer_size = kmer_size
        self.paralogs = OrderedDict()
        self.kmers = set()
        self.source_kmers = set()
        if derived is False:
            # these are edges that were pruned from the source graph and should NOT be re-introduced by the individual
            self.bad_source_edges = set()
            # these are masked kmers used to prevent repeat masked sequence being added
            self.masked_kmers = set()

    def _determine_orientation(self, prev, prev_canonical, kmer, kmer_canonical):
        if prev == prev_canonical:
            # exiting right side of previous kmer
            if kmer == kmer_canonical:
                # entering left side of next kmer
                l, r = prev + "_R", kmer + "_L"
            else:
                # entering right side of next kmer
                l, r = prev + "_R", reverse_complement(kmer) + "_R"
        else:
            # exiting left side of previous kmer
            if kmer == kmer_canonical:
                # entering left side of next kmer
                l, r = reverse_complement(prev) + "_L", kmer + "_L"
            else:
                # entering right side of next kmer
                l, r = reverse_complement(prev) + "_L", reverse_complement(kmer) + "_R"
        return l, r

    def _add_source_nodes(self, pos, name, a_canonical, b_canonical):
        """
        Produces two sets of sources nodes given a k+1mer, if they do not exist.
        """
        # keep track of the graph as it grows and make sure its always an even # of nodes
        prev_size = len(self)
        for i, k in enumerate([a_canonical, b_canonical]):
            l, r = labels_from_kmer(k)
            if self.has_node(l) is not True and self.has_node(r) is not True:
                # should not be possible to have just left or just right
                assert not (self.has_node(l) or self.has_node(r))
                self.add_node(l)
                self.add_node(r)
                self.add_edge(l, r, positions=defaultdict(set))
                self.edge[l][r]['positions'][name].add(pos + i)
                assert prev_size + 2 == len(self), (a_canonical, b_canonical)
                prev_size = len(self)
            else:
                self.edge[l][r]['positions'][name].add(pos + i)

    def _add_source_adjacency(self, a, a_canonical, b, b_canonical):
        self.kmers.update([a_canonical, b_canonical])
        # make sure we aren't adding nodes - they should already exist now
        prev_kmer_size = len(self)
        l, r = self._determine_orientation(a, a_canonical, b, b_canonical)
        self.add_edge(l, r)
        self.edge[l][r]['source'] = True
        assert prev_kmer_size == len(self), (a, b)

    def add_masked_kmers(self, masked_seq, unmasked_seq):
        """
        Builds the set of masked kmers. Should be ran before add_source_sequences()
        """
        for i in xrange(len(masked_seq) - self.kmer_size + 1):
            if "N" in masked_seq[i:i + self.kmer_size]:
                self.masked_kmers.add(canonical(unmasked_seq[i:i + self.kmer_size]))

    def add_source_sequence(self, name, offset, unmasked_seq):
        """
        masked kmers should already have been added to the graph via add_masked_kmers()
        """
        assert name not in self.paralogs
        self.paralogs[name] = [offset, offset + len(unmasked_seq)]
        self.paralogs = OrderedDict(sorted(self.paralogs.iteritems(), key=lambda x: x[0]))
        for i in xrange(len(unmasked_seq) - self.kmer_size):
            k1mer = unmasked_seq[i:i + self.kmer_size + 1]
            a, b = k1mer[:-1], k1mer[1:]
            a_canonical, b_canonical = canonical(a), canonical(b)
            if a_canonical not in self.masked_kmers and b_canonical not in self.masked_kmers:
                self._add_source_nodes(i, name, a_canonical, b_canonical)
                self._add_source_adjacency(a, a_canonical, b, b_canonical)

    def add_individual_sequence(self, seq):
        """
        Adds individual sequences from jellyfish kmer+1 counts.
        """
        assert len(seq) == self.kmer_size + 1
        prev, kmer = seq[:-1], seq[1:]
        prev_canonical = canonical(prev)
        kmer_canonical = canonical(kmer)
        # do not add individual sequence that were originally masked in the kmer masking procedure
        l, r = self._determine_orientation(prev, prev_canonical, kmer, kmer_canonical)
        if frozenset(sorted([l, r])) not in self.bad_source_edges:
            for k in [prev_canonical, kmer_canonical]:
                a, b = labels_from_kmer(k)
                if self.has_node(a) is not True and self.has_node(b) is not True:
                    # should not be possible to have just left or just right
                    assert not (self.has_node(a) or self.has_node(b))
                    self.add_node(a)
                    self.add_node(b)
                    # add sequence edge without any tags
                    self.add_edge(a, b)
            # build adjacency edge
            self.add_edge(l, r)
            self.kmers.update([prev_canonical, kmer_canonical])

    def prune_source_edges(self):
        """
        Prunes the source edges to create unitigs. This is defined as contiguous blocks of sequence that are the same
        in one or more paralogs. Thus, we remove all source adjacency edges that:
        1) connect two sequence edges that have different source sequences
        2) connect two sequence edges that have different counts (e.g. was in A twice then A once)
        """
        for a, b in self.edges_iter():
            if a == b:
                # ignore self loop edges
                continue
            if 'source' in self.edge[a][b]:
                # this is a source sequence adjacency edge
                # find the sequence edges this edge connect and check source sequences
                a_l, a_r = labels_from_node(a)
                b_l, b_r = labels_from_node(b)
                a_seqs = [para * len(positions) for para, positions in self.edge[a_l][a_r]['positions'].iteritems()]
                b_seqs = [para * len(positions) for para, positions in self.edge[b_l][b_r]['positions'].iteritems()]
                if sorted(a_seqs) != sorted(b_seqs):
                    self.bad_source_edges.add(frozenset(sorted([a, b])))
        for a, b in self.bad_source_edges:
            self.remove_edge(a, b)

    def _prune_individual_subgraph(self, parent, nodes_to_remove, edges_to_remove):
        these_nodes_to_remove = set()
        these_edges_to_remove = set()
        for subgraph in parent.connected_component_iter():
            source_sequences = set()
            for a, b in subgraph.edges_iter():
                if 'positions' in subgraph.edge[a][b]:
                    source_sequences.add(frozenset(subgraph.edge[a][b]['positions'].iterkeys()))
            if len(source_sequences) == 1:
                # no bubbles to resolve here
                continue
            elif len(source_sequences) == 0:
                # this subgraph has no anchors - can't know which paralog, if any, this came from. Remove this subgraph.
                these_nodes_to_remove.update(subgraph.nodes())
                nodes_to_remove.update(subgraph.nodes())
                continue
            # this subgraph needs to be resolved. Find the common source paralogs
            common_paralogs = frozenset.intersection(*source_sequences)
            # if there are no common paralogs, remove all individual sequence
            if len(common_paralogs) == 0:
                for n in subgraph.nodes_iter():
                    a, b = labels_from_node(n)
                    if 'positions' not in subgraph.edge[a][b]:
                        these_nodes_to_remove.update([a, b])
                        nodes_to_remove.update([a, b])
                for a, b in subgraph.edges_iter():
                    if 'positions' not in subgraph.edge[a][b] and 'source' not in subgraph.edge[a][b]:
                        these_edges_to_remove.add(frozenset(sorted([a, b])))
                        edges_to_remove.add(frozenset(sorted([a, b])))
            # remove any individual edge attached to a node whose source sequence is not common_paralogs
            else:
                for n in subgraph.nodes_iter():
                    l, r = labels_from_node(n)
                    if 'positions' in subgraph.edge[l][r]:
                        these_paralogs = frozenset(subgraph.edge[l][r]['positions'].iterkeys())
                        if these_paralogs != common_paralogs:
                            for next_node in subgraph.adj[n]:
                                if 'source' in subgraph.edge[n][next_node]:
                                    continue
                                elif remove_label(n) == remove_label(next_node):
                                    continue
                                these_edges_to_remove.add(frozenset(sorted([n, next_node])))
                                edges_to_remove.add(frozenset(sorted([n, next_node])))
            for a, b in these_edges_to_remove:
                subgraph.remove_edge(a, b)
            for n in these_nodes_to_remove:
                subgraph.remove_node(n)
            for new_subgraph in subgraph.connected_component_iter():
                source_sequences = {tuple(new_subgraph.edge[a][b]['positions'].keys()) for a, b in
                                    new_subgraph.edges_iter() if 'positions' in new_subgraph.edge[a][b]}
                if len(source_sequences) != 1:
                    self._prune_individual_subgraph(new_subgraph, nodes_to_remove, edges_to_remove)

    def prune_individual_edges(self):
        """
        For each remaining subgraph, determine if it is a unitig or if read data are joining unitigs.
        If no unitig information is present in this subgraph, remove it - we can't know which paralog(s), if any,
        these sequences came from. Otherwise, look for bubble nodes and remove the edge connecting the unitigs such that
        the unitig with the most common source sequences keeps the bubble. If it is a tie, discard these kmers.
        e.g. {A,B,C} - {A,B} -> {A,B}
             {A,B} - {A} -> {A}
             {A,B} - {A,C} -> {A}
             {A,A} - {A} -> {A}
             {A,C} - {B, D} -> discard
             {A} - {C} -> discard
        """
        nodes_to_remove = set()
        edges_to_remove = set()
        for subgraph in self.connected_component_iter():
            self._prune_individual_subgraph(subgraph, nodes_to_remove, edges_to_remove)
        for a, b in edges_to_remove:
            self.remove_edge(a, b)
        for n in nodes_to_remove:
            self.remove_node(n)
            k = remove_label(n)
            if k in self.kmers:
                self.kmers.remove(k)
        # make sure this worked...
        for i, subgraph in enumerate(self.connected_component_iter()):
            source_sequences = {tuple(subgraph.edge[a][b]['positions'].keys()) for a, b in
                                    subgraph.edges_iter() if 'positions' in subgraph.edge[a][b]}
            assert len(source_sequences) == 1

    def make_graphviz_labels(self):
        """
        Adds a label tag to each sequence edge to improve understandability in graphviz. Generally for debugging
        purposes.
        """
        self_loops = set(self.selfloop_edges())
        for edge in self.edges_iter():
            if edge in self_loops:
                continue
            l, r = edge
            self.node[l]['fontsize'] = self.node[r]['fontsize'] = 10
            self.edge[l][r]['penwidth'] = 3
            if remove_label(l) == remove_label(r):
                self.edge[l][r]['fontsize'] = 8
                # make a fancy label for this sequence edge if its from a source sequence
                if 'positions' in self.edge[l][r]:
                    self.edge[l][r]["label"] = "\\n".join(sorted(
                        [": ".join([y, ", ".join([str(x) for x in self.edge[l][r]['positions'][y]])]) for y in
                         self.edge[l][r]['positions']]))
                self.edge[l][r]["color"] = "purple"
            elif 'source' in self.edge[l][r]:
                self.edge[l][r]['color'] = "blue"
            else:
                self.edge[l][r]['color'] = "green"

    def connected_component_iter(self):
        """
        Yields connected components. Each subgraph is rebuilt to represent the kmers present within it.
        If a subgraph represents more than one source positions, the paralogs property will be updated accordingly.
        """
        for subgraph in nx.connected_component_subgraphs(self):
            # rebuild the kmer set based on this subgraph
            subgraph.kmers = {remove_label(n) for n in subgraph.nodes_iter()}
            # find which of these kmers are source kmers (used to normalize this unitig)
            subgraph.source_kmers = self.source_kmers & subgraph.kmers
            # build a flat list of all ranges for each paralog
            tmp_map = defaultdict(list)
            for k in subgraph.source_kmers:
                l, r = labels_from_kmer(k)
                if 'positions' not in subgraph.edge[l][r]:
                    continue
                for para, vals in subgraph.edge[l][r]['positions'].iteritems():
                    tmp_map[para].extend(vals)
            for para, data in tmp_map.iteritems():
                data = sorted(data)
                ranges = [map(itemgetter(1), g) for k, g in groupby(enumerate(data), lambda (i, x):i - x)]
                subgraph.paralogs[para] = [[x[0], x[-1]] for x in ranges]
            subgraph.paralogs = OrderedDict(sorted(subgraph.paralogs.iteritems(), key=lambda x: x[0]))
            yield subgraph

    def flag_nodes(self, kmer_iter):
        """
        Iterates over a kmer_iter and flags nodes as being bad.
        This is used to flag nodes whose kmer is represented elsewhere
        in the genome, so that we won't count it later.
        """
        for k in kmer_iter:
            k = k.rstrip()
            if k in self.kmers:
                l, r = labels_from_kmer(k)
                self.edge[l][r]['bad'] = True