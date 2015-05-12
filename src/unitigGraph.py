import networkx as nx
import cPickle as pickle
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
        self.kmers = set()
        self.source_kmers = set()
        # this dict stores the source paralogs and their original genome positions
        self.paralogs = OrderedDict()
        if derived is False:
            # masked kmers stores kmers that were repeat masked and should be ignored in individual data
            self.masked_kmers = set()
            # these are edges that were pruned from the source graph and should NOT be re-introduced by the individual
            self.bad_source_edges = set()

    def _build_source_nodes(self, kmer, pos, name):
        """
        Takes a kmer, and adds it to the graph, constructing a sequence edge for the kmer.
        Returns the canonical representation of this kmer.
        """
        kmer_canonical = canonical(kmer)
        l, r = labels_from_kmer(kmer_canonical)
        # keep track of the graph as it grows and make sure its always an even # of nodes
        prev_size = len(self)
        if self.has_node(l) is not True and self.has_node(r) is not True:
            # should not be possible to have just left or just right
            assert not (self.has_node(l) or self.has_node(r))
            self.add_node(l)
            self.add_node(r)
            self.add_edge(l, r, positions={})
            assert name not in self.edge[l][r]['positions'], "Error: this implementation cannot handle multiple" \
                                                             "instances of a kmer in a block. Try a larger size kmer" \
                                                             "or more repeat masking."
            self.edge[l][r]['positions'][name] = pos
            assert prev_size + 2 == len(self)
        else:
            assert name not in self.edge[l][r]['positions'], "Error: this implementation cannot handle multiple" \
                                                             "instances of a kmer in a block. Try a larger size kmer" \
                                                             "or more repeat masking."
            self.edge[l][r]['positions'][name] = pos
            assert prev_size == len(self)
        self.kmers.add(kmer_canonical)
        self.source_kmers.add(kmer_canonical)

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

    def _add_source_adjacency(self, prev_kmer, kmer):
        """
        Adds L/R nodes, sequence edge and adjacency edge to the graph
        """
        prev_kmer_canonical = canonical(prev_kmer)
        kmer_canonical = canonical(kmer)
        # make sure we aren't adding edges - they should already exist now
        prev_kmer_size = len(self)
        l, r = self._determine_orientation(prev_kmer, prev_kmer_canonical, kmer, kmer_canonical)
        self.add_edge(l, r)
        self.edge[l][r]['source'] = True
        # no new nodes should be created in this process
        assert prev_kmer_size == len(self), (prev_kmer, kmer)

    def add_masked_kmers(self, masked_seq, unmasked_seq):
        for i in xrange(len(masked_seq) - self.kmer_size + 1):
            if "N" in masked_seq[i:i + self.kmer_size]:
                self.masked_kmers.add(canonical(unmasked_seq[i:i + self.kmer_size]))

    def add_source_sequence(self, name, offset, masked_seq, unmasked_seq):
        """
        masked_seq should be the same length as unmasked_seq and represent the repeat masked version.
        """
        assert name not in self.paralogs
        self.paralogs[name] = [offset, offset + len(unmasked_seq)]
        self.paralogs = OrderedDict(sorted(self.paralogs.iteritems(), key=lambda x: x[0]))
        prev_kmer = masked_seq[:self.kmer_size]
        prev_kmer_unmasked = unmasked_seq[:self.kmer_size]
        prev_pos = 0
        start_flag = False
        for i in xrange(1, len(masked_seq) - self.kmer_size + 1):
            kmer = masked_seq[i:i + self.kmer_size]
            kmer_unmasked = unmasked_seq[i:i + self.kmer_size]
            if canonical(prev_kmer_unmasked) in self.masked_kmers:
                prev_kmer = kmer
                prev_kmer_unmasked = kmer_unmasked
                prev_pos = i
                continue
            elif canonical(kmer_unmasked) in self.masked_kmers:
                continue
            else:
                if start_flag is False:
                    self._build_source_nodes(prev_kmer, prev_pos, name)
                    start_flag = True
                prev_kmer = masked_seq[prev_pos:prev_pos + self.kmer_size]
                assert "N" not in prev_kmer and "N" not in kmer, (prev_pos, prev_kmer, i, kmer)
                self._build_source_nodes(kmer, i, name)
                self._add_source_adjacency(prev_kmer, kmer)
                prev_kmer = kmer
                prev_kmer_unmasked = kmer_unmasked
                prev_pos = i
        assert len(self.masked_kmers & self.kmers) == 0

    def add_individual_sequence(self, seq):
        """
        Adds individual sequences from jellyfish kmer+1 counts.
        """
        assert len(seq) == self.kmer_size + 1
        prev, kmer = seq[:-1], seq[1:]
        prev_canonical = canonical(prev)
        kmer_canonical = canonical(kmer)
        if prev_canonical not in self.masked_kmers and kmer_canonical not in self.masked_kmers:
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
                a_seqs = self.edge[a_l][a_r]['positions'].keys()
                b_seqs = self.edge[b_l][b_r]['positions'].keys()
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
                    self.edge[l][r]['label'] = "\\n".join(
                        sorted([": ".join(map(str, x)) for x in self.edge[l][r]['positions'].iteritems()]))
                self.edge[l][r]["color"] = "purple"
            elif 'source' in self.edge[l][r]:
                self.edge[l][r]['color'] = "blue"
            else:
                self.edge[l][r]['color'] = "green"

    def connected_component_iter(self):
        """
        Yields connected components. Each subgraph is rebuilt to represent the kmers present within it.
        """
        for subgraph in nx.connected_component_subgraphs(self):
            # rebuild the kmer set based on this subgraph
            subgraph.kmers = {remove_label(n) for n in subgraph.nodes_iter()}
            # rebuild the paralog dict, but only for paralogs in this subgraph
            # first we find any sequence edge and find the associated paralogs - they should all be the same
            for a, b in subgraph.edges_iter():
                if 'positions' in subgraph.edge[a][b]:
                    paralogs = sorted(subgraph.edge[a][b]['positions'].iterkeys())
                    break
            for para in paralogs:
                subgraph.paralogs[para] = self.paralogs[para]
            # find which of these kmers are source kmers (used to normalize this unitig)
            # TODO: does this introduce a bug where repeated calls in individual pruning leads to incorrectly small
            # numbers of source_kmers?
            subgraph.source_kmers = self.source_kmers & subgraph.kmers
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