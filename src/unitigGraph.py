import networkx as nx
from collections import defaultdict
from src.helperFunctions import remove_label, labels_from_kmer, labels_from_node, labeled_kmer_iter, canonical
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
        # masked kmers stores kmers that were repeat masked and should be ignored in individual data
        self.masked_kmers = set()
        self.normalizing_kmers = set()
        self.source_sequence_sizes = {}

    def add_normalizing(self, seq):
        """
        Adds normalizing kmers to the graph. These kmers are from a region of notch2
        that is after the duplication breakpoint, and so has exactly two copies in everyone.

        These kmers are stored as whichever strand comes first lexicographically.
        This is how jellyfish in the -C mode will report the kmer.
        """
        for i in xrange(len(seq) - self.kmer_size + 1):
            s = canonical(seq[i:i + self.kmer_size].upper())
            if "N" in s:
                continue
            self.normalizing_kmers.add(s)
        self.has_normalizing = True

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
            self.add_edge(l, r, positions=defaultdict(list))
            self.edge[l][r]['positions'][name].append(pos)
            assert prev_size + 2 == len(self)
        else:
            self.edge[l][r]['positions'][name].append(pos)
            assert prev_size == len(self)
            self.kmers.add(kmer_canonical)

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
        self.paralogs.append([name, offset])
        self.source_sequence_sizes[name] = len(masked_seq)
        prev_kmer = masked_seq[:self.kmer_size]
        prev_kmer_unmasked = unmasked_seq[:self.kmer_size]
        prev_pos = 0
        start_flag = False
        for i in xrange(1, len(masked_seq) - self.kmer_size + 1):
            kmer = masked_seq[i:i + self.kmer_size]
            kmer_unmasked = unmasked_seq[i:i + self.kmer_size]
            #assert len(self.masked_kmers & self.kmers) == 0, (i, kmer, prev_pos, prev_kmer)
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
            for k in [prev_canonical, kmer_canonical]:
                l, r = labels_from_kmer(k)
                if self.has_node(l) is not True and self.has_node(r) is not True:
                    # should not be possible to have just left or just right
                    assert not (self.has_node(l) or self.has_node(r))
                    self.add_node(l)
                    self.add_node(r)
                    self.add_edge(l, r)
            # build adjacency edge
            l, r = self._determine_orientation(prev, prev_canonical, kmer, kmer_canonical)
            self.add_edge(l, r)
            self.kmers.update([prev_canonical, kmer_canonical])
            # assert len(self.kmers & self.masked_kmers) == 0, ('source', prev, prev_canonical, kmer, kmer_canonical)

    def prune_source_edges(self):
        """
        Prunes the source edges to create unitigs. This is defined as contiguous blocks of sequence that are the same
        in one or more paralogs. Thus, we remove all source adjacency edges that:
        1) connect two sequence edges that have different source sequences
        2) connect two sequence edges that have different counts (e.g. was in A twice then A once)
        """
        to_remove = []
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
                    to_remove.append([a, b])
        for a, b in to_remove:
            self.remove_edge(a, b)

    def prune_individual_edges(self):
        """
        For each remaining subgraph, determine if it is a unitig or if read data are joining unitigs.
        If no unitig information is present in this subgraph, remove it - we can't know which paralog(s), if any,
        these sequences came from. Otherwise, look for bubble nodes and remove the edge connecting the unitigs such that
        the unitig with the most common source sequences keeps the bubble. If it is a tie, discard these kmers.
        e.g. {A,B,C} - {A,B} -> {A,B}
             {A,B} - {A} -> {A}
             {A,B} - {A,C} -> {A}
             {A,A - {A} -> {A}
             {A} - {C} -> discard
        """
        nodes_to_remove = set()
        edges_to_remove = set()
        for subgraph in nx.connected_component_subgraphs(self):
            source_sequences = set()
            for a, b in subgraph.edges_iter():
                if 'positions' in subgraph.edge[a][b]:
                    source_sequences.add(frozenset(self.edge[a][b]['positions'].iterkeys()))
            if len(source_sequences) == 1:
                # no bubbles to resolve here
                continue
            elif len(source_sequences) == 0:
                # this subgraph has no anchors - can't know which paralog, if any, this came from. Remove this subgraph.
                nodes_to_remove.update(subgraph.nodes())
                continue
            # this subgraph needs to be resolved. Find the common source paralogs
            common_paralogs = frozenset.intersection(*source_sequences)
            # if there are no common paralogs, remove all individual sequence
            if len(common_paralogs) == 0:
                for n in subgraph.nodes_iter():
                    a, b = labels_from_node(n)
                    if 'positions' not in subgraph.edge[a][b]:
                        nodes_to_remove.update([a, b])
            # remove any individual edge attached to a node whose source sequence is not common_paralogs
            else:
                for n in subgraph.nodes_iter():
                    l, r = labels_from_node(n)
                    if 'positions' in self.edge[l][r]:
                        these_paralogs = frozenset(self.edge[l][r]['positions'].iterkeys())
                        if these_paralogs != common_paralogs:
                            for next_node in subgraph.adj[n]:
                                if 'source' in subgraph.edge[n][next_node]:
                                    continue
                                elif remove_label(n) == remove_label(next_node):
                                    continue
                                edges_to_remove.add(frozenset(sorted([n, next_node])))
        for a, b in edges_to_remove:
            self.remove_edge(a, b)
        for n in nodes_to_remove:
            self.remove_node(n)
            k = remove_label(n)
            if k in self.kmers:
                self.kmers.remove(k)
        # debugging: make sure this worked
        for new_subgraph in nx.connected_component_subgraphs(self):
            source_sequences = {tuple(new_subgraph.edge[a][b]['positions'].keys()) for a, b in
                                new_subgraph.edges_iter() if 'positions' in new_subgraph.edge[a][b]}
            assert len(source_sequences) == 1, (source_sequences, len(new_subgraph))

    def finish_build(self, graphviz=False):
        """
        Finishes building the graph.

        If graphviz is True, adds a label tag to each sequence edge to improve understandability in graphviz
        """
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
                        self.edge[l][r]["label"] = " - ".join(sorted(
                            [": ".join([y, ", ".join([str(x) for x in self.edge[l][r]['positions'][y]])]) for y in
                             self.edge[l][r]['positions']]))
                    self.edge[l][r]["color"] = "purple"
                elif 'source' in self.edge[l][r]:
                    self.edge[l][r]['color'] = "blue"
                else:
                    self.edge[l][r]['color'] = "green"
        self.paralogs = sorted(self.paralogs, key=lambda x: x[0])
        self.is_built = True

    def connected_component_iter(self):
        """
        Yields connected components. Internal is used for pruning edges within this object and should not be set to
        True outside of this setup.
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
