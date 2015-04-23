from itertools import izip
from src.kmerIlpModel import KmerIlpModel
from src.unitigGraph import UnitigGraph
from src.helperFunctions import count_reader
from jobTree.src.bioio import fastaRead as fasta_read
from jobTree.scriptTree.target import Target


class KmerModel(Target):
    def __init__(self, paths, inferred_c, inferred_d):
        Target.__init__(self)
        self.paths = paths
        self.inferred_c = inferred_c
        self.inferred_d = inferred_d

    def run(self):
        paths = self.paths
        graph = UnitigGraph(paths.kmer_size)
        graph.add_normalizing(paths.normalizing)
        add_mole_to_graph(graph, paths.mole_seq)
        add_individual_to_graph(graph, paths.k1mer_counts)
        graph.prune_edges()


def add_individual_to_graph(graph, k1mer_counts):
    for count, seq in count_reader(k1mer_counts):
        graph.add_individual_sequence(seq)


def add_mole_to_graph(graph, unmasked_mole, masked_mole):
    """
    Loops over the input sequences twice, first time builds a list of masked kmers, second time adds unmasked kmers
    to the graph.
    """
    seq_map = {}
    for (unmasked_name, unmasked_seq), (masked_name, masked_seq) in izip(fasta_read(unmasked_mole),
                                                                         fasta_read(masked_mole)):
        assert unmasked_name == masked_name, (unmasked_name, masked_name)
        assert len(unmasked_seq) == len(masked_seq)
        assert "N" not in unmasked_seq
        assert "_" in unmasked_name and len(unmasked_name.split("_")) == 2
        name, offset = unmasked_name.split("_")
        try:
            offset = int(offset)
        except TypeError:
            raise TypeError("Naming convention for input reference fasta is not right. >Sequence_Offset")
        masked_seq = masked_seq.upper()
        unmasked_seq = unmasked_seq.upper()
        seq_map[(name, offset)] = [masked_seq, unmasked_seq]
        graph.add_masked_kmers(masked_seq, unmasked_seq)
        for (name, offset), (masked_seq, unmasked_seq) in seq_map.iteritems():
            graph.add_source_sequence(name, offset, masked_seq)