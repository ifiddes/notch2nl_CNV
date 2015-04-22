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
        graph.prune_source_edges()
        graph.prune_individual_edges()


def add_individual_to_graph(graph, k1mer_counts):
    for count, seq in count_reader(k1mer_counts):
        #if count > 1:
        graph.add_individual_sequences(seq)


def add_mole_to_graph(graph, unmasked_mole, masked_mole):
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
        graph.add_source_sequences(name, offset, masked_seq, unmasked_seq)