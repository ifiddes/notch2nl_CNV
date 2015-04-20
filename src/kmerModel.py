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
        graph.construct_individual_nodes(seq)
        graph.construct_adjacencies(seq, source_seq=False)


def add_mole_to_graph(graph, mole_seq):
    for name, seq in fasta_read(mole_seq):
        name, offset = name.split("_")
        graph.construct_ref_nodes(name, offset, seq)
        graph.construct_adjacencies(seq, source_seq=True)