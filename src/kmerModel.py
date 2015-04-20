from src.kmerIlpModel import KmerIlpModel
from src.unitigGraph import UnitigGraph
from jobTree.src.bioio import fastaRead, fastqRead
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
        add_individual_to_graph(graph, paths.fastq)
        graph.resolve_bubbles()



def add_individual_to_graph(graph, fastq):
    for name, seq, qual in fastqRead(fastq):
        graph.construct_individual_nodes(name, seq)
        graph.construct_adjacencies(seq, source_seq=False)


def add_mole_to_graph(graph, mole_seq):
    for name, seq in fastaRead(mole_seq):
        name, offset = name.split("_")
        graph.construct_ref_nodes(name, offset, seq)
        graph.construct_adjacencies(seq, source_seq=True)