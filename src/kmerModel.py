from itertools import izip
import cPickle as pickle
from src.kmerIlpModel import KmerIlpModel
from src.unitigGraph import UnitigGraph
from src.helperFunctions import count_reader, canonical
from jobTree.src.bioio import fastaRead as fasta_read
from jobTree.scriptTree.target import Target


class KmerModel(Target):
    def __init__(self, paths, ilp_config, inferred_c, inferred_d):
        Target.__init__(self)
        self.paths = paths
        self.ilp_config = ilp_config
        self.inferred_c = inferred_c
        self.inferred_d = inferred_d

    def run(self):
        graph = UnitigGraph(kmer_size=self.ilp_config.kmer_size, derived=False)
        add_mole_to_graph(graph, self.paths.unmasked_mole, self.paths.masked_mole)
        add_individual_to_graph(graph, self.paths.k_plus1_mer_counts)
        if self.paths.bad_kmers is not None:
            graph.flag_nodes(open(self.paths.bad_kmers))
        if self.paths.graph_out_path is not None:
            pickle.dump(graph, self.paths.graph_out_path)
        normalizing_kmers = get_normalizing_kmers(self.paths.normalizing)
        assert len(graph.kmers & graph.normalizing_kmers) == 0
        data_counts = {}
        normalizing = 0
        for count, seq in count_reader(self.paths.kmer_counts):
            if seq in graph.kmers:
                data_counts[seq] = count
            elif seq in normalizing_kmers:
                normalizing += count
        normalizing /= (1.0 * len(normalizing_kmers))
        ilp_model = KmerIlpModel(graph, normalizing, breakpoint_penalty=self.ilp_config.breakpoint_penalty,
                                 data_penalty=self.ilp_config.data_penalty,
                                 trash_penalty=self.ilp_config.trash_penalty,
                                 expected_value_penalty=self.ilp_config.expected_value_penalty,
                                 infer_c=self.inferred_c, infer_d=self.inferred_d, default_ploidy=2)


def get_normalizing_kmers(normalizing_path):
    return {canonical(seq) for name, seq in fasta_read(normalizing_path)}


def add_individual_to_graph(graph, k1mer_counts):
    for count, seq in count_reader(k1mer_counts):
        graph.add_individual_sequence(seq)
    graph.prune_individual_edges()


def add_mole_to_graph(graph, unmasked_mole, masked_mole):
    """
    Loops over the input sequences twice, first time builds a list of masked kmers, second time adds unmasked kmers
    to the graph.
    """
    seq_map = {}
    masked_seqs = {name: seq for name, seq in fasta_read(masked_mole)}
    unmasked_seqs = {name: seq for name, seq in fasta_read(unmasked_mole)}
    for name in masked_seqs:
        masked_seq = masked_seqs[name]
        unmasked_seq = unmasked_seqs[name]
        assert "N" not in unmasked_seq
        assert "_" in name and len(name.split("_")) == 2
        name, offset = name.split("_")
        try:
            offset = int(offset)
        except TypeError:
            raise RuntimeError("Naming convention for input reference fasta is not right. >Sequence_Offset")
        masked_seq = masked_seq.upper()
        unmasked_seq = unmasked_seq.upper()
        seq_map[(name, offset)] = [masked_seq, unmasked_seq]
        graph.add_masked_kmers(masked_seq, unmasked_seq)
    for (name, offset), (masked_seq, unmasked_seq) in seq_map.iteritems():
         graph.add_source_sequence(name, offset, masked_seq, unmasked_seq)
    graph.prune_source_edges()