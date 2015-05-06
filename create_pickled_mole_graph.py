from src.kmerModel import *
import cPickle as pickle
fastq_path = "/hive/users/ifiddes//notch_mike_snyder/snyder_notch.50mer.Counts.fa"
masked_ref_path = "data/kmer_model_data/notch2nl_masked_hg38.fa"
unmasked_ref_path = "data/kmer_model_data/notch2nl_unmasked_hg38.fa"
graph = UnitigGraph(49)
add_mole_to_graph(graph, unmasked_ref_path, masked_ref_path)
graph.prune_source_edges()
pickle.dump(graph, open("data/kmer_model_data/mole_graph_pruned.pickle","w"))
