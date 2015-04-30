import cPickle as pickle
from src.kmerModel import *
from src.kmerIlpModel import *
graph = pickle.load(open("b88_graph.pickle"))
ilp_model = KmerIlpModel(graph)
normalizing_kmers = get_normalizing_kmers("data/kmer_model_data/normalizing.fa", 49)
data_counts, normalizing = get_kmer_counts(graph, normalizing_kmers, "output/b88da23c/b88da23c.49mer.fa")
ilp_model.introduce_data(data_counts, normalizing)
ilp_model.solve()
sun_results = pickle.load(open("output/b88da23c/sun_results.pickle"))
result_dict = ilp_model.report_copy_map()
raw_counts = ilp_model.report_normalized_raw_data_map()
combined_plot(result_dict, raw_counts, graph.paralogs, sun_results, "b88", "./")