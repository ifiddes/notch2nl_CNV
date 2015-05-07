from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
from src.kmerIlpModel import *
import cPickle as pickle
#fastq_path = "/hive/users/ifiddes/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
#fastq_path = "/Users/ifiddes/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
fastq_path = "/home/ifiddes/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
unmasked_ref_path = "data/kmer_model_data/unmasked_last_2500bp.fa"
masked_ref_path = "data/kmer_model_data/masked_last_2500bp.fa"
graph = UnitigGraph(kmer_size=49, derived=False)
add_mole_to_graph(graph, unmasked_ref_path, masked_ref_path)
#add_individual_to_graph(graph, fastq_path)
ilp_model = KmerIlpModel(graph)
normalizing_kmers = get_normalizing_kmers("data/kmer_model_data/normalizing.fa", 49)
data_counts, normalizing = get_kmer_counts(graph, normalizing_kmers, "output/mike_sny/mike_sny.49mer.fa")
ilp_model.introduce_data(data_counts, normalizing)
ilp_model.solve()
result_dict = ilp_model.report_copy_map()
raw_dict = ilp_model.report_normalized_raw_data_map()
uuid = "test"
out_raw_path = "test.wiggle"
with open(out_raw_path, "w") as outf:
    for para in raw_dict:
        for start, stop, val in raw_dict[para]:
            outf.write("variableStep chrom=chr1 span={}\n".format(stop - start))
            outf.write("{} {}\n".format(start, val))