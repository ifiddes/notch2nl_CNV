#fastq_path = "/cluster/home/ifiddes/ifiddes_hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
#fastq_path = "/home/ifiddes/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
fastq_path = "/Users/ifiddes/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
masked_ref_path = "/Users/ifiddes/notch2nl_CNV/data/kmer_model_data/notch2nl_masked_hg38.fa"
unmasked_ref_path = "/Users/ifiddes/notch2nl_CNV/data/kmer_model_data/notch2nl_unmasked_hg38.fa"

from src.kmerModel import *


def main():
    graph = UnitigGraph(49)
    add_mole_to_graph(graph, unmasked_ref_path, masked_ref_path)

if __name__ == "__main__":
    main()