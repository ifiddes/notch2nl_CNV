import os
import math
import cPickle as pickle
from collections import defaultdict
from itertools import izip

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from src.kmerIlpModel import KmerIlpModel
from src.unitigGraph import UnitigGraph
from src.helperFunctions import count_reader, canonical
from jobTree.src.bioio import fastaRead as fasta_read
from jobTree.scriptTree.target import Target


class KmerModel(Target):
    """
    This class acts as a wrapper for the KmerIlpModel. Given the results from the SUN model, and the jellyfish counts
    containing both k sized kmer counts and k+1 sized kmer counts, this target creates a unitig graph, individualizes
    it, then passes this graph off to KmerIlpModel for copy number inference.
    """

    def __init__(self, paths, uuid, ilp_config, sun_results, fastq_path, kmer_counts_path, k_plus1_mer_counts_path,
                 inferred_c, inferred_d, add_individual):
        Target.__init__(self)
        self.paths = paths
        self.uuid = uuid
        self.ilp_config = ilp_config
        self.sun_results = sun_results
        self.fastq_path = fastq_path
        self.kmer_counts_path = kmer_counts_path
        self.k_plus1_mer_counts_path = k_plus1_mer_counts_path
        self.inferred_c = inferred_c
        self.inferred_d = inferred_d
        self.add_individual = add_individual

    def run(self):
        graph = UnitigGraph(kmer_size=self.ilp_config.kmer_size, derived=False)
        add_mole_to_graph(graph, self.paths.ilp_ref)
        if self.add_individual is True:
            add_individual_to_graph(graph, self.k_plus1_mer_counts_path)
            pickle.dump(graph, open(
                os.path.join(self.paths.out_dir, self.uuid, "{}_individual_graph.pickle".format(self.uuid)), "w"))
        else:
            pickle.dump(graph,
                        open(os.path.join(self.paths.out_dir, self.uuid, "{}_graph.pickle".format(self.uuid)), "w"))
        # graph.flag_nodes(open(self.paths.bad_kmers))
        normalizing_kmers = get_normalizing_kmers(self.paths.normalizing, self.ilp_config.kmer_size)
        try:
            assert len(graph.kmers & normalizing_kmers) == 0
        except AssertionError:
            raise RuntimeError("Normalizing kmer file contains kmers in the input reference files.")
        data_counts, normalizing = get_kmer_counts(graph, normalizing_kmers, self.kmer_counts_path)
        ilp_model = KmerIlpModel(graph, breakpoint_penalty=self.ilp_config.breakpoint_penalty,
                                 data_penalty=self.ilp_config.data_penalty,
                                 trash_penalty=self.ilp_config.trash_penalty,
                                 expected_value_penalty=self.ilp_config.expected_value_penalty,
                                 infer_c=self.inferred_c, infer_d=self.inferred_d, default_ploidy=2)
        ilp_model.introduce_data(data_counts, normalizing)
        ilp_model.solve()
        result_dict = ilp_model.report_copy_map()
        raw_dict = ilp_model.report_normalized_raw_data_map()
        if self.add_individual is True:
            out_png_path = os.path.join(self.paths.out_dir, self.uuid, self.uuid + ".Individual.png")
            out_raw_path = os.path.join(self.paths.out_dir, self.uuid, "tracks",
                                        "{}.Individual.RawCounts.hg38.wiggle".format(self.uuid))
            out_ilp_path = os.path.join(self.paths.out_dir, self.uuid, "tracks",
                                        "{}.Individual.ILP.hg38.wiggle".format(self.uuid))
            # debugging code
            r_d_path = os.path.join(self.paths.out_dir, self.uuid, self.uuid + ".RawData.pickle")
            pickle.dump(raw_dict, open(r_d_path, "w"))
        else:
            out_png_path = os.path.join(self.paths.out_dir, self.uuid, self.uuid + ".png")
            out_raw_path = os.path.join(self.paths.out_dir, self.uuid, "tracks",
                                        "{}.RawCounts.hg38.wiggle".format(self.uuid))
            out_ilp_path = os.path.join(self.paths.out_dir, self.uuid, "tracks",
                                        "{}.ILP.hg38.wiggle".format(self.uuid))
            # debugging code
            r_d_path = os.path.join(self.paths.out_dir, self.uuid, self.uuid + ".RawData.pickle")
            pickle.dump(raw_dict, open(r_d_path, "w"))
        combined_plot(result_dict, raw_dict, graph, self.sun_results, self.uuid, out_png_path)
        generate_wiggle_plots(result_dict, raw_dict, self.uuid, out_raw_path, out_ilp_path)


def get_normalizing_kmers(normalizing_path, kmer_size):
    """
    Takes a fasta file containing sequence that will be used to normalize the copy number inference. Generates a set of
    kmers of size kmer_size for this.
    """
    normalizing_kmers = set()
    for _, seq in fasta_read(normalizing_path):
        seq = seq.upper()
        for i in xrange(len(seq) - kmer_size + 1):
            k = seq[i:i + kmer_size]
            if "N" not in k:
                normalizing_kmers.add(canonical(k))
    return normalizing_kmers


def get_kmer_counts(graph, normalizing_kmers, kmer_counts_file):
    """
    Parses a kmer_count file from jellyfish.
    """
    data_counts = {}
    normalizing = 0
    for count, seq in count_reader(kmer_counts_file):
        # just in case...
        assert seq == canonical(seq)
        if seq in graph.kmers:
            data_counts[seq] = count
        elif seq in normalizing_kmers:
            normalizing += count
    normalizing /= (1.0 * len(normalizing_kmers))
    return data_counts, normalizing


def add_individual_to_graph(graph, k1mer_counts):
    """
    Adds individual sequence to the graph using k+1 kmers.
    """
    for count, seq in count_reader(k1mer_counts):
        graph.add_individual_sequence(seq)
    graph.prune_individual_edges()


def add_mole_to_graph(graph, mole):
    """
    Loops over the input sequences twice, first time builds a list of masked kmers, second time adds unmasked kmers
    to the graph.
    """
    for name, seq in fasta_read(mole):
        assert "N" not in seq
        assert "_" in name and len(name.split("_")) == 2
        name, offset = name.split("_")
        try:
            offset = int(offset)
        except TypeError:
            raise RuntimeError("Naming convention for input reference fasta is not right. >Sequence_Offset")
        seq = seq.upper()
        graph.add_source_sequence(name, offset, seq)
    graph.prune_source_edges()


def combined_plot(result_dict, raw_dict, graph, sun_results, uuid, out_path):
    """
    Generates a final combined plot overlaying both ILP and SUN results.
    """
    colors = ["#9b59b6", "#3498db", "#e74c3c", "#34495e", "#2ecc71"]
    # used because the SUN model uses single letter labels
    para_map = {"Notch2NL-A": "A", "Notch2NL-B": "B", "Notch2NL-C": "C", "Notch2NL-D": "D", "Notch2": "N"}
    fig, plots = plt.subplots(len(graph.paralogs), sharey=True)
    plt.yticks((0, 1, 2, 3, 4))
    plt.suptitle("Unitig ILP and SUN results for {}".format(uuid))
    max_gap = max(stop - start for start, stop in graph.paralogs.itervalues())
    for i, (p, para) in enumerate(izip(plots, graph.paralogs.iterkeys())):
        data = explode_result(result_dict[para], graph.paralogs[para])
        raw_data = explode_result(raw_dict[para], graph.paralogs[para])
        windowed_raw_data = [1.0 * sum(raw_data[k:k + 300]) / 300 for k in xrange(0, len(raw_data), 300)]
        start, stop = graph.paralogs[para]
        rounded_max_gap = int(math.ceil(1.0 * max_gap / 10000) * 10000)
        p.axis([start, start + rounded_max_gap, 0, 4])
        x_ticks = [start] + range(start + 20000, start + rounded_max_gap + 20000, 20000)
        p.axes.set_xticks(x_ticks)
        p.axes.set_xticklabels(["{:.3e}".format(start)] + [str(20000 * x) for x in xrange(1, len(x_ticks))])
        p.fill_between(range(start, stop), data, color=colors[i], alpha=0.8)
        p.plot(range(start, stop, 300), windowed_raw_data, alpha=0.8, linewidth=1.5)
        if len(sun_results[para_map[para]]) > 0:
            sun_pos, sun_vals = zip(*sun_results[para_map[para]])
            p.vlines(np.asarray(sun_pos), np.zeros(len(sun_pos)), sun_vals, color="#E83535")
        for i in range(1, 4):
            p.axhline(y=i, ls="--", lw=0.7)
        p.set_title("{}".format(para))
    fig.subplots_adjust(hspace=0.8)
    plt.savefig(out_path, format="png")
    plt.close()


def generate_wiggle_plots(result_dict, raw_dict, uuid, out_raw_path, out_ilp_path):
    """
    Generates wiggle plots for the UCSC browser on hg38
    """
    with open(out_ilp_path, "w") as outf:
        #outf.write(
        #    "track type=wiggle_0 name={} color=35,125,191 autoScale=off visibility=full alwaysZero=on "
        #    "yLineMark=2 viewLimits=0:4 yLineOnOff=on maxHeightPixels=100:75:50\n".format(uuid))
        for para in result_dict:
            for start, stop, val in result_dict[para]:
                outf.write("variableStep chrom=chr1 span={}\n".format(stop - start))
                outf.write("{} {}\n".format(start, val))
    with open(out_raw_path, "w") as outf:
        #outf.write(
        #    "track type=wiggle_0 name={} color=35,125,191 autoScale=off visibility=full alwaysZero=on "
        #    "yLineMark=2 viewLimits=0:4 yLineOnOff=on maxHeightPixels=100:75:50\n".format(uuid))
        for para in raw_dict:
            for start, stop, val in raw_dict[para]:
                outf.write("variableStep chrom=chr1 span={}\n".format(stop - start))
                outf.write("{} {}\n".format(start, val))


def explode_result(result, graph_positions):
    r = []
    graph_start, graph_stop = graph_positions
    prev_start = result[0][0]
    r.extend([result[0][2]] * (prev_start - graph_start))
    for start, stop, val in result:
        r.extend([val] * (start - prev_start))
        prev_start = start
    r.extend([val] * (stop - start))
    r.extend([val] * (graph_stop - stop))
    return r