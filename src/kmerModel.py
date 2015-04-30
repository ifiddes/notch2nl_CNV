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
        add_mole_to_graph(graph, self.paths.unmasked_ref, self.paths.masked_ref)
        if self.add_individual is True:
            add_individual_to_graph(graph, self.k_plus1_mer_counts_path)
            pickle.dump(graph, "{}_test_graph.pickle".format(self.uuid))
        #graph.flag_nodes(open(self.paths.bad_kmers))
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
        raw_counts = ilp_model.report_normalized_raw_data_map()
        if self.add_individual is True:
            out_path = os.path.join(self.paths.out_dir, self.uuid, self.uuid + ".Individual.png")
        else:
            out_path = os.path.join(self.paths.out_dir, self.uuid, self.uuid + ".png")
        combined_plot(result_dict, raw_counts, graph.paralogs, self.sun_results, self.uuid, out_path)


def get_normalizing_kmers(normalizing_path, kmer_size):
    normalizing_kmers = set()
    for _, seq in fasta_read(normalizing_path):
        seq = seq.upper()
        for i in xrange(len(seq) - kmer_size + 1):
            k = seq[i:i + kmer_size]
            if "N" not in k:
                normalizing_kmers.add(canonical(k))
    return normalizing_kmers


def get_kmer_counts(graph, normalizing_kmers, kmer_counts_file):
        data_counts = {}
        normalizing = 0
        for count, seq in count_reader(kmer_counts_file):
            if seq in graph.kmers:
                data_counts[seq] = count
            elif seq in normalizing_kmers:
                normalizing += count
        normalizing /= (1.0 * len(normalizing_kmers))
        return data_counts, normalizing


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


def combined_plot(ilpDict, rawCounts, offsetMap, unfilteredSunDict, uuid, outPath):
    """
    Generates a final combined plot overlaying both ILP and SUN results.
    """
    colors = ["#9b59b6", "#3498db", "#e74c3c", "#34495e", "#2ecc71"]
    explodedRawCounts = explodeResultDict(rawCounts, offsetMap)
    explodedData = explodeResultDict(ilpDict, offsetMap)
    # used because the SUN model uses single letter labels
    paraMap = {"Notch2NL-A": "A", "Notch2NL-B": "B", "Notch2NL-C": "C", "Notch2NL-D": "D", "Notch2": "N"}
    sortedParalogs = ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C", "Notch2NL-D", "Notch2"]
    fig, plots = plt.subplots(len(sortedParalogs), sharey=True)
    plt.yticks((0, 1, 2, 3, 4))
    plt.suptitle("kmer-DeBruijn ILP and SUN results")
    maxGap = max(len(x) for x in explodedRawCounts.itervalues())
    for i, (p, para) in enumerate(izip(plots, sortedParalogs)):
        data = explodedData[para]
        rawData = explodedRawCounts[para]
        windowedRawData = [sum(rawData[k:k+300])/300 for k in xrange(0, len(rawData), 300)]
        start = offsetMap[para]
        stop = start + len(data)
        rounded_max_gap = int(math.ceil(1.0 * maxGap / 10000) * 10000)
        p.axis([start, start + rounded_max_gap, 0, 4])
        x_ticks = [start] + range(start + 20000, start + rounded_max_gap, 20000) + [stop]
        p.axes.set_xticks(x_ticks)
        p.axes.set_xticklabels([str(start)] + [str(20000 * x) for x in xrange(1, len(x_ticks) - 1)] + [stop])
        p.fill_between(range(start, stop), data, color=colors[i], alpha=0.8)
        p.plot(range(start, stop, 300), windowedRawData, alpha=0.8, linewidth=1.5)
        if len(unfilteredSunDict[paraMap[para]]) > 0:
            sunPos, sunVals = zip(*unfilteredSunDict[paraMap[para]])
            p.vlines(np.asarray(sunPos), np.zeros(len(sunPos)), sunVals, color="#E83535")
        for i in range(1, 4):
            p.axhline(y=i, ls="--", lw=0.8)
        p.set_title("{}".format(para))
    fig.subplots_adjust(hspace=0.8)
    plt.savefig(outPath, format="png")
    plt.close()


def explodeResultDict(rd, offsetMap):
    r = defaultdict(list)
    for para in rd:
        prevStart = offsetMap[para]
        for start, span, val in rd[para]:
            r[para].extend([val] * (start - prevStart))
            prevStart = start
    return r