import os
import math
import argparse
import cPickle as pickle
from collections import defaultdict
from itertools import izip

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from src.unitigGraph import UnitigGraph


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--individual_graph", help="Pickled individual graph", required=True)
    parser.add_argument("--mole_graph", help="Pickled mole graph", required=True)
    parser.add_argument("--individual_raw_data", help="Individual raw data", required=True)
    parser.add_argument("--mole_raw_data", help="Mole raw data", required=True)
    parser.add_argument("--sun_results", help="Pickled SUN results", required=True)
    parser.add_argument("--uuid", required=True)
    parser.add_argument("--out_dir", required=True)
    return parser.parse_args()


def overlay_raw_data(raw_dict, individual_raw_dict, mole_graph, individual_graph, sun_results, uuid, out_path):
    # used because the SUN model uses single letter labels
    para_map = {"Notch2NL-A": "A", "Notch2NL-B": "B", "Notch2NL-C": "C", "Notch2NL-D": "D", "Notch2": "N"}
    fig, plots = plt.subplots(len(mole_graph.paralogs), sharey=True)
    plt.yticks((0, 1, 2, 3, 4))
    plt.suptitle("Overlaid Raw Data for {} individualized and mole mole_graphs".format(uuid))
    max_gap = max(stop - start for start, stop in mole_graph.paralogs.itervalues())
    for i, (p, para) in enumerate(izip(plots, mole_graph.paralogs.iterkeys())):
        raw_data = explode_result(raw_dict[para], mole_graph.paralogs[para])
        windowed_raw_data = [1.0 * sum(raw_data[k:k + 300]) / 300 for k in xrange(0, len(raw_data), 300)]
        start, stop = mole_graph.paralogs[para]
        rounded_max_gap = int(math.ceil(1.0 * max_gap / 10000) * 10000)
        p.axis([start, start + rounded_max_gap, 0, 4])
        x_ticks = [start] + range(start + 20000, start + rounded_max_gap + 20000, 20000)
        p.axes.set_xticks(x_ticks)
        p.axes.set_xticklabels(["{:.3e}".format(start)] + [str(20000 * x) for x in xrange(1, len(x_ticks))])
        p.plot(range(start, stop, 300), windowed_raw_data, alpha=0.8, linewidth=1.5)
        if len(sun_results[para_map[para]]) > 0:
            sun_pos, sun_vals = zip(*sun_results[para_map[para]])
            p.vlines(np.asarray(sun_pos), np.zeros(len(sun_pos)), sun_vals, color="#E83535")
        for i in range(1, 4):
            p.axhline(y=i, ls="--", lw=0.7)
        p.set_title("{}".format(para))
    for i, (p, para) in enumerate(izip(plots, individual_graph.paralogs.iterkeys())):
        raw_data = explode_result(raw_dict[para], individual_graph.paralogs[para])
        windowed_raw_data = [1.0 * sum(raw_data[k:k + 300]) / 300 for k in xrange(0, len(raw_data), 300)]
        start, stop = individual_graph.paralogs[para]
        p.plot(range(start, stop, 300), windowed_raw_data, alpha=0.8, linewidth=1.5)
        p.set_title("{}".format(para))    
    fig.subplots_adjust(hspace=0.8)
    plt.savefig(out_path, format="png")
    plt.close()    


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


def main():
    args = parse_args()
    individual_graph = pickle.load(open(args.individual_graph))
    mole_graph = pickle.load(open(args.mole_graph))
    individual_raw_data = pickle.load(open(args.individual_raw_data))
    mole_raw_data = pickle.load(open(args.mole_raw_data))
    sun_results = pickle.load(open(args.sun_results))
    path = os.path.join(args.out_dir, args.uuid, args.uuid + ".OverlaidRawData.png")
    overlay_raw_data(mole_raw_data, individual_raw_data, mole_graph, individual_graph, sun_results, args.uuid, path)


if __name__ == "__main__":
    main()