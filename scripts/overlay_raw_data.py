import os
import math
import argparse
import cPickle as pickle
from itertools import izip

import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


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


def overlay_raw_data(raw_dict, individual_raw_dict, mole_graph, individual_graph, sun_results, out_path):
    color_palette = sns.color_palette()
    # used because the SUN model uses single letter labels
    para_map = {"Notch2NL-A": "A", "Notch2NL-B": "B", "Notch2NL-C": "C", "Notch2NL-D": "D", "Notch2": "N"}
    fig, plots = plt.subplots(len(mole_graph.paralogs), sharey=True, figsize=(10.0, 5.0))
    individual_patch = matplotlib.patches.Patch(color=color_palette[0], fill="true")
    mole_patch = matplotlib.patches.Patch(color=color_palette[1], fill="true")
    plt.figlegend((individual_patch, mole_patch), ("Individual", "Mole"), loc='upper right', ncol=2)
    max_gap = max(stop - start for start, stop in mole_graph.paralogs.itervalues())
    for i, (p, para) in enumerate(izip(plots, mole_graph.paralogs.iterkeys())):
        positions, vals, raw_vals = zip(*raw_dict[para])
        windowed_raw_data = [1.0 * sum(raw_vals[k:k + 200]) / 200 for k in xrange(0, len(raw_vals) - 200, 200)]
        start, stop = mole_graph.paralogs[para]
        windowed_positions = [int(round(start + 1.0 * sum(positions[k:k + 200]) / 200))
                              for k in xrange(0, len(positions) - 200, 200)]
        rounded_max_gap = int(math.ceil(1.0 * max_gap / 10000) * 10000)
        p.axis([start, start + rounded_max_gap, 0, 4])
        x_ticks = [start] + range(start + 20000, start + rounded_max_gap + 20000, 20000)
        p.axes.set_yticks(range(5))
        p.axes.set_yticklabels(map(str, range(5)), fontsize=9)
        p.axes.set_xticks(x_ticks)
        p.axes.set_xticklabels(["{:.3e}".format(start)] + [str(20000 * x) for x in xrange(1, len(x_ticks))])
        #p.plot(windowed_positions, windowed_raw_data, alpha=0.8, color=color_palette[1], linewidth=1)
        if len(sun_results[para_map[para]]) > 0:
            sun_pos, sun_vals = zip(*sun_results[para_map[para]])
            p.vlines(np.asarray(sun_pos), np.zeros(len(sun_pos)), sun_vals, color="#E83535", linewidth=0.9, alpha=0.5)
        p.set_title("{}".format(para))
    for i, (p, para) in enumerate(izip(plots, individual_graph.paralogs.iterkeys())):
        positions, vals, raw_vals = zip(*individual_raw_dict[para])
        individual_windowed_raw_data = [1.0 * sum(raw_vals[k:k + 200]) / 200 for k in xrange(0, len(raw_vals) - 200,
                                                                                             200)]
        start, stop = individual_graph.paralogs[para]
        windowed_positions = [int(round(start + 1.0 * sum(positions[k:k + 200]) / 200))
                              for k in xrange(0, len(positions) - 200, 200)]
        p.plot(windowed_positions, individual_windowed_raw_data, alpha=0.8, color=color_palette[0], linewidth=1)
        p.set_title("{}".format(para))
        if para in ["Notch2NL-C", "Notch2NL-D"]:
            p.invert_xaxis()
    fig.subplots_adjust(hspace=0.8)
    plt.savefig(out_path, format="png", dpi=300)
    plt.close()


def main():
    args = parse_args()
    individual_graph = pickle.load(open(args.individual_graph))
    mole_graph = pickle.load(open(args.mole_graph))
    individual_raw_data = pickle.load(open(args.individual_raw_data))
    mole_raw_data = pickle.load(open(args.mole_raw_data))
    sun_results = pickle.load(open(args.sun_results))
    path = os.path.join(args.out_dir, args.uuid, args.uuid + ".OverlaidRawData.png")
    overlay_raw_data(mole_raw_data, individual_raw_data, mole_graph, individual_graph, sun_results, path)


if __name__ == "__main__":
    main()