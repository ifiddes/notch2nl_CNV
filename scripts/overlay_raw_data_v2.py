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

color_palette = [( 93, 165, 218),  # m blue
                 (250, 164,  58),  # m orange
                 ( 96, 189, 104),  # m green
                 (241, 124, 167),  # m red
                 (178, 145,  47),  # m brown
                 (178, 118, 178),  # m purple
                 (241,  88,  84),  # m magenta
                 ( 77,  77,  77),  # m grey
                 (222, 207,  63)   # m yellow
                ]  

# put on a 0-1 scale
color_palette = np.asarray(color_palette) / 255.0


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
    fig, plots = plt.subplots(len(mole_graph.paralogs), sharey=True, figsize=(12.0, 5.0))
    individual_patch = matplotlib.patches.Patch(color=color_palette[0], fill="true")
    mole_patch = matplotlib.patches.Patch(color=color_palette[1], fill="true")
    plt.figlegend((individual_patch, mole_patch), ("Individual", "Mole"), loc='upper right', ncol=2)
    max_gap = max(stop - start for start, stop in mole_graph.paralogs.itervalues())
    for i, (p, para) in enumerate(izip(plots, mole_graph.paralogs.iterkeys())):
        start, stop = mole_graph.paralogs[para]
        rounded_max_gap = int(math.ceil(1.0 * max_gap / 10000) * 10000)
        p.axis([start, start + rounded_max_gap, 0, 4])
        x_ticks = [start] + range(start + 20000, start + rounded_max_gap + 20000, 20000)
        p.axes.set_yticks(range(5))
        p.axes.set_yticklabels(map(str, range(5)), fontsize=9)
        p.axes.set_xticks(x_ticks)
        p.axes.set_xticklabels(["{:.3e}".format(start)] + [str(20000 * x) for x in xrange(1, len(x_ticks))])
        starts, stops, vals = zip(*raw_dict[para])
        p.hlines(vals, starts, stops, color=color_palette[1], alpha=0.8, linewidth=2.0)
        if len(sun_results[para_map[para]]) > 0:
            sun_pos, sun_vals = zip(*sun_results[para_map[para]])
            p.vlines(np.asarray(sun_pos), np.zeros(len(sun_pos)), sun_vals, color="#E83535", linewidth=0.5, alpha=0.5)
        for x in range(1, 4):
            p.axhline(y=x, ls="--", lw=0.5)
        p.set_title("{}".format(para))
    for i, (p, para) in enumerate(izip(plots, individual_graph.paralogs.iterkeys())):
        starts, stops, vals = zip(*individual_raw_dict[para])
        p.hlines(vals, starts, stops, color=color_palette[0], alpha=0.8, linewidth=2.0)
        p.set_title("{}".format(para))    
    fig.subplots_adjust(hspace=0.8)
    plt.savefig(out_path, format="png", dpi=500)
    plt.close()


def main():
    args = parse_args()
    individual_graph = pickle.load(open(args.individual_graph))
    mole_graph = pickle.load(open(args.mole_graph))
    individual_raw_data = pickle.load(open(args.individual_raw_data))
    mole_raw_data = pickle.load(open(args.mole_raw_data))
    sun_results = pickle.load(open(args.sun_results))
    path = os.path.join(args.out_dir, args.uuid, args.uuid + ".OverlaidRawData.V2.png")
    overlay_raw_data(mole_raw_data, individual_raw_data, mole_graph, individual_graph, sun_results, args.uuid, path)


if __name__ == "__main__":
    main()