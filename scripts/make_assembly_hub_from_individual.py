#!/usr/bin/env python2.7
"""
Takes the output folder for one individual and creates a trackHub out of that individual.
"""

import sys
import os
import argparse
import subprocess

colors = {"ILP": "242,148,176", "SUN": "191,191,191", "RawCounts": "137,137,137", "RawIndividual": "164,239,245"}
paths = namedtuple("paths", "individual_ilp, individual_raw, raw, sun")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", "-d", required=True, help="Dir to dig through. If does not end in tracks, will append")
    parser.add_argument("--out_dir", "-o", required=True, help="Dir to put assembly hub into.")
    parser.add_argument("--name", required=True, help="assembly hub name")
    parser.add_argument("--uuid", required=True, help='UUID')
    parser.add_argument("--hg38_chrom_sizes", help="hg38 chrom.sizes file", default="data/hg38.chrom.sizes")
    return parser.parse_args()


def start_hub(dir_path, name):
    with open(os.path.join(dir_path, "hub.txt"), "w") as outf:
        hub_str = "hub {}\nshortLabel {}\nlongLabel {}\ngenomesFile genomes.txt\nemail ian.t.fiddes@gmail.com\n"
        outf.write(hub_str.format(name, name + " Notch2NL", name + " Notch2NL"))
    with open(os.path.join(dir_path, "genomes.txt"), "w") as outf:
        outf.write("genome hg38\ntrackDb hg38/trackDb.txt\ndefaultPos chr1:145987622-149723055\n\n")
    if not os.path.exists(os.path.join(dir_path, "hg38")):
        os.mkdir(os.path.join(dir_path, "hg38"))
    with open(os.path.join(dir_path, "hg38", "Notch2NL.html"), "w") as outf:
        outf.write("Notch2NL {}\n".format(name))


def track_db(out_dir, uuid):
    """
    This generates three types of tracks.
    1) all: overlays individual raw data, raw data and SUN
    2) raw: overlays individual raw and raw
    3) ILP: overlays individual ILP result with individual raw data and SUN information
    """
    container = ("track {0}_{1}\ncontainer multiWig\nshortLabel {0} all\nlongLabel Individual Raw/Raw/SUN Combined\n"
                 "type bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\n"
                  "yLineOnOff on\nmaxHeightPixels 100:75:50\n")
    all = ("track {}_ind_raw_1\nshortLabel {}_ind_raw\nlongLabel {}_ind_raw\ntype bigWig 0 4\nparent {}_all\n"
            "bigDataUrl {}.Individual.RawCounts.hg38.bw\n")
    all += ("track {}_raw_1\nshortLabel {}_raw\nlongLabel {}_raw\ntype bigWig 0 4\nparent {}_all"
            "bigDataUrl {}.Individual.RawCounts.hg38.bw\n")
    all += ("track {}_sun_1\nshortLabel {}_sun\nlongLabel {}_sun\ntype bigWig 0 4\nparent {}_all"
            "bigDataUrl {}.sun_model.hg38.bedGraph\n")

    raw = ("track {}_ind_raw_2\nshortLabel {}_ind_raw\nlongLabel {}_ind_raw\ntype bigWig 0 4\nparent {}_raw\n"
            "bigDataUrl {}.Individual.RawCounts.hg38.bw\n")
    raw += ("track {}_raw_2\nshortLabel {}_raw\nlongLabel {}_raw\ntype bigWig 0 4\nparent {}_all"
            "bigDataUrl {}.RawCounts.hg38.bw\n")

    ilp = ("track {}_ind_raw_3\nshortLabel {}_ind_raw\nlongLabel {}_ind_raw\ntype bigWig 0 4\nparent {}_raw\n"
           "bigDataUrl {}.Individual.RawCounts.hg38.bw")
    ilp += ("track {}_ilp\nshortLabel {}_raw\nlongLabel {}_raw\ntype bigWig 0 4\nparent {}_all"
            "bigDataUrl {}.ILP.hg38.bw")
    ilp += ("track {}_sun_1\nshortLabel {}_sun\nlongLabel {}_sun\ntype bigWig 0 4\nparent {}_all"
            "bigDataUrl {}.sun_model.hg38.bedGraph")
    result = []
    result.append("".join([container.format(uuid, "all"), all.format(uuid)]))
    result.append("".join([container.format(uuid, "raw"), raw.format(uuid)]))
    result.append("".join([container.format(uuid, "ilp"), ilp.format(uuid)]))
    return "\n".join(result)


def build_tracks(individual_ilp, individual_raw, raw, sun, out_dir, uuid, chrom_sizes):
    for p in [individual_ilp, individual_raw, raw, sun]:
        subprocess.call("wigToBigWig {} {} {}".format(p, chrom_sizes, os.path.join(out_dir, "hg38",
                                                                                   os.path.basename(p))))

def main():
    args = parse_args()
    individual_ilp = os.path.join(args.dir, args.uuid + ".Individual.ILP.hg38.wiggle")
    individual_raw = os.path.join(args.dir, args.uuid + ".Individual.RawCounts.hg38.wiggle")
    raw = os.path.join(args.dir, args.uuid + ".RawCounts.hg38.wiggle")
    sun = os.path.join(args.dir, args.uuid + ".sun_model.hg38.bedGraph")
    start_hub(args.out_dir, args.uuid)
    track_db(args.out_dir, args.name)
    build_tracks(individual_ilp, individual_raw, raw, sun, args.out_dir, args.uuid, args.hg38_chrom_sizes)


if __name__ == "__main__":
    main()