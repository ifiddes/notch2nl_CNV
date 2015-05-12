#!/usr/bin/env python2.7
"""
Takes the output folder for one individual and creates a trackHub out of that individual.
"""

import sys
import os
import argparse
import subprocess

colors = {"ILP": "242,148,176", "SUN": "180,0,0", "Raw": "194,194,143", "RawIndividual": "143,162,194"}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", "-d", required=True, help="Dir to dig through. If does not end in tracks, will append")
    parser.add_argument("--out_dir", "-o", required=True, help="Dir to put assembly hub into.")
    parser.add_argument("--name", required=True, help="assembly hub name")
    parser.add_argument("--uuid", required=True, help='UUID')
    parser.add_argument("--hg38_chrom_sizes", help="hg38 chrom.sizes file", default="data/hg38.chrom.sizes")
    return parser.parse_args()


def start_hub(dir_path, name):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    with open(os.path.join(dir_path, "hub.txt"), "w") as outf:
        hub_str = "hub {0}\nshortLabel {0}\nlongLabel {0}\ngenomesFile genomes.txt\nemail ian.t.fiddes@gmail.com\n"
        outf.write(hub_str.format(name, name + " Notch2NL", name + " Notch2NL"))
    with open(os.path.join(dir_path, "genomes.txt"), "w") as outf:
        outf.write("genome hg38\ntrackDb hg38/trackDb.txt\ndefaultPos chr1:145987622-149723055\n\n")
    if not os.path.exists(os.path.join(dir_path, "hg38")):
        os.mkdir(os.path.join(dir_path, "hg38"))
    with open(os.path.join(dir_path, "hg38", "Notch2NL.html"), "w") as outf:
        outf.write("Notch2NL {0}\n".format(name))


def track_db(out_dir, uuid):
    """
    This generates three types of tracks.
    1) all: overlays individual raw data, raw data and SUN
    2) raw: overlays individual raw and raw
    3) ILP: overlays individual ILP result with individual raw data and SUN information
    """
    container = ("track {0}_{1}\ncontainer multiWig\nshortLabel {0}_{1}\nlongLabel {2} Combined\n"
                 "type bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\n"
                  "yLineOnOff on\nmaxHeightPixels 100:75:50\n\n")
    all = ("track {0}_ind_raw_1\nshortLabel {0}_ind_raw\nlongLabel {0}_ind_raw\ntype bigWig 0 4\nparent {0}_all\n"
           "bigDataUrl {0}.Individual.RawCounts.hg38.wiggle\ncolor {RawIndividual}\n\n"
           "track {0}_raw_1\nshortLabel {0}_raw\nlongLabel {0}_raw\ntype bigWig 0 4\nparent {0}_all\n"
           "bigDataUrl {0}.RawCounts.hg38.wiggle\ncolor {Raw}\n\n"
           "track {0}_sun_1\nshortLabel {0}_sun\nlongLabel {0}_sun\ntype bigWig 0 4\nparent {0}_all\n"
           "bigDataUrl {0}.sun_model.hg38.bedGraph\ncolor {SUN}\n\n")

    raw = ("track {0}_ind_raw_2\nshortLabel {0}_ind_raw\nlongLabel {0}_ind_raw\ntype bigWig 0 4\nparent {0}_raw\n"
           "bigDataUrl {0}.Individual.RawCounts.hg38.wiggle\ncolor {RawIndividual}\n\n"
           "track {0}_raw_2\nshortLabel {0}_raw\nlongLabel {0}_raw\ntype bigWig 0 4\nparent {0}_all\n"
           "bigDataUrl {0}.RawCounts.hg38.wiggle\ncolor {Raw}\n\n")

    ilp = ("track {0}_ind_raw_3\nshortLabel {0}_ind_raw\nlongLabel {0}_ind_raw\ntype bigWig 0 4\nparent {0}_raw\n"
           "bigDataUrl {0}.Individual.RawCounts.hg38.wiggle\ncolor {RawIndividual}\n\n"
           "track {0}_ilp_1\nshortLabel {0}_raw\nlongLabel {0}_raw\ntype bigWig 0 4\nparent {0}_all\n"
           "bigDataUrl {0}.ILP.hg38.wiggle\ncolor {ILP}\n\n"
           "track {0}_sun_1\nshortLabel {0}_sun\nlongLabel {0}_sun\ntype bigWig 0 4\nparent {0}_all\n"
           "bigDataUrl {0}.sun_model.hg38.bedGraph\ncolor {SUN}\n\n")
    result = []
    result.append("".join([container.format(uuid, "all", "IndividualRaw/Raw/SUN"), all.format(uuid, **colors)]))
    #result.append("".join([container.format(uuid, "raw", "IndividualRaw/Raw"), raw.format(uuid, **colors)]))
    #result.append("".join([container.format(uuid, "ilp"), ilp.format(uuid, **colors)]))
    with open(os.path.join(out_dir, "hg38", "trackDb.txt"), "w") as outf:
        outf.write("\n".join(result))


def build_tracks(individual_ilp, individual_raw, raw, sun, out_dir, uuid, chrom_sizes):
    for p in [individual_ilp, individual_raw, raw]:
        subprocess.call("wigToBigWig {} {} {}".format(p, chrom_sizes, os.path.join(out_dir, "hg38",
                                                                                   os.path.basename(p))),
                        shell=True)
    subprocess.call("bedGraphToBigWig {} {} {}".format(sun, chrom_sizes, os.path.join(out_dir, "hg38",
                                                                                      os.path.basename(sun))),
                    shell=True)


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