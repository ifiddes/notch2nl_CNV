#!/usr/bin/env python2.7
"""
Takes the output folder for one individual and creates a trackHub out of that individual.
"""

import sys
import os
import argparse
import subprocess


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", "-d", required=True, help="Dir to dig through. If does not end in tracks, will append")
    parser.add_argument("--out_dir", "-o", required=True, help="Dir to put assembly hub into.")
    parser.add_argument("--name", required=True, help="assembly hub name")
    parser.add_argument("--hg38_chrom_sizes", help="hg38 chrom.sizes file", default="hg38.chrom.sizes")
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


def trackdb(dir_path, track_paths):



def main():
    #do stuff



if __name__ == "__main__":
    main()