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
    parser.add_argument("--dir", "-d", required=True, help="Dir to dig through. If does not end in tracks, will append"
    parser.add_argument("--out_dir", "-o", required=True, help="Dir to put assembly hub into.")
    parser.add_argument("--name", required=True, help="assembly hub name")
    parser.add_argument("--hg38_chrom_sizes", help="hg38 chrom.sizes file", default="hg38.chrom.sizes")
    return parser.parse_args()