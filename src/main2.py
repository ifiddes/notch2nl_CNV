import argparse
import os
import cPickle as pickle
from collections import namedtuple
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import logger, setLoggingFromOptions
from src.prepareData import PrepareData

ilp_tuple = namedtuple("ilp_tuple", "breakpoint_penalty, data_penalty, expected_value_penalty, trash_penalty, "
                                  "kmer_size")
paths_tuple = namedtuple("paths_tuple", "out_dir, aln_index, whitelist, masked_ref, unmasked_ref, bad_kmers, "
                         "normalizing, key_file")


def build_analyses(target, paths, ilp_config, cgquery_dict):
    for uuid, query_string in cgquery_dict.iteritems():
        # unlikely to have uuid collisions at 8 characters, easier to read
        uuid = uuid[:8]
        target.addChildTarget(PrepareData(paths, ilp_config, uuid, query_string))


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    ilp_config = ilp_tuple(args.breakpoint_penalty, args.data_penalty, args.expected_value_penalty, args.trash_penalty,
                           args.kmer_size)
    paths = paths_tuple(args.out_dir, args.aln_index, args.whitelist, args.masked_ref, args.unmasked_ref,
                        args.bad_kmers, args.normalizing, args.key_file)
    try:
        cgquery_dict = pickle.load(open(args.cgquery_file))
    except IOError:
        raise IOError("Cgquery dict does not exist.")

    if not os.path.exists(paths.out_dir):
        os.makedirs(paths.out_dir)

    i = Stack(Target.makeTargetFn(build_analyses, args=(paths, ilp_config, cgquery_dict))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == "__main__":
    from src.main import *
    main()