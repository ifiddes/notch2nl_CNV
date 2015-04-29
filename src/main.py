import argparse
import cPickle as pickle
from collections import namedtuple
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import logger, setLoggingFromOptions
from src.prepareData import PrepareData


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cgquery_file", required=True,
                        help="output file from cgqueryHandler.py.")
    ####################################################################################################################
    # paths
    ####################################################################################################################
    parser.add_argument("--out_dir", default="./output/",
                        help="base output directory that results will be written to. Default is ./output/")
    parser.add_argument("--key_file", default="/inside/home/cwilks/haussl_new.key",
                        help="text file containing the cgHub security key.")
    parser.add_argument("--aln_index", default="data/sun_model_data/hs_n2.masked.fa",
                        help="index for SUN model.")
    parser.add_argument("--whitelist", default="data/sun_model_data/hg38_unfiltered_whitelist.txt",
                        help="whitelist for SUN model.")
    parser.add_argument("--masked_ref_path", default="data/kmer_model_data/notch2nl_masked_hg38.fa",
                        help="repeat masked (to N) reference for ILP model.")
    parser.add_argument("--unmasked_ref_path", default="data/kmer_model_data/notch2nl_unmasked_hg38.fa",
                        help="unmasked reference for ILP model.")
    parser.add_argument("--bad_kmers", default="data/kmer_model_data/bad_kmers.txt",
                        help="bad kmers. These are kmers present elsewhere in genome that are in the reference files.")
    parser.add_argument("--normalizing", default="data/kmer_model_data/normalizing.fa",
                        help="fasta file containing a normalizing region with expected copy number 2.")
    ####################################################################################################################
    # ilp_config
    ####################################################################################################################
    parser.add_argument("--breakpoint_penalty", type=float, default=75.0,
                        help="breakpoint penalty used for ILP model.")
    parser.add_argument("--data_penalty", type=float, default=1.0,
                        help="data penalty used for ILP model.")
    parser.add_argument("--expected_value_penalty", type=float, default=0.25,
                        help="How closely should a copy number of 2 be enforced?")
    parser.add_argument("--trash_penalty", type=float, default=0.1,
                        help="How closely should we keep the trash to 0? (How many kmers do we expect to have excess"
                             "counts?)")
    parser.add_argument("--kmer_size", type=int, default=49, help="kmer size")
    return parser


def build_analyses(target, paths, ilp_config, cgquery_dict):
    for uuid, query_string in cgquery_dict.iteritems():
        # unlikely to have uuid collisions at 8 characters, easier to read
        uuid = uuid[:8]
        target.addChildTarget(PrepareData(paths, ilp_config, uuid, query_string))


def main():
    ilp_tuple = namedtuple("ilp_config", "breakpoint_penalty, data_penalty, expected_value_penalty, trash_penalty, "
                                      "kmer_size")
    paths_tuple = namedtuple("paths", "out_dir, aln_index, whitelist, masked_ref_path, unmasked_ref_path, bad_kmers, "
                             "normalizing, key_file")
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    ilp_config = ilp_tuple(args.breakpoint_penalty, args.data_penalty, args.expected_value_penalty, args.trash_penalty,
                           args.kmer_size)
    paths = paths_tuple(args.out_dir, args.aln_index, args.whitelist, args.masked_ref_path, args.unmasked_ref_path,
                        args.bad_kmers, args.normalizing, args.key_file)
    try:
        cgquery_dict = pickle.load(open(args.cgquery_file))
    except IOError:
        raise IOError "Cgquery dict does not exist."

    i = Stack(Target.makeTargetFn(build_analyses, args=(paths, ilp_config, cgquery_dict)))

    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == "__main__":
    from src.main import *
    main()