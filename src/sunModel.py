import os
import pysam
import cPickle as pickle
from collections import defaultdict, Counter
from src.helperFunctions import format_ratio
from src.kmerModel import KmerModel
from jobTree.scriptTree.target import Target


class SunModel(Target):
    """
    Runs the SUN model, where alelle fraction at unique sites is used to infer copy number
    """
    def __init__(self, paths, ilp_config, uuid, fastq_path, bam_path, kmer_counts_path, k_plus1_mer_counts_path):
        Target.__init__(self)
        self.paths = paths
        self.ilp_config = ilp_config
        self.uuid = uuid
        self.bam_path = bam_path
        self.fastq_path = fastq_path
        self.kmer_counts_path = kmer_counts_path
        self.k_plus1_mer_counts_path = k_plus1_mer_counts_path
        # whitelist is a text file of whitelisted SUN positions - ignore commented out lines
        # format is paralog, hg19_pos, ref_base, alt_base, hg38_pos
        wl_list = [x.split() for x in open(paths.whitelist) if not x.startswith("#")]
        # dict mapping genome positions to which paralog has a SUN at that position
        self.whitelist = {int(hg19_pos): (para, ref, alt, hg38_pos) for para, hg19_pos, ref, alt, hg38_pos in wl_list}

    def find_site_coverages(self, bam_in, min_depth=20):
        """
        Runs mpileup on each SUN position and finds allele fractions
        """
        bases = {"A", "T", "G", "C", "a", "t", "g", "c"}
        sun_results = defaultdict(list)
        n_vals = []
        for pos, (para, ref, alt, hg38_pos) in self.whitelist.iteritems():
            pos_str = "chr1:{0}-{0}".format(pos)
            pile_up = pysam.mpileup("-q", "20", "-Q", "20", "-r", pos_str, bam_in)
            if len(pile_up) == 0:
                continue
            pile_up_str = pile_up[0].split()
            if len(pile_up_str) != 6:
                continue
            pile_up_result = Counter(x.upper() for x in pile_up_str[4] if x in bases)
            if ref not in pile_up_result or alt not in pile_up_result:
                continue
            if sum(pile_up_result.itervalues()) < min_depth:
                continue
            frac = format_ratio(pile_up_result[alt], sum(pile_up_result.values()))
            # invert fraction for Notch2 paralogs
            if para == "N":
                frac = 1 - frac
                n_vals.append(frac)
            else:
                sun_results[para].append([hg38_pos, frac])
        # calculate normalizing factor from mean value of Notch2 SUNs, then normalize other SUNs
        normalizing_factor = 1.0 * sum(n_vals) / len(n_vals)
        for para, result in sun_results.iteritems():
            sun_results[para] = [[pos, 2.0 * val / normalizing_factor] for pos, val in result]
        return sun_results

    def make_bedgraphs(self, sun_results):
        """
        Generates hg38 bedGraphs for the SUN results.
        """
        if not os.path.exists(os.path.join(self.paths.out_dir, self.uuid, "tracks")):
            os.mkdir(os.path.join(self.paths.out_dir, self.uuid, "tracks"))
        path = os.path.join(self.paths.out_dir, self.uuid, "tracks", "{}.sun_model.hg38.bedGraph".format(self.uuid))
        merged_results = [[pos, val] for result in sun_results.itervalues() for pos, val in result]
        sorted_merged_results = sorted(merged_results, key=lambda x: x[0])
        with open(path, "w") as outf:
            bed_header = ("track type=bedGraph name={} autoScale=off visibility=full alwaysZero=on yLineMark=2 "
                          "viewLimits=0:4 yLineOnOff=on maxHeightPixels=100:75:50\n")
            outf.write(bed_header.format(self.uuid))
            for pos, val in sorted_merged_results:
                outf.write("\t".join(map(str, ["chr1", pos, int(pos) + 1, val])) + "\n")

    def run(self):
        sun_results = self.find_site_coverages(self.bam_path)
        inferred_c, inferred_d = infer_copy_number(sun_results)
        self.make_bedgraphs(sun_results)
        pickle.dump(sun_results, open(os.path.join(self.paths.out_dir, self.uuid, "sun_results.pickle"), "w"))
        self.setFollowOnTargetFn(kmerModelWrapperFn, args=(self.paths, self.uuid, self.ilp_config, sun_results,
                                                           self.fastq_path, self.kmer_counts_path,
                                                           self.k_plus1_mer_counts_path, inferred_c, inferred_d))


def kmerModelWrapperFn(target, paths, uuid, ilp_config, sun_results, fastq_path, kmer_counts_path, 
                       k_plus1_mer_counts_path, inferred_c, inferred_d):
    #target.addChildTarget(KmerModel(paths, uuid, ilp_config, sun_results, fastq_path, kmer_counts_path, 
    #                                k_plus1_mer_counts_path, inferred_c, inferred_d, add_individual=True))
    target.addChildTarget(KmerModel(paths, uuid, ilp_config, sun_results, fastq_path, kmer_counts_path, 
                                    k_plus1_mer_counts_path, inferred_c, inferred_d, add_individual=False))


def infer_copy_number(sun_results):
    """
    Infer copy number of Notch2NL-C and Notch2NL-D
    """
    if len(sun_results['C']) == 0:
        c = 0
    else:
        c_vals = zip(*sun_results["C"])[1]
        c = int(round(sum(c_vals) / len(c_vals), 0))
    if len(sun_results['D']) == 0:
        d = 0
    else:
        d_vals = zip(*sun_results["D"])[1]
        d = int(round(sum(d_vals) / len(d_vals), 0))
    return c, d