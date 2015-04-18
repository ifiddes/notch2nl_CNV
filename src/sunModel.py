import os
import pysam
from collections import defaultdict, Counter
from src.helperFunctions import format_ratio
from src.kmerModel import KmerModel
from jobTree.scriptTree.target import Target


class SunModel(Target):
    """
    Runs the SUN model, where alelle fraction at unique sites is used to infer copy number
    """
    def __init__(self, paths):
        Target.__init__(self)
        self.paths = paths
        self.bam = paths.bam
        self.uuid = paths.uuid
        self.out_dir = paths.out_dir
        # whitelist is a text file of whitelisted SUN positions - ignore commented out lines
        # format is paralog, hg19_pos, ref_base, alt_base, hg38_pos
        wl_list = [x.split() for x in open(paths.whitelist_path) if not x.startswith("#")]
        # dict mapping genome positions to which paralog has a SUN at that position
        self.whitelist = {int(hg19_pos): (para, ref, alt, hg38_pos) for para, hg19_pos, ref, alt, hg38_pos in wl_list}
        
    def find_site_coverages(self, bam_in, min_depth=40):
        """
        Runs mpileup on each SUN position and finds allele fractions
        """
        bases = {"A", "T", "G", "C", "a", "t", "g", "c"}
        results_dict = defaultdict(list)
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
                results_dict[para].append([hg38_pos, frac])
        # calculate normalizing factor from mean value of Notch2 SUNs, then normalize other SUNs
        normalizing_factor = 1.0 * sum(n_vals) / len(n_vals)
        for para, result in results_dict.iteritems():
            results_dict[para] = [[pos, 2.0 * val / normalizing_factor] for pos, val in result]
        return results_dict

    def make_bedgraphs(self, results_dict):
        """
        Generates hg38 bedGraphs for the SUN results.
        """
        if not os.path.exists(os.path.join(self.out_dir, "tracks")):
            os.mkdir(os.path.join(self.out_dir, "tracks"))
        path = os.path.join(self.out_dir, "tracks", "{}.sun_model.hg38.bedGraph".format(self.uuid))
        merged_results = [[pos, val] for result in results_dict.itervalues() for pos, val in result]
        sorted_merged_results = sorted(merged_results, key=lambda x: x[0])
        with open(path, "w") as outf:
            bed_header = ("track type=bedGraph name={} autoScale=off visibility=full alwaysZero=on yLineMark=2 "
                          "viewLimits=0:4 yLineOnOff=on maxHeightPixels=100:75:50\n")
            outf.write(bed_header.format(self.uuid))
            for pos, val in sorted_merged_results:
                outf.write("\t".join(map(str, ["chr1", pos, pos + 1, val])) + "\n")

    def run(self):
        results_dict = self.find_site_coverages(self.paths.bam)
        inferred_c, inferred_d = infer_copy_number(results_dict)
        self.make_bedgraphs(results_dict)
        self.setFollowOnTarget(KmerModel(self.paths, inferred_c, inferred_d))


def infer_copy_number(results_dict):
    """
    Infer copy number of Notch2NL-C and Notch2NL-D
    """
    if len(results_dict['C']) == 0:
        c = 0
    else:
        c_vals = zip(*results_dict["C"])[1]
        c = int(round(sum(c_vals) / len(c_vals), 0))
    if len(results_dict['D']) == 0:
        d = 0
    else:
        d_vals = zip(*results_dict["D"])[1]
        d = int(round(sum(d_vals) / len(d_vals), 0))
    return c, d