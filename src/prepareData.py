import os
import pysam
from src.sunModel import SunModel
from jobTree.scriptTree.target import Target
from jobTree.src.bioio import system



class PrepareData(Target):
    """Takes the information from a paths namedtuple and runs all data preparation steps before running models"""
    def __init__(self, paths, ilp_config):
        Target.__init__(self)
        self.paths = paths
        self.ilp_config = ilp_config

    def run(self):
        paths = self.paths
        if not os.path.exists(paths.fastq):
            download_query(paths.fastq, paths.out_dir, paths.cghub_key, paths.query_string, paths.uuid)
        if not os.path.exists(paths.bam):
            Target.addChildTargetFn(align_query, args=(paths.fastq, paths.bam, paths.out_dir, paths.uuid,
                                                       paths.aln_index))
        if not os.path.exists(paths.jf_counts) or not os.path.exist(paths.jf_kplus1_counts):
            Target.addChildTargetFn(run_jellyfish, args=(paths.out_dir, paths.jf_counts, paths.jf_kplus1_counts,
                                                         paths.fastq, paths.uuid, paths.kmer_size))
        Target.setFollowOnTarget(SunModel(paths, ilp_config))


def download_query(fastq, out_dir, cghub_key, query_string, uuid):
    """
    Downloads data from CGHub BAM Slicer
    """
    system("""curl --silent "{}" -u "{}" | samtools bamshuf -Ou - {} | samtools bam2fq - > {}""".format(
           query_string, "haussler:" + cghub_key, os.path.join(out_dir, "tmp"), fastq))
    if os.path.getsize(fastq) < 513:
        raise RuntimeError("curl did not download a BAM for {}. exiting.".format(uuid))
    os.remove(os.path.join(out_dir, "tmp"))


def run_jellyfish(target, out_dir, jf_counts, jf_kplus1_counts, fastq, uuid, kmer_size):
    """
    Runs jellyfish twice: the first time counts kmers with -C at kmer_size kmers. This is the raw data for the ILP
    model. Runs jellyfish a second time, with kmer_size +1 and the bloom filter which removes most kmers with counts of
    oe. This will be used to add individual nodes to the graph.
    """
    jf_file = os.path.join(out_dir, uuid + ".jf")
    system("jellyfish count -C -m {} -s 200M -o {} {}".format(kmer_size, jf_file, fastq))
    system("jellyfish dump {} > {}".format(jf_file, jf_counts))
    system("jellyfish count -C -m {} --bf-size 1G -s 200M -o {} {}".format(kmer_size + 1, jf_file, fastq))
    system("jellyfish dump {} > {}".format(jf_file, jf_kplus1_counts))
    os.remove(jf_file)

def align_query(target, fastq, bam, out_dir, uuid, index):
    """
    Aligns the extracted reads to the notch locus, filtering for unmapped reads and creating a custom reheadered BAM.
    """
    # align the extracted reads to the index
    sorted_bam = os.path.join(out_dir, "{}.sorted.bam".format(uuid))
    system("bwa mem -v 1 {} {} | samtools view -F 4 -bS - | samtools sort - > {}".format(index, fastq, sorted_bam))
    header = {"HD": {"VN": "1.3"}, "SQ": [{"LN": 248956422, "SN": "chr1"}]}
    outfile = pysam.Samfile(bam, "wb", header=header)
    bamfile = pysam.Samfile(sorted_bam, "rb")
    for record in bamfile:
        chrom, span = bamfile.getrname(record.tid).split(":")
        start, end = map(int, span.split("-"))
        record.pos = record.pos + start - 1
        outfile.write(record)
    outfile.close()
    system("samtools index {}".format(bam))
