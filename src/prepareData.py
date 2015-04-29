import os
import pysam
from src.sunModel import SunModel
from jobTree.scriptTree.target import Target
from jobTree.src.bioio import system


class PrepareData(Target):
    def __init__(self, paths, ilp_config, uuid, query_string):
        Target.__init__(self)
        self.paths = paths
        self.ilp_config = ilp_config
        self.uuid = uuid
        self.query_string = query_string

    def run(self):
        out_dir = os.path.join(self.paths.out_dir, self.uuid)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        fastq_path = os.path.join(out_dir, self.uuid + ".fastq")
        bam_path = os.path.join(out_dir, self.uuid + ".bam")
        kmer_counts_path = os.path.join(out_dir, self.uuid + ".{}mer.fa".format(self.ilp_config.kmer_size))
        k_plus1_mer_counts_path = os.path.join(out_dir, self.uuid + ".{}mer.fa".format(self.ilp_config.kmer_size + 1))
        if not os.path.exists(fastq_path):
            download_query(fastq_path, self.getLocalTempDir(), self.paths.key_file, self.query_string, self.uuid)
        if not os.path.exists(bam_path):
            self.addChildTargetFn(align_query, args=(fastq_path, bam_path, self.uuid, self.paths.aln_index))
        if not os.path.exists(kmer_counts_path) or not os.path.exists(k_plus1_mer_counts_path):
            self.addChildTargetFn(run_jellyfish, args=(kmer_counts_path, k_plus1_mer_counts_path, fastq_path, 
                                                         self.uuid, self.ilp_config.kmer_size))
        self.setFollowOnTarget(SunModel(self.paths, self.ilp_config, self.uuid, fastq_path, bam_path,
                                          kmer_counts_path, k_plus1_mer_counts_path))


def download_query(fastq_tmp_path, tmp_dir, key_file, query_string, uuid):
    """
    Downloads data from CGHub BAM Slicer
    """
    key = open(key_file).readline().rstrip()
    system("""curl --silent "{}" -u "{}" | samtools bamshuf -Ou - {} | samtools bam2fq - > {}""".format(
           query_string, "haussler:" + key, os.path.join(tmp_dir, "tmp"), fastq_tmp_path))
    os.remove(os.path.join(tmp_dir, "tmp"))
    if os.path.getsize(fastq_tmp_path) < 513:
        raise RuntimeError("curl did not download a BAM for {}. exiting.".format(uuid))


def run_jellyfish(target, jf_counts, k_plus1_mer_counts, fastq, uuid, kmer_size):
    """
    Runs jellyfish twice: the first time counts kmers with -C at kmer_size kmers. This is the raw data for the ILP
    model. Runs jellyfish a second time, with kmer_size +1 and the bloom filter which removes most kmers with counts of
    oe. This will be used to add individual nodes to the graph.
    """
    jf_file = os.path.join(target.getLocalTempDir(), uuid + ".jf")
    system("jellyfish count -C -m {} -s 200M -o {} {}".format(kmer_size, jf_file, fastq))
    system("jellyfish dump {} > {}".format(jf_file, jf_counts))
    os.remove(jf_file)
    system("jellyfish count -C -m {} --bf-size 1G -s 200M -o {} {}".format(kmer_size + 1, jf_file, fastq))
    system("jellyfish dump {} > {}".format(jf_file, k_plus1_mer_counts))
    os.remove(jf_file)


def align_query(target, fastq, bam, uuid, index):
    """
    Aligns the extracted reads to the notch locus, filtering for unmapped reads and creating a custom reheadered BAM.
    """
    # align the extracted reads to the index
    tmp = os.path.join(target.getLocalTempDir(), "tmp")
    sorted_bam = os.path.join(target.getLocalTempDir(), "{}.sorted.bam".format(uuid))
    system("bwa mem -v 1 {} {} | samtools view -F 4 -bS - | samtools sort -O bam -T {} - > {}".format(index, fastq, tmp, sorted_bam))
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
