from itertools import izip
from sonLib.bioio import _getFileHandle as get_file_handle
from jobTree.src.bioio import reverseComplement as reverse_complement


def remove_label(n):
    """
    removes the left/right label from a node, returning the actual kmer
    """
    return n[:-2]


def labels_from_node(n):
    """
    given one node, returns both nodes
    """
    k = remove_label(n)
    return k + "_L", k + "_R"


def labels_from_kmer(k):
    """
    given a kmer, returns a tuple of the two nodes
    """
    return k + "_L", k + "_R"


def labeled_kmer_iter(kmers):
    """
    Given a set of kmers, yields each node individually
    """
    for k in kmers:
        for n in labels_from_kmer(k):
            yield n


def canonical(k):
    """
    Returns the canonical version of this kmer. This is defined as whichever comes first, the kmer or the
    reverse complement of the kmer lexicographically.
    """
    return sorted([k, reverse_complement(k)])[0]


def fastq_reader(path_or_file):
    """
    This is a super fast fastq reader that expects things to be correct, I.E. 4 lines per entry.
    """
    handle = get_file_handle(path_or_file)
    for name, seq, _, qual in izip(*[handle] * 4):
        yield seq.rstrip()


def count_reader(path_or_file):
    """
    This is a super fast jellyfish count reader that expects things to be correct, I.E. 2 lines per entry.
    """
    handle = get_file_handle(path_or_file)
    rm = ">\n"
    for count, seq in izip(*[handle] * 2):
        yield int(count.translate(None, rm)), seq.translate(None, rm)


def format_ratio(numerator, denominator):
    if denominator == 0:
        return float("nan")
    return float(numerator)/denominator