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


def strandless(k):
    """
    Returns the strandless version of this kmer. This is defined as whichever comes first, the kmer or the
    reverse complement of the kmer lexicographically.
    """
    return sorted([k, reverse_complement(k)])[0]


def format_ratio(numerator, denominator):
    if denominator == 0:
        return float("nan")
    return float(numerator)/denominator