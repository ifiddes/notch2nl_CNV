from jobTree.src.bioio import reverseComplement as reverse_complement


def remove_label(edge):
    """
    removes the left/right label from a edge, returning the actual kmer
    """
    return edge[:-2]


def create_labels(edge):
    """
    returns the node names between a sequence edge given the edge (kmer)
    """
    k = remove_label(edge)
    return k + "_L", k + "_R"


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