from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
from src.kmerIlpModel import *
very_mini_graph = UnitigGraph(5, derived=False)
add_mole_to_graph(very_mini_graph, "very_short_unmasked.fa", "very_short_masked.fa")
subgraphs = list(very_mini_graph.connected_component_iter())
print very_mini_graph.paralogs
import cPickle as pickle
#fastq_path = "/hive/users/ifiddes/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
#fastq_path = "/Users/ifiddes/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
fastq_path = "output/mike_sny/mike_sny.50mer.fa"
masked_ref_path = "data/kmer_model_data/notch2nl_masked_hg38.fa"
unmasked_ref_path = "data/kmer_model_data/notch2nl_unmasked_hg38.fa"
graph = UnitigGraph(kmer_size=49, derived=False)
add_mole_to_graph(graph, unmasked_ref_path, masked_ref_path)
graph.make_graphviz_labels()
subgraphs = list(graph.connected_component_iter())
from collections import defaultdict
r = defaultdict(list)

def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

for s in subgraphs:
    for para in s.paralogs:
        r[para].append(s.paralogs[para])

for para in r:
    r[para] = sorted(r[para])

bad = []
for para in r:
    for i in xrange(1, len(r[para])):
        a, b = r[para][i - 1], r[para][i]
        if get_overlap(a, b) != 0:
            bad.append([a, b])

#first hit
b = [[120752460, 120753187], [120752869, 120752870]]
for s in subgraphs:
    for para in s.paralogs:
        if s.paralogs[para] == b[0]:
            left = s
            print para
        elif s.paralogs[para] == b[1]:
            right = s
            print para
#45094, 45093
#encased by: [44684, 45411]
#           CAATTTTTATTTATTTTATTCTCCCTTGCTTTGTTTCTTCAGACATACA, AAGTGCCTCAAAAGGTCCTAAACCCTTGACTCATAGCCAGGATATTTCT
# found this region in the graph, extracted it
mini_graph = UnitigGraph(49, derived=False)
add_mole_to_graph(mini_graph, "test_unmasked.fasta", "test_masked.fasta")
subgraphs = list(mini_graph.connected_component_iter())


very_mini_graph = UnitigGraph(5, derived=False)
add_mole_to_graph(very_mini_graph, "very_short_unmasked.fa", "very_short_masked.fa")
subgraphs = list(very_mini_graph.connected_component_iter())
for s in subgraphs:
    print s.paralogs

add_individual_to_graph(graph, fastq_path)

#unpruned = graph.copy()
#graph.prune_individual_edges()
subgraphs = list(graph.connected_component_iter())
test = []
for x in subgraphs:
    for a, b in x.edges_iter():
        if 'positions' not in x.edge[a][b] and 'source' not in x.edge[a][b]:
            test.append(x)
            break

ilp_model = KmerIlpModel(graph, set())

# remove all individual edges to determine which subgraphs were lost
tmp = graph.copy()
t = []
subgraphs = []
for subgraph in tmp.connected_component_iter():
    q = subgraph.copy()
    nodes_to_remove = set()
    edges_to_remove = set()
    for n in subgraph.nodes_iter():
        a, b = labels_from_node(n)
        if 'positions' not in subgraph.edge[a][b]:
            nodes_to_remove.update([a, b])
    for a, b in subgraph.edges_iter():
        if 'positions' not in subgraph.edge[a][b] and 'source' not in subgraph.edge[a][b]:
            edges_to_remove.add(frozenset(sorted([a, b])))
    for a, b in edges_to_remove:
        subgraph.remove_edge(a, b)
    for n in nodes_to_remove:
        subgraph.remove_node(n)
    try:
        assert len(list(subgraph.connected_component_iter())) == 1
    except AssertionError:
        t.append(q)
    subgraphs.append(subgraph)

test = []
for subgraph in tmp.connected_component_iter():
    for a, b in subgraph.edges_iter():
        if 'positions' not in subgraph.edge[a][b] and 'source' not in subgraph.edge[a][b]:
            test.append(subgraph)
            break


from src.kmerModel import *
from src.kmerIlpModel import *
import cPickle as pickle
graph = UnitigGraph(5, derived=False)
add_mole_to_graph(graph, "test_ref_short.fa", "test_ref_short.fa")
add_individual_to_graph(graph, "test_counts_short.fa")
ilp_model = KmerIlpModel(graph, set())

from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
masked_seq = "GCATTNTAAGACTNNGCTGTTATAA"
unmasked_seq = "GCATTTTAAGACTGTGCTGTTATAA"
name = "A"
offset = 0
graph = UnitigGraph(5)
graph.add_masked_kmers(masked_seq, unmasked_seq)
graph.add_source_sequence(name, offset, masked_seq, unmasked_seq)

import cPickle as pickle
from src.kmerModel import *
from src.kmerIlpModel import *
graph = pickle.load(open("output/b88da23c/b88da23c_graph.pickle"))
ilp_model = KmerIlpModel(graph)
normalizing_kmers = get_normalizing_kmers("data/kmer_model_data/normalizing.fa", 49)
data_counts, normalizing = get_kmer_counts(graph, normalizing_kmers, "output/b88da23c/b88da23c.49mer.fa")
ilp_model.introduce_data(data_counts, normalizing)
ilp_model.solve()
sun_results = pickle.load(open("output/b88da23c/sun_results.pickle"))
result_dict = ilp_model.report_copy_map()
raw_counts = ilp_model.report_normalized_raw_data_map()
combined_plot(result_dict, raw_counts, graph.paralogs, sun_results, "b88", "./")