from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
from src.kmerIlpModel import *
import cPickle as pickle
#fastq_path = "/hive/users/ifiddes/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
#fastq_path = "/Users/ifiddes/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
fastq_path = "/home/ifiddes/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
masked_ref_path = "data/kmer_model_data/notch2nl_masked_hg38.fa"
unmasked_ref_path = "data/kmer_model_data/notch2nl_unmasked_hg38.fa"
graph = UnitigGraph(kmer_size=49, derived=False)
add_mole_to_graph(graph, unmasked_ref_path, masked_ref_path)
#source_graph = graph.copy()
add_individual_to_graph(graph, fastq_path)
graph.make_graphviz_labels()
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

