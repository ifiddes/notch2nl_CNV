#fastq_path = "/cluster/home/ifiddes/ifiddes_hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
#fastq_path = "/home/ifiddes/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
fastq_path = "~/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
masked_ref_path = "data/kmer_model_data/notch2nl_masked_hg38.fa"
unmasked_ref_path = "data/kmer_model_data/notch2nl_unmasked_hg38.fa"
from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
import cPickle as pickle
graph = UnitigGraph(49)
add_mole_to_graph(graph, unmasked_ref_path, masked_ref_path)
len(graph)
add_individual_to_graph(graph, fastq_path)
len(graph)
unpruned = graph.copy()
#pickle.dump(unpruned, open("unpruned_graph.pickle", "w"))
graph.prune_source_edges()
graph.prune_individual_edges()
graph.finish_build(graphviz=True)
#pickle.dump(graph, open("pruned_graph.pickle", "w"))

import cPickle as pickle
from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
graph = pickle.load(open("data/kmer_model_data/mole_graph_pruned.pickle"))
fastq_path = "/hive/users/ifiddes/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
add_individual_to_graph(graph, fastq_path)

"""
combos = set()
for subgraph in graph.connected_component_iter(internal=True):
    self_loop_nodes = [a for a, b in subgraph.selfloop_edges() if 'source' in graph.edge[a][b]]
    bubble_nodes = [n for n in subgraph.nodes_iter() if subgraph.degree(n) > 2 and n not in self_loop_nodes]
    # we count the number of source sequences each bubble node is attached to
    source_sequences = frozenset(tuple(subgraph.edge[a][b]['positions'].keys()) for a, b in subgraph.edges_iter() if
                     'positions' in subgraph.edge[a][b])
    combos.add(source_sequences)
"""

subgraphs = []
for subgraph in graph.connected_component_iter(internal=True):
    if len(frozenset(tuple(subgraph.edge[a][b]['positions'].keys()) for a, b in subgraph.edges_iter() if
                     'positions' in subgraph.edge[a][b])) > 1:
            subgraphs.append(subgraph)

pickle.dump(subgraph, open("subgraph.pickle", "w"))


ends = [x for x in subgraph.nodes_iter() if subgraph.degree(x) == 1]
dist = 30
tmp = subgraph.copy()
for n in ends:
    for i in xrange(dist):
        adj = unpruned.adj[n]
        for inner_node in adj:
            tmp.add_node(inner_node)
            tmp.node[inner_node] = graph.node[inner_node]
            for a, b in unpruned.edges_iter(inner_node):
                tmp.add_edge(a, b)
                tmp.edge[a][b] = graph.edge[a][b]

positions = defaultdict(list)
for a, b in unpruned.edges_iter():
    if 'positions' in unpruned.edge[a][b]:
        for para in unpruned.edge[a][b]['positions']:
            for val in unpruned.edge[a][b]['positions'][para]:
                positions[para].append((val, remove_label(a)))

for x in positions:
    positions[x] = sorted(positions[x], key = lambda x: x[0])

subgraphs = list(graph.connected_component_iter(internal=True))
subgraphs = sorted(subgraphs, key=lambda x: len(x))


from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
masked_seq = "NATGCACAANNAAGAGAG"
unmasked_seq = "AATGCACAACGAAGAGAG"
name = "A"
offset = 0
graph = UnitigGraph(5)
graph.add_source_sequence(name, offset, masked_seq, unmasked_seq)


tmp = []
for l, r in graph.edges_iter():
    if 'positions' in graph.edge[l][r] and 'Notch2NL-A' in graph.edge[l][r]['positions']:
        for x in graph.edge[l][r]['positions']['Notch2NL-A']:
            if x == 1224:
                print l[:-2]

import random
unmasked_seq = "GCATTTTAAGACTGTGCTGTTATAAAACCTCGGGCCACTTAACTGATTAATCATGGCAATGAGGGCAGGGACCAGAAAAGAGTCTTTAGAACCTGTCATCCCCACACAGAAGAGCAACTTTCAGGGAAACACCCTTATCTTTCCATTTTCAGACCCCGGGAGGTGTGAGGGTGGAAAGGCTAGGTAGAGAAGAGAGCAGAAAGGAGATGAGATGACACAACCAGGATTCTCCGAAGCTGGGCTTGAAGTCCTCAAGAAAAACTCCCATGAACAAGGAAGGAAGAGTGAAGAAAAAAACAGGGATACCTGGAACTGGACAAAAGTAAAAAGATAGAAGGATACTTTTTTTCCCCCAGAAGAAGTCTGTCACAAAAGCAAACCTGCAAATATACGATCAGTATAACACCCAAGAAAATGACACATGCGGC"
masked_seq = list(unmasked_seq)
for i in xrange(30):
    p = int(len(masked_seq) * random.random())
    masked_seq[p] = "N"

masked_seq = "".join(masked_seq)
graph = UnitigGraph(49)
graph.add_source_sequence("A", 0, masked_seq, unmasked_seq)


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




from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
import cPickle as pickle
#fastq_path = "/hive/users/ifiddes/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
fastq_path = "/Users/ifiddes/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
masked_ref_path = "data/kmer_model_data/notch2nl_masked_hg38.fa"
unmasked_ref_path = "data/kmer_model_data/notch2nl_unmasked_hg38.fa"
graph = UnitigGraph(49)
add_mole_to_graph(graph, unmasked_ref_path, masked_ref_path)
graph.prune_source_edges()
source_graph = graph.copy()
add_individual_to_graph(graph, fastq_path)
graph.make_graphviz_labels()
unpruned = graph.copy()
graph.prune_individual_edges()
subgraphs = list(graph.connected_component_iter())
test = []
for x in subgraphs:
    for a, b in x.edges_iter():
        if 'positions' not in x.edge[a][b] and 'source' not in x.edge[a][b]:
            test.append(x)
            break

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
import cPickle as pickle
import networkx as nx
from jobTree.src.bioio import system
graph = UnitigGraph(5)
add_mole_to_graph(graph, "test_ref_short.fa", "test_ref_short.fa")
nx.write_dot(graph, "unpruned_test.dot")
graph.prune_source_edges()
nx.write_dot(graph, "pruned_test.dot")
system("dot -Tpdf pruned_test.dot > pruned_test.pdf")
system("dot -Tpdf unpruned_test.dot > unpruned_test.pdf")
add_individual_to_graph(graph, "test_counts_short.fa")
nx.write_dot(graph, "with_individual.dot")
system("dot -Tpdf with_individual.dot > with_individual.pdf")
graph.prune_individual_edges()
nx.write_dot(graph, "with_individual_pruned.dot")
system("dot -Tpdf with_individual_pruned.dot > with_individual_pruned.pdf")