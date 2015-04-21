#fastq_path = "/cluster/home/ifiddes/ifiddes_hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
fastq_path = "/home/ifiddes/hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
masked_ref_path = "data/kmer_model_data/notch2nl_masked_hg38.fa"
unmasked_ref_path = "data/kmer_model_data/notch2nl_unmasked_hg38.fa"
from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
import cPickle as pickle
graph = UnitigGraph(49)
add_mole_to_graph(graph, unmasked_ref_path, masked_ref_path)
add_individual_to_graph(graph, fastq_path)
unpruned = graph.copy()
#pickle.dump(unpruned, open("unpruned_graph.pickle", "w"))
graph.prune_source_edges()
graph.prune_individual_edges()
graph.finish_build(graphviz=True)
#pickle.dump(graph, open("pruned_graph.pickle", "w"))
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
subgraphs = sorted(subgraphs, key = lambda x: len(x))