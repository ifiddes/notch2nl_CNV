fastq_path = "/cluster/home/ifiddes/ifiddes_hive/notch_mike_snyder/snyder_notch.50mer.Counts.fa"
ref_path = "data/kmer_model_data/notch2nl_graph_hg38.fa"
from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
import cPickle as pickle
graph = UnitigGraph(49)
add_mole_to_graph(graph, ref_path)
add_individual_to_graph(graph, fastq_path)
unpruned = graph.copy()
pickle.dump(unpruned, open("unpruned_graph.pickle", "w"))
graph.prune_source_edges()
graph.prune_individual_edges()
graph.finish_build(graphviz=True)
pickle.dump(graph, open("pruned_graph.pickle", "w"))
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

t = None
for subgraph in graph.connected_component_iter(internal=True):
    if len(frozenset(tuple(subgraph.edge[a][b]['positions'].keys()) for a, b in subgraph.edges_iter() if
                     'positions' in subgraph.edge[a][b])) > 1:
        if t is None:
            t = subgraph
        if len(t) > len(subgraph):
            t = subgraph

pickle.dump(subgraph, open("subgraph.pickle", "w"))