fastq_path = "/inside/home/ifiddes/depot4/notch2nl_copy_number_variation/output_weighted/b88da23c/b88da23c.fastq"
ref_path = "data/kmer_model_data/notch2nl_graph_hg38.fa"
from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
graph = UnitigGraph(49)
add_mole_to_graph(graph, ref_path)
add_individual_to_graph(graph, fastq_path)
graph.prune_source_edges()

combos = set()
combos = set()
for subgraph in graph.connected_component_iter(internal=True):
    self_loop_nodes = [a for a, b in subgraph.selfloop_edges() if 'source' in graph.edge[a][b]]
    bubble_nodes = [n for n in subgraph.nodes_iter() if subgraph.degree(n) > 2 and n not in self_loop_nodes]
    # we count the number of source sequences each bubble node is attached to
    source_sequences = frozenset(tuple(subgraph.edge[a][b]['positions'].keys()) for a, b in subgraph.edges_iter() if
                     'positions' in subgraph.edge[a][b])
    combos.add(source_sequences)