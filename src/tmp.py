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




from src.unitigGraph import *
from src.helperFunctions import *
from src.kmerModel import *
from src.kmerIlpModel import *
graph = UnitigGraph(7, derived=False)
masked_mole = unmasked_mole = "test/test_ref_short.fa"
seq_map = {}
masked_seqs = {name: seq for name, seq in fasta_read(masked_mole)}
unmasked_seqs = {name: seq for name, seq in fasta_read(unmasked_mole)}
for name in masked_seqs:
    masked_seq = masked_seqs[name]
    unmasked_seq = unmasked_seqs[name]
    assert "N" not in unmasked_seq
    assert "_" in name and len(name.split("_")) == 2
    name, offset = name.split("_")
    try:
        offset = int(offset)
    except TypeError:
        raise RuntimeError("Naming convention for input reference fasta is not right. >Sequence_Offset")
    masked_seq = masked_seq.upper()
    unmasked_seq = unmasked_seq.upper()
    seq_map[(name, offset)] = [masked_seq, unmasked_seq]
    graph.add_masked_kmers(masked_seq, unmasked_seq)

for (name, offset), (masked_seq, unmasked_seq) in seq_map.iteritems():
    graph.add_source_sequence(name, offset, masked_seq, unmasked_seq)


for subgraph in graph.connected_component_iter():
    sequence_edges = [e for e in subgraph.edges_iter() if 'positions' in subgraph.edge[e[0]][e[1]]]
    ordering = defaultdict(list)
    for a, b in sequence_edges:
        for para in subgraph.edge[a][b]['positions']:
            pos = subgraph.edge[a][b]['positions'][para][0]
            ordering[para].append([pos, a[:-2]])


for para, vals in ordering.iteritems():
    ordering[para] = sorted(vals, key=lambda x: x[0])


new_graph = nx.DiGraph()
for subgraph in graph.connected_component_iter():
    sequence_edges = [e for e in subgraph.edges_iter() if 'positions' in subgraph.edge[e[0]][e[1]]]
    adjacency_edges = [e for e in subgraph.edges_iter() if 'source' in subgraph.edge[e[0]][e[1]]]
    new_edges = [e for e in subgraph.edges_iter() if 'source' not in subgraph.edge[e[0]][e[1]] and 'positions' not in subgraph.edge[e[0]][e[1]]]
    for a, b in adjacency_edges:
        l, r = a[:-2], b[:-2]
        new_graph.add_edge(l, r)
        new_graph.node[l]['positions'] = subgraph.edge[l + "_R"][l + "_L"]
        new_graph.node[r]['positions'] = subgraph.edge[r + "_R"][r + "_L"]
        new_graph.node[l]['source'] = new_graph.node[r]['source'] = True
    for a, b in new_edges:
        l, r = a[:-2], b[:-2]
        new_graph.add_edge(l, r)

for a, b in new_graph.edges_iter():
    new_graph.node[a]['fontsize'] = new_graph.node[b]['fontsize'] = 10
    new_graph.edge[a][b]['penwidth'] = 3
    if 'positions' in new_graph.node[a]:
        new_graph.node[a]["label"] = "\\n".join(sorted(
            [": ".join([y, ", ".join([str(x) for x in new_graph.node[a]['positions'][y]])]) for y in
             new_graph.node[a]['positions']]))
        new_graph.node[a]["color"] = "purple"
    if 'positions' in new_graph.node[b]:
        new_graph.node[b]["label"] = "\\n".join(sorted(
            [": ".join([y, ", ".join([str(x) for x in new_graph.node[b]['positions'][y]])]) for y in
             new_graph.node[b]['positions']]))
        new_graph.node[b]["color"] = "purple"
    else:
        new_graph