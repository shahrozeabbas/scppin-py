import pandas as pd
from scppin import scPPIN
import matplotlib.pyplot as plt
import networkx as nx


network_df = pd.read_csv('tests/denoised_mg_ad_fava_network.tsv', sep='\t', index_col=0)

pvals_df = pd.read_csv('tests/magma_gene_symbol_results.tsv', sep='\t')
pvalues = dict(zip(pvals_df['GENE'], pvals_df['P']))

analyzer = scPPIN()
analyzer.load_network(network_df, weight_column='Score')

analyzer.set_node_weights(pvalues)
analyzer.detect_module(fdr=0.01, edge_weight_attr='weight', normalization='minmax')

# Compute a custom layout before plotting
layout = nx.fruchterman_reingold_layout(analyzer.module)

# Pass the layout to plot_module
analyzer.plot_module(layout=layout, figsize=(10, 8), node_size=500)
plt.savefig('tests/module_plot.png', dpi=300, bbox_inches='tight')
plt.close()


print(f"Module nodes: {analyzer.module.number_of_nodes()}")
print(f"Module edges: {analyzer.module.number_of_edges()}")
print(f"Nodes (sorted): {sorted(list(analyzer.module.nodes()))[:10]}")

# Calculate node metrics
print("\n" + "="*60)
print("NODE METRICS")
print("="*60)

# Degree centrality (how many neighbors each node has)
degree_centrality = nx.degree_centrality(analyzer.module)
print("\nDegree Centrality (normalized by number of nodes):")
for node, centrality in sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f"  {node}: {centrality:.4f} (degree: {analyzer.module.degree(node)})")

# Betweenness centrality (how often a node lies on shortest paths)
betweenness_centrality = nx.betweenness_centrality(analyzer.module)
print("\nBetweenness Centrality (importance in paths):")
for node, centrality in sorted(betweenness_centrality.items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f"  {node}: {centrality:.4f}")

# Closeness centrality (average distance to all other nodes)
closeness_centrality = nx.closeness_centrality(analyzer.module)
print("\nCloseness Centrality (proximity to all nodes):")
for node, centrality in sorted(closeness_centrality.items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f"  {node}: {centrality:.4f}")

# Eigenvector centrality (influence based on connections)
try:
    eigenvector_centrality = nx.eigenvector_centrality(analyzer.module, max_iter=100)
    print("\nEigenvector Centrality (influence via network):")
    for node, centrality in sorted(eigenvector_centrality.items(), key=lambda x: x[1], reverse=True)[:5]:
        print(f"  {node}: {centrality:.4f}")
except:
    print("\nEigenvector Centrality: Could not compute (disconnected graph)")

# Clustering coefficient (how clique-like neighborhoods are)
clustering_coeff = nx.clustering(analyzer.module)
print("\nClustering Coefficient (local cliqueness):")
for node, coeff in sorted(clustering_coeff.items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f"  {node}: {coeff:.4f}")

# Average shortest path length (if connected)
if nx.is_connected(analyzer.module):
    avg_shortest_path = nx.average_shortest_path_length(analyzer.module)
    print(f"\nAverage Shortest Path Length: {avg_shortest_path:.4f}")
else:
    print(f"\nGraph is disconnected ({nx.number_connected_components(analyzer.module)} components)")
    for i, component in enumerate(nx.connected_components(analyzer.module)):
        subgraph = analyzer.module.subgraph(component)
        if subgraph.number_of_nodes() > 1:
            avg_path = nx.average_shortest_path_length(subgraph)
            print(f"  Component {i}: {len(component)} nodes, avg path length: {avg_path:.4f}")

# Density (proportion of possible edges present)
density = nx.density(analyzer.module)
print(f"\nNetwork Density: {density:.4f}")

# Assortativity (tendency of nodes with high degree to connect)
try:
    assortativity = nx.degree_assortativity_coefficient(analyzer.module)
    print(f"Degree Assortativity: {assortativity:.4f}")
except:
    print("Degree Assortativity: Could not compute")

# Node scores summary
print("\n" + "="*60)
print("NODE SCORES SUMMARY")
print("="*60)
scores_list = [analyzer.module.nodes[node].get('score', 0.0) for node in analyzer.module.nodes()]
print(f"Score range: [{min(scores_list):.6f}, {max(scores_list):.6f}]")
print(f"Score mean: {sum(scores_list) / len(scores_list):.6f}")
print(f"Score median: {sorted(scores_list)[len(scores_list)//2]:.6f}")