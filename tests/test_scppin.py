import pandas as pd
from scppin import scPPIN
import matplotlib.pyplot as plt
import igraph as ig
import numpy as np
from collections import Counter


network_df = pd.read_csv('tests/raw_mg_ad_fava_string_network.tsv', sep='\t', index_col=0)
# network_df['merge_score'] = (network_df['fava_score'] + network_df['string_score']) / 2
network_df['merge_score'] = (network_df['fava_score'] * network_df['string_score'])
# network_df['merge_score'] = (network_df['fava_score'] * network_df['string_score']) ** 0.5


pvals_df = pd.read_csv('tests/magma_gene_symbol_results.tsv', sep='\t')
pvalues = dict(zip(pvals_df['GENE'], pvals_df['P']))

analyzer = scPPIN()
analyzer.load_network(network_df, weight_column='merge_score')
analyzer.set_node_weights(pvalues)
analyzer.detect_module(fdr=0.01, edge_weight_attr='weight', normalization=None)

# Define different layout algorithms to compare
layouts_to_try = {
    'fruchterman_reingold': analyzer.module.layout_fruchterman_reingold(),
    'kamada_kawai': analyzer.module.layout_kamada_kawai(),
    'drl': analyzer.module.layout_drl(),
    'circle': analyzer.module.layout_circle(),
    'mds': analyzer.module.layout_mds(),
}

# Plot and save each layout
for layout_name, layout_coords in layouts_to_try.items():
    layout = {analyzer.module.vs[i]['name']: layout_coords[i] for i in range(analyzer.module.vcount())}
    ax = analyzer.plot_module(layout=layout, figsize=(10, 8), node_size=500)
    
    # Get the figure from the axes and save it
    fig = ax.figure
    fig.savefig(f'tests/module_plot_{layout_name}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved layout: module_plot_{layout_name}.png")

print(f"Module connected: {analyzer.module.is_connected()}")
print(f"Components: {len(analyzer.module.components())}")

# Get prizes from module and min-max normalize them for Leiden
# node_names = analyzer.module.vs['name']
# prizes = [v['prize'] if v['prize'] is not None else 0.0 for v in analyzer.module.vs]

# # Min-max normalize prizes
# prizes_array = np.array(prizes)
# if prizes_array.max() > prizes_array.min():
#     prizes_norm = (prizes_array - prizes_array.min()) / (prizes_array.max() - prizes_array.min())
# else:
#     prizes_norm = np.ones_like(prizes_array) * 0.5  # All same value
# analyzer.module.vs['prize_norm'] = prizes_norm.tolist()

# # Check if edge weights exist (should be merge_score)
# has_edge_weights = 'weight' in analyzer.module.es.attributes()
# print(f"\nEdge weights (merge_score) available: {has_edge_weights}")

# # Run Leiden with both node and edge weights
# if has_edge_weights:
#     # Leiden with both edge weights (merge_score) and node weights (normalized prizes)
#     communities = analyzer.module.community_leiden(
#         weights='weight',
#         node_weights='prize_norm',
#         resolution_parameter=0.1
#     )
# else:
#     # Leiden with only node weights
#     communities = analyzer.module.community_leiden(node_weights='prize_norm')

# membership = communities.membership

# # Display results
# print(f"\nLeiden Community Detection Results:")
# print(f"Number of communities: {len(communities)}")
# print(f"\nCommunity sizes:")
# community_counts = Counter(membership)
# for comm_id, count in sorted(community_counts.items()):
#     print(f"  Community {comm_id}: {count} nodes")

# # Show nodes in each community
# print(f"\nNodes per community:")
# for comm_id in sorted(community_counts.keys()):
#     comm_nodes = [node_names[i] for i, m in enumerate(membership) if m == comm_id]
#     print(f"  Community {comm_id}: {comm_nodes[:10]}{'...' if len(comm_nodes) > 10 else ''}")

print(f"\nModule nodes: {analyzer.module.vcount()}")
print(f"Module edges: {analyzer.module.ecount()}")
