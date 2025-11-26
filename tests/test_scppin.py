import pandas as pd
from scppin import scPPIN
import matplotlib.pyplot as plt
import networkx as nx


network_df = pd.read_csv('tests/raw_mg_ad_fava_network.tsv', sep='\t', index_col=0)

pvals_df = pd.read_csv('tests/magma_gene_symbol_results.tsv', sep='\t')
pvalues = dict(zip(pvals_df['GENE'], pvals_df['P']))

analyzer = scPPIN()
analyzer.load_network(network_df, weight_column='Score')

analyzer.set_node_weights(pvalues)
analyzer.detect_module(fdr=0.01, edge_weight_attr='weight', normalization='power')

# Compute a custom layout before plotting
layout = nx.fruchterman_reingold_layout(analyzer.module)

# Pass the layout to plot_module
ax = analyzer.plot_module(layout=layout, figsize=(10, 8), node_size=500)

# Get the figure from the axes and save it
fig = ax.figure
fig.savefig('tests/module_plot.png', dpi=300, bbox_inches='tight')
plt.close(fig)


print(f"Module nodes: {analyzer.module.number_of_nodes()}")
print(f"Module edges: {analyzer.module.number_of_edges()}")
print(f"Nodes (sorted): {sorted(list(analyzer.module.nodes()))[:10]}")
