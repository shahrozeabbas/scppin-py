import scppin
import pandas as pd


network_df = pd.read_csv('tests/raw_mg_ad_fava_string_network.tsv', sep='\t', index_col=0)
# network_df['merge_score'] = (network_df['fava_score'] + network_df['string_score']) / 2
network_df['merge_score'] = (network_df['fava_score'] * network_df['string_score'])
# network_df['merge_score'] = (network_df['fava_score'] * network_df['string_score']) ** 0.5


pvals_df = pd.read_csv('tests/magma_gene_symbol_results.tsv', sep='\t')
pvalues = dict(zip(pvals_df['GENE'], pvals_df['P']))

analyzer = scppin.scPPIN()
analyzer.load_network(network_df, weight_column='merge_score')

analyzer.set_node_weights(pvalues)
analyzer.detect_module(fdr=0.01, edge_weight_attr='weight', normalization=None, use_max_prize_root=True, missing_data_score)


print(f"Module connected: {analyzer.module.is_connected()}")
print(f"Components: {len(analyzer.module.components())}")


print(f"\nModule nodes: {analyzer.module.vcount()}")
print(f"Module edges: {analyzer.module.ecount()}")
