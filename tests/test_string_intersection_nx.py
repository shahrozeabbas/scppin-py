"""Test STRING PPI intersection using NetworkX graphs."""

import pandas as pd
import networkx as nx
import stringdb

print("=" * 70)
print("NetworkX-based STRING PPI Intersection")
print("=" * 70)

# 1. Load network graph
network_df = pd.read_csv('tests/raw_mg_ad_fava_network.tsv', sep='\t', index_col=0)
net_g = nx.from_pandas_edgelist(network_df, 'Protein_1', 'Protein_2', 'Score', create_using=nx.Graph())
print(f"Network loaded: {net_g.number_of_nodes()} nodes, {net_g.number_of_edges()} edges")

# Load genes from magma file and filter to network genes
pvals_df = pd.read_csv('tests/magma_gene_symbol_results.tsv', sep='\t')
magma_genes = pvals_df['GENE'].tolist()
network_genes_set = set(net_g.nodes())
genes_to_query = [g for g in magma_genes if g in network_genes_set]
print(f"Magma file: {len(magma_genes)} genes, {len(genes_to_query)} in network")

# 2. Batch query STRING for genes from magma file (human: species=9606)
batch_size = 500
species = 9606  # Human (Homo sapiens)
all_string_interactions = []
all_string_ids = []

print(f"Querying STRING for {len(genes_to_query)} genes from magma file (species={species}) in batches of {batch_size}...")
for i in range(0, len(genes_to_query), batch_size):
    batch = genes_to_query[i:i+batch_size]
    print(f"  Batch {i//batch_size + 1}/{(len(genes_to_query)-1)//batch_size + 1} ({len(batch)} genes)...", end=' ')
    try:
        string_ids = stringdb.get_string_ids(batch, species=species)
        string_interactions = stringdb.get_network(string_ids.queryItem, species=species)
        all_string_interactions.append(string_interactions)
        all_string_ids.append(string_ids)
        print(f"{len(string_interactions)} interactions")
    except Exception as e:
        print(f"Error: {e}")
        exit(1)

# Combine batch results
string_interactions_df = pd.concat(all_string_interactions, ignore_index=True)
string_ids_df = pd.concat(all_string_ids, ignore_index=True)
print(f"Total STRING interactions: {len(string_interactions_df)}")

# 3. Build STRING graph with relabeling
string_to_symbol = dict(zip(string_ids_df['stringId'], string_ids_df['queryItem']))
attr_cols = ['score']
if 'confidence_score' in string_interactions_df.columns:
    attr_cols.append('confidence_score')
string_g = nx.from_pandas_edgelist(string_interactions_df, 'stringId_A', 'stringId_B', 
                                    attr_cols, create_using=nx.Graph())
string_g = nx.relabel_nodes(string_g, string_to_symbol, copy=True)

# Verify gene symbol matching
print(f"\nVerifying gene symbol matching:")
print(f"  Network uses gene symbols: {list(net_g.nodes())[:5]}")
print(f"  STRING after relabeling uses gene symbols: {list(string_g.nodes())[:5]}")
print(f"  Mapping verification (first 3):")
for ensp_id, symbol in list(string_to_symbol.items())[:3]:
    print(f"    {ensp_id} -> {symbol}")

# Rename edge attributes for consistency
for u, v, d in string_g.edges(data=True):
    d['string_score'] = d.pop('score', 0)
    d['confidence'] = d.pop('confidence_score', 0)

# 4. Compute intersection using NetworkX built-in function
intersect_g = nx.intersection(net_g, string_g)

# Copy attributes from both graphs to intersection edges
for gene1, gene2 in intersect_g.edges():
    # Network attributes (from first graph - already preserved as 'Score')
    if 'Score' in net_g[gene1][gene2]:
        intersect_g[gene1][gene2]['network_score'] = net_g[gene1][gene2]['Score']
    
    # STRING attributes (need to copy from second graph)
    if string_g.has_edge(gene1, gene2):
        intersect_g[gene1][gene2]['string_score'] = string_g[gene1][gene2].get('string_score', 0)
        intersect_g[gene1][gene2]['confidence'] = string_g[gene1][gene2].get('confidence', 0)

# Verify edge format matches
print(f"\nEdge format verification:")
print(f"  Network edge sample: {list(net_g.edges())[:3] if net_g.number_of_edges() > 0 else 'No edges'}")
print(f"  STRING edge sample: {list(string_g.edges())[:3] if string_g.number_of_edges() > 0 else 'No edges'}")
print(f"  Intersection edges: {intersect_g.number_of_edges()}")

# 5. Summary
common_nodes = set(net_g.nodes()) & set(string_g.nodes())
print("\n" + "=" * 70)
print("Summary:")
print(f"  Common nodes: {len(common_nodes)}")
print(f"  Network nodes: {net_g.number_of_nodes()}, edges: {net_g.number_of_edges()}")
print(f"  STRING nodes: {string_g.number_of_nodes()}, edges: {string_g.number_of_edges()}")
print(f"  Intersection edges: {intersect_g.number_of_edges()}")
if net_g.number_of_edges() > 0:
    print(f"  Overlap: {intersect_g.number_of_edges() / net_g.number_of_edges() * 100:.1f}%")
print("=" * 70)
