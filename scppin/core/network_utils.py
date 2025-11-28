"""Network loading, filtering, and manipulation utilities."""

import igraph as ig
import numpy as np
from typing import Dict, List, Optional, Set, Tuple
import warnings
from pathlib import Path


def load_ppin(
    filepath: Optional[str] = None,
    format: str = 'graphml'
) -> ig.Graph:
    """
    Load protein-protein interaction network.
    
    Parameters
    ----------
    filepath : str, optional
        Path to network file. If None, loads default BioGRID human PPIN.
    format : str, optional
        Network file format: 'graphml', 'gml', 'edgelist' (default: 'graphml')
        
    Returns
    -------
    ig.Graph
        Protein-protein interaction network
    """
    if filepath is None:
        # Try to load default network from package data
        try:
            from importlib import resources
            # Try to find the default network in package data
            data_path = Path(__file__).parent.parent / 'data' / 'networks'
            default_file = data_path / 'biogridHomoSapiens3.5.166.graphml'
            
            if default_file.exists():
                filepath = str(default_file)
            else:
                raise FileNotFoundError(
                    "Default PPIN not found. Please provide filepath or copy "
                    "network files to scppin/data/networks/"
                )
        except Exception as e:
            raise FileNotFoundError(
                f"Could not load default PPIN: {e}\n"
                "Please provide filepath explicitly."
            )
    
    # Load network based on format
    if format == 'graphml':
        network = ig.Graph.Read_GraphML(filepath)
    elif format == 'gml':
        network = ig.Graph.Read_GML(filepath)
    elif format == 'edgelist':
        # For edgelist, read as simple edge list
        # Read file and create graph from tuples
        edge_tuples = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        edge_tuples.append((parts[0], parts[1]))
        network = ig.Graph.TupleList(edge_tuples, directed=False, vertex_name_attr='name')
    else:
        raise ValueError(f"Unknown format: {format}. "
                        "Use 'graphml', 'gml', or 'edgelist'.")
    
    # Ensure graph is undirected (igraph graphs are undirected by default, but check)
    if network.is_directed():
        network.to_undirected(mode='collapse')
    
    # Ensure node names are stored in 'name' attribute
    # GraphML/GML files may have 'id' or 'label' attributes
    if 'name' in network.vs.attributes():
        # Already has name attribute, ensure they're strings
        network.vs['name'] = [str(n) for n in network.vs['name']]
    elif 'id' in network.vs.attributes():
        # Use 'id' attribute as name
        network.vs['name'] = [str(n) for n in network.vs['id']]
    elif 'label' in network.vs.attributes():
        # Use 'label' attribute as name
        network.vs['name'] = [str(n) for n in network.vs['label']]
    else:
        # If no name/id/label attribute, use vertex indices as names
        network.vs['name'] = [str(i) for i in range(network.vcount())]
    
    return network


def simplify_network(network: ig.Graph, copy: bool = True) -> ig.Graph:
    """
    Simplify network by removing self-loops and parallel edges.
    
    Parameters
    ----------
    network : ig.Graph
        Network to simplify
    copy : bool, optional
        If True, work on a copy. If False, modify in place. (default: True)
    
    Returns
    -------
    ig.Graph
        Simplified network
    """
    if copy:
        network_simple = network.copy()
    else:
        network_simple = network
    
    # Simplify removes self-loops and parallel edges
    # Preserve edge attributes by taking the maximum weight for parallel edges
    edge_attrs = network_simple.es.attributes()
    if edge_attrs:
        # For each attribute, use 'max' to combine parallel edges (keeps highest weight)
        combine_edges = {attr: 'max' for attr in edge_attrs}
        network_simple.simplify(multiple=True, loops=True, combine_edges=combine_edges)
    else:
        network_simple.simplify(multiple=True, loops=True)
    
    return network_simple


def filter_network(
    network: ig.Graph,
    genes: Set[str],
    keep_only_genes: bool = True
) -> ig.Graph:
    """
    Filter network to include only specified genes.
    
    Parameters
    ----------
    network : ig.Graph
        Network to filter
    genes : Set[str]
        Set of gene names to keep
    keep_only_genes : bool, optional
        If True, keep only nodes in genes. If False, remove nodes in genes.
        (default: True)
        
    Returns
    -------
    ig.Graph
        Filtered network (creates a copy)
    """
    if keep_only_genes:
        # Keep only nodes that are in genes
        vertices_to_keep = [v.index for v in network.vs if v['name'] in genes]
    else:
        # Remove nodes that are in genes
        vertices_to_keep = [v.index for v in network.vs if v['name'] not in genes]
    
    if not vertices_to_keep:
        warnings.warn("No nodes remain after filtering", UserWarning, stacklevel=2)
        return ig.Graph()
    
    return network.subgraph(vertices_to_keep)


def filter_network_by_pvalues(
    network: ig.Graph,
    pvalues: Dict[str, float],
    missing_data_score: bool = False
) -> ig.Graph:
    """
    Filter network to genes with p-values (or keep all if missing_data_score=True).
    
    Parameters
    ----------
    network : ig.Graph
        Network to filter
    pvalues : Dict[str, float]
        Dictionary of gene p-values
    missing_data_score : bool, optional
        If True, keep all nodes even without p-values (default: False)
        
    Returns
    -------
    ig.Graph
        Filtered network
    """
    if missing_data_score:
        # Keep all nodes - no filtering needed, but return copy for safety
        return network.copy()
    else:
        # Filter to genes with p-values
        genes_with_pvalues = set(pvalues.keys())
        vertices_to_keep = [
            v.index for v in network.vs 
            if v['name'] in genes_with_pvalues
        ]
        
        if not vertices_to_keep:
            warnings.warn("No nodes remain after filtering", UserWarning, stacklevel=2)
            return ig.Graph()
        
        return network.subgraph(vertices_to_keep)


def get_largest_connected_component(network: ig.Graph) -> ig.Graph:
    """
    Extract the largest connected component from the network.
    
    Parameters
    ----------
    network : ig.Graph
        Network (possibly disconnected)
        
    Returns
    -------
    ig.Graph
        Largest connected component
    """
    if network.vcount() == 0:
        return network
    
    # Get connected components
    components = network.components()
    
    if len(components) == 1:
        return network
    
    # Find largest component
    largest_component = max(components, key=len)
    
    warnings.warn(
        f"Extracted largest component: {len(largest_component)} nodes "
        f"(from {len(components)} components)",
        UserWarning,
        stacklevel=2
    )
    
    return network.subgraph(largest_component)


def network_statistics(network: ig.Graph) -> Dict:
    """
    Compute comprehensive network statistics.
    
    Parameters
    ----------
    network : ig.Graph
        Network to analyze
        
    Returns
    -------
    Dict
        Dictionary with network statistics including:
        - Basic stats: num_nodes, num_edges, density, num_components
        - Degree statistics: avg_degree, max_degree, min_degree
        - Clustering: avg_clustering_coefficient
        - Path metrics: avg_shortest_path_length, diameter (if connected)
        - Centrality: avg_degree_centrality, avg_betweenness_centrality
    """
    stats = {
        'num_nodes': network.vcount(),
        'num_edges': network.ecount(),
        'density': network.density(),
        'num_components': len(network.components()),
    }
    
    if stats['num_nodes'] == 0:
        return stats
    
    # Degree statistics
    degrees = network.degree()
    stats['avg_degree'] = np.mean(degrees) if degrees else 0.0
    stats['max_degree'] = max(degrees) if degrees else 0
    stats['min_degree'] = min(degrees) if degrees else 0
    
    # Clustering coefficient (local clustering)
    clustering = network.transitivity_local_undirected(mode='zero')
    stats['avg_clustering_coefficient'] = np.mean(clustering) if clustering else 0.0
    
    # Centrality measures
    n = stats['num_nodes']
    if n > 1:
        degree_centrality = [d / (n - 1) for d in degrees]
        stats['avg_degree_centrality'] = np.mean(degree_centrality) if degree_centrality else 0.0
    else:
        stats['avg_degree_centrality'] = 0.0
    
    betweenness = network.betweenness()
    stats['avg_betweenness_centrality'] = np.mean(betweenness) if betweenness else 0.0
    
    # Path metrics (only if connected)
    if network.is_connected():
        stats['avg_shortest_path_length'] = network.average_path_length()
        stats['diameter'] = network.diameter()
    else:
        stats['avg_shortest_path_length'] = None
        stats['diameter'] = None
    
    return stats


def validate_network(network: ig.Graph) -> bool:
    """
    Validate that network is suitable for scPPIN analysis.
    
    Parameters
    ----------
    network : ig.Graph
        Network to validate
        
    Returns
    -------
    bool
        True if valid
        
    Raises
    ------
    ValueError
        If network is invalid
    """
    if network.vcount() == 0:
        raise ValueError("Network has no nodes")
    
    if network.ecount() == 0:
        raise ValueError("Network has no edges")
    
    num_components = len(network.components())
    if num_components > 1:
        warnings.warn(
            f"Network has {num_components} connected components. "
            f"Consider using get_largest_connected_component() to extract the main component.",
            UserWarning,
            stacklevel=2
        )
    
    # Check for self-loops
    num_selfloops = sum(1 for e in network.es if e.source == e.target)
    if num_selfloops > 0:
        warnings.warn(
            f"Network has {num_selfloops} self-loops. Consider using simplify_network().",
            UserWarning,
            stacklevel=2
        )
    
    # Check for multiple edges
    if network.has_multiple():
        warnings.warn(
            "Network has multiple edges. Consider using simplify_network().",
            UserWarning,
            stacklevel=2
        )
    
    return True


def add_pvalues_to_network(
    network: ig.Graph,
    pvalues: Dict[str, float],
    attr_name: str = 'pvalue'
) -> ig.Graph:
    """
    Add p-values as node attributes to the network.
    
    Parameters
    ----------
    network : ig.Graph
        Network to modify
    pvalues : Dict[str, float]
        Dictionary mapping gene names to p-values
    attr_name : str, optional
        Name of the node attribute (default: 'pvalue')
        
    Returns
    -------
    ig.Graph
        Network with p-values added (modifies in place and returns)
    """
    # Initialize attribute with None
    pvalue_list = [None] * network.vcount()
    
    # Set p-values for nodes that have them
    for v in network.vs:
        node_name = v['name']
        if node_name in pvalues:
            pvalue_list[v.index] = pvalues[node_name]
    
    # Batch assign
    network.vs[attr_name] = pvalue_list
    
    return network

