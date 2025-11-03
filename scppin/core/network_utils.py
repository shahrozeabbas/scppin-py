"""Network loading, filtering, and manipulation utilities."""

import networkx as nx
from typing import Dict, List, Optional, Set
import warnings
from pathlib import Path


def load_ppin(
    filepath: Optional[str] = None,
    format: str = 'graphml'
) -> nx.Graph:
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
    nx.Graph
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
        network = nx.read_graphml(filepath)
    elif format == 'gml':
        network = nx.read_gml(filepath)
    elif format == 'edgelist':
        network = nx.read_edgelist(filepath)
    else:
        raise ValueError(f"Unknown format: {format}. "
                        "Use 'graphml', 'gml', or 'edgelist'.")
    
    # Convert to undirected if needed
    if isinstance(network, nx.DiGraph):
        network = network.to_undirected()
    
    return network


def simplify_network(network: nx.Graph, copy: bool = True) -> nx.Graph:
    """
    Simplify network by removing self-loops and parallel edges.
    
    Parameters
    ----------
    network : nx.Graph
        Network to simplify
    copy : bool, optional
        If True, work on a copy. If False, modify in place. (default: True)
        Note: MultiGraph conversion always creates a new object regardless of copy flag.
    
    Returns
    -------
    nx.Graph
        Simplified network
    """
    # If it's a MultiGraph, we need to convert anyway (creates new object)
    if isinstance(network, nx.MultiGraph):
        network_simple = nx.Graph(network)  # Conversion creates new Graph
        # Remove self-loops from the new graph
        network_simple.remove_edges_from(list(nx.selfloop_edges(network_simple)))
        return network_simple
    
    # For regular Graph, we can work in-place if copy=False
    if copy:
        network_simple = network.copy()
    else:
        network_simple = network
    
    # Remove self-loops (in-place operation)
    network_simple.remove_edges_from(list(nx.selfloop_edges(network_simple)))
    
    return network_simple


def filter_network(
    network: nx.Graph,
    genes: Set[str],
    keep_only_genes: bool = True
) -> nx.Graph:
    """
    Filter network to include only specified genes.
    
    Parameters
    ----------
    network : nx.Graph
        Network to filter
    genes : Set[str]
        Set of gene names to keep
    keep_only_genes : bool, optional
        If True, keep only nodes in genes. If False, remove nodes in genes.
        (default: True)
        
    Returns
    -------
    nx.Graph
        Filtered network (creates a copy)
    """
    if keep_only_genes:
        # Keep only nodes that are in genes
        nodes_to_keep = [node for node in network.nodes() if str(node) in genes]
    else:
        # Remove nodes that are in genes
        nodes_to_keep = [node for node in network.nodes() if str(node) not in genes]
    
    if not nodes_to_keep:
        warnings.warn("No nodes remain after filtering")
    
    return network.subgraph(nodes_to_keep).copy()


def filter_network_by_pvalues(
    network: nx.Graph,
    pvalues: Dict[str, float],
    missing_data_score: bool = False
) -> nx.Graph:
    """
    Filter network to genes with p-values (or keep all if missing_data_score=True).
    
    Parameters
    ----------
    network : nx.Graph
        Network to filter
    pvalues : Dict[str, float]
        Dictionary of gene p-values
    missing_data_score : bool, optional
        If True, keep all nodes even without p-values (default: False)
        
    Returns
    -------
    nx.Graph
        Filtered network
    """
    if missing_data_score:
        # Keep all nodes - no filtering needed, but return copy for safety
        return network.copy()
    else:
        # Optimize: directly filter without calling filter_network to avoid double copy
        genes_with_pvalues = set(pvalues.keys())
        nodes_to_keep = [
            node for node in network.nodes() 
            if str(node) in genes_with_pvalues
        ]
        
        if not nodes_to_keep:
            warnings.warn("No nodes remain after filtering")
        
        return network.subgraph(nodes_to_keep).copy()


def get_largest_connected_component(network: nx.Graph) -> nx.Graph:
    """
    Extract the largest connected component from the network.
    
    Parameters
    ----------
    network : nx.Graph
        Network (possibly disconnected)
        
    Returns
    -------
    nx.Graph
        Largest connected component
    """
    if not network.nodes():
        return network
    
    # Get connected components
    components = list(nx.connected_components(network))
    
    if len(components) == 1:
        return network
    
    # Find largest component
    largest_component = max(components, key=len)
    
    warnings.warn(f"Network has {len(components)} connected components. "
                 f"Returning largest with {len(largest_component)} nodes.")
    
    return network.subgraph(largest_component).copy()


def network_statistics(network: nx.Graph) -> Dict:
    """
    Compute basic network statistics.
    
    Parameters
    ----------
    network : nx.Graph
        Network to analyze
        
    Returns
    -------
    Dict
        Dictionary with network statistics
    """
    stats = {
        'num_nodes': network.number_of_nodes(),
        'num_edges': network.number_of_edges(),
        'density': nx.density(network),
        'num_components': nx.number_connected_components(network),
    }
    
    if stats['num_nodes'] > 0:
        degrees = [d for n, d in network.degree()]
        stats['avg_degree'] = sum(degrees) / len(degrees)
        stats['max_degree'] = max(degrees)
        stats['min_degree'] = min(degrees)
    
    return stats


def validate_network(network: nx.Graph) -> bool:
    """
    Validate that network is suitable for scPPIN analysis.
    
    Parameters
    ----------
    network : nx.Graph
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
    if network.number_of_nodes() == 0:
        raise ValueError("Network has no nodes")
    
    if network.number_of_edges() == 0:
        raise ValueError("Network has no edges")
    
    if nx.number_connected_components(network) > 1:
        warnings.warn(f"Network has {nx.number_connected_components(network)} "
                     "connected components. Consider using get_largest_connected_component()")
    
    # Check for self-loops
    num_selfloops = nx.number_of_selfloops(network)
    if num_selfloops > 0:
        warnings.warn(f"Network has {num_selfloops} self-loops. "
                     "Consider using simplify_network()")
    
    return True


def add_pvalues_to_network(
    network: nx.Graph,
    pvalues: Dict[str, float],
    attr_name: str = 'pvalue'
) -> nx.Graph:
    """
    Add p-values as node attributes to the network.
    
    Parameters
    ----------
    network : nx.Graph
        Network to modify
    pvalues : Dict[str, float]
        Dictionary mapping gene names to p-values
    attr_name : str, optional
        Name of the node attribute (default: 'pvalue')
        
    Returns
    -------
    nx.Graph
        Network with p-values added (modifies in place and returns)
    """
    # Add p-values to nodes
    for node in network.nodes():
        node_str = str(node)
        if node_str in pvalues:
            network.nodes[node][attr_name] = pvalues[node_str]
    
    return network

