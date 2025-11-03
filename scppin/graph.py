"""Graph building from edge lists."""

import networkx as nx
import pandas as pd
from typing import Union, List, Tuple, Optional, Dict
from pathlib import Path


def _build_graph(
    edges: Union[str, List[Tuple], pd.DataFrame],
    weights: Optional[Union[str, Dict, List]] = None,
    directed: bool = False
) -> nx.Graph:
    """
    Build network from edge list with optional weights.
    
    Parameters
    ----------
    edges : str, list of tuples, or DataFrame
        Edge list. Can be:
        - Path to CSV/TXT file with columns: source, target, [weight]
        - List of tuples: [('A', 'B'), ('B', 'C'), ...]
        - DataFrame with columns: ['source', 'target'] or ['source', 'target', 'weight']
    weights : str, dict, or list, optional
        Edge weights. Can be:
        - Column name in file/DataFrame (str)
        - Dict mapping (source, target) -> weight
        - List of weights (same order as edges)
        If None, creates unweighted graph
    directed : bool, optional
        Create directed graph (default: False)
        
    Returns
    -------
    nx.Graph or nx.DiGraph
        Network with optional edge weights as 'weight' attribute
        
    Examples
    --------
    >>> # From file
    >>> network = scppin.build_graph('edges.csv')
    
    >>> # From file with weights column
    >>> network = scppin.build_graph('edges.csv', weights='confidence')
    
    >>> # From list of tuples
    >>> edges = [('TP53', 'MDM2'), ('TP53', 'CDKN1A')]
    >>> network = scppin.build_graph(edges)
    
    >>> # From DataFrame with weights
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'source': ['TP53', 'TP53'],
    ...     'target': ['MDM2', 'CDKN1A'],
    ...     'confidence': [0.95, 0.90]
    ... })
    >>> network = scppin.build_graph(df, weights='confidence')
    
    >>> # With weight dictionary
    >>> edges = [('A', 'B'), ('B', 'C')]
    >>> weights = {('A', 'B'): 0.9, ('B', 'C'): 0.8}
    >>> network = scppin.build_graph(edges, weights=weights)
    """
    # Create appropriate graph type
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    
    # Handle different edge input types
    if isinstance(edges, (str, Path)):
        # Load from file
        edges_df = _load_edges_from_file(edges)
    elif isinstance(edges, pd.DataFrame):
        edges_df = edges
    elif isinstance(edges, list):
        # Convert list of tuples to DataFrame
        edges_df = pd.DataFrame(edges, columns=['source', 'target'])
    else:
        raise TypeError(f"edges must be str, list, or DataFrame, got {type(edges)}")
    
    # Validate DataFrame has required columns
    if 'source' not in edges_df.columns or 'target' not in edges_df.columns:
        # Try first two columns
        if len(edges_df.columns) >= 2:
            edges_df = edges_df.rename(columns={
                edges_df.columns[0]: 'source',
                edges_df.columns[1]: 'target'
            })
        else:
            raise ValueError("Edge list must have at least 2 columns (source, target)")
    
    # Add edges
    for idx, row in edges_df.iterrows():
        source = str(row['source'])
        target = str(row['target'])
        G.add_edge(source, target)
    
    # Add weights if provided
    if weights is not None:
        _add_weights_to_graph(G, edges_df, weights)
    
    return G


def _load_edges_from_file(filepath: Union[str, Path]) -> pd.DataFrame:
    """Load edge list from file."""
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    # Determine delimiter
    if filepath.suffix in ['.csv']:
        delimiter = ','
    elif filepath.suffix in ['.tsv']:
        delimiter = '\t'
    elif filepath.suffix in ['.txt']:
        # Try to detect delimiter
        with open(filepath, 'r') as f:
            first_line = f.readline()
            if '\t' in first_line:
                delimiter = '\t'
            elif ',' in first_line:
                delimiter = ','
            else:
                delimiter = r'\s+'  # whitespace
    else:
        # Default to comma
        delimiter = ','
    
    # Load file
    try:
        df = pd.read_csv(filepath, sep=delimiter)
    except Exception as e:
        raise ValueError(f"Failed to load edge list from {filepath}: {e}")
    
    return df


def _add_weights_to_graph(
    G: nx.Graph,
    edges_df: pd.DataFrame,
    weights: Union[str, Dict, List]
) -> None:
    """Add weights to graph edges."""
    if isinstance(weights, str):
        # Weight column name
        if weights not in edges_df.columns:
            raise ValueError(f"Weight column '{weights}' not found in edge list")
        
        for idx, row in edges_df.iterrows():
            source = str(row['source'])
            target = str(row['target'])
            weight = float(row[weights])
            if G.has_edge(source, target):
                G[source][target]['weight'] = weight
    
    elif isinstance(weights, dict):
        # Weight dictionary
        for (source, target), weight in weights.items():
            source = str(source)
            target = str(target)
            if G.has_edge(source, target):
                G[source][target]['weight'] = float(weight)
            elif G.has_edge(target, source):  # Undirected
                G[target][source]['weight'] = float(weight)
    
    elif isinstance(weights, list):
        # Weight list (same order as edges)
        if len(weights) != len(edges_df):
            raise ValueError(f"Weight list length ({len(weights)}) must match "
                           f"number of edges ({len(edges_df)})")
        
        for idx, row in edges_df.iterrows():
            source = str(row['source'])
            target = str(row['target'])
            weight = float(weights[idx])
            if G.has_edge(source, target):
                G[source][target]['weight'] = weight
    
    else:
        raise TypeError(f"weights must be str, dict, or list, got {type(weights)}")

