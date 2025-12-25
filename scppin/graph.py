"""Graph building from edge lists."""

import igraph as ig
import pandas as pd
import numpy as np
from typing import Union, List, Tuple, Optional, Dict
from pathlib import Path


def _build_graph(
    edges: Union[str, List[Tuple], pd.DataFrame],
    weights: Optional[Union[str, Dict, List]] = None,
    directed: bool = False
) -> ig.Graph:
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
    ig.Graph
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
    # Step 1: Normalize all inputs to DataFrame
    edges_df = _normalize_to_dataframe(edges)
    
    # Step 2: Prepare edge list and weights for igraph
    edge_tuples, edge_weights = _prepare_igraph_edges(edges_df, weights)
    
    # Step 3: Create graph using igraph
    G = ig.Graph.TupleList(edge_tuples, directed=directed, vertex_name_attr='name')
    
    # Step 4: Add edge weights if provided (batch operation)
    if edge_weights is not None:
        G.es['weight'] = edge_weights
    
    return G


def _normalize_to_dataframe(
    edges: Union[str, List[Tuple], pd.DataFrame]
) -> pd.DataFrame:
    """
    Normalize all input types to DataFrame format.
    
    Parameters
    ----------
    edges : str, list of tuples, or DataFrame
        Edge list input
        
    Returns
    -------
    pd.DataFrame
        DataFrame with 'source' and 'target' columns
    """
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
    
    return edges_df


def _prepare_igraph_edges(
    edges_df: pd.DataFrame,
    weights: Optional[Union[str, Dict, List]] = None
) -> Tuple[List[Tuple], Optional[List[float]]]:
    """
    Convert DataFrame to edge list format for igraph and extract weights.
    
    Parameters
    ----------
    edges_df : pd.DataFrame
        DataFrame with 'source' and 'target' columns
    weights : str, dict, or list, optional
        Edge weights specification
        
    Returns
    -------
    Tuple[List[Tuple], Optional[List[float]]]
        Edge tuples list and optional weights list
    """
    n_edges = len(edges_df)
    
    # Pre-extract weight data based on type (validate once, not per row)
    weight_values = None
    if weights is not None:
        if isinstance(weights, str):
            if weights not in edges_df.columns:
                raise ValueError(f"Weight column '{weights}' not found in edge list")
            weight_values = edges_df[weights].values  # Extract as numpy array
        elif isinstance(weights, list):
            if len(weights) != n_edges:
                raise ValueError(f"Weight list length ({len(weights)}) must match "
                               f"number of edges ({n_edges})")
            weight_values = weights  # Use directly
    
    edge_tuples = []
    edge_weights_list = []
    
    # Convert DataFrame rows to edge tuples
    for idx, row in enumerate(edges_df.itertuples(index=False)):
        source = str(row.source)
        target = str(row.target)
        edge_tuples.append((source, target))
        
        # Determine weight for this edge
        weight = None
        if weights is not None:
            if isinstance(weights, (str, list)):
                weight = float(weight_values[idx])
            elif isinstance(weights, dict):
                # Weight dictionary - lookup by edge tuple
                edge_key = (source, target)
                weight = weights.get(edge_key) or weights.get((target, source))
                if weight is not None:
                    weight = float(weight)
        
        edge_weights_list.append(weight)
    
    # Return weights only if at least one weight was found
    if weights is not None and any(w is not None for w in edge_weights_list):
        return edge_tuples, edge_weights_list
    else:
        return edge_tuples, None


def _load_edges_from_file(filepath: Union[str, Path]) -> pd.DataFrame:
    """Load edge list from file."""
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    # Determine delimiter
    if filepath.suffix == '.csv':
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
                delimiter = r'\s+'  # whitespace (raw string for regex)
    else:
        # Default to comma
        delimiter = ','
    
    # Load file
    try:
        df = pd.read_csv(filepath, sep=delimiter)
    except Exception as e:
        raise ValueError(f"Failed to load edge list from {filepath}: {e}")
    
    return df
