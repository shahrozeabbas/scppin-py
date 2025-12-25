"""scPPIN analyzer class for stateful module detection."""

import igraph as ig
import numpy as np
import pandas as pd
from typing import Dict, Optional, Tuple, Union, List
import warnings

# Import internal helpers
from .graph import _build_graph, _load_edges_from_file
from .pvalues import _extract_pvalues
from .module import _detect_module, _validate_pvalues
from .core.network_utils import load_ppin as _load_ppin


def _normalize_gene_name(name: str, case_sensitive: bool = False) -> str:
    """
    Normalize gene name for consistent matching.
    
    Parameters
    ----------
    name : str
        Gene name to normalize
    case_sensitive : bool, optional
        If False, convert to uppercase (default: False)
        
    Returns
    -------
    str
        Normalized gene name
    """
    name = str(name).strip()
    if not case_sensitive:
        name = name.upper()
    return name


class scPPIN:
    """
    scPPIN analyzer for detecting functional modules in protein-protein interaction networks.
    
    This class manages network, node weights (p-values), edge weights, and analysis results
    as object attributes, ensuring data consistency and providing a clean workflow.
    
    Attributes
    ----------
    network : Optional[ig.Graph]
        Protein-protein interaction network filtered to genes with weights
    node_weights : Optional[Dict[str, float]]
        Node weights (p-values from differential expression)
    edge_weights : Dict[Tuple[str, str], float]
        Edge weights dictionary
    module : Optional[ig.Graph]
        Detected functional module
    
    Examples
    --------
    >>> import scppin
    >>> 
    >>> # Create analyzer
    >>> analyzer = scppin.scPPIN()
    >>> 
    >>> # Load network with weights from CSV column
    >>> analyzer.load_network('edges.csv', weight_column='confidence')
    >>> 
    >>> # Set node weights (p-values)
    >>> analyzer.set_node_weights({'TP53': 0.0001, 'MDM2': 0.001})
    >>> 
    >>> # Detect module
    >>> module = analyzer.detect_module(fdr=0.01)
    """
    
    def __init__(
        self,
        network: Optional[ig.Graph] = None,
        node_weights: Optional[Dict[str, float]] = None,
        edge_weights: Optional[Dict[Tuple[str, str], float]] = None
    ):
        """
        Initialize scPPIN analyzer.
        
        Parameters
        ----------
        network : Optional[ig.Graph]
            Initial network (optional)
        node_weights : Optional[Dict[str, float]]
            Initial node weights (optional)
        edge_weights : Optional[Dict[Tuple[str, str], float]]
            Initial edge weights dictionary (optional)
        """
        self.network = network
        self.node_weights = None  # Will be set by set_node_weights if provided
        self.edge_weights: Dict[Tuple[str, str], float] = {}
        self.module: Optional[ig.Graph] = None
        
        # Normalize network nodes if network provided
        if self.network is not None:
            self._normalize_network_nodes()
        
        # Set node weights (which will normalize and filter network)
        if node_weights is not None:
            self.set_node_weights(node_weights)
        
        # Set edge weights (must be after network normalization and node weight filtering)
        if edge_weights is not None:
            if self.network is None:
                warnings.warn(
                    'Edge_weights provided but no network. Edge weights will be ignored. '
                    'Load network first, then call set_edge_weights().'
                )
            else:
                self.set_edge_weights(weights=edge_weights)
    
    def _normalize_node_weights(self) -> None:
        """Normalize gene names in node_weights."""
        if self.node_weights is None:
            return
        
        normalized = {}
        for gene, weight in self.node_weights.items():
            norm_gene = _normalize_gene_name(gene)
            normalized[norm_gene] = weight
        
        self.node_weights = normalized
    
    def _filter_network_to_node_weights(self) -> None:
        """Filter network to genes with node weights."""
        if self.network is None or self.node_weights is None:
            return
        
        # Nodes are already normalized, so we can use direct lookup
        genes_with_weights = set(self.node_weights.keys())
        vertices_to_keep = [
            v.index for v in self.network.vs
            if v['name'] in genes_with_weights
        ]
        
        if vertices_to_keep:
            self.network = self.network.subgraph(vertices_to_keep)
        else:
            warnings.warn("No nodes in network match node_weights after normalization")
    
    def load_network(
        self,
        source: Union[str, List[Tuple], pd.DataFrame, ig.Graph],
        weight_column: Optional[str] = None,
        fmt: str = 'auto'
    ) -> 'scPPIN':
        """
        Load network from file, list, DataFrame, or igraph graph.
        
        Parameters
        ----------
        source : Union[str, List[Tuple], pd.DataFrame, ig.Graph]
            Network source:
            - String: Path to CSV/TXT/GraphML file
            - List: List of edge tuples
            - DataFrame: Edge list DataFrame
            - ig.Graph: Existing igraph graph
        weight_column : Optional[str]
            Column name in CSV/DataFrame to use as edge weights.
            If provided and source is file/DataFrame, loads weights from that column.
            Sets edge weights as 'weight' attribute on network edges.
        fmt : str, optional
            File format hint ('auto', 'csv', 'graphml', 'gml') (default: 'auto')
            
        Returns
        -------
        scPPIN
            self (for method chaining)
            
        Examples
        --------
        >>> analyzer = scppin.scPPIN()
        >>> analyzer.load_network('edges.csv')
        >>> analyzer.load_network('edges.csv', weight_column='confidence')
        >>> analyzer.load_network([('A', 'B'), ('B', 'C')])
        """
        # Handle igraph graph directly
        if isinstance(source, ig.Graph):
            self.network = source.copy()
        # Handle GraphML/GML files
        elif isinstance(source, str) and (fmt in ['graphml', 'gml'] or \
             source.endswith(('.graphml', '.gml'))):
            if fmt == 'auto':
                # Auto-detect format from extension
                if source.endswith('.graphml'):
                    file_fmt = 'graphml'
                elif source.endswith('.gml'):
                    file_fmt = 'gml'
                else:
                    file_fmt = 'graphml'  # default
            else:
                file_fmt = fmt
            self.network = _load_ppin(source, fmt=file_fmt)
        # Handle CSV/TXT/list/DataFrame using build_graph
        else:
            # If weight_column specified and source is file/DataFrame, use it
            weights_param = weight_column if weight_column else None
            self.network = _build_graph(source, weights=weights_param, directed=False)
        
        # Normalize gene names in network
        if self.network is not None:
            self._normalize_network_nodes()
        
        # Filter to node_weights if already set (after normalization)
        if self.node_weights:
            self._filter_network_to_node_weights()
        
        # Extract edge weights if weight_column was used
        if weight_column and self.network is not None:
            self._extract_edge_weights_from_network(attr_name='weight')
        
        return self
    
    def _normalize_network_nodes(self) -> None:
        """Normalize node names in network."""
        if self.network is None:
            return
        
        # Batch update node names
        new_names = [_normalize_gene_name(v['name']) for v in self.network.vs]
        self.network.vs['name'] = new_names
    
    def _extract_edge_weights_from_network(self, attr_name: str = 'weight') -> None:
        """Extract edge weights from network attributes to self.edge_weights dict."""
        if self.network is None:
            return
        
        # Nodes are already normalized, so we can use them directly
        self.edge_weights = {}
        node_names = self.network.vs['name']
        
        try:
            weights = self.network.es[attr_name]
            for e in self.network.es:
                u_name = node_names[e.source]
                v_name = node_names[e.target]
                weight = weights[e.index]
                if weight is not None:
                    self.edge_weights[(u_name, v_name)] = float(weight)
        except (KeyError, TypeError):
            # No weights attribute, leave edge_weights empty
            pass
    
    def set_node_weights(
        self,
        weights: Union[Dict[str, float], object],
        groupby: Optional[str] = None,
        group: Optional[str] = None
    ) -> 'scPPIN':
        """
        Set node weights (p-values) and filter network to genes with weights.
        
        Parameters
        ----------
        weights : Union[Dict[str, float], AnnData]
            Node weights:
            - Dict: Dictionary mapping gene names to weights (p-values)
            - AnnData: Extract from rank_genes_groups (requires groupby/group)
        groupby : Optional[str]
            Key in adata.obs for grouping labels (required if weights is AnnData)
        group : Optional[str]
            Specific group to extract (required if weights is AnnData)
            
        Returns
        -------
        scPPIN
            self (for method chaining)
            
        Examples
        --------
        >>> analyzer.set_node_weights({'TP53': 0.0001, 'MDM2': 0.001})
        >>> analyzer.set_node_weights(adata, groupby='louvain', group='0')
        
        Note
        ----
        This method automatically:
        - Normalizes gene names for matching
        - Filters network to only include genes with weights
        """
        # Handle AnnData input
        if hasattr(weights, 'var_names'):  # AnnData object
            if groupby is None or group is None:
                raise ValueError("groupby and group required when weights is AnnData")
            weights_dict = _extract_pvalues(weights, groupby, group)
        elif isinstance(weights, dict):
            weights_dict = weights
        else:
            raise TypeError(f"weights must be dict or AnnData, got {type(weights)}")
        
        # Validate p-values early (fail fast)
        _validate_pvalues(weights_dict)
        
        # Store node weights
        self.node_weights = weights_dict
        
        # Always normalize gene names
        self._normalize_node_weights()
        
        # Always filter network to genes with weights
        if self.network is not None:
            self._filter_network_to_node_weights()
        
        return self
    
    def set_edge_weights(
        self,
        weights: Dict[Tuple[str, str], float],
        attr_name: str = 'weight'
    ) -> 'scPPIN':
        """
        Set edge weights from user-provided dictionary.
        
        Parameters
        ----------
        weights : Dict[Tuple[str, str], float]
            User-provided edge weights dictionary.
            Sets these weights on network edges.
            Only edges in network are set (automatically filtered).
        attr_name : str, optional
            Edge attribute name to store weights (default: 'weight')
            
        Returns
        -------
        scPPIN
            self (for method chaining)
            
        Examples
        --------
        >>> # From dictionary
        >>> weights = {('TP53', 'MDM2'): 0.9, ('TP53', 'CDKN1A'): 0.8}
        >>> analyzer.set_edge_weights(weights=weights)
        """
        if self.network is None:
            raise ValueError('Network must be loaded before setting edge weights. '
                           'Call load_network() first.')
        
        if weights is None:
            raise ValueError('weights dictionary must be provided')
        
        filtered_weights = {}
        vertex_names = set(self.network.vs['name'])
        
        for (u, v), weight in weights.items():
            # Normalize user input to match network node names
            u_norm = _normalize_gene_name(str(u))
            v_norm = _normalize_gene_name(str(v))
            
            # Only proceed if both vertices exist in network
            if (u_norm in vertex_names and v_norm in vertex_names):
                eid = self.network.get_eid(u_norm, v_norm, directed=False, error=False)
                if eid != -1:
                    self.network.es[eid][attr_name] = float(weight)
                    filtered_weights[(u_norm, v_norm)] = float(weight)
        
        self.edge_weights = filtered_weights
        
        if not filtered_weights:
            warnings.warn('No edges from weights dictionary matched network edges')
        
        return self
    
    def detect_module(
        self,
        fdr: float = 0.01,
        edge_weight_attr: Optional[str] = None,
        c0: float = 0.01,
        normalization: Optional[str] = 'minmax',
        simplify: bool = True,
        validate: bool = True,
        use_max_prize_root: bool = False
    ) -> ig.Graph:
        """
        Detect functional module using PCST optimization.
        
        Parameters
        ----------
        fdr : float, optional
            False discovery rate threshold (default: 0.01)
        edge_weight_attr : Optional[str], optional
            Edge attribute name for weights (default: None = uniform costs).
            If None, uses uniform edge costs matching R implementation.
        c0 : Optional[float], optional
            Minimum edge cost (default: 0.01)
        normalization : Optional[str], optional
            Normalization method for edge weights: 'minmax', 'log1p', 'power', or None
            (default: 'minmax'). If None, uses weights directly without normalization
            (assumes weights are already in [0, 1] range). Only used when edge_weight_attr is provided.
        simplify : bool, optional
            Simplify network (default: True)
        validate : bool, optional
            Validate network (default: True)
        use_max_prize_root : bool, optional
            If True, use the node with highest prize as root (default: False)
            
        Returns
        -------
        ig.Graph
            Detected functional module (also stored on ``self.module``)
            
        Examples
        --------
        >>> # Default PCST
        >>> module = analyzer.detect_module(fdr=0.01)
        >>> 
        >>> # PCST with edge weights
        >>> module = analyzer.detect_module(fdr=0.01, edge_weight_attr='weight')
        >>> 
        >>> # PCST with max prize root (more deterministic)
        >>> module = analyzer.detect_module(
        ...     fdr=0.01, use_max_prize_root=True
        ... )
        """
        if self.network is None:
            raise ValueError('Network must be loaded. Call load_network() first.')
        
        if self.node_weights is None:
            raise ValueError('Node weights must be set. Call set_node_weights() first.')
        
        # Use internal detect_module function
        self.module = _detect_module(
            self.network,
            self.node_weights,
            fdr=fdr,
            edge_weight_attr=edge_weight_attr,
            c0=c0,
            normalization=normalization,
            simplify=simplify,
            validate=validate,
            use_max_prize_root=use_max_prize_root
        )
        
        return self.module
    
    def plot_module(self, fdr: float = 0.01, **kwargs):
        """
        Visualize detected module.
        
        Parameters
        ----------
        fdr : float, optional
            FDR threshold for visualization (default: 0.01)
        **kwargs
            Additional plotting arguments passed to plot_functional_module
            
        Returns
        -------
        matplotlib.figure.Figure
            Figure object
        """
        if self.module is None:
            raise ValueError("No module detected. Call detect_module() first.")
        
        from .visualization.plotting import _plot_functional_module
        return _plot_functional_module(self.module, fdr=fdr, **kwargs)
    
    def network_statistics(self, graph: Optional[ig.Graph] = None) -> Dict:
        """
        Compute comprehensive network statistics.
        
        Parameters
        ----------
        graph : Optional[ig.Graph], optional
            Network to analyze. If None, uses self.module if available,
            otherwise uses self.network (default: None)
            
        Returns
        -------
        Dict
            Dictionary with network statistics including:
            - Basic stats: num_nodes, num_edges, density, num_components
            - Degree statistics: avg_degree, max_degree, min_degree
            - Clustering: avg_clustering_coefficient
            - Path metrics: avg_shortest_path_length, diameter (if connected)
            - Centrality: avg_degree_centrality, avg_betweenness_centrality
            
        Examples
        --------
        >>> analyzer.detect_module(fdr=0.01)
        >>> stats = analyzer.network_statistics()  # Statistics for module
        >>> print(f"Density: {stats['density']:.4f}")
        >>> print(f"Avg clustering: {stats['avg_clustering_coefficient']:.4f}")
        """
        from .core.network_utils import network_statistics
        
        if graph is None:
            if self.module is not None:
                graph = self.module
            elif self.network is not None:
                graph = self.network
            else:
                raise ValueError("No network available. Load a network or detect a module first.")
        
        return network_statistics(graph)
