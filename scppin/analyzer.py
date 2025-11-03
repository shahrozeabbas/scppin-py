"""scPPIN analyzer class for stateful module detection."""

import networkx as nx
import numpy as np
import pandas as pd
from typing import Dict, Optional, Tuple, Union, List
import warnings

# Import internal helpers
from .graph import _build_graph, _load_edges_from_file
from .pvalues import _extract_pvalues
from .weights import _normalize_weights
from .module import _detect_module
from .core import (
    fit_bum,
    compute_node_scores,
    filter_network_by_pvalues,
    simplify_network,
    validate_network
)
from .core.network_utils import load_ppin as _load_ppin
from .core.edge_weights import (
    compute_edge_weights_from_expression,
    compute_edge_weights_from_scanpy
)


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
    network : Optional[nx.Graph]
        Protein-protein interaction network
    node_weights : Optional[Dict[str, float]]
        Node weights (p-values from differential expression)
    edge_weights : Dict[Tuple[str, str], float]
        Edge weights dictionary
    adata : Optional[AnnData]
        Expression data (if provided)
    module : Optional[nx.Graph]
        Detected functional module
    bum_params : Optional[Tuple[float, float]]
        Cached BUM parameters (lambda, alpha)
    node_scores : Optional[Dict[str, float]]
        Cached node scores
    _gene_normalization : Dict[str, str]
        Gene name normalization mapping (original -> normalized)
    
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
        network: Optional[nx.Graph] = None,
        node_weights: Optional[Dict[str, float]] = None,
        edge_weights: Optional[Dict[Tuple[str, str], float]] = None,
        adata = None
    ):
        """
        Initialize scPPIN analyzer.
        
        Parameters
        ----------
        network : Optional[nx.Graph]
            Initial network (optional)
        node_weights : Optional[Dict[str, float]]
            Initial node weights (optional)
        edge_weights : Optional[Dict[Tuple[str, str], float]]
            Initial edge weights dictionary (optional)
        adata : Optional[AnnData]
            Expression data (optional)
        """
        self.network = network
        self.node_weights = None  # Will be set by set_node_weights if provided
        self.adata = adata
        self.edge_weights: Dict[Tuple[str, str], float] = {}
        self.module: Optional[nx.Graph] = None
        self.bum_params: Optional[Tuple[float, float]] = None
        self.node_scores: Optional[Dict[str, float]] = None
        self._gene_normalization: Dict[str, str] = {}
        
        # Normalize network nodes if network provided
        if self.network is not None:
            self._normalize_network_nodes()
        
        # Set node weights (which will normalize and filter network)
        if node_weights is not None:
            self.set_node_weights(node_weights)
        
        # Set edge weights (must be after network normalization and node weight filtering)
        if edge_weights is not None:
            if self.network is None:
                warnings.warn("edge_weights provided but no network. Edge weights will be ignored. "
                             "Load network first, then call set_edge_weights().")
            else:
                self.set_edge_weights(weights=edge_weights)
    
    def _normalize_node_weights(self) -> None:
        """Normalize gene names in node_weights and create mapping."""
        if self.node_weights is None:
            return
        
        normalized = {}
        for gene, weight in self.node_weights.items():
            norm_gene = _normalize_gene_name(gene)
            normalized[norm_gene] = weight
            self._gene_normalization[norm_gene] = gene  # Store original -> normalized
        
        self.node_weights = normalized
    
    def _filter_network_to_node_weights(self) -> None:
        """Filter network to genes with node weights."""
        if self.network is None or self.node_weights is None:
            return
        
        # Normalize network node names and filter
        genes_with_weights = set(self.node_weights.keys())
        nodes_to_keep = []
        
        for node in self.network.nodes():
            node_str = str(node)
            norm_node = _normalize_gene_name(node_str)
            if norm_node in genes_with_weights:
                nodes_to_keep.append(node)
        
        if nodes_to_keep:
            self.network = self.network.subgraph(nodes_to_keep).copy()
        else:
            warnings.warn("No nodes in network match node_weights after normalization")
    
    def load_network(
        self,
        source: Union[str, List[Tuple], pd.DataFrame, nx.Graph],
        weight_column: Optional[str] = None,
        format: str = 'auto'
    ) -> 'scPPIN':
        """
        Load network from file, list, DataFrame, or NetworkX graph.
        
        Parameters
        ----------
        source : Union[str, List[Tuple], pd.DataFrame, nx.Graph]
            Network source:
            - String: Path to CSV/TXT/GraphML file
            - List: List of edge tuples
            - DataFrame: Edge list DataFrame
            - nx.Graph: Existing NetworkX graph
        weight_column : Optional[str]
            Column name in CSV/DataFrame to use as edge weights.
            If provided and source is file/DataFrame, loads weights from that column.
            Sets edge weights as 'weight' attribute on network edges.
        format : str, optional
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
        # Handle NetworkX graph directly
        if isinstance(source, nx.Graph):
            self.network = source.copy()
        # Handle GraphML/GML files
        elif isinstance(source, str) and (format in ['graphml', 'gml'] or \
             source.endswith(('.graphml', '.gml'))):
            if format == 'auto':
                # Auto-detect format from extension
                if source.endswith('.graphml'):
                    fmt = 'graphml'
                elif source.endswith('.gml'):
                    fmt = 'gml'
                else:
                    fmt = 'graphml'  # default
            else:
                fmt = format
            self.network = _load_ppin(source, format=fmt)
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
            self._normalize_node_weights()
            self._filter_network_to_node_weights()
        
        # Extract edge weights if weight_column was used
        if weight_column and self.network is not None:
            self._extract_edge_weights_from_network(attr_name='weight')
        
        return self
    
    def _normalize_network_nodes(self) -> None:
        """Normalize node names in network."""
        if self.network is None:
            return
        
        # Create mapping from old to new names
        node_mapping = {}
        for node in list(self.network.nodes()):
            old_name = str(node)
            new_name = _normalize_gene_name(old_name)
            if old_name != new_name:
                node_mapping[old_name] = new_name
                self._gene_normalization[new_name] = old_name
        
        # Relabel nodes if needed
        if node_mapping:
            # Only relabel nodes that actually changed
            mapping_dict = {k: v for k, v in node_mapping.items() if k != v}
            if mapping_dict:
                self.network = nx.relabel_nodes(self.network, mapping_dict, copy=False)
    
    def _extract_edge_weights_from_network(self, attr_name: str = 'weight') -> None:
        """Extract edge weights from network attributes to self.edge_weights dict."""
        if self.network is None:
            return
        
        self.edge_weights = {}
        for u, v in self.network.edges():
            if attr_name in self.network[u][v]:
                weight = self.network[u][v][attr_name]
                # Normalize node names
                u_norm = _normalize_gene_name(str(u))
                v_norm = _normalize_gene_name(str(v))
                self.edge_weights[(u_norm, v_norm)] = float(weight)
    
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
        
        # Store node weights
        self.node_weights = weights_dict
        
        # Always normalize gene names
        self._normalize_node_weights()
        
        # Always filter network to genes with weights
        if self.network is not None:
            self._filter_network_to_node_weights()
        
        # Clear caches (new weights = new BUM fit needed)
        self.bum_params = None
        self.node_scores = None
        
        return self
    
    def set_edge_weights(
        self,
        weights: Optional[Dict[Tuple[str, str], float]] = None,
        adata = None,
        groupby: Optional[str] = None,
        group: Optional[str] = None,
        method: str = 'pearson',
        attr_name: str = 'weight'
    ) -> 'scPPIN':
        """
        Set edge weights either from expression data or user-provided dictionary.
        
        Parameters
        ----------
        weights : Optional[Dict[Tuple[str, str], float]]
            User-provided edge weights dictionary.
            If provided, sets these weights on network edges.
            Only edges in network are set (automatically filtered).
        adata : Optional[AnnData]
            Expression data for computing correlations.
            If provided, computes edge weights from expression automatically.
        groupby : Optional[str]
            Grouping key for group-specific correlations
        group : Optional[str]
            Specific group ID
        method : str, optional
            Correlation method: 'pearson' or 'spearman' (default: 'pearson')
        attr_name : str, optional
            Edge attribute name to store weights (default: 'weight')
            
        Returns
        -------
        scPPIN
            self (for method chaining)
            
        Examples
        --------
        >>> # From expression data
        >>> analyzer.set_edge_weights(adata=adata, groupby='louvain', group='0')
        >>> 
        >>> # From dictionary
        >>> weights = {('TP53', 'MDM2'): 0.9, ('TP53', 'CDKN1A'): 0.8}
        >>> analyzer.set_edge_weights(weights=weights)
        """
        if self.network is None:
            raise ValueError("Network must be loaded before setting edge weights. "
                           "Call load_network() first.")
        
        # Handle AnnData input - compute from expression
        if adata is not None:
            # Use compute_edge_weights_from_scanpy which only computes for network edges
            network_copy = self.network.copy()
            compute_edge_weights_from_scanpy(
                network_copy,
                adata,
                groupby=groupby,
                group=group,
                use_abs=True,
                min_correlation=0.0,
                output_attr=attr_name
            )
            
            # Extract computed weights
            computed_weights = {}
            for u, v in network_copy.edges():
                if attr_name in network_copy[u][v]:
                    u_norm = _normalize_gene_name(str(u))
                    v_norm = _normalize_gene_name(str(v))
                    computed_weights[(u_norm, v_norm)] = network_copy[u][v][attr_name]
            
            # Normalize to [0, 1] using min-max
            if computed_weights:
                computed_weights = _normalize_weights(computed_weights)
            
            # Set on network and store
            for (u, v), weight in computed_weights.items():
                # Find original node names in network
                u_orig = self._find_node_in_network(u)
                v_orig = self._find_node_in_network(v)
                if u_orig is not None and v_orig is not None:
                    if self.network.has_edge(u_orig, v_orig):
                        self.network[u_orig][v_orig][attr_name] = weight
                    elif self.network.has_edge(v_orig, u_orig):  # Undirected
                        self.network[v_orig][u_orig][attr_name] = weight
            
            self.edge_weights = computed_weights
        
        # Handle user-provided weights dictionary
        elif weights is not None:
            # Filter to only edges in network
            filtered_weights = {}
            for (u, v), weight in weights.items():
                u_norm = _normalize_gene_name(str(u))
                v_norm = _normalize_gene_name(str(v))
                
                # Check if edge exists in network
                u_orig = self._find_node_in_network(u_norm)
                v_orig = self._find_node_in_network(v_norm)
                
                if u_orig is not None and v_orig is not None:
                    if self.network.has_edge(u_orig, v_orig):
                        filtered_weights[(u_norm, v_norm)] = float(weight)
                        self.network[u_orig][v_orig][attr_name] = float(weight)
                    elif self.network.has_edge(v_orig, u_orig):  # Undirected
                        filtered_weights[(v_norm, u_norm)] = float(weight)
                        self.network[v_orig][u_orig][attr_name] = float(weight)
            
            self.edge_weights = filtered_weights
            
            if not filtered_weights:
                warnings.warn("No edges from weights dictionary matched network edges")
        
        else:
            raise ValueError("Either 'weights' dict or 'adata' must be provided")
        
        return self
    
    def _find_node_in_network(self, normalized_name: str) -> Optional[str]:
        """Find original node name in network given normalized name."""
        if self.network is None:
            return None
        
        # First check if normalized name matches directly
        for node in self.network.nodes():
            if _normalize_gene_name(str(node)) == normalized_name:
                return str(node)
        
        # Check reverse mapping
        if normalized_name in self._gene_normalization:
            orig_name = self._gene_normalization[normalized_name]
            if orig_name in self.network.nodes():
                return orig_name
        
        return None
    
    def detect_module(
        self,
        fdr: float = 0.01,
        edge_weight_attr: Optional[str] = None,
        edge_weight_scale: float = 1.0,
        c0: Optional[float] = None,
        missing_data_score: bool = False,
        simplify: bool = True,
        validate: bool = True
    ) -> nx.Graph:
        """
        Detect functional module using current network, node weights, and edge weights.
        
        Parameters
        ----------
        fdr : float, optional
            False discovery rate threshold (default: 0.01)
        edge_weight_attr : Optional[str], optional
            Edge attribute name for weights (default: None = uniform costs).
            If None, uses uniform edge costs matching R implementation.
        edge_weight_scale : float, optional
            Scaling factor for edge weights (default: 1.0)
        c0 : Optional[float], optional
            Minimum edge cost (default: None)
        missing_data_score : bool, optional
            Include genes without weights (default: False)
        simplify : bool, optional
            Simplify network (default: True)
        validate : bool, optional
            Validate network (default: True)
            
        Returns
        -------
        nx.Graph
            Detected functional module
            
        Examples
        --------
        >>> module = analyzer.detect_module(fdr=0.01)
        >>> module = analyzer.detect_module(fdr=0.01, edge_weight_attr='weight')
        """
        if self.network is None:
            raise ValueError("Network must be loaded. Call load_network() first.")
        
        if self.node_weights is None:
            raise ValueError("Node weights must be set. Call set_node_weights() first.")
        
        # Use internal detect_module function
        self.module = _detect_module(
            self.network,
            self.node_weights,
            fdr=fdr,
            edge_weight_attr=edge_weight_attr,
            edge_weight_scale=edge_weight_scale,
            c0=c0,
            missing_data_score=missing_data_score,
            simplify=simplify,
            validate=validate
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

