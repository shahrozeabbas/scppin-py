"""Plotting functions for functional modules."""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import ScalarMappable
import networkx as nx
import numpy as np
from typing import Optional, Dict, Tuple, List
import warnings


def _plot_functional_module(
    module: nx.Graph,
    fdr: float = 0.0,
    layout: Optional[Dict] = None,
    node_size: float = 300,
    figsize: Tuple[float, float] = (12, 8),
    cmap: str = 'Purples',
    show_colorbar: bool = True,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    save_path: Optional[str] = None,
    **kwargs
) -> plt.Axes:
    """
    Plot functional module with nodes colored by p-values.
    
    Parameters
    ----------
    module : nx.Graph
        Functional module to plot (should have 'score' node attribute)
    fdr : float, optional
        False discovery rate threshold for highlighting (default: 0.0)
    layout : Dict, optional
        Node positions. If None, uses spring layout.
    node_size : float, optional
        Size of nodes (default: 300)
    figsize : Tuple[float, float], optional
        Figure size (default: (12, 8))
    cmap : str, optional
        Colormap name (default: 'Purples')
    show_colorbar : bool, optional
        Whether to show colorbar (default: True)
    title : str, optional
        Plot title
    ax : plt.Axes, optional
        Axes to plot on. If None, creates new figure.
    save_path : str, optional
        Path to save figure
    **kwargs
        Additional arguments passed to nx.draw_networkx
        
    Returns
    -------
    plt.Axes
        Axes object with the plot
        
    Examples
    --------
    >>> import scppin
    >>> module = scppin.detect_functional_module(network, pvalues, fdr=0.01)
    >>> scppin.plot_functional_module(module, fdr=0.01, title='Cluster 0')
    """
    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    
    if module.number_of_nodes() == 0:
        ax.text(0.5, 0.5, 'Empty module', ha='center', va='center')
        return ax
    
    # Get node scores and convert to p-values for coloring
    node_scores = {}
    for node in module.nodes():
        if 'score' in module.nodes[node]:
            node_scores[node] = module.nodes[node]['score']
        else:
            node_scores[node] = 0.0
    
    # Compute layout if not provided
    if layout is None:
        try:
            layout = nx.spring_layout(module, k=1/np.sqrt(module.number_of_nodes()), 
                                     iterations=50, seed=42)
        except:
            layout = nx.spring_layout(module, seed=42)
    
    # Prepare node colors based on scores
    # Higher scores (more significant) = darker colors
    node_list = list(module.nodes())
    scores = [node_scores.get(node, 0.0) for node in node_list]
    
    # Normalize scores for coloring
    if len(scores) > 0 and max(scores) > min(scores):
        scores_norm = (np.array(scores) - min(scores)) / (max(scores) - min(scores))
    else:
        scores_norm = np.ones(len(scores)) * 0.5
    
    # Create colormap
    cmap_obj = plt.get_cmap(cmap)
    node_colors = [cmap_obj(score) for score in scores_norm]
    
    # Determine node shapes based on FDR
    # Nodes above FDR threshold get different marker
    if fdr > 0:
        # This is simplified - NetworkX doesn't easily support different node shapes
        # We'll use edge colors to indicate FDR threshold
        node_edge_colors = []
        node_edge_widths = []
        for node in node_list:
            score = node_scores.get(node, 0.0)
            # Higher score = more significant = below FDR
            # This is a simplification; ideally we'd recompute p-value
            if score > np.median(scores):  # Simplified FDR check
                node_edge_colors.append('black')
                node_edge_widths.append(1.0)
            else:
                node_edge_colors.append('red')
                node_edge_widths.append(2.0)
    else:
        node_edge_colors = 'black'
        node_edge_widths = 1.0
    
    # Draw network
    nx.draw_networkx_nodes(
        module, layout,
        node_list,
        node_color=node_colors,
        node_size=node_size,
        edgecolors=node_edge_colors,
        linewidths=node_edge_widths,
        ax=ax
    )
    
    # Draw edges
    nx.draw_networkx_edges(
        module, layout,
        width=1.5,
        alpha=0.5,
        edge_color='gray',
        ax=ax
    )
    
    # Draw labels
    nx.draw_networkx_labels(
        module, layout,
        font_size=8,
        font_weight='bold',
        ax=ax
    )
    
    # Add colorbar
    if show_colorbar:
        sm = ScalarMappable(cmap=cmap_obj, 
                           norm=mcolors.Normalize(vmin=min(scores), vmax=max(scores)))
        sm.set_array([])
        
        # Create colorbar
        cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Node Score', rotation=270, labelpad=15)
    
    # Set title
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title(f'Functional Module ({module.number_of_nodes()} nodes, '
                    f'{module.number_of_edges()} edges)', fontsize=12)
    
    ax.axis('off')
    
    # Save if requested
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return ax


def plot_multiple_modules(
    modules: Dict[str, nx.Graph],
    fdr: float = 0.0,
    ncols: int = 3,
    figsize: Optional[Tuple[float, float]] = None,
    **kwargs
) -> plt.Figure:
    """
    Plot multiple functional modules in a grid.
    
    Parameters
    ----------
    modules : Dict[str, nx.Graph]
        Dictionary mapping cluster IDs to modules
    fdr : float, optional
        False discovery rate threshold (default: 0.0)
    ncols : int, optional
        Number of columns in grid (default: 3)
    figsize : Tuple[float, float], optional
        Figure size. If None, computed automatically.
    **kwargs
        Additional arguments passed to plot_functional_module
        
    Returns
    -------
    plt.Figure
        Figure object with all plots
        
    Examples
    --------
    >>> modules = scppin.detect_modules_for_all_clusters(adata, network)
    >>> scppin.plot_multiple_modules(modules, fdr=0.01)
    """
    n_modules = len(modules)
    nrows = int(np.ceil(n_modules / ncols))
    
    if figsize is None:
        figsize = (5 * ncols, 4 * nrows)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    
    if n_modules == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for idx, (cluster_id, module) in enumerate(modules.items()):
        ax = axes[idx]
        plot_functional_module(
            module,
            fdr=fdr,
            title=f'Cluster {cluster_id}',
            ax=ax,
            show_colorbar=False,
            **kwargs
        )
    
    # Hide empty subplots
    for idx in range(n_modules, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    
    return fig


def plot_module_statistics(
    modules: Dict[str, nx.Graph],
    figsize: Tuple[float, float] = (12, 4)
) -> plt.Figure:
    """
    Plot statistics for multiple modules.
    
    Parameters
    ----------
    modules : Dict[str, nx.Graph]
        Dictionary mapping cluster IDs to modules
    figsize : Tuple[float, float], optional
        Figure size (default: (12, 4))
        
    Returns
    -------
    plt.Figure
        Figure with statistics plots
    """
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    cluster_ids = list(modules.keys())
    num_nodes = [modules[cid].number_of_nodes() for cid in cluster_ids]
    num_edges = [modules[cid].number_of_edges() for cid in cluster_ids]
    densities = [nx.density(modules[cid]) if modules[cid].number_of_nodes() > 0 else 0 
                for cid in cluster_ids]
    
    # Number of nodes
    axes[0].bar(range(len(cluster_ids)), num_nodes, color='steelblue')
    axes[0].set_xlabel('Cluster')
    axes[0].set_ylabel('Number of Nodes')
    axes[0].set_title('Module Size')
    axes[0].set_xticks(range(len(cluster_ids)))
    axes[0].set_xticklabels(cluster_ids, rotation=45)
    
    # Number of edges
    axes[1].bar(range(len(cluster_ids)), num_edges, color='coral')
    axes[1].set_xlabel('Cluster')
    axes[1].set_ylabel('Number of Edges')
    axes[1].set_title('Module Connectivity')
    axes[1].set_xticks(range(len(cluster_ids)))
    axes[1].set_xticklabels(cluster_ids, rotation=45)
    
    # Density
    axes[2].bar(range(len(cluster_ids)), densities, color='mediumseagreen')
    axes[2].set_xlabel('Cluster')
    axes[2].set_ylabel('Density')
    axes[2].set_title('Module Density')
    axes[2].set_xticks(range(len(cluster_ids)))
    axes[2].set_xticklabels(cluster_ids, rotation=45)
    
    plt.tight_layout()
    
    return fig

