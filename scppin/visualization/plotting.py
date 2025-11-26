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
    node_size: float = 1000,
    figsize: Tuple[float, float] = (16, 12),
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
        Size of nodes (default: 1000)
    figsize : Tuple[float, float], optional
        Figure size (default: (16, 12))
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
    >>> from scppin import scPPIN
    >>> analyzer = scPPIN().load_network('edges.csv')
    >>> analyzer.set_node_weights(pvalues)
    >>> module = analyzer.detect_module(fdr=0.01)
    >>> analyzer.plot_module(fdr=0.01, title='Cluster 0')
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
            layout = nx.spring_layout(module, k=2.0, iterations=100, seed=42)
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
    
    # Simple node edge styling
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
        width=1.0,
        alpha=0.3,
        edge_color='lightgray',
        ax=ax
    )
    
    # Draw labels with white background for readability
    nx.draw_networkx_labels(
        module, layout,
        font_size=5,
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
