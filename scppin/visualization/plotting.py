"""Plotting functions for functional modules."""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.cm import ScalarMappable
import igraph as ig
import numpy as np
from typing import Optional, Dict, Tuple, List
import warnings


def _plot_functional_module(
    module: ig.Graph,
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
    module : ig.Graph
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
        Additional arguments (currently unused, reserved for future use)
        
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
    
    if module.vcount() == 0:
        ax.text(0.5, 0.5, 'Empty module', ha='center', va='center')
        return ax
    
    # Get node scores
    node_scores = {}
    node_names = module.vs['name']
    for i, v in enumerate(module.vs):
        try:
            score = v['score']
        except KeyError:
            score = 0.0
        node_scores[node_names[i]] = score if score is not None else 0.0
    
    # Compute layout if not provided
    if layout is None:
        try:
            layout_coords = module.layout_spring(seed=42, niter=100)
        except:
            layout_coords = module.layout_spring(seed=42)
        # Convert to dict format
        layout = {node_names[i]: layout_coords[i] for i in range(module.vcount())}
    
    # Prepare node colors based on scores
    # Higher scores (more significant) = darker colors
    scores = [node_scores.get(node_names[i], 0.0) for i in range(module.vcount())]
    
    # Normalize scores for coloring
    if len(scores) > 0 and max(scores) > min(scores):
        scores_norm = (np.array(scores) - min(scores)) / (max(scores) - min(scores))
    else:
        scores_norm = np.ones(len(scores)) * 0.5
    
    # Create colormap
    cmap_obj = plt.get_cmap(cmap)
    node_colors = [cmap_obj(score) for score in scores_norm]
    
    # Get node positions
    node_positions = np.array([layout[node_names[i]] for i in range(module.vcount())])
    
    # Draw edges first (so nodes appear on top)
    if module.ecount() > 0:
        edge_lines = []
        for e in module.es:
            u_pos = layout[node_names[e.source]]
            v_pos = layout[node_names[e.target]]
            edge_lines.append([u_pos, v_pos])
        
        edge_collection = LineCollection(
            edge_lines,
            colors='lightgray',
            linewidths=1.0,
            alpha=0.3,
            zorder=1
        )
        ax.add_collection(edge_collection)
    
    # Draw nodes
    ax.scatter(
        node_positions[:, 0],
        node_positions[:, 1],
        s=node_size,
        c=node_colors,
        edgecolors='black',
        linewidths=1.0,
        zorder=2
    )
    
    # Draw labels
    for i, name in enumerate(node_names):
        x, y = layout[name]
        ax.text(
            x, y,
            name,
            fontsize=5,
            fontweight='bold',
            ha='center',
            va='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7),
            zorder=3
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
        ax.set_title(f'Functional Module ({module.vcount()} nodes, '
                    f'{module.ecount()} edges)', fontsize=12)
    
    ax.axis('off')
    
    # Save if requested
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return ax
