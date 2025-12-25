"""P-value extraction from scanpy results."""

from typing import Dict


def _extract_pvalues(
    adata,
    groupby: str,
    group: str
) -> Dict[str, float]:
    """
    Extract p-values from scanpy rank_genes_groups results.
    
    Parameters
    ----------
    adata : AnnData
        Scanpy AnnData object.
        Must have run sc.tl.rank_genes_groups() first.
    groupby : str
        Key in adata.obs for grouping labels
    group : str
        Specific group to extract p-values for
        
    Returns
    -------
    Dict[str, float]
        Dictionary mapping gene names to p-values
        
    Raises
    ------
    ValueError
        If rank_genes_groups results not found or cluster not found
        
    Examples
    --------
    >>> import scanpy as sc
    >>> import scppin
    >>> 
    >>> # Run differential expression
    >>> sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
    >>> 
    >>> # Extract p-values for cluster 0
    >>> pvalues = scppin.extract_pvalues(adata, 'louvain', '0')
    >>> analyzer = scppin.scPPIN()
    >>> analyzer.load_network('edges.csv')
    >>> analyzer.set_node_weights(pvalues)
    >>> module = analyzer.detect_module(fdr=0.01)
    """
    # Check if rank_genes_groups has been run
    if 'rank_genes_groups' not in adata.uns:
        raise ValueError(
            "rank_genes_groups results not found in adata.uns. "
            "Please run sc.tl.rank_genes_groups() first."
        )
    
    # Get results
    rank_genes = adata.uns['rank_genes_groups']
    
    # Validate results
    if 'names' not in rank_genes or 'pvals' not in rank_genes:
        raise ValueError("rank_genes_groups results are incomplete")
    
    # Get group names
    group_names = rank_genes['names'].dtype.names
    if group not in group_names:
        raise ValueError(
            f"Group '{group}' not found. "
            f"Available groups: {group_names}"
        )
    
    # Extract gene names and p-values for this group
    genes = rank_genes['names'][group]
    pvals = rank_genes['pvals'][group]
    
    # Create p-value dictionary (only valid p-values)
    pvalues = {gene: float(pval) for gene, pval in zip(genes, pvals) if 0 < pval <= 1}
    
    if not pvalues:
        raise ValueError(f"No valid p-values found for group '{group}'")
    
    return pvalues
