# scPPIN-py: Single-cell Protein-Protein Interaction Network Analysis

**Pure Python implementation of scPPIN with class-based API**

This Python package detects functional modules in protein-protein interaction networks by integrating single-cell RNA sequencing data. Reimplementation of the original R package with a clean, object-oriented API.

## Original Method

The method is described in:

> Functional module detection through integration of single-cell RNA sequencing data with protein–protein interaction networks.  
> Florian Klimm, Enrique M. Toledo, Thomas Monfeuga, Fang Zhang, Charlotte M. Deane, and Gesine Reinert  
> BMC Genomics 21, Article number: 756 (2020) https://doi.org/10.1186/s12864-020-07144-2

## Key Features

- **Class-Based API**: Object-oriented design with method chaining
- **Pure Python**: No external C++ binaries (uses `pcst_fast`)
- **Edge Weights**: Supports confidence scores or computed correlations
- **Author's Formula**: Implements recommended edge cost formula from [GitHub Issue #10](https://github.com/floklimm/scPPIN/issues/10)
- **Scanpy Integration**: Works seamlessly with AnnData objects
- **Fast**: Vectorized NumPy operations, ~5-10x faster than R
- **Easy Installation**: Single `pip install` command
- **Automatic Filtering**: Network automatically filtered to genes with p-values

## Installation

### Quick Install

```bash
pip install scppin
```

### Development Install

1. **Clone the repository:**
```bash
git clone https://github.com/floklimm/scPPIN.git
cd scPPIN/scppin-py
```

2. **Create a virtual environment (recommended):**
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install in development mode:**
```bash
pip install -e ".[all]"
```

### Dependencies

#### Core Dependencies (required)
- numpy >= 1.20.0
- scipy >= 1.7.0
- networkx >= 2.6.0
- matplotlib >= 3.3.0
- pandas >= 1.3.0
- pcst-fast >= 1.0.0

#### Optional Dependencies
- scanpy >= 1.8.0 (for scanpy integration)
- anndata >= 0.8.0 (for AnnData support)

#### Development Dependencies
- pytest >= 7.0.0
- pytest-cov >= 3.0.0
- black >= 22.0.0
- flake8 >= 4.0.0

### Installing pcst-fast

The `pcst-fast` package is required for solving the Prize-Collecting Steiner Tree problem. It should install automatically, but if you encounter issues:

```bash
pip install pcst-fast
```

If compilation fails, you may need to install build tools:

**On Ubuntu/Debian:**
```bash
sudo apt-get install build-essential python3-dev
```

**On macOS:**
```bash
xcode-select --install
```

**On Windows:**
Install Visual Studio Build Tools from https://visualstudio.microsoft.com/downloads/

### Verifying Installation

```python
from scppin import scPPIN

# Run a simple test
analyzer = scPPIN()
analyzer.load_network([('A', 'B'), ('B', 'C')])
analyzer.set_node_weights({'A': 0.001, 'B': 0.005, 'C': 0.01})
analyzer.detect_module(fdr=0.01)

print(f"Module detected: {analyzer.module.number_of_nodes()} nodes")
```

### Platform-Specific Notes

**macOS (Apple Silicon M1/M2)**
- May need to install Rosetta for some dependencies
- Use native ARM Python if possible

**Windows**
- May need Visual Studio Build Tools
- Use Anaconda for easier dependency management

**Linux**
- Should work out of the box on most distributions
- May need `python3-dev` package

## Quick Start

### Basic Workflow

```python
from scppin import scPPIN

# Create analyzer instance
analyzer = scPPIN()

# 1. Load network from edge list
analyzer.load_network('edges.csv')

# 2. Set node weights (p-values)
pvalues = {'TP53': 0.0001, 'MDM2': 0.001, 'CDKN1A': 0.005}
analyzer.set_node_weights(pvalues)

# 3. Detect functional module
analyzer.detect_module(fdr=0.01)

# 4. Visualize
analyzer.plot_module(fdr=0.01)
```

### Method Chaining

```python
from scppin import scPPIN

analyzer = (scPPIN()
           .load_network('edges.csv')
           .set_node_weights({'TP53': 0.0001, 'MDM2': 0.001})
           .detect_module(fdr=0.01))

# Access results
print(f"Module has {analyzer.module.number_of_nodes()} nodes")
```

### With Scanpy

```python
import scanpy as sc
from scppin import scPPIN

# Process single-cell data
adata = sc.read_h5ad('data.h5ad')
sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')

# Create analyzer and set data
analyzer = scPPIN()
analyzer.load_network('ppi_edges.csv')
analyzer.set_node_weights(adata, groupby='louvain', group='0')

# Detect module
analyzer.detect_module(fdr=0.01)
```

### With Edge Weights

#### Option 1: Pre-existing Weights from File

```python
from scppin import scPPIN

analyzer = scPPIN()

# Load network with weights from CSV column
analyzer.load_network('edges.csv', weight_column='confidence')
analyzer.set_node_weights(pvalues)
analyzer.detect_module(fdr=0.01, edge_weight_scale=0.5)
```

#### Option 2: Pre-existing Weights from Dictionary

```python
from scppin import scPPIN

analyzer = scPPIN()
analyzer.load_network('edges.csv')

weights = {('TP53', 'MDM2'): 0.95, ('TP53', 'CDKN1A'): 0.90}
analyzer.set_edge_weights(weights=weights)
analyzer.set_node_weights(pvalues)
analyzer.detect_module(fdr=0.01, edge_weight_scale=0.5)
```

#### Option 3: Compute from Expression

```python
from scppin import scPPIN

analyzer = scPPIN()
analyzer.load_network('edges.csv')
analyzer.set_node_weights(adata, groupby='louvain', group='0')

# Compute edge weights from expression data
analyzer.set_edge_weights(
    adata=adata,
    groupby='louvain',
    group='0'
)

analyzer.detect_module(fdr=0.01, edge_weight_scale=0.5)
```

### API Summary

| Method | Purpose |
|--------|---------|
| `load_network(source, weight_column)` | Load network from file/list/DataFrame |
| `set_node_weights(weights, ...)` | Set p-values from dict or AnnData |
| `set_edge_weights(weights/adata, ...)` | Set edge weights from dict or compute from expression |
| `detect_module(fdr, ...)` | Detect functional module |
| `plot_module(fdr, ...)` | Visualize module |

## API Reference

### `scPPIN` Class

Main class for scPPIN analysis. All analysis is performed through an instance of this class.

#### `__init__(network=None, adata=None)`

Create a new scPPIN analyzer instance.

**Parameters:**
- `network`: Optional NetworkX graph (can be loaded later)
- `adata`: Optional AnnData object (can be set later)

**Attributes:**
- `network`: NetworkX graph (PPIN)
- `node_weights`: Dict mapping gene names to p-values
- `edge_weights`: Dict mapping (gene1, gene2) tuples to weights
- `module`: Detected functional module (NetworkX subgraph)
- `bum_params`: BUM model parameters (lambda, alpha)
- `fdr`: False discovery rate threshold

#### `load_network(source, weight_column=None, format='auto')`

Load network from various sources.

**Parameters:**
- `source`: Network source:
  - Path to CSV/TXT file with edge list
  - List of (gene1, gene2) tuples
  - pandas DataFrame with 'source'/'target' columns
  - NetworkX graph
  - Path to GraphML/GML file
- `weight_column`: Optional column name in CSV/DataFrame for edge weights
- `format`: File format hint ('auto', 'csv', 'graphml', 'gml')

**Returns:** `self` (for method chaining)

**Examples:**
```python
# From file
analyzer.load_network('edges.csv')

# From file with weights
analyzer.load_network('edges.csv', weight_column='confidence')

# From list
edges = [('TP53', 'MDM2'), ('TP53', 'CDKN1A')]
analyzer.load_network(edges)

# From DataFrame
import pandas as pd
df = pd.DataFrame({'source': ['TP53'], 'target': ['MDM2']})
analyzer.load_network(df)
```

#### `set_node_weights(weights, groupby=None, group=None)`

Set node weights (p-values) and filter network to genes with weights.

**Parameters:**
- `weights`: Node weights:
  - Dict mapping gene names to p-values
  - AnnData object (requires groupby/group)
- `groupby`: Key in adata.obs for grouping labels (if weights is AnnData)
- `group`: Specific group to extract (if weights is AnnData)

**Returns:** `self` (for method chaining)

**Note:** This method automatically normalizes gene names and filters the network to only include genes with weights.

**Examples:**
```python
# From dictionary
analyzer.set_node_weights({'TP53': 0.0001, 'MDM2': 0.001})

# From AnnData
analyzer.set_node_weights(adata, groupby='louvain', group='0')
```

#### `set_edge_weights(weights=None, adata=None, groupby=None, group=None, method='pearson', attr_name='weight')`

Set edge weights either from expression data or user-provided dictionary.

**Parameters:**
- `weights`: Optional dict mapping (gene1, gene2) to weights
- `adata`: Optional AnnData object for computing correlations
- `groupby`: Grouping key for group-specific correlations
- `group`: Specific group ID
- `method`: Correlation method ('pearson' or 'spearman')
- `attr_name`: Edge attribute name to store weights (default: 'weight')

**Returns:** `self` (for method chaining)

**Examples:**
```python
# From dictionary
weights = {('TP53', 'MDM2'): 0.95}
analyzer.set_edge_weights(weights=weights)

# From expression data
analyzer.set_edge_weights(
    adata=adata,
    groupby='louvain',
    group='0'
)
```

#### `detect_module(fdr=0.01, edge_weight_scale=1.0, c0=None, missing_data_score=False, simplify=True, validate=True)`

Detect functional module using loaded network, node weights, and edge weights.

**Parameters:**
- `fdr`: False discovery rate threshold (default: 0.01)
- `edge_weight_scale`: Scaling factor for edge weights (default: 1.0)
- `c0`: Minimum edge cost (default: 0.1 * base_cost)
- `missing_data_score`: Include genes without p-values (default: False)
- `simplify`: Remove self-loops and parallel edges (default: True)
- `validate`: Validate inputs (default: True)

**Returns:** `self` (for method chaining)

**Sets:** `self.module` to detected functional module (NetworkX subgraph)

**Examples:**
```python
# Basic detection
analyzer.detect_module(fdr=0.01)

# With edge weights
analyzer.detect_module(fdr=0.01, edge_weight_scale=0.5)
```

#### `plot_module(fdr=None, **kwargs)`

Plot the detected functional module.

**Parameters:**
- `fdr`: FDR threshold for visualization (default: uses analyzer.fdr)
- `**kwargs`: Additional arguments passed to plotting function

**Returns:** matplotlib Figure object

**Examples:**
```python
# Basic plot
analyzer.plot_module()

# Customized plot
analyzer.plot_module(fdr=0.01, figsize=(10, 8), node_size=500)
```

## Complete Workflow Examples

### Example 1: Basic (Unweighted)

```python
from scppin import scPPIN

analyzer = scPPIN()
analyzer.load_network('ppi_edges.csv')
analyzer.set_node_weights({'TP53': 0.0001, 'MDM2': 0.001, 'CDKN1A': 0.005})
analyzer.detect_module(fdr=0.01)
analyzer.plot_module()
```

### Example 2: With Pre-computed Weights

```python
from scppin import scPPIN

analyzer = scPPIN()
analyzer.load_network('ppi_edges.csv', weight_column='confidence')
analyzer.set_node_weights(pvalues)
analyzer.detect_module(fdr=0.01, edge_weight_scale=0.5)
```

### Example 3: Complete Scanpy Workflow

```python
import scanpy as sc
from scppin import scPPIN

# Load and process data
adata = sc.read_h5ad('data.h5ad')
sc.tl.rank_genes_groups(adata, 'louvain')

# For each group
for group_id in ['0', '1', '2']:
    analyzer = scPPIN()
    analyzer.load_network('ppi_edges.csv')
    analyzer.set_node_weights(adata, groupby='louvain', group=group_id)
    analyzer.set_edge_weights(
        adata=adata,
        groupby='louvain',
        group=group_id
    )
    analyzer.detect_module(fdr=0.01)
    
    print(f"Group {group_id}: {analyzer.module.number_of_nodes()} nodes")
```

## Edge Weight Formula

This implementation uses the edge cost formula recommended by the author (Florian Klimm) in response to [GitHub Issue #10](https://github.com/floklimm/scPPIN/issues/10):

```
cost = base_cost * (1 - weight * scale) + c0
```

Where:
- `base_cost`: Typically `-minimumScore` (positive value)
- `weight`: Edge weight in [0, 1]
- `scale`: Scaling factor for weight influence
- `c0`: Minimum cost to prevent zeros (default: 0.1 * base_cost)

**Key properties:**
- Higher edge weights → Lower costs → More likely to include edge
- `c0` prevents zero-cost edges (important for numerical stability)
- FDR guarantees are approximate when using edge weights

## Algorithm Overview

The pipeline:

1. **Fit BUM Model**: Distinguish true signals from noise in p-values
2. **Compute Node Scores**: Transform p-values to node weights
3. **Prepare Edge Costs**: Apply author's formula with edge weights
4. **Solve PCST**: Find optimal subgraph maximizing prize - cost
5. **Return Module**: Connected subgraph of significant genes

## Running Tests

```bash
pytest tests/ -v
```

With coverage:
```bash
pytest tests/ -v --cov=scppin --cov-report=html
```

## Running Examples

```bash
cd examples
python basic_tutorial.py
python new_api_example.py
python edge_weights_example.py
```

## Troubleshooting

**Issue: pcst-fast installation fails**
- Solution: Install build tools (see Installation section above)
- Alternative: Try installing from conda-forge

**Issue: "No module named 'scppin'"**
- Solution: Make sure you're in the correct virtual environment
- Solution: Try `pip install -e .` from the scppin-py directory

**Issue: Network file not found**
- Solution: Copy network files from R package:
  ```bash
  cp ../R/inst/extdata/*.graphml scppin/data/networks/
  ```

**Issue: Import errors with scanpy**
- Solution: Install optional dependencies:
  ```bash
  pip install "scppin[scanpy]"
  ```

**Issue: pcst-fast not found**
```bash
pip install pcst-fast
```

**Issue: Import errors**
```bash
pip install -e .
```

## Requirements

### Core (Required)
- numpy >= 1.20.0
- scipy >= 1.7.0
- networkx >= 2.6.0
- matplotlib >= 3.3.0
- pandas >= 1.3.0
- pcst-fast >= 1.0.0

### Optional
- scanpy >= 1.8.0 (for scanpy integration)
- anndata >= 0.8.0 (for AnnData support)

## Citation

If you use this package, please cite the original paper:

```bibtex
@article{klimm2020functional,
  title={Functional module detection through integration of single-cell RNA sequencing data with protein--protein interaction networks},
  author={Klimm, Florian and Toledo, Enrique M and Monfeuga, Thomas and Zhang, Fang and Deane, Charlotte M and Reinert, Gesine},
  journal={BMC genomics},
  volume={21},
  number={1},
  pages={1--10},
  year={2020},
  publisher={BioMed Central}
}
```

## License

GPL-3.0 (same as original R package)

## Acknowledgments

- Original R implementation by Florian Klimm et al.
- Edge weight formula from author's response to [GitHub Issue #10](https://github.com/floklimm/scPPIN/issues/10)
- Python implementation uses `pcst_fast` library

## Support

For issues and questions:
- Original R package: https://github.com/floklimm/scPPIN
- Python implementation: Open an issue on GitHub
