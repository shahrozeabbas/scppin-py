# scPPIN-py

**Pure Python implementation of scPPIN for single-cell protein-protein interaction network analysis**

[![Documentation](https://readthedocs.org/projects/scppin-py/badge/?version=latest)](https://scppin-py.readthedocs.io/en/latest/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Paper](https://img.shields.io/badge/BMC%20Genomics-2020-blue)](https://doi.org/10.1186/s12864-020-07144-2)

scPPIN-py detects functional modules in protein-protein interaction networks by integrating single-cell RNA sequencing data. This is a reimplementation of the [original R package](https://github.com/floklimm/scPPIN) with a clean, object-oriented Python API.

**Original method**: [Klimm et al. (2020)](https://doi.org/10.1186/s12864-020-07144-2), *BMC Genomics*

## Key Features

- 🐍 **Pure Python** — No external C++ binaries (uses `pcst_fast`)
- ⚡ **Fast** — Vectorized NumPy operations, ~5-10x faster than R
- 🔗 **Class-Based API** — Object-oriented design with setup method chaining
- 🎯 **Edge Weights** — Supports pre-computed edge weights from dictionaries
- 📊 **Scanpy Integration** — Works seamlessly with AnnData objects for p-value extraction
- 📦 **Easy Installation** — Single `pip install` command
- 🎯 **Standalone PCST** — Direct PCST implementation without dependencies on expression data

## Quick Start

### Installation

```bash
pip install scppin
```

### Basic Usage

```python
from scppin import scPPIN

# Create analyzer and load network
analyzer = scPPIN()
analyzer.load_network('edges.csv')

# Set node weights (p-values from differential expression)
pvalues = {'TP53': 0.0001, 'MDM2': 0.001, 'CDKN1A': 0.005}
analyzer.set_node_weights(pvalues)

# Optionally set edge weights (from pre-computed dictionary)
edge_weights = {('TP53', 'MDM2'): 0.9, ('TP53', 'CDKN1A'): 0.8}
analyzer.set_edge_weights(weights=edge_weights)

# Detect functional module using PCST
module = analyzer.detect_module(fdr=0.01, edge_weight_attr='weight')

# Visualize
analyzer.plot_module()
```

## Documentation

📚 **Full documentation**: https://scppin-py.readthedocs.io

- [Installation Guide](https://scppin-py.readthedocs.io/en/latest/installation.html)
- [Quick Start Tutorial](https://scppin-py.readthedocs.io/en/latest/quickstart.html)
- [API Reference](https://scppin-py.readthedocs.io/en/latest/api/index.html)
- [Examples](https://scppin-py.readthedocs.io/en/latest/examples/index.html)
- [Algorithm Overview](https://scppin-py.readthedocs.io/en/latest/algorithm/overview.html)

## Citation

If you use scPPIN-py in your research, please cite the original paper:

```bibtex
@article{klimm2020functional,
  title={Functional module detection through integration of single-cell RNA sequencing data with protein--protein interaction networks},
  author={Klimm, Florian and Toledo, Enrique M and Monfeuga, Thomas and Zhang, Fang and Deane, Charlotte M and Reinert, Gesine},
  journal={BMC Genomics},
  volume={21},
  number={1},
  pages={756},
  year={2020},
  publisher={BioMed Central},
  doi={10.1186/s12864-020-07144-2}
}
```

## License

GPL-3.0 (same as original R package)

## Links

- [Documentation](https://scppin-py.readthedocs.io)
- [Original R Package](https://github.com/floklimm/scPPIN)
- [Original Paper](https://doi.org/10.1186/s12864-020-07144-2)
