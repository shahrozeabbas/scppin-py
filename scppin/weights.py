"""Edge weight normalization utilities."""

from typing import Dict, Tuple


def _normalize_weights(weights: Dict[Tuple[str, str], float]) -> Dict[Tuple[str, str], float]:
    """Normalize weights to [0, 1] using min-max scaling."""
    if not weights:
        return weights
    
    values = list(weights.values())
    min_val = min(values)
    max_val = max(values)
    
    if max_val > min_val:
        # Min-max normalization
        normalized = {
            edge: (weight - min_val) / (max_val - min_val)
            for edge, weight in weights.items()
        }
    else:
        # All weights are the same - set to 0.5
        normalized = {edge: 0.5 for edge in weights.keys()}
    
    return normalized

