"""
Define some common types used throughout the project
"""

from typing import List, Tuple, Callable
from src.utils.Protein import Protein

ProteinDataset = List[Tuple[Tuple[Protein, Protein], Tuple[Protein, Protein]]]
# BindingSiteMapping = Tuple[]
Solution = Tuple[Protein, Protein, float]
Solutions = List[Solution]