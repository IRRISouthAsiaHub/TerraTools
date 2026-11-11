"""GenoSee (https://github.com/hashimotoshumpei/GenoSee), MIT licence, Copyright (c) 2024 Shumpei Hashimoto.

Hashimoto (2024) GenoSee: A Novel Visualization Tool for Graphical Genotypes, Breeding Science 74(5): 454-461,
https://doi.org/10.1270/jsbbs.24041

plotting.py, utilities.py, color_set.json, chromosome_length_database.json and add_chromosome_length_to_database.py
are GenoSee's own files (commit 95090df7). Four changes in plotting.py: it imports utilities as a sibling module
("from .utilities import ...") so that it works inside this package, and the legend texts REF / ALT / HET / MISSING are read
from the dictionary LEGEND_LABELS (default: the same four words) so that they can be renamed, and the horizontal legend of the zoomed view spaces its entries by the length of those words (legend_step()), and the dashed marker lines of the compare and zoomed views can be switched off (GUIDE_LINES), and GRAPHIC (default True) restyles the finished figure: rounded, shaded chromosomes without the black lines inside, a light shadow and grid (_beautify).
"""
