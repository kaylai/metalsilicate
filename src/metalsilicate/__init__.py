"""
metalsilicate

Python library for calculating thermodynamic parameters of geochemical systems with both silicate
melt and Fe-rich metal in equilibrium.

Originally developed by Iacovino et al. to process a partitioning database and data from completed
gas levitation laser melting experiments to develop a thermodynamic model for Si partitioning
between silicate melt and Fe-rich metal in support of assessing the cataclysmic reduction hypothesis
for Mercury’s origin
"""

__version__ = "0.1.0"
__author__ = "Kayla Iacovino"

# ----------------- IMPORTS ----------------- #
import warnings as w
import pandas as pd

import metalsilicate.core