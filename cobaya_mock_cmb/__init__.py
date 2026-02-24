"""
Mock CMB likelihoods for Cobaya
"""

# Base class
from .mock_cmb_likelihood import MockCMBLikelihood

# Particular examples
from .mock_Planck import MockPlanck
from .mock_Planck_highl import MockPlanckHighL
from .mock_Planck_lowl import MockPlanckLowL
from .mock_SO import MockSO
from .mock_SO_baseline import MockSOBaseline
from .mock_SO_goal import MockSOGoal
from .mock_CMBS4 import MockCMBS4
from .mock_CMBS4sens0 import MockCMBS4sens0
from .mock_SO_clumping import MockSOClumping
from .mock_CMBS4_clumping import MockCMBS4Clumping
from .mock_SO_ax28 import MockSOAx28
from .mock_SO_ax27 import MockSOAx27
from .mock_SO_ax26 import MockSOAx26
from .mock_SO_ax25 import MockSOAx25
from .mock_SO_ax24 import MockSOAx24
from .mock_SO_ax23 import MockSOAx23

# Metadata
__author__ = "Michael Rashkovetskyi, Julian B. Muñoz, Daniel J. Eisenstein and Cora Dvorkin"
__version__ = "0.1.0"
__obsolete__ = False
__year__ = "2021"
__url__ = "https://github.com/misharash/cobaya_mock_cmb"
