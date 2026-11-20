# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
#
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Copyright (C) 2026 Pedro Sordo Martínez
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with this program. If not, see
# <https://www.gnu.org/licenses/>.

"""
statphys-viz: Interactive statistical physics visualizations for education.

Provides Monte Carlo simulations of the 2D Ising model with real-time
visualization, Jupyter widgets, and validation against the Onsager
exact solution.
"""

from .benchmark import CriticalExponents, OnsagerBenchmark
from .ising import IsingConfig, IsingModel, MetropolisSampler, WolffClusterSampler

__version__ = "0.2.1"
__author__ = "Pedro Sordo Martínez"
__email__ = "amurlaniakea@gmail.com"
__license__ = "AGPL-3.0-or-later"

__all__ = [
    "IsingConfig",
    "IsingModel",
    "MetropolisSampler",
    "WolffClusterSampler",
    "OnsagerBenchmark",
    "CriticalExponents",
]
