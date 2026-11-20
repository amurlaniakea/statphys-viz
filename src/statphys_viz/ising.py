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
Ising model implementation with Monte Carlo samplers.

Implements the 2D Ising model on a square lattice with periodic boundary
conditions. Provides both Metropolis single-spin-flip and Wolff cluster
algorithms for sampling.

Reference:
- Metropolis et al. (1953) J. Chem. Phys. 21, 1087
- Wolff (1989) Phys. Rev. Lett. 62, 361
- Onsager (1944) Phys. Rev. 65, 117 (exact solution)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numba import njit


@njit(cache=True, fastmath=True)
def _periodic_idx(i: int, L: int) -> int:
    """Periodic boundary condition index."""
    if i < 0:
        return i + L
    if i >= L:
        return i - L
    return i


@njit(cache=True, fastmath=True)
def _delta_energy(spins: np.ndarray, i: int, j: int, J: float, h: float) -> float:
    """Energy change for flipping spin at (i, j)."""
    L = spins.shape[0]
    s = spins[i, j]
    nb_sum = (
        spins[_periodic_idx(i + 1, L), j]
        + spins[_periodic_idx(i - 1, L), j]
        + spins[i, _periodic_idx(j + 1, L)]
        + spins[i, _periodic_idx(j - 1, L)]
    )
    return 2.0 * J * s * nb_sum + 2.0 * h * s


@njit(cache=True, fastmath=True)
def _seed_numba(seed: int) -> None:
    """Siembra el RNG global de numba (np.random dentro de funciones njit)."""
    np.random.seed(seed)


@njit(cache=True, fastmath=True)
def _metropolis_sweep(
    spins: np.ndarray,
    J: float,
    h: float,
    beta: float,
    rng_state: np.ndarray,
) -> int:
    """Single Metropolis sweep (L*L spin flip attempts). Returns accepted flips."""
    L = spins.shape[0]
    accepted = 0

    # Barrido secuencial: prange sobre vecinos acoplados es una carrera de datos
    # y rompe el balance detallado.
    for idx in range(L * L):
        i = idx // L
        j = idx % L

        dE = _delta_energy(spins, i, j, J, h)

        if dE <= 0.0 or np.random.random() < np.exp(-beta * dE):
            spins[i, j] *= -1
            accepted += 1

    return accepted


@njit(cache=True, fastmath=True)
def _wolff_cluster(
    spins: np.ndarray,
    i: int,
    j: int,
    J: float,
    beta: float,
    cluster: np.ndarray,
    p_add: float,
    rng_state: np.ndarray,
) -> int:
    """Grow Wolff cluster from seed (i, j). Returns cluster size."""
    L = spins.shape[0]
    spin_val = spins[i, j]
    cluster[i, j] = 1
    stack = [(i, j)]
    size = 1

    while stack:
        ci, cj = stack.pop()
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ni = _periodic_idx(ci + di, L)
            nj = _periodic_idx(cj + dj, L)
            if cluster[ni, nj] == 0 and spins[ni, nj] == spin_val and np.random.random() < p_add:
                cluster[ni, nj] = 1
                stack.append((ni, nj))
                size += 1
    return size


@njit(cache=True, fastmath=True)
def _wolff_sweep(
    spins: np.ndarray,
    J: float,
    beta: float,
    rng_state: np.ndarray,
) -> int:
    """Single Wolff cluster flip. Returns cluster size."""
    L = spins.shape[0]
    p_add = 1.0 - np.exp(-2.0 * beta * J)

    # Random seed site
    i = np.random.randint(0, L)
    j = np.random.randint(0, L)

    cluster = np.zeros((L, L), dtype=np.int8)
    size = _wolff_cluster(spins, i, j, J, beta, cluster, p_add, rng_state)

    # Flip entire cluster
    for ci in range(L):
        for cj in range(L):
            if cluster[ci, cj]:
                spins[ci, cj] *= -1

    return size


@dataclass(frozen=True)
class IsingConfig:
    """Configuration for Ising model simulation."""

    L: int = 32
    J: float = 1.0
    h: float = 0.0
    T: float = 2.269  # Near Tc for 2D Ising
    algorithm: Literal["metropolis", "wolff"] = "wolff"
    seed: int | None = None


class IsingModel:
    """
    2D Ising model on square lattice with periodic boundaries.

    Attributes:
        L: Linear lattice size (L x L spins)
        J: Coupling constant (ferromagnetic if > 0)
        h: External magnetic field
        T: Temperature
        beta: Inverse temperature (1/T)
    """

    def __init__(self, config: IsingConfig):
        self.config = config
        self.L = config.L
        self.J = config.J
        self.h = config.h
        self.T = config.T
        self.beta = 1.0 / config.T if config.T > 0 else np.inf

        if config.algorithm == "wolff" and config.h != 0.0:
            raise ValueError("Wolff no es válido con campo externo h != 0; usa 'metropolis'.")
        if config.seed is not None:
            _seed_numba(config.seed)

        # Initialize spins
        rng = np.random.default_rng(config.seed)
        self.spins = rng.choice([-1, 1], size=(self.L, self.L), p=[0.5, 0.5])

        # Observables storage
        self._energy = 0.0
        self._magnetization = 0.0
        self._update_observables()

    def _update_observables(self) -> None:
        """Compute energy and magnetization from current spin config."""
        J = self.J
        h = self.h
        s = self.spins

        # Nearest neighbor sum with periodic boundaries
        nb_sum = (
            np.roll(s, 1, axis=0)
            + np.roll(s, -1, axis=0)
            + np.roll(s, 1, axis=1)
            + np.roll(s, -1, axis=1)
        )
        self._energy = -0.5 * J * np.sum(s * nb_sum) - h * np.sum(s)
        self._magnetization = np.sum(s)

    @property
    def energy(self) -> float:
        return float(self._energy)

    @property
    def energy_per_spin(self) -> float:
        return float(self._energy / (self.L * self.L))

    @property
    def magnetization(self) -> float:
        return float(self._magnetization)

    @property
    def magnetization_per_spin(self) -> float:
        return float(self._magnetization / (self.L * self.L))

    @property
    def abs_magnetization_per_spin(self) -> float:
        return float(abs(self._magnetization) / (self.L * self.L))

    def metropolis_step(self, n_sweeps: int = 1) -> int:
        """Perform n_sweeps Metropolis updates. Returns total accepted."""
        total_accepted = 0
        for _ in range(n_sweeps):
            total_accepted += _metropolis_sweep(
                self.spins, self.J, self.h, self.beta, np.array([0])
            )
            self._update_observables()
        return total_accepted

    def wolff_step(self, n_clusters: int = 1) -> int:
        """Perform n_clusters Wolff cluster flips. Returns total cluster size."""
        total_size = 0
        for _ in range(n_clusters):
            total_size += _wolff_sweep(self.spins, self.J, self.beta, np.array([0]))
            self._update_observables()
        return total_size

    def step(self, n: int = 1) -> int:
        """Perform n steps of the configured algorithm."""
        if self.config.algorithm == "metropolis":
            return self.metropolis_step(n)
        return self.wolff_step(n)

    def thermalize(self, n_sweeps: int = 1000) -> None:
        """Thermalize the system (discard observables)."""
        for _ in range(n_sweeps):
            self.step(1)

    def measure(
        self,
        n_samples: int,
        sample_every: int = 1,
    ) -> dict[str, np.ndarray]:
        """
        Collect measurements.

        Returns dict with keys: energy, magnetization, energy2, magnetization2,
        susceptibility, specific_heat (per spin).
        """
        energies = np.zeros(n_samples)
        mags = np.zeros(n_samples)

        for i in range(n_samples):
            self.step(sample_every)
            energies[i] = self.energy_per_spin
            mags[i] = self.abs_magnetization_per_spin

        # Fluctuation observables
        e_mean = np.mean(energies)
        m_mean = np.mean(mags)
        e2_mean = np.mean(energies**2)
        m2_mean = np.mean(mags**2)

        chi = self.L * self.L * self.beta * (m2_mean - m_mean**2)  # susceptibility
        cv = self.L * self.L * self.beta**2 * (e2_mean - e_mean**2)  # specific heat

        return {
            "energy": energies,
            "magnetization": mags,
            "energy_mean": e_mean,
            "magnetization_mean": m_mean,
            "susceptibility": chi,
            "specific_heat": cv,
            "energy_var": e2_mean - e_mean**2,
            "magnetization_var": m2_mean - m_mean**2,
        }

    def copy(self) -> IsingModel:
        """Return a copy of the current state."""
        new_model = IsingModel(self.config)
        new_model.spins = self.spins.copy()
        new_model._update_observables()
        return new_model


class MetropolisSampler:
    """Metropolis single-spin-flip sampler."""

    def __init__(self, model: IsingModel):
        self.model = model

    def run(self, n_sweeps: int, thermalize: int = 1000) -> dict:
        self.model.thermalize(thermalize)
        return self.model.measure(n_sweeps, sample_every=1)


class WolffClusterSampler:
    """Wolff cluster sampler."""

    def __init__(self, model: IsingModel):
        self.model = model

    def run(self, n_clusters: int, thermalize: int = 100) -> dict:
        self.model.thermalize(thermalize)
        return self.model.measure(n_clusters, sample_every=1)
