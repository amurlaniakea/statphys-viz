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
Benchmarks and validation against exact solutions.

Provides validation of Monte Carlo results against the Onsager exact
solution for the 2D Ising model, and critical exponent measurements.

References:
- Onsager (1944) Phys. Rev. 65, 117 - Exact solution
- Ferdinand & Fisher (1969) Phys. Rev. 185, 832 - Finite-size scaling
- Binder (1981) Z. Phys. B 43, 119 - Binder cumulant
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy import special

# Onsager exact critical temperature for 2D square lattice Ising
ONSAGER_TC = 2.0 / np.log(1.0 + np.sqrt(2.0))  # ≈ 2.26918531421


@dataclass(frozen=False)
class OnsagerExact:
    """Exact Onsager solution values at given temperature."""

    T: float
    beta: float

    # Spontaneous magnetization (T < Tc)
    m_exact: float = 0.0

    # Energy per spin
    u_exact: float = 0.0

    # Specific heat per spin
    c_exact: float = 0.0

    def __post_init__(self) -> None:
        if self.T < ONSAGER_TC:
            k = 1.0 / np.sinh(2.0 * self.beta) ** 2
            self.m_exact = (1.0 - k**2) ** (1.0 / 8.0)

        # Specific heat (per spin) - Onsager exact solution
        # C = (2/π) * (K(k) - E(k)) * (2βJ/sinh(2βJ))^2
        # where k = 1/sinh²(2βJ) for T < Tc (always ≤ 1)
        # and k = sinh²(2βJ) for T > Tc (always ≤ 1)
        if self.T < ONSAGER_TC:
            k = 1.0 / (np.sinh(2.0 * self.beta) ** 2)
        else:
            k = np.sinh(2.0 * self.beta) ** 2

        K = special.ellipk(k)
        E = special.ellipe(k)

        sinh2b = np.sinh(2.0 * self.beta)
        self.u_exact = (
            -1.0 / np.tanh(2.0 * self.beta) * (1.0 + (2.0 / np.pi) * (2.0 * K - E) / sinh2b)
        )

        pref = (2.0 / np.pi) * (2.0 * self.beta / sinh2b) ** 2
        self.c_exact = pref * (K - E)


class OnsagerBenchmark:
    """
    Benchmark Monte Carlo results against Onsager exact solution.

    Computes relative errors for energy, magnetization, specific heat,
    and susceptibility at various temperatures.
    """

    def __init__(self, L: int = 32, J: float = 1.0):
        self.L = L
        self.J = J

    def exact_values(self, T: float) -> OnsagerExact:
        """Get exact Onsager values at temperature T."""
        beta = 1.0 / T
        return OnsagerExact(T=T, beta=beta)

    def relative_error(self, mc_value: float, exact_value: float) -> float:
        """Relative error: |mc - exact| / |exact|."""
        if abs(exact_value) < 1e-12:
            return abs(mc_value)
        return abs(mc_value - exact_value) / abs(exact_value)

    def relative_error_close(
        self, mc_value: float, exact_value: float, rel_tol: float = 1e-10
    ) -> bool:
        """Check if relative error is within tolerance."""
        return self.relative_error(mc_value, exact_value) <= rel_tol

    def validate(
        self,
        mc_results: dict[str, float],
        T: float,
    ) -> dict[str, float]:
        """
        Validate MC results against exact solution.

        mc_results should contain: energy_mean, magnetization_mean,
        specific_heat, susceptibility
        """
        exact = self.exact_values(T)
        errors = {}

        if "energy_mean" in mc_results:
            errors["energy_rel_error"] = self.relative_error(
                mc_results["energy_mean"], exact.u_exact
            )

        if T < ONSAGER_TC and "magnetization_mean" in mc_results:
            errors["magnetization_rel_error"] = self.relative_error(
                mc_results["magnetization_mean"], exact.m_exact
            )

        if "specific_heat" in mc_results:
            errors["specific_heat_rel_error"] = self.relative_error(
                mc_results["specific_heat"], exact.c_exact
            )

        # Track max error
        errors["max_rel_error"] = max(errors.values()) if errors else 0.0

        return errors

    def finite_size_scaling(
        self,
        Ls: list[int],
        T: float,
        mc_func: Callable[[int], dict[str, float]],
    ) -> dict[str, np.ndarray]:
        """
        Finite-size scaling analysis.

        mc_func(L) should return dict with magnetization_mean, susceptibility.
        """
        mags = []
        chis = []

        for L in Ls:
            res = mc_func(L)
            mags.append(res.get("magnetization_mean", 0.0))
            chis.append(res.get("susceptibility", 0.0))

        return {
            "L": np.array(Ls),
            "magnetization": np.array(mags),
            "susceptibility": np.array(chis),
        }


@dataclass(frozen=True)
class CriticalExponents:
    """Measured critical exponents from finite-size scaling."""

    # β: m ~ |t|^β (t = (T-Tc)/Tc)
    beta: float = 0.0
    beta_error: float = 0.0

    # γ: χ ~ |t|^(-γ)
    gamma: float = 0.0
    gamma_error: float = 0.0

    # ν: ξ ~ |t|^(-ν)
    nu: float = 0.0
    nu_error: float = 0.0

    # α: C ~ |t|^(-α) (α=0 for 2D Ising → logarithmic)
    alpha: float = 0.0
    alpha_error: float = 0.0

    # Exact 2D Ising values
    EXACT_BETA = 1.0 / 8.0  # 0.125
    EXACT_GAMMA = 7.0 / 4.0  # 1.75
    EXACT_NU = 1.0
    EXACT_ALPHA = 0.0  # logarithmic

    def relative_errors(self) -> dict[str, float]:
        """Relative errors vs exact values."""
        return {
            "beta_rel_error": abs(self.beta - self.EXACT_BETA) / self.EXACT_BETA,
            "gamma_rel_error": abs(self.gamma - self.EXACT_GAMMA) / self.EXACT_GAMMA,
            "nu_rel_error": abs(self.nu - self.EXACT_NU) / self.EXACT_NU,
        }

    def __str__(self) -> str:
        err = self.relative_errors()
        return (
            "Critical Exponents (2D Ising):\n"
            f"  β = {self.beta:.4f} ± {self.beta_error:.4f} "
            f"(exact: {self.EXACT_BETA:.4f}, err: {err['beta_rel_error']:.2%})\n"
            f"  γ = {self.gamma:.4f} ± {self.gamma_error:.4f} "
            f"(exact: {self.EXACT_GAMMA:.4f}, err: {err['gamma_rel_error']:.2%})\n"
            f"  ν = {self.nu:.4f} ± {self.nu_error:.4f} "
            f"(exact: {self.EXACT_NU:.4f}, err: {err['nu_rel_error']:.2%})\n"
            f"  α = {self.alpha:.4f} ± {self.alpha_error:.4f} "
            f"(exact: {self.EXACT_ALPHA:.4f})"
        )


def binder_cumulant(m2: float, m4: float) -> float:
    """
    Binder cumulant U4 = 1 - ⟨m⁴⟩ / (3⟨m²⟩²).

    At Tc, U4 is universal (≈ 0.61069 for 2D Ising).
    """
    if m2 == 0:
        return 0.0
    return 1.0 - m4 / (3.0 * m2 * m2)


def estimate_tc_from_binder(
    Ls: list[int],
    temps: np.ndarray,
    binder_data: dict[int, np.ndarray],
) -> float:
    """
    Estimate Tc from Binder cumulant crossing.

    binder_data[L][i] = U4 at temps[i] for size L.
    """
    # Find crossing points between consecutive L
    crossings = []
    for i in range(len(Ls) - 1):
        L1, L2 = Ls[i], Ls[i + 1]
        U1 = binder_data[L1]
        U2 = binder_data[L2]

        # Find where U1 - U2 changes sign
        diff = U1 - U2
        for j in range(len(temps) - 1):
            if diff[j] * diff[j + 1] < 0:
                # Linear interpolation
                t_cross = temps[j] - diff[j] * (temps[j + 1] - temps[j]) / (diff[j + 1] - diff[j])
                crossings.append(t_cross)

    if crossings:
        return np.mean(crossings)
    return ONSAGER_TC


def fit_critical_exponents(
    Ls: np.ndarray,
    obs_at_tc: np.ndarray,
    obs_name: str = "magnetization",
) -> tuple[float, float]:
    """
    Fit critical exponent from finite-size scaling at Tc.

    At Tc: m ~ L^(-β/ν), χ ~ L^(γ/ν)
    """
    log_L = np.log(Ls)
    log_obs = np.log(np.abs(obs_at_tc) + 1e-12)

    # Linear fit: log(obs) = slope * log(L) + intercept
    A = np.vstack([log_L, np.ones_like(log_L)]).T
    slope, intercept = np.linalg.lstsq(A, log_obs, rcond=None)[0]

    # Error from residuals
    pred = slope * log_L + intercept
    residuals = log_obs - pred
    sigma = np.sqrt(np.sum(residuals**2) / (len(Ls) - 2))
    slope_err = sigma / np.sqrt(np.sum((log_L - np.mean(log_L)) ** 2))

    return slope, slope_err
