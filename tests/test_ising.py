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
Tests for statphys-viz package.

Validates Ising model implementation against Onsager exact solution
and tests Monte Carlo samplers.
"""

from __future__ import annotations

import numpy as np
import pytest

from statphys_viz.benchmark import (
    ONSAGER_TC,
    CriticalExponents,
    OnsagerBenchmark,
    OnsagerExact,
    binder_cumulant,
    fit_critical_exponents,
)
from statphys_viz.ising import (
    IsingConfig,
    IsingModel,
    MetropolisSampler,
    WolffClusterSampler,
    _delta_energy,
    _periodic_idx,
)
from statphys_viz.visualize import SPIN_CMAP


class TestPeriodicIndex:
    """Test periodic boundary conditions."""

    def test_positive_index(self):
        assert _periodic_idx(5, 10) == 5

    def test_negative_index(self):
        assert _periodic_idx(-1, 10) == 9
        assert _periodic_idx(-5, 10) == 5

    def test_overflow_index(self):
        assert _periodic_idx(10, 10) == 0
        assert _periodic_idx(15, 10) == 5


class TestDeltaEnergy:
    """Test energy change calculation."""

    def test_all_up_spins(self):
        """All spins up: flipping one should cost energy."""
        spins = np.ones((4, 4), dtype=np.int8)
        # Interior spin: 4 neighbors up
        dE = _delta_energy(spins, 1, 1, J=1.0, h=0.0)
        assert dE == 8.0  # 2*J*1*4 = 8

    def test_all_down_spins(self):
        """All spins down: flipping one should cost energy."""
        spins = -np.ones((4, 4), dtype=np.int8)
        dE = _delta_energy(spins, 1, 1, J=1.0, h=0.0)
        assert dE == 8.0

    def test_antiferromagnetic_neighbors(self):
        """Spin surrounded by opposite spins: flipping lowers energy."""
        spins = np.ones((4, 4), dtype=np.int8)
        spins[1, 1] = -1  # center down, neighbors up
        dE = _delta_energy(spins, 1, 1, J=1.0, h=0.0)
        assert dE == -8.0  # 2*J*(-1)*(-4) = -8

    def test_with_field(self):
        """External field contributes to energy change."""
        spins = np.ones((4, 4), dtype=np.int8)
        dE = _delta_energy(spins, 1, 1, J=1.0, h=0.5)
        # Neighbors: 4 up, field: 2*h*s = 2*0.5*1 = 1
        assert dE == 9.0


class TestIsingModel:
    """Test Ising model basic functionality."""

    def setup_method(self):
        self.config = IsingConfig(L=16, T=2.269, J=1.0, h=0.0, seed=42)
        self.model = IsingModel(self.config)

    def test_initialization(self):
        assert self.model.L == 16
        assert self.model.J == 1.0
        assert self.model.h == 0.0
        assert abs(self.model.T - 2.269) < 1e-10
        assert self.model.spins.shape == (16, 16)
        assert np.all(np.isin(self.model.spins, [-1, 1]))

    def test_observables_computed(self):
        """Energy and magnetization should be computed."""
        assert isinstance(self.model.energy, float)
        assert isinstance(self.model.magnetization, float)
        assert isinstance(self.model.energy_per_spin, float)
        assert isinstance(self.model.magnetization_per_spin, float)

    def test_metropolis_step(self):
        """Metropolis step should run and update observables."""
        _ = self.model.energy_per_spin
        _ = self.model.abs_magnetization_per_spin
        accepted = self.model.metropolis_step(1)
        assert accepted >= 0
        assert accepted <= 256  # L*L
        # Observables should have changed (or stayed same)
        assert isinstance(self.model.energy_per_spin, float)

    def test_wolff_step(self):
        """Wolff step should run and update observables."""
        _ = self.model.energy_per_spin
        size = self.model.wolff_step(1)
        assert size >= 1
        assert size <= 256
        assert isinstance(self.model.energy_per_spin, float)

    def test_thermalize(self):
        """Thermalization should run without error."""
        self.model.thermalize(100)
        assert isinstance(self.model.energy_per_spin, float)

    def test_measure(self):
        """Measure should return all expected observables."""
        results = self.model.measure(100, sample_every=1)
        expected_keys = {
            "energy",
            "magnetization",
            "energy_mean",
            "magnetization_mean",
            "susceptibility",
            "specific_heat",
            "energy_var",
            "magnetization_var",
        }
        assert set(results.keys()) == expected_keys
        assert len(results["energy"]) == 100
        assert len(results["magnetization"]) == 100

    def test_copy(self):
        """Copy should create independent model with same state."""
        copy = self.model.copy()
        assert np.array_equal(copy.spins, self.model.spins)
        assert copy.energy_per_spin == self.model.energy_per_spin
        # Modify original
        self.model.step(1)
        assert not np.array_equal(copy.spins, self.model.spins)


class TestSamplers:
    """Test Monte Carlo samplers."""

    def test_metropolis_sampler(self):
        config = IsingConfig(L=16, T=3.0, seed=42)
        model = IsingModel(config)
        sampler = MetropolisSampler(model)
        results = sampler.run(100, thermalize=50)
        assert "energy" in results
        assert len(results["energy"]) == 100

    def test_wolff_sampler(self):
        config = IsingConfig(L=16, T=2.269, seed=42)
        model = IsingModel(config)
        sampler = WolffClusterSampler(model)
        results = sampler.run(100, thermalize=50)
        assert "energy" in results
        assert len(results["energy"]) == 100


class TestOnsagerExact:
    """Test Onsager exact solution values."""

    def test_Tc_value(self):
        assert abs(ONSAGER_TC - 2.26918531421) < 1e-10

    def test_above_Tc_magnetization_zero(self):
        exact = OnsagerExact(T=3.0, beta=1 / 3.0)
        assert exact.m_exact == 0.0

    def test_below_Tc_magnetization_positive(self):
        exact = OnsagerExact(T=2.0, beta=0.5)
        assert exact.m_exact > 0.0
        assert exact.m_exact < 1.0

    def test_energy_negative(self):
        """Energy should be negative for ferromagnetic coupling."""
        for T in [1.5, 2.0, 2.269, 3.0, 5.0]:
            exact = OnsagerExact(T=T, beta=1 / T)
            assert exact.u_exact < 0.0

    def test_specific_heat_positive(self):
        """Specific heat should be positive."""
        for T in [1.5, 2.0, 2.269, 3.0, 5.0]:
            exact = OnsagerExact(T=T, beta=1 / T)
            assert exact.c_exact >= 0.0


class TestOnsagerBenchmark:
    """Test benchmark against exact solution."""

    def setup_method(self):
        self.benchmark = OnsagerBenchmark(L=32)

    def test_exact_values(self):
        exact = self.benchmark.exact_values(2.269)
        assert exact.T == 2.269
        assert abs(exact.beta - 1 / 2.269) < 1e-10

    def test_relative_error(self):
        assert abs(self.benchmark.relative_error(1.0, 1.0) - 0.0) < 1e-10
        assert abs(self.benchmark.relative_error(1.1, 1.0) - 0.1) < 1e-10
        assert abs(self.benchmark.relative_error(0.9, 1.0) - 0.1) < 1e-10
        assert abs(self.benchmark.relative_error(0.01, 0.0) - 0.01) < 1e-10

    def test_validate_above_Tc(self):
        """Validate MC results above Tc (no magnetization check)."""
        mc_results = {
            "energy_mean": -1.0,
            "specific_heat": 0.5,
            "susceptibility": 10.0,
        }
        errors = self.benchmark.validate(mc_results, T=3.0)
        assert "energy_rel_error" in errors
        assert "specific_heat_rel_error" in errors
        assert "magnetization_rel_error" not in errors

    def test_validate_below_Tc(self):
        """Validate MC results below Tc (with magnetization)."""
        mc_results = {
            "energy_mean": -1.5,
            "magnetization_mean": 0.8,
            "specific_heat": 0.5,
            "susceptibility": 50.0,
        }
        errors = self.benchmark.validate(mc_results, T=2.0)
        assert "magnetization_rel_error" in errors


class TestCriticalExponents:
    """Test critical exponent measurements."""

    def test_exact_values(self):
        exp = CriticalExponents()
        assert exp.EXACT_BETA == 1 / 8
        assert exp.EXACT_GAMMA == 7 / 4
        assert exp.EXACT_NU == 1.0
        assert exp.EXACT_ALPHA == 0.0

    def test_relative_errors(self):
        exp = CriticalExponents(
            beta=0.13,
            beta_error=0.01,
            gamma=1.7,
            gamma_error=0.1,
            nu=1.05,
            nu_error=0.05,
        )
        errors = exp.relative_errors()
        assert "beta_rel_error" in errors
        assert "gamma_rel_error" in errors
        assert "nu_rel_error" in errors
        assert errors["beta_rel_error"] > 0
        assert errors["gamma_rel_error"] > 0


class TestBinderCumulant:
    """Test Binder cumulant calculation."""

    def test_disordered(self):
        """For disordered state, m2 > 0, m4 ≈ 3*m2^2 -> U4 ≈ 0."""
        U4 = binder_cumulant(m2=1.0, m4=3.0)
        assert abs(U4) < 0.01

    def test_ordered(self):
        """For fully ordered, m2 = m4 = 1 -> U4 = 2/3."""
        U4 = binder_cumulant(m2=1.0, m4=1.0)
        assert abs(U4 - 2 / 3) < 0.01

    def test_zero_m2(self):
        """Edge case: m2 = 0."""
        U4 = binder_cumulant(m2=0.0, m4=0.0)
        assert U4 == 0.0


class TestFitCriticalExponents:
    """Test finite-size scaling exponent fitting."""

    def test_power_law_fit(self):
        """Test fit on perfect power law: y = L^(-0.125)."""
        Ls = np.array([8, 16, 32, 64, 128])
        obs = Ls ** (-1 / 8)  # β/ν = 1/8 for 2D Ising
        slope, err = fit_critical_exponents(Ls, obs, "magnetization")
        assert abs(slope - (-1 / 8)) < 0.01
        assert err < 0.01

    def test_susceptibility_exponent(self):
        """Test fit for susceptibility: χ ~ L^(γ/ν) = L^(7/4)."""
        Ls = np.array([8, 16, 32, 64, 128])
        obs = Ls ** (7 / 4)
        slope, err = fit_critical_exponents(Ls, obs, "susceptibility")
        assert abs(slope - 7 / 4) < 0.01


class TestVisualization:
    """Test visualization components."""

    def test_spin_cmap(self):
        assert SPIN_CMAP.N == 2
        # Color 0 should be blue-ish, color 1 orange-ish
        c0 = SPIN_CMAP(0)[:3]
        c1 = SPIN_CMAP(1)[:3]
        assert c0[2] > c0[0]  # blue dominant
        assert c1[0] > c1[2]  # red dominant


# Integration test: full simulation cycle
@pytest.mark.slow
class TestIntegration:
    """Integration tests (marked slow)."""

    def test_metropolis_equilibration(self):
        """Metropolis should equilibrate at high T."""
        config = IsingConfig(L=16, T=5.0, seed=42)
        model = IsingModel(config)
        model.thermalize(500)
        results = model.measure(200)
        # At high T, energy should be near 0, magnetization near 0
        assert abs(results["energy_mean"]) < 0.5
        assert abs(results["magnetization_mean"]) < 0.2

    def test_wolff_at_Tc(self):
        """Wolff should sample efficiently at Tc."""
        config = IsingConfig(L=32, T=ONSAGER_TC, seed=42)
        model = IsingModel(config)
        model.thermalize(200)
        results = model.measure(500)
        # At Tc, susceptibility grows with L; for L=32 it's moderate
        assert results["susceptibility"] > 10
        assert results["specific_heat"] > 0.5
