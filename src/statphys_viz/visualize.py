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
Visualization components for Ising model.

Provides matplotlib-based static plots and ipycanvas-based
interactive animations for Jupyter notebooks.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap

# Color maps for spin visualization
SPIN_CMAP = ListedColormap(["#1f77b4", "#ff7f0e"])  # blue for -1, orange for +1
ENERGY_CMAP = "RdBu_r"
MAG_CMAP = "RdBu"


def plot_spin_config(
    spins: np.ndarray,
    ax: plt.Axes | None = None,
    title: str = "",
    cmap: ListedColormap | None = None,
    show_colorbar: bool = False,
) -> plt.Axes:
    """
    Plot a single spin configuration.

    Args:
        spins: 2D array of spins (±1)
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        cmap: Colormap (uses SPIN_CMAP if None)
        show_colorbar: Whether to show colorbar

    Returns:
        The axes object
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    if cmap is None:
        cmap = SPIN_CMAP

    im = ax.imshow(spins, cmap=cmap, vmin=-1, vmax=1, interpolation="nearest")
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])

    if show_colorbar:
        cbar = plt.colorbar(im, ax=ax, ticks=[-1, 1], shrink=0.8)
        cbar.ax.set_yticklabels(["−1", "+1"])

    return ax


def plot_observables(
    results: dict[str, np.ndarray],
    temps: np.ndarray,
    exact_func: callable | None = None,
    figsize: tuple = (12, 8),
) -> plt.Figure:
    """
    Plot Monte Carlo observables vs temperature with exact solution.

    Args:
        results: Dict with keys 'energy_mean', 'magnetization_mean',
                'specific_heat', 'susceptibility' (arrays matching temps)
        temps: Temperature array
        exact_func: Function(T) returning exact values (optional)
        figsize: Figure size

    Returns:
        Matplotlib figure
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.ravel()

    # Energy
    ax = axes[0]
    ax.plot(temps, results["energy_mean"], "o-", label="MC", markersize=4)
    if exact_func:
        exact_e = [exact_func(T).u_exact for T in temps]
        ax.plot(temps, exact_e, "--", label="Onsager exact", color="red")
    ax.set_xlabel("Temperature (T)")
    ax.set_ylabel("Energy per spin")
    ax.axvline(ONSAGER_TC, color="gray", linestyle=":", alpha=0.5)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Magnetization
    ax = axes[1]
    ax.plot(temps, results["magnetization_mean"], "o-", label="MC", markersize=4)
    if exact_func:
        exact_m = [exact_func(T).m_exact for T in temps]
        ax.plot(temps, exact_m, "--", label="Onsager exact", color="red")
    ax.set_xlabel("Temperature (T)")
    ax.set_ylabel("|Magnetization| per spin")
    ax.axvline(ONSAGER_TC, color="gray", linestyle=":", alpha=0.5)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Specific heat
    ax = axes[2]
    ax.plot(temps, results["specific_heat"], "o-", label="MC", markersize=4)
    if exact_func:
        exact_c = [exact_func(T).c_exact for T in temps]
        ax.plot(temps, exact_c, "--", label="Onsager exact", color="red")
    ax.set_xlabel("Temperature (T)")
    ax.set_ylabel("Specific heat per spin")
    ax.axvline(ONSAGER_TC, color="gray", linestyle=":", alpha=0.5)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Susceptibility
    ax = axes[3]
    ax.plot(temps, results["susceptibility"], "o-", label="MC", markersize=4)
    ax.set_xlabel("Temperature (T)")
    ax.set_ylabel("Susceptibility")
    ax.axvline(ONSAGER_TC, color="gray", linestyle=":", alpha=0.5)
    ax.set_yscale("log")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def animate_spins(
    spin_history: list[np.ndarray],
    interval: int = 100,
    figsize: tuple = (6, 6),
    title_prefix: str = "Step",
) -> FuncAnimation:
    """
    Create animation of spin configurations.

    Args:
        spin_history: List of spin configurations (each 2D array)
        interval: Milliseconds between frames
        figsize: Figure size
        title_prefix: Prefix for frame title

    Returns:
        FuncAnimation object (call .save() or display in Jupyter)
    """
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(spin_history[0], cmap=SPIN_CMAP, vmin=-1, vmax=1, interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    title = ax.set_title(f"{title_prefix} 0")

    def update(frame):
        im.set_array(spin_history[frame])
        title.set_text(f"{title_prefix} {frame}")
        return im, title

    anim = FuncAnimation(fig, update, frames=len(spin_history), interval=interval, blit=True)
    return anim


def plot_binder_cumulant(
    binder_data: dict[int, np.ndarray],
    temps: np.ndarray,
    figsize: tuple = (8, 6),
) -> plt.Figure:
    """
    Plot Binder cumulant for different lattice sizes.

    Args:
        binder_data: Dict[L] = array of U4 values at each temperature
        temps: Temperature array
        figsize: Figure size

    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    for L, U4 in binder_data.items():
        ax.plot(temps, U4, "o-", label=f"L={L}", markersize=4)

    ax.axhline(0.61069, color="red", linestyle="--", label="Universal U* ≈ 0.61069")
    ax.axvline(ONSAGER_TC, color="gray", linestyle=":", alpha=0.5, label=f"Tc = {ONSAGER_TC:.4f}")

    ax.set_xlabel("Temperature (T)")
    ax.set_ylabel("Binder cumulant U₄")
    ax.set_title("Binder Cumulant Crossing (Tc estimation)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    return fig


def plot_finite_size_scaling(
    Ls: np.ndarray,
    obs: np.ndarray,
    obs_name: str = "magnetization",
    figsize: tuple = (8, 6),
) -> tuple[plt.Figure, tuple[float, float]]:
    """
    Plot finite-size scaling and fit exponent.

    At Tc: obs ~ L^(-β/ν) for magnetization, L^(γ/ν) for susceptibility

    Returns:
        (figure, (slope, slope_error))
    """
    fig, ax = plt.subplots(figsize=figsize)

    log_L = np.log(Ls)
    log_obs = np.log(np.abs(obs) + 1e-12)

    ax.loglog(Ls, np.abs(obs) + 1e-12, "o-", markersize=8)

    # Linear fit in log-log
    A = np.vstack([log_L, np.ones_like(log_L)]).T
    slope, intercept = np.linalg.lstsq(A, log_obs, rcond=None)[0]

    # Plot fit
    fit_line = np.exp(slope * log_L + intercept)
    ax.loglog(Ls, fit_line, "--", color="red", label=f"Fit: slope = {slope:.4f}")

    # Error
    pred = slope * log_L + intercept
    residuals = log_obs - pred
    sigma = np.sqrt(np.sum(residuals**2) / (len(Ls) - 2))
    slope_err = sigma / np.sqrt(np.sum((log_L - np.mean(log_L)) ** 2))

    ax.set_xlabel("L")
    ax.set_ylabel(f"|{obs_name}| at Tc")
    ax.set_title(f"Finite-Size Scaling at Tc ({obs_name})")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")

    return fig, (slope, slope_err)


# Constant for Onsager Tc
ONSAGER_TC = 2.0 / np.log(1.0 + np.sqrt(2.0))
