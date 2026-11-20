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
Jupyter widgets for interactive Ising model exploration.

Provides ipywidgets-based controls for real-time parameter adjustment,
live spin configuration display, and observable plotting.
"""

from __future__ import annotations

import threading
import time
from typing import Literal

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import clear_output, display
from ipywidgets import (
    Button,
    Checkbox,
    Dropdown,
    FloatSlider,
    HBox,
    IntSlider,
    Layout,
    Output,
    VBox,
)

from .ising import IsingConfig, IsingModel
from .visualize import ONSAGER_TC, plot_spin_config


class IsingWidget:
    """
    Interactive Ising model widget for Jupyter notebooks.

    Features:
    - Real-time spin lattice visualization (ipycanvas)
    - Parameter sliders (T, L, J, h, algorithm)
    - Live observable plots (energy, magnetization)
    - Thermalization and measurement controls
    - Export capabilities
    """

    def __init__(
        self,
        L: int = 32,
        T: float = 2.269,
        J: float = 1.0,
        h: float = 0.0,
        algorithm: Literal["metropolis", "wolff"] = "wolff",
        canvas_size: int = 300,
    ):
        self.L = L
        self.config = IsingConfig(L=L, T=T, J=J, h=h, algorithm=algorithm)
        self.model = IsingModel(self.config)
        self.canvas_size = canvas_size

        # State
        self._running = False
        self._thread: threading.Thread | None = None
        self._spin_history: list[np.ndarray] = []
        self._obs_history: dict[str, list] = {
            "energy": [],
            "magnetization": [],
            "step": [],
        }

        # Create UI components
        self._create_widgets()
        self._create_canvas()
        self._create_plots()
        self._layout = self._build_layout()

    def _create_widgets(self) -> None:
        """Create all control widgets."""
        style = {"description_width": "initial"}
        slider_layout = Layout(width="300px")

        # Temperature (highlight Tc)
        self.temp_slider = FloatSlider(
            value=self.config.T,
            min=0.5,
            max=5.0,
            step=0.01,
            description="Temperature (T):",
            style=style,
            layout=slider_layout,
            continuous_update=False,
            readout_format=".3f",
        )
        # Mark Tc
        self.temp_slider.description = f"Temperature (Tc≈{ONSAGER_TC:.3f}):"

        # Lattice size
        self.L_slider = IntSlider(
            value=self.L,
            min=8,
            max=128,
            step=4,
            description="Lattice size (L):",
            style=style,
            layout=slider_layout,
            continuous_update=False,
        )

        # Coupling
        self.J_slider = FloatSlider(
            value=self.config.J,
            min=0.1,
            max=3.0,
            step=0.1,
            description="Coupling (J):",
            style=style,
            layout=slider_layout,
            continuous_update=False,
        )

        # External field
        self.h_slider = FloatSlider(
            value=self.config.h,
            min=-2.0,
            max=2.0,
            step=0.05,
            description="Field (h):",
            style=style,
            layout=slider_layout,
            continuous_update=False,
        )

        # Algorithm
        self.algo_dropdown = Dropdown(
            options=["metropolis", "wolff"],
            value=self.config.algorithm,
            description="Algorithm:",
            style=style,
            layout=Layout(width="200px"),
        )

        # Control buttons
        self.run_button = Button(
            description="▶ Run",
            button_style="success",
            layout=Layout(width="100px"),
        )
        self.stop_button = Button(
            description="■ Stop",
            button_style="danger",
            layout=Layout(width="100px"),
            disabled=True,
        )
        self.reset_button = Button(
            description="↻ Reset",
            button_style="warning",
            layout=Layout(width="100px"),
        )
        self.step_button = Button(
            description="Step ×100",
            layout=Layout(width="120px"),
        )

        # Measurement controls
        self.measure_checkbox = Checkbox(
            value=True,
            description="Measure observables",
            style=style,
        )
        self.thermalize_slider = IntSlider(
            value=100,
            min=0,
            max=5000,
            step=50,
            description="Thermalize:",
            style=style,
            layout=slider_layout,
        )

        # Output areas
        self.canvas_output = Output()
        self.plot_output = Output()
        self.info_output = Output()

        # Link callbacks
        self.temp_slider.observe(self._on_param_change, "value")
        self.L_slider.observe(self._on_L_change, "value")
        self.J_slider.observe(self._on_param_change, "value")
        self.h_slider.observe(self._on_param_change, "value")
        self.algo_dropdown.observe(self._on_algo_change, "value")
        self.run_button.on_click(self._on_run)
        self.stop_button.on_click(self._on_stop)
        self.reset_button.on_click(self._on_reset)
        self.step_button.on_click(self._on_step)

    def _create_canvas(self) -> None:
        """Create ipycanvas for real-time spin display."""
        try:
            from ipycanvas import Canvas

            self.canvas = Canvas(
                width=self.canvas_size,
                height=self.canvas_size,
                sync_image_data=True,
            )
            self._draw_spins()
        except ImportError:
            self.canvas = None
            with self.canvas_output:
                msg = "ipycanvas not installed. Install with: pip install ipycanvas\n"
                self.canvas_output.append_stdout(msg)

    def _create_plots(self) -> None:
        """Create matplotlib figure for observables."""
        self.fig, self.axes = plt.subplots(2, 2, figsize=(10, 7))
        self.fig.tight_layout(pad=3.0)

        self.ax_energy = self.axes[0, 0]
        self.ax_mag = self.axes[0, 1]
        self.ax_spins = self.axes[1, 0]
        self.ax_hist = self.axes[1, 1]

        for ax in self.axes.flat:
            ax.grid(True, alpha=0.3)

        self.ax_energy.set_xlabel("Step")
        self.ax_energy.set_ylabel("Energy/spin")
        self.ax_mag.set_xlabel("Step")
        self.ax_mag.set_ylabel("|M|/spin")
        self.ax_spins.set_title("Current Configuration")
        self.ax_spins.set_xticks([])
        self.ax_spins.set_yticks([])
        self.ax_hist.set_xlabel("Magnetization")
        self.ax_hist.set_ylabel("Frequency")

    def _build_layout(self) -> VBox:
        """Build the complete widget layout."""
        controls = VBox(
            [
                widgets.HTML("<h3>Ising Model Controls</h3>"),
                self.temp_slider,
                self.L_slider,
                self.J_slider,
                self.h_slider,
                self.algo_dropdown,
                widgets.HTML("<hr>"),
                self.thermalize_slider,
                self.measure_checkbox,
                HBox([self.run_button, self.stop_button, self.reset_button, self.step_button]),
                widgets.HTML("<hr>"),
                self.info_output,
            ],
            layout=Layout(width="350px", border="1px solid #ddd", padding="10px"),
        )

        display_area = VBox(
            [
                self.canvas_output,
                self.plot_output,
            ]
        )

        return HBox([controls, display_area])

    def _draw_spins(self) -> None:
        """Draw current spin configuration on canvas."""
        if self.canvas is None:
            return

        spins = self.model.spins
        L = self.L
        cell_size = self.canvas_size / L

        self.canvas.fill_style = "#ffffff"
        self.canvas.fill_rect(0, 0, self.canvas_size, self.canvas_size)

        for i in range(L):
            for j in range(L):
                x = j * cell_size
                y = i * cell_size
                if spins[i, j] == 1:
                    self.canvas.fill_style = "#ff7f0e"  # orange for +1
                else:
                    self.canvas.fill_style = "#1f77b4"  # blue for -1
                self.canvas.fill_rect(x, y, cell_size, cell_size)

    def _update_plots(self) -> None:
        """Update observable plots."""
        if not self._obs_history["step"]:
            return

        steps = self._obs_history["step"]
        energies = self._obs_history["energy"]
        mags = self._obs_history["magnetization"]

        # Energy
        self.ax_energy.clear()
        self.ax_energy.plot(steps, energies, "b-", linewidth=0.8)
        self.ax_energy.set_xlabel("Step")
        self.ax_energy.set_ylabel("Energy/spin")
        self.ax_energy.grid(True, alpha=0.3)

        # Magnetization
        self.ax_mag.clear()
        self.ax_mag.plot(steps, mags, "r-", linewidth=0.8)
        self.ax_mag.set_xlabel("Step")
        self.ax_mag.set_ylabel("|M|/spin")
        self.ax_mag.grid(True, alpha=0.3)

        # Spin config
        self.ax_spins.clear()
        plot_spin_config(self.model.spins, ax=self.ax_spins, title="Current State")
        self.ax_spins.set_xticks([])
        self.ax_spins.set_yticks([])

        # Magnetization histogram
        self.ax_hist.clear()
        self.ax_hist.hist(mags, bins=30, density=True, alpha=0.7, color="purple")
        self.ax_hist.set_xlabel("Magnetization")
        self.ax_hist.set_ylabel("Frequency")
        self.ax_hist.grid(True, alpha=0.3)

        self.fig.tight_layout()

        with self.plot_output:
            clear_output(wait=True)
            display(self.fig)

    def _on_param_change(self, change) -> None:
        """Handle parameter changes (T, J, h)."""
        if change["name"] == "value":
            # Wolff no admite campo externo: si se mueve h con Wolff, pasar a Metropolis
            # (el cambio del desplegable reconstruye el modelo).
            if (
                change.get("owner") is self.h_slider
                and self.h_slider.value != 0
                and self.algo_dropdown.value == "wolff"
            ):
                self.algo_dropdown.value = "metropolis"
                return
            self.config = IsingConfig(
                L=self.L,
                T=self.temp_slider.value,
                J=self.J_slider.value,
                h=self.h_slider.value,
                algorithm=self.algo_dropdown.value,
            )
            self.model = IsingModel(self.config)
            self._draw_spins()
            self._update_info()

    def _on_L_change(self, change) -> None:
        """Handle lattice size change (requires full reset)."""
        if change["name"] == "value":
            self.L = self.L_slider.value
            self._on_param_change(change)

    def _on_algo_change(self, change) -> None:
        """Handle algorithm change."""
        if change["name"] == "value":
            # Wolff solo es válido con h = 0: poner h a 0 (esto reconstruye el modelo).
            if self.algo_dropdown.value == "wolff" and self.h_slider.value != 0:
                self.h_slider.value = 0
                return
            self.config = IsingConfig(
                L=self.L,
                T=self.temp_slider.value,
                J=self.J_slider.value,
                h=self.h_slider.value,
                algorithm=self.algo_dropdown.value,
            )
            self.model = IsingModel(self.config)
            self._draw_spins()
            self._update_info()

    def _on_run(self, btn) -> None:
        """Start simulation."""
        if self._running:
            return

        self._running = True
        self.run_button.disabled = True
        self.stop_button.disabled = False
        self.reset_button.disabled = True
        self.temp_slider.disabled = True
        self.L_slider.disabled = True
        self.J_slider.disabled = True
        self.h_slider.disabled = True
        self.algo_dropdown.disabled = True

        # Thermalize first
        thermalize_steps = self.thermalize_slider.value
        if thermalize_steps > 0:
            self.model.thermalize(thermalize_steps)

        self._spin_history.clear()
        self._obs_history = {"energy": [], "magnetization": [], "step": []}

        self._thread = threading.Thread(target=self._run_loop, daemon=True)
        self._thread.start()

    def _run_loop(self) -> None:
        """Main simulation loop (runs in thread)."""
        step = 0
        measure = self.measure_checkbox.value

        while self._running:
            # Run one step
            self.model.step(1)

            if measure:
                self._obs_history["step"].append(step)
                self._obs_history["energy"].append(self.model.energy_per_spin)
                self._obs_history["magnetization"].append(self.model.abs_magnetization_per_spin)

            # Update canvas (must be in main thread for ipycanvas)
            # We'll update via the canvas_output context
            with self.canvas_output:
                self._draw_spins()

            step += 1

            # Update plots every 10 steps
            if step % 10 == 0:
                self._update_plots()

            time.sleep(0.01)  # ~100 steps/sec max

    def _on_stop(self, btn) -> None:
        """Stop simulation."""
        self._running = False
        if self._thread:
            self._thread.join(timeout=1.0)
            self._thread = None

        self.run_button.disabled = False
        self.stop_button.disabled = True
        self.reset_button.disabled = False
        self.temp_slider.disabled = False
        self.L_slider.disabled = False
        self.J_slider.disabled = False
        self.h_slider.disabled = False
        self.algo_dropdown.disabled = False

        self._update_plots()
        self._update_info()

    def _on_reset(self, btn) -> None:
        """Reset to initial state."""
        self._on_stop(btn)
        self.model = IsingModel(self.config)
        self._spin_history.clear()
        self._obs_history = {"energy": [], "magnetization": [], "step": []}
        self._draw_spins()

        # Clear plots
        for ax in self.axes.flat:
            ax.clear()
            ax.grid(True, alpha=0.3)

        self.ax_energy.set_xlabel("Step")
        self.ax_energy.set_ylabel("Energy/spin")
        self.ax_mag.set_xlabel("Step")
        self.ax_mag.set_ylabel("|M|/spin")
        self.ax_spins.set_title("Current Configuration")
        self.ax_spins.set_xticks([])
        self.ax_spins.set_yticks([])
        self.ax_hist.set_xlabel("Magnetization")
        self.ax_hist.set_ylabel("Frequency")

        with self.plot_output:
            clear_output(wait=True)
            display(self.fig)

        self._update_info()

    def _on_step(self, btn) -> None:
        """Run 100 steps and measure."""
        for _ in range(100):
            self.model.step(1)
            self._obs_history["step"].append(len(self._obs_history["step"]))
            self._obs_history["energy"].append(self.model.energy_per_spin)
            self._obs_history["magnetization"].append(self.model.abs_magnetization_per_spin)
        self._draw_spins()
        self._update_plots()
        self._update_info()

    def _update_info(self) -> None:
        """Update info panel."""
        with self.info_output:
            clear_output(wait=True)
            info_lines = [
                f"L = {self.L}  |  T = {self.config.T:.3f}  |  "
                f"J = {self.config.J:.1f}  |  h = {self.config.h:.2f}",
                f"Algorithm: {self.config.algorithm}",
                f"Tc (Onsager) = {ONSAGER_TC:.4f}",
                f"Energy/spin = {self.model.energy_per_spin:.6f}",
                f"|M|/spin = {self.model.abs_magnetization_per_spin:.6f}",
            ]
            if self._obs_history["step"]:
                info_lines.append(f"Steps measured: {len(self._obs_history['step'])}")
            self.info_output.append_stdout("\n".join(info_lines) + "\n")

    def display(self) -> None:
        """Display the widget."""
        display(self._layout)
        self._update_info()

    def get_results(self) -> dict[str, np.ndarray]:
        """Return measured observables as numpy arrays."""
        return {k: np.array(v) for k, v in self._obs_history.items()}

    def save_state(self, filepath: str) -> None:
        """Save current spin configuration."""
        np.savez(
            filepath,
            spins=self.model.spins,
            config=np.array([self.L, self.config.T, self.config.J, self.config.h]),
            history=self._obs_history,
        )

    @classmethod
    def load_state(cls, filepath: str) -> IsingWidget:
        """Load spin configuration from file."""
        data = np.load(filepath, allow_pickle=True)
        L, T, J, h = data["config"]
        widget = cls(L=int(L), T=float(T), J=float(J), h=float(h))
        widget.model.spins = data["spins"]
        widget.model._update_observables()
        widget._obs_history = data["history"].item()
        widget._draw_spins()
        widget._update_plots()
        return widget


def create_temperature_scan_widget(
    L: int = 32,
    T_min: float = 1.0,
    T_max: float = 4.0,
    n_points: int = 30,
    n_samples: int = 1000,
    thermalize: int = 500,
) -> widgets.VBox:
    """
    Create widget for automated temperature scan with exact comparison.

    Returns a VBox with run button and output area.
    """
    run_btn = Button(description="Run Temperature Scan", button_style="success")
    progress = widgets.IntProgress(min=0, max=n_points, description="Progress:")
    output = Output()

    def on_run(btn):
        btn.disabled = True
        progress.value = 0

        temps = np.linspace(T_min, T_max, n_points)
        results = {
            "energy_mean": [],
            "magnetization_mean": [],
            "specific_heat": [],
            "susceptibility": [],
        }

        from .benchmark import OnsagerBenchmark

        benchmark = OnsagerBenchmark(L=L)

        with output:
            clear_output(wait=True)
            output.append_stdout(f"Scanning T ∈ [{T_min:.2f}, {T_max:.2f}] with L={L}\n")
            display(progress)

        for i, T in enumerate(temps):
            config = IsingConfig(L=L, T=T, algorithm="wolff")
            model = IsingModel(config)
            model.thermalize(thermalize)
            meas = model.measure(n_samples)

            results["energy_mean"].append(meas["energy_mean"])
            results["magnetization_mean"].append(meas["magnetization_mean"])
            results["specific_heat"].append(meas["specific_heat"])
            results["susceptibility"].append(meas["susceptibility"])

            progress.value = i + 1

        # Plot results
        from .visualize import plot_observables

        fig = plot_observables(results, temps, benchmark.exact_values)

        with output:
            clear_output(wait=True)
            display(fig)

            # Print errors
            output.append_stdout("\nRelative Errors vs Onsager Exact:\n")
            for i, T in enumerate(temps):
                mc_res = {k: v[i] for k, v in results.items()}
                errors = benchmark.validate(mc_res, T)
                if errors:
                    msg = f"  T={T:.3f}: max_rel_error={errors['max_rel_error']:.2%}\n"
                    output.append_stdout(msg)

        btn.disabled = False

    run_btn.on_click(on_run)
    return VBox([run_btn, progress, output])
