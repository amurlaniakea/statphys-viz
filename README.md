<!-- SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com> -->
<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# 🧲 statphys-viz

> Interactive statistical physics visualizations for education: Ising model Monte Carlo with Jupyter widgets and live web demo.

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0.txt)
[![CI](https://github.com/amurlaniakea/statphys-viz/actions/workflows/ci.yml/badge.svg)](https://github.com/amurlaniakea/statphys-viz/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24389055.svg)](https://doi.org/10.5281/zenodo.24389055)

🌐 **Demo en vivo:** <https://amurlaniakea.github.io/statphys-viz/>

📚 **Documentación:** <https://amurlaniakea.github.io/statphys-viz/docs/>

**statphys-viz** provides a complete, verified implementation of the 2D Ising model for statistical physics education. It includes both Metropolis and Wolff cluster Monte Carlo algorithms, real-time visualization via Jupyter widgets, and quantitative validation against the Onsager exact solution.

---

## Problema

Teaching statistical physics requires bridging abstract theory (partition functions, phase transitions, critical exponents) with concrete computational experiments. Existing tools are either:
- **Too simple**: Fixed demos without parameter exploration
- **Too complex**: Full research codes (ALPS, LAMMPS) with steep learning curves
- **Unverified**: No quantitative benchmarks against exact solutions

Students need an interactive environment where they can *see* the phase transition, *measure* critical exponents, and *compare* against the Onsager exact solution — all in a browser or notebook.

## Solución

statphys-viz delivers:

| Component | Function | Validation |
|-----------|----------|------------|
| **IsingModel** | 2D square lattice with periodic boundaries | Energy/magnetization conservation |
| **MetropolisSampler** | Single-spin-flip MCMC | Detailed balance test |
| **WolffClusterSampler** | Cluster algorithm (critical slowing down ↓) | Autocorrelation time measurement |
| **OnsagerBenchmark** | Exact solution comparison | Relative error tracking |
| **IsingWidget** | Jupyter ipywidgets + ipycanvas live demo | Real-time parameter sweeps |
| **Web Demo** | Static GitHub Pages deployment | Zero-install access |

All observables are validated against the Onsager exact solution (1944) with relative errors reported quantitatively.

## Arquitectura

```
statphys-viz/
├── src/statphys_viz/
│   ├── __init__.py          # Public API
│   ├── ising.py             # Core Ising model + MCMC samplers (numba JIT)
│   ├── benchmark.py         # Onsager exact solution + finite-size scaling
│   ├── visualize.py         # Matplotlib plots + animations
│   └── widgets.py           # Jupyter interactive widgets
├── tests/
│   └── test_ising.py        # Unit + integration tests (validated vs Onsager)
├── examples/
│   └── temperature_scan.ipynb  # Notebook demo
├── docs/                    # Sphinx documentation
└── .github/workflows/ci.yml # CI: pytest + ruff + bandit on 3.10/3.11/3.12
```

**Key technical choices:**
- **numba JIT** on MCMC loops → ~100× speedup vs pure Python
- **Wolff cluster** → eliminates critical slowing down at Tc
- **Onsager elliptic integrals** (scipy.special) → exact reference
- **ipycanvas** → 60 FPS spin lattice rendering in notebooks
- **GitHub Pages** → zero-install web demo for classrooms

## Instalación

```bash
# From PyPI (not yet published — install from source below)
# pip install statphys-viz

# From source (development)
git clone https://github.com/amurlaniakea/statphys-viz.git
cd statphys-viz
pip install -e ".[dev]"   # includes pytest, ruff, bandit
```

**Requirements:** Python ≥ 3.10, NumPy, Numba, Matplotlib, SciPy, ipywidgets.

## Uso

### Quick Start (Python)

```python
from statphys_viz import IsingConfig, IsingModel, OnsagerBenchmark

# Configure the model (T=2.0: ordered phase; near Tc the error vs Onsager is large due to finite size)
config = IsingConfig(L=64, T=2.0, J=1.0, algorithm="wolff")
model = IsingModel(config)

# Thermalize and measure
model.thermalize(500)
results = model.measure(2000)

print(f"Energy/spin: {results['energy_mean']:.6f}")
print(f"|M|/spin:    {results['magnetization_mean']:.6f}")
print(f"χ:           {results['susceptibility']:.2f}")
print(f"C:           {results['specific_heat']:.4f}")

# Validate against Onsager exact
benchmark = OnsagerBenchmark(L=64)
errors = benchmark.validate(results, T=2.0)
print(f"Max relative error: {errors['max_rel_error']:.2%}")
```

### Jupyter Interactive Widget

```python
from statphys_viz.widgets import IsingWidget

widget = IsingWidget(L=64, T=2.269, algorithm="wolff")
widget.display()
```

Features:
- **Live spin lattice** (ipycanvas, 60 FPS)
- **Temperature slider** with Tc marker
- **Real-time energy/magnetization plots**
- **Algorithm switch** (Metropolis ↔ Wolff)
- **One-click thermalization & measurement**

### Temperature Scan (Automated Benchmark)

```python
from statphys_viz.widgets import create_temperature_scan_widget

scan = create_temperature_scan_widget(L=32, T_min=1.5, T_max=3.5, n_points=25)
scan  # displays run button + progress + auto-plot with exact comparison
```

Outputs a 2×2 plot grid (E, |M|, C, χ) with Onsager exact curves and relative error table.

### Command Line

```bash
# Run benchmark suite
python -m statphys_viz.benchmark --L 32 --T-min 1.5 --T-max 3.5 --points 20

# Generate validation plot
python -m statphys_viz.visualize --output validation.png
```

## Tests

```bash
# All tests
python -m pytest tests/ -v

# Fast tests only (skip slow integration)
python -m pytest tests/ -v -m "not slow"

# With coverage
python -m pytest tests/ --cov=statphys_viz --cov-report=term-missing
```

**Test categories:**
- **Unit**: Periodic boundaries, ΔE calculation, sampler interfaces
- **Validation**: Onsager exact solution (energy, M, C, χ)
- **Critical exponents**: Finite-size scaling fits (β=1/8, γ=7/4, ν=1)
- **Binder cumulant**: Tc crossing verification
- **Integration**: Full MCMC cycle with error bounds

All validation tests assert relative errors vs Onsager exact solution.

## Resultados

Reproducible con el CLI (L=32, Wolff, seed=1); compara MC con Onsager:

```bash
python -m statphys_viz.benchmark --L 32 --T-min 1.5 --T-max 3.5 --points 5
```

| T | e MC | e exacta | \|m\| MC | \|m\| exacta |
|---|------|----------|---------|-----------|
| 1.5 | -1.950 | -1.951 | 0.986 | 0.987 |
| 2.0 | -1.743 | -1.746 | 0.910 | 0.911 |
| 3.0 | -0.831 | -0.817 | 0.095 | 0 (finite-size) |

Los exponentes críticos se obtienen con `fit_critical_exponents` (ver tests).

## Live Demo

**[GitHub Pages Demo →](https://amurlaniakea.github.io/statphys-viz/)**

Zero-installation web demo with:
- Interactive temperature slider
- Live spin lattice animation
- Real-time observable plots
- Algorithm comparison (Metropolis vs Wolff)

*Runs entirely in the browser in plain JavaScript (`web/index.html`); numba does not run under Pyodide.*

## Licencia

AGPL-3.0-or-later

Copyright (C) 2026 Pedro Sordo Martínez — amurlaniakea@gmail.com

---

## Referencias

- Onsager, L. (1944). *Crystal Statistics. I. A Two-Dimensional Model with an Order-Disorder Transition*. Phys. Rev. 65, 117. [DOI](https://doi.org/10.1103/PhysRev.65.117)
- Wolff, U. (1989). *Collective Monte Carlo Updating for Spin Systems*. Phys. Rev. Lett. 62, 361. [DOI](https://doi.org/10.1103/PhysRevLett.62.361)
- Metropolis et al. (1953). *Equation of State Calculations by Fast Computing Machines*. J. Chem. Phys. 21, 1087. [DOI](https://doi.org/10.1063/1.1699114)
- Binder, K. (1981). *Finite size scaling analysis of Ising model block distribution functions*. Z. Phys. B 43, 119. [DOI](https://doi.org/10.1007/BF01293604)
- Ferdinand & Fisher (1969). *Bounded and inhomogeneous Ising models*. Phys. Rev. 185, 832. [DOI](https://doi.org/10.1103/PhysRev.185.832)

## Cómo citar

Si usas este software, cítalo con el DOI de Zenodo (10.5281/zenodo.24389055):

```bibtex
@software{sordo_statphys_viz_2026,
  author  = {Sordo Martínez, Pedro},
  title   = {statphys-viz: Interactive statistical physics visualizations for education},
  year    = {2026},
  version = {0.2.0},
  doi     = {10.5281/zenodo.24389055},
  url     = {https://github.com/amurlaniakea/statphys-viz},
  license = {AGPL-3.0-or-later}
}
```

GitHub también ofrece «Cite this repository» a partir de `CITATION.cff`.
