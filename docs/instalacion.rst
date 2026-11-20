Instalación
===========

El paquete todavía no está en PyPI. Se instala desde el código fuente:

.. code-block:: bash

   git clone https://github.com/amurlaniakea/statphys-viz.git
   cd statphys-viz
   python3 -m venv .venv && . .venv/bin/activate
   pip install -e ".[dev]"      # pytest, ruff, bandit, hypothesis

.. note::
   En distribuciones como Debian o Ubuntu recientes ``pip`` global está bloqueado (PEP 668),
   por eso se recomienda el entorno virtual.

**Requisitos:** Python ≥ 3.10, NumPy, Numba, SciPy, Matplotlib, ipywidgets e ipycanvas.

Comprobar la instalación
------------------------

.. code-block:: bash

   python -m pytest tests -q

Construir esta documentación
----------------------------

.. code-block:: bash

   pip install -e ".[docs]"
   make -C docs html            # resultado en docs/_build/html
