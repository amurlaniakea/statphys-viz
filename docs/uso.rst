Uso
===

Inicio rápido en Python
-----------------------

.. code-block:: python

   from statphys_viz import IsingConfig, IsingModel, OnsagerBenchmark

   config = IsingConfig(L=64, T=2.0, J=1.0, algorithm="wolff", seed=1)
   model = IsingModel(config)
   model.thermalize(500)
   results = model.measure(2000)

   print(results["energy_mean"], results["magnetization_mean"])
   print(results["specific_heat"], results["susceptibility"])

   errors = OnsagerBenchmark(L=64).validate(results, T=2.0)
   print(f"Error relativo máximo: {errors['max_rel_error']:.2%}")   # ≈ 0.6 %

``measure`` devuelve un diccionario con las series (``energy``, ``magnetization``) y los promedios
por espín (``energy_mean``, ``magnetization_mean`` = ⟨\|m\|⟩, ``specific_heat``, ``susceptibility``,
``energy_var``, ``magnetization_var``).

.. warning::
   **Cerca de Tc el error es grande y es normal.** A T = 2.269 con L = 64 el error máximo ronda el 60 %:
   la solución exacta es para L → ∞ y en la criticalidad el tamaño finito domina (el calor específico
   diverge y la magnetización exacta cae a 0). Para validar el código elige T lejos de Tc, y para
   estudiar la región crítica usa escalado de tamaño finito (:class:`statphys_viz.benchmark.OnsagerBenchmark`).

Reglas del algoritmo
--------------------

* ``algorithm="wolff"`` (por defecto): voltea clusters, sin ralentización crítica. **Solo es válido con h = 0**;
  con ``h != 0`` lanza ``ValueError``.
* ``algorithm="metropolis"``: voltea espines sueltos; admite campo externo, pero cerca de Tc converge muy despacio.
* En el widget de Jupyter, mover ``h`` con Wolff seleccionado cambia a Metropolis, y elegir Wolff pone ``h`` a 0.
* ``seed`` siembra la red inicial **y** el generador de numba. Es un estado global del proceso: dos modelos
  creados en el mismo proceso comparten ese generador.

Widgets de Jupyter
------------------

.. code-block:: python

   from statphys_viz.widgets import IsingWidget, create_temperature_scan_widget

   IsingWidget(L=64, T=2.269, algorithm="wolff").display()
   create_temperature_scan_widget(L=32, T_min=1.5, T_max=3.5, n_points=25)

Línea de comandos
-----------------

.. code-block:: bash

   # Barrido de temperatura (Wolff) comparado con Onsager
   python -m statphys_viz.benchmark --L 32 --T-min 1.5 --T-max 3.5 --points 20

   # Instantánea de la red cerca de Tc
   python -m statphys_viz.visualize --output spins.png
