Historial de cambios
====================

0.2.2
-----

* Widget de Jupyter: ya no lanza ``ValueError`` al combinar Wolff con un campo externo; ahora cambia
  de algoritmo o pone ``h`` a 0 (con tests).
* README: tabla de validación, árbol, tests y mediciones corregidos para reflejar lo que el código hace.

0.2.1
-----

* Demo web: gráfica de e y \|m\| frente a T con la curva exacta de Onsager, DOI en el pie, secciones desplegables,
  fórmulas en LaTeX y tres idiomas (ES/EN/FR).
* Documentación completa publicada en ``/docs/`` junto a la demo.
* README: icono, enlaces a la demo y a la documentación, sección «Cómo citar» y ejemplo de inicio rápido
  fuera de la región crítica.

0.2.0
-----

* Corregidas las referencias exactas de Onsager (energía y calor específico) y añadidos tests cuantitativos.
* Metropolis secuencial (el ``prange`` anterior era una carrera de datos), semilla de numba reproducible y
  Wolff rechaza ``h != 0``.
* CLI funcional: ``python -m statphys_viz.benchmark`` y ``python -m statphys_viz.visualize``.
* Demo web didáctica en JavaScript puro, desplegada con GitHub Pages.
* DOI de Zenodo.

0.1.0
-----

* Primera versión: modelo de Ising 2D con Metropolis y Wolff, widgets de Jupyter y referencias de Onsager
  (con errores en energía y calor específico, corregidos en 0.2.0).
