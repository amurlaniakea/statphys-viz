Demo web
========

La demo está en https://amurlaniakea.github.io/statphys-viz/ y funciona en el navegador, sin instalar nada.
Está escrita en JavaScript puro (``web/index.html``) porque Numba no se ejecuta en Pyodide. Es la misma
física (Metropolis y Wolff sobre una red 64×64), pero **otro generador de números aleatorios**: los resultados
coinciden estadísticamente con el paquete Python, no dato a dato.

Qué incluye
-----------

* Red de espines animada, control de temperatura, selector de algoritmo y lectura en vivo de e/espín y \|m\|
  junto a los valores exactos de Onsager.
* Gráfica de e y \|m\| frente a T con la curva exacta, los puntos del barrido y el estado actual.
* Teoría en secciones desplegables (fórmulas en LaTeX) en español, inglés y francés.

Datos exportables (CSV)
-----------------------

Las líneas que empiezan por ``#`` son metadatos.

``ising_time_series.csv``
   ``step,T,e_per_spin,abs_m_per_spin``. Un paso = 5 clusters de Wolff o 1 barrido de Metropolis. Máximo 20 000 filas.

``ising_lattice.csv``
   Matriz 64×64 de espines ±1 en el estado actual.

``ising_sweep_vs_onsager.csv``
   ``T,e_MC,e_exact,absm_MC,absm_exact,C_MC,C_exact,chi_MC``. L = 32, Wolff, 300 clusters de termalización y 1500 de medida
   por punto, 21 temperaturas de 1.5 a 3.5, semilla fija (``mulberry32``, semilla 1): es reproducible.
   ``C_exact`` es ``inf`` en :math:`T_c`.

Ejemplo de análisis:

.. code-block:: python

   import pandas as pd

   df = pd.read_csv("ising_sweep_vs_onsager.csv", comment="#")
   df["err_e"] = (df.e_MC - df.e_exact).abs() / df.e_exact.abs()
   print(df[["T", "e_MC", "e_exact", "err_e"]])
