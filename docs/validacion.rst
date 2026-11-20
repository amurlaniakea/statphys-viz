Validación
==========

Qué se comprueba
----------------

La suite (``python -m pytest tests -q``) incluye comparaciones **cuantitativas** con la solución exacta:

* :math:`u(T_c) = -\sqrt{2}` y :math:`u(2.0) \approx -1.7456`.
* El calor específico exacto coincide con la derivada numérica :math:`du/dT` en cinco temperaturas.
* La energía de Monte Carlo (Wolff y Metropolis) concuerda con la exacta dentro de 0.02 en T = 2.0.
* Wolff rechaza ``h != 0`` y la semilla hace reproducible Metropolis.

Resultados de referencia
------------------------

Generados con ``python -m statphys_viz.benchmark --L 32 --T-min 1.5 --T-max 3.5 --points 5 --samples 1500``
(Wolff, semilla 1):

.. list-table::
   :header-rows: 1

   * - T
     - e MC
     - e exacta
     - \|m\| MC
     - \|m\| exacta
     - C MC
     - C exacta
   * - 1.5
     - -1.9503
     - -1.9511
     - 0.9863
     - 0.9865
     - 0.199
     - 0.197
   * - 2.0
     - -1.7431
     - -1.7456
     - 0.9095
     - 0.9113
     - 0.740
     - 0.725
   * - 2.5
     - -1.0990
     - -1.1061
     - 0.1935
     - 0
     - 0.984
     - 0.862
   * - 3.0
     - -0.8308
     - -0.8173
     - 0.0948
     - 0
     - 0.339
     - 0.401
   * - 3.5
     - -0.6410
     - -0.6601
     - 0.0659
     - 0
     - 0.206
     - 0.248

Las diferencias por encima de :math:`T_c` en :math:`|m|` (y la desviación de :math:`C` cerca de :math:`T_c`) son
efectos de tamaño finito (L = 32), no errores del código.

Historial
---------

La versión 0.1.0 tenía **incorrectas** las referencias exactas de energía y calor específico y los tests solo
comprobaban signos. La versión 0.2.0 las corrige y añade las comprobaciones anteriores.
