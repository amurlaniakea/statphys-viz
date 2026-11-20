Teoría
======

Modelo
------

Espines :math:`s_i = \pm 1` en una red cuadrada :math:`L \times L` con bordes periódicos:

.. math::

   H = -J \sum_{\langle ij \rangle} s_i s_j - h \sum_i s_i, \qquad P(\mathbf{s}) \propto e^{-\beta H}, \qquad \beta = \frac{1}{T}

con :math:`J = k_B = 1`.

Algoritmos
----------

**Metropolis** (1953). Se voltea un espín con probabilidad

.. math::

   p = \min\left(1,\, e^{-\beta \Delta E}\right), \qquad \Delta E = 2 s_i \Big( J \sum_{j \in \mathrm{nn}(i)} s_j + h \Big)

**Wolff** (1989). Se crece un cluster de espines iguales añadiendo cada vecino con probabilidad

.. math::

   p_{\mathrm{add}} = 1 - e^{-2\beta J}

y se voltea siempre entero. Es un algoritmo sin rechazo que evita la ralentización crítica, pero solo es correcto con :math:`h = 0`.

Observables (por espín, con :math:`N = L^2`)
--------------------------------------------

.. math::

   e = \frac{\langle E \rangle}{N}, \quad m = \frac{\langle |M| \rangle}{N}, \quad
   C = N \beta^2 \left( \langle e^2 \rangle - \langle e \rangle^2 \right), \quad
   \chi = N \beta \left( \langle m^2 \rangle - \langle |m| \rangle^2 \right)

El cumulante de Binder :math:`U_4 = 1 - \langle m^4 \rangle / (3 \langle m^2 \rangle^2)` se cruza en :math:`T_c`
para distintos :math:`L`.

Solución exacta de Onsager (L → ∞)
----------------------------------

.. math::

   T_c = \frac{2J}{\ln(1 + \sqrt{2})} \approx 2.26919

.. math::

   m = \left( 1 - \sinh^{-4}(2\beta J) \right)^{1/8} \quad (T < T_c), \qquad m = 0 \quad (T > T_c)

.. math::

   u = -J \coth(2\beta J) \left[ 1 + \frac{2}{\pi} \left( 2 \tanh^2(2\beta J) - 1 \right) K(\kappa) \right],
   \qquad \kappa = \frac{2 \sinh(2\beta J)}{\cosh^2(2\beta J)}

con :math:`K` la integral elíptica completa de primera especie (de módulo :math:`\kappa`; en SciPy,
``scipy.special.ellipk`` toma el parámetro :math:`\kappa^2`). En :math:`T_c`, :math:`u = -\sqrt{2}`. El calor
específico es :math:`C = du/dT`, con divergencia logarítmica en :math:`T_c`.

Exponentes críticos y tamaño finito
-----------------------------------

.. math::

   m \sim |T - T_c|^{1/8}, \qquad \chi \sim |T - T_c|^{-7/4}, \qquad \xi \sim |T - T_c|^{-1}, \qquad C \sim \ln|T - T_c|

Exactamente en :math:`T_c` y con :math:`L` finito:

.. math::

   m \sim L^{-1/8}, \qquad \chi \sim L^{7/4}

Por eso con :math:`L` finito los picos se redondean y :math:`|m|` no cae a 0 por encima de :math:`T_c`.

Referencias
-----------

* L\. Onsager (1944), *Crystal Statistics. I.*, Phys. Rev. 65, 117. https://doi.org/10.1103/PhysRev.65.117
* N\. Metropolis et al. (1953), J. Chem. Phys. 21, 1087. https://doi.org/10.1063/1.1699114
* U\. Wolff (1989), Phys. Rev. Lett. 62, 361. https://doi.org/10.1103/PhysRevLett.62.361
* K\. Binder (1981), Z. Phys. B 43, 119. https://doi.org/10.1007/BF01293604
