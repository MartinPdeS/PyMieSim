Theoretical background
======================




Lorenz-Mie Theory (LMT)
-----------------------

Here some usefull relations:

Calculates S\ :sub:`1` and S\ :sub:`2` at μ=cos(θ), where θ is the scattering angle.

 S\ :sub:`1` and S\ :sub:`2` are calculated by:

             :math:`{\displaystyle S_1=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n\pi_n+b_n\tau_n)}`

             :math:`{\displaystyle S_2=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n\tau_n+b_n\pi_n)}`


 Computes Mie efficencies *Q* and asymmetry parameter *g* of a single, homogeneous particle:

.. math::
  Q_{sca} &= \frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)(|a_n|^2+|b_n|^2)

  Q_{ext} &= \frac{2}{x^2} \sum_{n=1}^{n_{max}} \frac{(2n+1)}{\mathcal Re \{ a_n+b_n \}}

  Q_{abs} &= Q_{ext}-Q_{sca}


**For spherical scatterer**

.. math::

  a_n &= \frac{\mu_{sp} \Psi_n(\alpha) \Psi_n^\prime(\beta) - \mu M \Psi_n^\prime(\alpha) \Psi_n(\beta)}
  {\mu_{sp} \xi_n(\alpha) \Psi_n^\prime(\beta)- \mu M \xi_n^\prime (\alpha) \Psi_n(\beta)}

  b_n &= \frac{\mu M \Psi_n(\alpha) \Psi_n^\prime(\beta) -\mu_{sp} \Psi_n^\prime(\alpha) \Psi_n(\beta)}
  {\mu M \xi_n(\alpha) \Psi_n^\prime(\beta)-\mu_{sp} \xi_n^\prime (\alpha) \Psi_n(\beta)}

  c_n &= \frac{\mu_{sp} M \big[ \xi_n(\alpha) \Psi_n^\prime(\alpha) - \xi_n^\prime(\alpha) \Psi_n(\alpha) \big]}
  {\mu_{sp} \xi_n(\alpha) \Psi_n^\prime(\beta)- \mu M \xi_n^\prime (\alpha) \Psi_n(\beta)}

  d_n &= \frac{ \mu M^2 \big[ \xi_n(\alpha) \Psi_n^\prime(\alpha) -\xi_n^\prime(\alpha) \Psi_n(\alpha) \big]}
  {\mu M \xi_n(\alpha) \Psi_n^\prime(\beta)-\mu_{sp} M \xi_n^\prime (\alpha) \Psi_n(\beta)}

| With:
|   :math:`M = \frac{k_{sp}}{k}`


**For cylindrical scatterer**

.. math::

    a_n & = \frac{ m_t J_n(m_t x) J_n^\prime (m x) - m J_n^\prime (m_t x) J_n(m x) }
    { m_t J_n(m_t x) H_n^\prime (m x) - m J_n^\prime (m_t x) H_n(m x) }

    b_n & = \frac{ m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x) }
    { m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x) }

| With:
|   :math:`m` being the refractive index of the medium and
|   :math:`m_t` being the refractive index of the index.

-----

Generalized Lorenz-Mie Theory (GLMT)
------------------------------------


**Coming soon**




-----

Coupling mechanism
-------------------


There is two main coupling mechanism, **coherent coupling** and non-coherent coupling.
For instance photodiode collect light via an **non-coherent mechanism**, on the other part
fiber optic LP mode collect light in a coherent way and as such they usually
collect a lot less light but they add additional information on the sample studied.


Mathematically they are defined as follows:

.. math::
    C_{coh.} &= \Big| \iint_{\\Omega}  \Phi_{det} \, . \, \Psi_{scat}^* \,  d \Omega \Big|^2

    C_{Noncoh.} &=  \iint_{\\Omega}  \Big| \Phi_{det} \Big|^2 \,.\, \Big| \Psi_{scat} \Big|^2 \,  d \Omega



It is to be noted that the **coherent coupling** definition is derived from the coupled mode theory
which remains true as long as the parallax approximation is also true.
Also this coupling are what we would call **centered coupling**. It means that the
scatterer is perfectly centered with the detector. As much as it doesn't affect
so much the **non-coherent coupling** coupling, it can largely affect **coherent coupling**.

In order to take account of the effect of transversal offset of the scatterer we define
the footprint of the scatterer.


.. math::
  \eta_{l,m}(\delta_x, \delta_y) = \Big| \mathcal{F}^{-1} \big\{ \Phi_{det} \, . \, \Psi_{scat} \big\}  \Big|^2

Thus we can compute the **mean coupling** as the mean value of :math:`\eta_{l,m}`

.. math::
  \widetilde{\eta}_{l,m} = \big< \eta_{l,m}(\delta_x, \delta_y) \big>
