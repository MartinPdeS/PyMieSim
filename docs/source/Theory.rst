Theoretical background
======================

Lorenz-Mie Theory (LMT)
-----------------------

.. note::

  The Lorenz-Mie Theory or (LMT for short) is a framework that can be used to find
  exact solution of the scattered field considering a plane wave incident to a
  scatterer with a certain geometry.
  The solution is usually written in the form of an infinite summation which, of
  course, has to be truncated. PyMieSim is a library which does solve the
  equations in order to retrieve plenty of important informations.
  It is to be noted that in all the library, the angles :math:`\theta` and
  :math:`\phi` are defined as in a spherical coordinate system as shown in the
  following figure.

  .. image:: ../images/opticalSetup.png
    :width: 600


.. note::

  Here are few of the most important relations governing the PyMieSim library.

  .. math::
     &S_1=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n\pi_n+b_n\tau_n)

     &.

     &S_2=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n\tau_n+b_n\pi_n)

     .&

     &\text{Fields} = E_{\theta}(\phi,\theta) \vec{\theta} +   E_{\phi}(\phi,\theta) \vec{\phi}

     .&

     &\text{SPF} = \sqrt{ E_{\parallel}(\phi,\theta)^2 + E_{\perp}(\phi,\theta)^2 }

  **Stokes parameters:**

   .. math::

     &I = \big| E_x \big|^2 + \big| E_y \big|^2

     .&

     &Q = \big| E_x \big|^2 - \big| E_y \big|^2

     .&

     &U = 2 \mathcal{Re} \big\{ E_x E_y^* \big\}

     .&

     &V = 2 \mathcal{Im} \big\{ E_x E_y^* \big\}


Scattering properties
---------------------

.. note::
  There are many properties of the scatterer that might be useful to know such as:
  - scattering efficiency
  - extinction efficiency
  - absorption efficiency
  - back-scattering efficiency
  - ratio of front and back scattering
  - optical pressure efficiency
  - anisotropy factor g

  Those parameters can be computed using PyMieSim according to those equations.

  .. math::
    &Q_{sca} = \frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)(|a_n|^2+|b_n|^2)

    &Q_{ext} = \frac{2}{x^2} \sum_{n=1}^{n_{max}} \frac{(2n+1)}{\mathcal Re \{ a_n+b_n \}}

    &Q_{abs} = Q_{ext}-Q_{sca}

    &Q_{back} = \frac{1}{x^2} \Big| \sum\limits_{n=1}^{n_{max}} (2n+1)(-1)^n (a_n - b_n) \Big|^2

    &Q_{ratio} = \frac{Q_{back}}{Q_{sca}}

    &Q_{pr} = Q_{ext} - g * Q_{sca}

    &g = \frac{4}{Q_{sca} x^2}
            \Big[ \sum\limits_{n=1}^{n_{max}} \frac{n(n+2)}{n+1} \text{Re} \left\{ a_n a_{n+1}^* + b_n b_{n+1}^*\right\} +
            \sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)} \text{Re}\left\{ a_n b_n^* \right\} \Big]

    &A_s = \pi r^2

    &\sigma_{i} = Q_{i} A

    &\mu_{sca} = \sigma_{sca} C

    &\mu_{ext} = \sigma_{ext} C

    &\mu_{abs} = \sigma_{abs} C

  | With:
  |   C: the scatterer concentration in the sample.

-----


An and Bn coefficients:
-----------------------



From the An and Bn coefficients, we can retrieve many useful properties of
the scatterer and scattered far-fields. Those are complementary to the
Cn and Dn coefficient (for near-field properties) which we do no compute
with PyMieSim at the moment.
Depending on the scatterer geometry, all those coefficient may vary. Here we
have three example which are available with the PyMieSim library.

.. note::
  **Sphere**

  .. math::

      a_n = \frac{
      \mu_{sp} \Psi_n(x) \Psi_n^\prime(M x) - \mu M \Psi_n^\prime(x) \Psi_n(M x)}
      {\mu_{sp} \xi_n(x) \Psi_n^\prime(M x)- \mu M \xi_n^\prime (x) \Psi_n(M x)}

  .. math::

      b_n = \frac{
       \mu M \Psi_n(x) \Psi_n^\prime(M x) - \mu_{sp} \Psi_n^\prime(x) \Psi_n(M x)}
      {\mu M \xi_n(x) \Psi_n^\prime(M x) - \mu_{sp} \xi_n^\prime (x) \Psi_n(M x)}


  |   With:
  |     :math:`\psi_n = x \psi^{(1)}_n (x) = \sqrt{x \pi/2} J_{n+1/2} (x)`.
  |     :math:`M = k_{sp}/k` is the relative complex refractive index.
  |     :math:`x = \pi d / \lambda`.
  |     :math:`\lambda` is the wavelength in the surrounding medium.
  |     :ref:`References` [1] Eq(III.88-91).
  |
  |
  |  **important**: It is to be noted that PyMieSim assume :math:`\mu_{sp} = \mu` at the moment.
       It might change in a future update.


----

.. note::
  **Cylinder**

  .. math::

      a_n = \frac{ M J_n(M x) J_n^\prime (m x) - m J_n^\prime (M x) J_n(m x) }
      { m_t J_n(M x) H_n^\prime (m x) - m J_n^\prime (M x) H_n(m x) }

  .. math::

      b_n = \frac{ m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x) }
      { m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x) }


  |   With:
  |     :math:`M` is the refractive index of the scatterer.
  |     :math:`m` is the refractive index of the medium.
  |     :math:`H_n` is the Hankel function of first kind of order n.
  |     :ref:`References` [5] Eq(8.30-32).

----

.. note::
  **Core/Shell sphere**


  .. math::

      a_n = \frac{
      \psi_n \left[ \psi_n' (m_2 y) - A_n \chi_n' (m_s) \right] -
      m_2 \psi_n'(y) \left[ \psi_n (m_2 y) - A_n \chi_n (m_2 y) \right]}
      {\xi_n (y) \left[ \psi_n' (m_2 y) -A_n \chi_n' (m_2 y) \right] -
      m_2 \xi_n' (y) \left[ \psi_n(m_2 y) - A_n \chi_n (m_2 y) \right]}

  .. math::

      b_n = \frac{
        m_2 \psi_n(y) \left[ \psi_n' (m_2 y) - B_n \chi_n' (m_2 y) \right] -
        \psi_n' (y) \left[ \psi_n (m_2 y) - B_n \chi_n (m_2 y) \right]}
        {m_2 \xi_n(y) \left[ \psi_n' (m_2 y) - B_n \xi_n' (m_2 y) \right] -
        \xi_n' \left[ \psi_n (m_2 y) -A_n \chi_n (m_2 y)  \right]

        }

  |   With:

  .. math::
    A_n = \frac{ m_2 \psi_n (m_2 x) \psi_n' (m_1 x) - m1 \psi_n'(m_2 x) \psi_n(m_1x)}
      {m_2 \xi_n (m_2x) \psi_n' (m_1 x) - m_1 \xi_n' (m_2 x) \psi_n (m_1 x)}

  .. math::
    B_n = \frac{m_2 \psi_n (m_1 x) \psi_n' (m_2 x) - m_1 \psi_n (m_2 x) \psi_n' (m_1 x)}
      {m_2 \chi_n' (m_2 x) \psi_n (m_1 x) - m_1 \psi_n' (m_1 x) \chi_n (m_2 x)}

  |   and:

  .. math:: x = \frac{2 \pi R_{core}}{\lambda}, \:
            y = \frac{2 \pi R_{shell}}{\lambda}, \:
            m_1 = \frac{n_{core}}{n_{medium}}, \:
            m_2 = \frac{n_{shell}}{n_{medium}}.

  .. math::
    \chi_n (x) = -x\sqrt{\frac{\pi}{2x}} N_{n+1/2} (x) \leftarrow \text{Neumann}

  |     :ref:`References` [8] Eq(4-5).



Generalized Lorenz-Mie Theory (GLMT)
------------------------------------

.. note::
  **Coming soon**




-----

Coupling mechanism
-------------------

.. note::

  There are two main coupling mechanisms, **coherent coupling** and non-coherent coupling.
  For instance, photodiode collect light via an **non-coherent mechanism**, on the other part
  fiber optic LP mode collect light in a coherent way and as such they usually
  collect a lot less light but they add additional information on the sample studied.


  Mathematically they are defined as follows:

  .. math::
      C_{coh.} &= \Big| \iint_{\Omega}  \Phi_{det} \, . \, \Psi_{scat}^* \,  d \Omega \Big|^2

      C_{Noncoh.} &=  \iint_{\Omega}  \Big| \Phi_{det} \Big|^2 \,.\, \Big| \Psi_{scat} \Big|^2 \,  d \Omega



  It is to be noted that the **coherent coupling** definition is derived from the coupled mode theory
  which remains true as long as the parallax approximation is also true.
  Also, this coupling are what we would call **centered coupling**. It means that the
  scatterer is perfectly centered with the detector. As much as it doesn't affect
  so much the **non-coherent coupling** coupling, it can largely affect **coherent coupling**.

  In order to take account of the effect of transversal offset of the scatterer we define
  the footprint of the scatterer.


  .. math::
    \eta_{l,m}(\delta_x, \delta_y) = \Big| \mathcal{F}^{-1} \big\{ \Phi_{det} \, . \, \Psi_{scat} \big\}  \Big|^2

  Thus we can compute the **mean coupling** as the mean value of :math:`\eta_{l,m}`

  .. math::
    \widetilde{\eta}_{l,m} = \big< \eta_{l,m}(\delta_x, \delta_y) \big>
