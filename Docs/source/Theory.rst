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

              :math:`{\displaystyle Q_{sca}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)(|a_n|^2+|b_n|^2)}`

              :math:`{\displaystyle Q_{ext} = \frac{2}{x^2} \sum_{n=1}^{n_{max}} \frac{(2n+1)}{\mathcal Re \{ a_n+b_n \}}}`

              :math:`Q_{abs}=Q_{ext}-Q_{sca}`


              :math:`{\displaystyle a_n = \frac{\mu_{sp} \Psi_n(\alpha) \Psi_n^\prime(\beta) - \mu M \Psi_n^\prime(\alpha) \Psi_n(\beta)}
              {\mu_{sp} \xi_n(\alpha) \Psi_n^\prime(\beta)- \mu M \xi_n^\prime (\alpha) \Psi_n(\beta)} }`

              :math:`{\displaystyle b_n = \frac{\mu M \Psi_n(\alpha) \Psi_n^\prime(\beta) -\mu_{sp} \Psi_n^\prime(\alpha) \Psi_n(\beta)}
              {\mu M \xi_n(\alpha) \Psi_n^\prime(\beta)-\mu_{sp} \xi_n^\prime (\alpha) \Psi_n(\beta)} }`

              :math:`{\displaystyle c_n = \frac{\mu_{sp} M \big[ \xi_n(\alpha) \Psi_n^\prime(\alpha) - \xi_n^\prime(\alpha) \Psi_n(\alpha) \big]}
              {\mu_{sp} \xi_n(\alpha) \Psi_n^\prime(\beta)- \mu M \xi_n^\prime (\alpha) \Psi_n(\beta)} }`

              :math:`{\displaystyle d_n = \frac{ \mu M^2 \big[ \xi_n(\alpha) \Psi_n^\prime(\alpha) -\xi_n^\prime(\alpha) \Psi_n(\alpha) \big]}
              {\mu M \xi_n(\alpha) \Psi_n^\prime(\beta)-\mu_{sp} M \xi_n^\prime (\alpha) \Psi_n(\beta)} }`



Generalized Lorenz-Mie Theory (GLMT)
------------------------------------
