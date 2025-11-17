InGaAsPElectroOpticalModel
==========================

An ``ElectroOpticalModel`` is a ``PhotonicPolygon``-specific object that will calculate the change in optical permitivity, :math:`\Delta \epsilon`, locally, due to an applied electric field to your structure. These models are not meant to be used as a stand-alone package, but to be fed into a ``PhotonicPolygon``, which in turn `ElectroOpticalSimulator` will make use of to calculate the full electro-optical response. 

InGaAsP-based modulators are one of the most mature technologies in integrated photonics, which translates into well established models that predict the electro-optic response of these alloys. All of the models below are specific to alloys lattice matched to InP, where the :math:`x` can be predicted as:

.. math::

    x = \frac{y}{2.2020 - 0.659y}

In order to use many of the models, we first need to be able to predict some physical properties for all alloys. We will go over them one my one.

Effective mass
--------------

The effective masses are calculated according to :cite:t:`fiedler_optical_1987`:

.. math::

    \begin{align}
        m_e &= \left( 0.07 - 0.0308y \right)m_0 \\
        m_{hh} &= \left( 0.6 - 0.218y + 0.07y^2 \right)m_0 \\
        m_{lh} &= \left( 0.12 - 0.078y + 0.002y^2 \right)m_0
    \end{align}

where :math:`m_0` is the vacuum electron mass.

Carrier effective density of states
-----------------------------------

The effective density of states are calculated from :cite:t:`bennett_carrier-induced_1990`. For electrons we have:

.. math::

    N_C = 2 \left( \frac{m_e k_b T}{2\pi \hbar^2}\right)^{3/2}

and for holes:

.. math::

    N_v = 2 \left( \frac{m_h k_b T}{2\pi \hbar^2}\right)^{3/2}

with:

.. math::
    m_h = \left(m_{hh}^{1.5} + m_{lh}^{1.5}\right)^{2/3}

where :math:`k_b` is the boltzman constant, :math:`T` is the temperature and :math:`\hbar` is the planck constant.

Spin orbit splitting
--------------------

The spin orbit splitting is taken from :cite:t:`fiedler_optical_1987`:

.. math::
    \Delta_{so} = 0.119 + 0.30y - 0.107y^2 \hspace{1cm} [eV]

Bandgap narrowing
-----------------

The bandgap narrowing is a well known effect in III-V semiconductors, where it is observed that the bandgap of heavily doped materials shrinks. Here we follow the work of :cite:t:`jain_bandgap_1990` where the band gap narrowing (BGN) is estimated via:

.. math::
    \Delta E_g^{BGN} = A \times N^{1/3} + B \times N^{1/4} + C \times N^{1/2}

where :math:`N` is the carrier concentration and the constants :math:`A`, :math:`B` and :math:`C` have been empirically found for GaAs, GaP, InP and InAs, allowing for interpolation for InGaAsP alloys via :cite:t:`fiedler_optical_1987`:

.. math::
    Q(x,y) = xyQ_{GaAs} + x(1-y)Q_{GaP} + y(1-x)Q_{InAs} + (1-x)(1-y)Q_{InP}

The constants are:

.. list-table:: P-doped material
    :widths: 12 10 10 10
    :align: center
    :header-rows: 1

    * - Material
      - A :math:`\times 10^{-9}`
      - B :math:`\times 10^{-7}`
      - C :math:`\times 10^{-12}`
    * - GaAs
      - 9.83
      - 3.90
      - 3.90
    * - InAs
      - 8.34
      - 2.91
      - 4.53
    * - InP
      - 10.3
      - 4.43
      - 3.38
    * - GaP
      - 12.7
      - 5.85
      - 3.90


.. list-table:: N-doped material
    :widths: 12 10 10 10
    :align: center
    :header-rows: 1

    * - Material
      - A :math:`\times 10^{-9}`
      - B :math:`\times 10^{-7}`
      - C :math:`\times 10^{-12}`
    * - GaAs
      - 16.5
      - 2.39
      - 91.4
    * - InAs
      - 14.0
      - 1.97
      - 57.9
    * - InP
      - 17.2
      - 2.62
      - 98.4
    * - GaP
      - 10.7
      - 3.45
      - 9.97

Charge carrier mobility
-----------------------

Despite most commercial software for the simulation of charge transport in semiconductors have a good material database, III-V (in particular InGaAsP) materials seldom have a correct value for mobility, and in particular doping dependent mobility values. For that reason we will follow the empirical model of :cite:t:`sotoodeh_empirical_2000` for low field mobility:

.. math::

  \mu (N, T) = \mu_{min} + \frac{\mu_{max}(300K)\left(\frac{300K}{T}\right)^{\theta_1}-\mu_{min}}{1+\left(\frac{N}{N_{ref}(300K)\left(\frac{T}{300K}\right)^{\theta_2}}\right)^\lambda}

Where the values of :math:`\lambda_n`, :math:`\theta_{n2}`, :math:`\log_10{N_{n,ref}(300K)}`, :math:`\lambda_p`, :math:`\mu_{p, min}`, :math:`\theta_{p2}` and :math:`\log_10{N_{p,ref}(300K)}` for :math:`In_{1-x}Ga_{x}As_{y}P_{1-y}` are found via:

.. math::
  Q(In_{1-x}Ga_{x}As_{y}P_{1-y}) = yQ(In_{1-x}Ga_{x}As) + (1-y)Q(In_{1-x}Ga_{x}P)

The :math:`\mu_{n,max}(300K)`, :math:`\mu_{n,min}(300K)`, :math:`\mu_{p,max}(300K)`, :math:`\theta_{n1}`, :math:`\theta_{p1}` are found via:

.. math::
  Q(In_{1-x}Ga_{x}As_{y}P_{1-y}) = \frac{yQ(In_{1-x}Ga_{x}As) + (1-y)Q(In_{1-x}Ga_{x}P)}{1+my(1-y)}

with:

.. list-table::
    :widths: 8 8 8 8 8 8
    :align: center
    :header-rows: 1

    * - 
      - :math:`\mu_{n,max}(300K)`
      - :math:`\mu_{n,min}(300K)`
      - :math:`\mu_{p,max}(300K)`
      - :math:`\theta_{n1}`
      - :math:`\theta_{p1}`
    * - m
      - 6
      - 6
      - 6
      - 1
      - 1

Therefore, we just need to find the values for InGaAs and InGaP. To do so, we do a quadratic interpolation between each of the parameters, unless there are only two data points, in which case we do a linear interpolation. The values used are layed down in the tables below:


.. list-table:: InGaP parameter list
    :widths: 8 8 8 8
    :align: center
    :header-rows: 1

    * - x
      - 0
      - 0.51
      - 1
    * - :math:`\mu_{n, max}`
      - 5200
      - 4300
      - 152
    * - :math:`\mu_{n, min}`
      - 400
      - 400
      - 10
    * - :math:`N_{n, ref}`
      - log10(3e17)
      - log10(2e16)
      - log10(4.4e18)
    * - :math:`\lambda_{n}`
      - 0.47
      - 0.70
      - 0.80
    * - :math:`\theta_{n,1}`
      - 2.0
      - 1.66
      - 1.60
    * - :math:`\theta_{n, 2}`
      - 3.25
      - -
      - 0.71
    * - :math:`\mu_{p, max}`
      - 170
      - 150
      - 147
    * - :math:`\mu_{p, min}`
      - 10
      - 15
      - 10
    * - :math:`N_{p, ref}`
      - log10(4.87e17)
      - log10(1.5e17)
      - log10(1.0e18)
    * - :math:`\lambda_{p}`
      - 0.62
      - 0.80
      - 0.85
    * - :math:`\theta_{p,1}`
      - 2.0
      - 2.0
      - 1.98
    * - :math:`\theta_{p, 2}`
      - 3.0
      - -
      - 0.0


.. list-table:: InGaAs parameter list
    :widths: 8 8 8 8
    :align: center
    :header-rows: 1

    * - x
      - 0
      - 0.47
      - 1
    * - :math:`\mu_{n, max}`
      - 34000
      - 14000
      - 9400
    * - :math:`\mu_{n, min}`
      - 1000
      - 300
      - 500
    * - :math:`N_{n, ref}`
      - log10(1.1e18)
      - log10(1.3e17)
      - log10(6.0e16)
    * - :math:`\lambda_{n}`
      - 0.32
      - 0.48
      - 0.394
    * - :math:`\theta_{n,1}`
      - 1.57
      - 1.59
      - 2.1
    * - :math:`\theta_{n, 2}`
      - 3.0
      - 3.68
      - 3.0
    * - :math:`\mu_{p, max}`
      - 530
      - 320
      - 491.5
    * - :math:`\mu_{p, min}`
      - 20
      - 10
      - 20
    * - :math:`N_{p, ref}`
      - log10(1.1e17)
      - log10(4.9e17)
      - log10(1.48e17)
    * - :math:`\lambda_{p}`
      - 0.46
      - 0.403
      - 0.38
    * - :math:`\theta_{p,1}`
      - 2.3
      - 1.59
      - 2.2
    * - :math:`\theta_{p, 2}`
      - 3.0
      - 3.0
      - 3.0


Refractive index
----------------

The optical refractive index above the bandgap absorption edge is calculated via the modified single oscillator model through :cite:t:`fiedler_optical_1987`:

.. math::
  n^2 = 1+\frac{E_d}{E_0} + E_d \frac{E_d (\hbar \omega)^2}{E_0^3} + \frac{E_d}{2E_0^3(E_0^2 - E_g^2)} (\hbar \omega)^4ln\left[\frac{2E_0^2 - E_g^2 - (\hbar \omega)^2}{E_g^2 - (\hbar \omega)^2}\right]

where :math:`E_0` and :math:`E_d` are given by:

.. math::
  E_0 = 3.391 - 1.652y + 0.863y^2 - 0.123y^3

.. math::
  E_d = 28.91 - 9.278y + 5.626y^2

Whereas if any information above the bandgap is needed, we resort to the following formula from :cite:t:`seifert_revised_2016`:

.. math::
  \begin{aligned}
    n^* - 1 &\approx \frac{a}{b - (E + i\Gamma)^*} + \frac{A\sqrt{R}}{(E + i\Gamma)^*} 
    \Bigg\{ \ln \frac{E_z^*}{E_z^* - (E + i\Gamma)^*} + \pi \Bigg[ 2 \cot \left( \pi \sqrt{\frac{R}{E_z^*}} \right) \\
    &\quad - \cot \left( \pi \sqrt{\frac{R}{E_z^* - (E + i\Gamma)^*}} \right) - \cot \left( \pi \sqrt{\frac{R}{E_z^* + (E + i\Gamma)}} \right) \Bigg] \Bigg\}
  \end{aligned}

where:

.. math::
  \begin{aligned}
    R &= -0.00115 + 0.0191E_g \\
    \Gamma &= -0.000691 + 0.00433 E_g \\
    A &= -0.0453 + 2.1103 E_g \\
    a &= 72.32 + 12.78 E_g \\
    b &= 4.84 + 4.66 E_g \\
    c &= -0.015 + 0.02 E_g \\
    d &= -0.178 + 1.042 E_g
  \end{aligned}

Charge carrier electro-optic effects
------------------------------------

Now that we have defined all the necessary physical properties needed to do our calculations, we can finally dive deeper into the electro-optic effects that take place. We will start with the effects that are governed by the electrons and holes. We will not consider in extreme details the physics behind each effect, for that we reccomend you follow the cited references. Instead, we will focus on the models that are employed and their possible shortcomings.

Band filling effect
~~~~~~~~~~~~~~~~~~~

When doping a semiconductor with additional donors, the fermi level gets closer to the conduction band, which ultimately can cause the occupation of energy levels above the minimum of the conduction band to be occupied. This means that the excitation of electrons will only occur at higher energy levels than the bandgap energy. This change in the absorption spectrum will translate into a change in the refractive index via the Kramers-Kronig relations. Here we follow the work of :cite:t:`bennett_carrier-induced_1990`, with the difference that the quasi-fermi levels will be taken from numerical calculations of the charge transport simulator.

The absorption due to the bandfilling effect is:

.. warning::
    In the models below, we do not consider the bandgap narrowing as a separate effect. Instead we consider it in-built into every model that is dependent on the bandgap value by considering :math:`E_g \to E_g - \Delta E_{BGN}`

.. math::
    \begin{aligned}
        \alpha(N,P,E) &= \frac{C_{hh}}{E} \sqrt{E - E_g - \Delta E_{BGN}} 
        \left[ f_v(E_{ah}) - f_c(E_{bh}) - 1 \right] \\
        &\quad + \frac{C_{lh}}{E} \sqrt{E - E_g - \Delta E_{BGN}} 
        \left[ f_v(E_{al}) - f_c(E_{bl}) - 1 \right]
    \end{aligned}

where :math:`f_v` and :math:`f_c` are the fermi-dirac distributions considering the fermi level as the quasi-fermi level for holes and electrons, respectively. :math:`E` is the energy of the incoming photon. The change in absorption is then given by:

.. math::
    \Delta \alpha(N,P,E) = \alpha(N,P,E) - \alpha_0

where:

.. math::
    \alpha_0(N,P,E) = \frac{C_{hh}}{E}\sqrt{E - E_g - \Delta E_{BGN}}  + \frac{C_{lh}}{E}\sqrt{E - E_g - \Delta E_{BGN}} 

The constants :math:`C_{hh}` and :math:`C_{lh}` are adapted from a constant :math:`C = 4.4e12 cm^{-1} s^{-0.5}` in :cite:t:`bennett_carrier-induced_1990`. 

.. warning::
    To determine the constants :math:`C_{hh}` and :math:`C_{lh}` for arbitrary concentrations of InGaAsP, we have adopted a new approach. We have noticed that in the work of :cite:t:`bennett_carrier-induced_1990`, the quaternary interpolation is done via the charge carrier masses alone. However, the parabolic absorption formula given by

    .. math::
        \alpha_0(N,P,E) = \frac{C}{E}\sqrt{E - E_g}

    tells us that the constant :math:`C` can be written as :cite:t:`moss_semiconductor_1973`

    .. math::
        C = \frac{2\pi e^2 (2m_r)^{3/2} |p_{m0}|^2}{3m_0^2 n_0 \epsilon_0 c h^3 \nu}

    we see that the constant is dependent on the refractive index as well, which will have an impact and is not accounted by Bennet. At the same time, :cite:t:`vinchant_inpgainasp_1992` states that a scaling factor is applied to the absorption spectrum so as to fit measurements at :math:`E_g + 0.2eV`, however, such scaling factor is not disclosed. For these reasons we have decided to employ a new model.

The reduced mass of the electron-heavy/light hole pairs is:

.. math::
    \begin{aligned}
        m_{r, ehh} &= \left(\frac{1}{m_e} + \frac{1}{m_hh}\right) \\
        m_{r, elh} &= \left(\frac{1}{m_e} + \frac{1}{m_lh}\right)
    \end{aligned}

The constants :math:`C_{hh}` and :math:`C_{lh}` are now calculated as:

.. math::
    C_{hh} = C \frac{m_{r,hh, InP}^{3/2}}{m_{r,hh, InP}^{3/2} + m_{r,lh, InP}^{3/2}} \left(\frac{m_{r, hh}}{m_{r,hh,InP}}\right)^{3/2} \frac{n_{0,InP}}{n_0}

.. math::
    C_{lh} = C \frac{m_{r,lh, InP}^{3/2}}{m_{r,hh, InP}^{3/2} + m_{r,lh, InP}^{3/2}} \left(\frac{m_{r, lh}}{m_{r,lh,InP}}\right)^{3/2} \frac{n_{0,InP}}{n_0}

The change in refractive index can now be calculated via the Kramers-Kronig integral. We have found, in accordance to the literature :cite:`vinchant_inpgainasp_1992`, that for all relevant alloys, the dependence with carrier concentration is linear. Therefore, we now emply the following slopes from :math:`\Delta n = m_n N + m_p P`:

.. list-table:: Our model with bandgap narrowing
    :widths: 8 8 8
    :align: center
    :header-rows: 1

    * - x
      - :math:`m_n` 
      - :math:`m_p`
    * - 0.0
      - -4.935e-21
      - -1.339e-21
    * - 0.1
      - -6.053e-21
      - -1.582e-21
    * - 0.2
      - -7.589e-21
      - -1.884e-21
    * - 0.3
      - -9.738e-21
      - -2.266e-21
    * - 0.4
      - -12.96e-21
      - -2.758e-21
    * - 0.53
      - -20.48e-21
      - -3.6560e-21
    * - 0.6
      - -28.27e-21
      - -4.340e-21
    * - 0.7
      - -58.43e-21
      - -5.789e-21
    
.. list-table:: Our model without bandgap narrowing
    :widths: 8 8 8
    :align: center
    :header-rows: 1

    * - x
      - :math:`m_n` 
      - :math:`m_p`
    * - 0.0
      - -3.520e-21
      - -1.258e-21
    * - 0.1
      - -4.284e-21
      - -1.479e-21
    * - 0.2
      - -5.276e-21
      - -1.752e-21
    * - 0.3
      - -6.597e-21
      - -2.092e-21
    * - 0.4
      - -8.410e-21
      - -2.523e-21
    * - 0.53
      - -12.01e-21
      - -3.289e-21
    * - 0.6
      - -21.67e-21
      - -4.982e-21
    * - 0.7
      - -35.52e-21
      - -6.936e-21

.. list-table:: Vinchant model
   :align: center
   :header-rows: 1

   * - x
     - :math:`m_n` 
     - :math:`m_p`
   * - -0.003
     - -5.625e-21
     - 0
   * - 0.047
     - -6.192e-21
     - 0
   * - 0.089
     - -6.637e-21
     - 0
   * - 0.127
     - -7.122e-21
     - 0
   * - 0.166
     - -7.446e-21
     - 0
   * - 0.2
     - -7.932e-21
     - 0
   * - 0.235
     - -8.255e-21
     - 0
   * - 0.27
     - -8.984e-21
     - 0
   * - 0.307
     - -9.469e-21
     - 0
   * - 0.352
     - -10.85e-21
     - 0
   * - 0.398
     - -12.06e-21
     - 0
   * - 0.44
     - -13.35e-21
     - 0
   * - 0.473
     - -14.65e-21
     - 0
   * - 0.515
     - -16.43e-21
     - 0
   * - 0.557
     - -18.86e-21
     - 0
   * - 0.603
     - -21.93e-21
     - 0
   * - 0.636
     - -25.33e-21
     - 0
   * - 0.669
     - -28.25e-21
     - 0
   * - 0.692
     - -31.48e-21
     - 0
   * - 0.72
     - -36.34e-21
     - 0
   * - 0.742
     - -40.87e-21
     - 0
   * - 0.762
     - -47.35e-21
     - 0
   * - 0.778
     - -52.85e-21
     - 0
   * - 0.793
     - -58.35e-21
     - 0

Plasma effect
~~~~~~~~~~~~~

In this model we only consider the plasma effect coming from n-dopants. The reason for this is that, in p-dopants, the inter-valence absorption mechanism is so much stronger that scattering effects in holes are negligible. We follow the model from :cite:t:`walukiewicz_electron_1980`, which is in accordance with :cite:t:`dumke_intra-_1970`. The change in absorption is given by:

.. math::
    \Delta \alpha(N, \lambda) = A_{\text{imp}}(N)\left(\frac{\lambda}{\lambda_0}\right)^{3.5} 
    + A_{\text{op}}(N)\left(\frac{\lambda}{\lambda_0}\right)^{2.5} 
    + A_{\text{ac}}(N)\left(\frac{\lambda}{\lambda_0}\right)^{1.5}


where :math:`\lambda_0 = 10\mu m`, and the constants :math:`A_{imp}`, :math:`A_{op}`, :math:`A_{ac}` are interpolated from:

.. list-table:: 
   :align: center
   :header-rows: 1

   * - :math:`N \ (\text{cm}^{-3})`
     - :math:`A_{\text{imp}}`
     - :math:`A_{\text{op}}`
     - :math:`A_{\text{ac}}`
   * - 1.0e16
     - 0.004
     - 0.623
     - 0.034
   * - 1.5e16
     - 0.008
     - 0.932
     - 0.052
   * - 2.0e16
     - 0.014
     - 1.239
     - 0.069
   * - 3.0e16
     - 0.031
     - 1.850
     - 0.104
   * - 4.0e16
     - 0.056
     - 2.456
     - 0.139
   * - 5.0e16
     - 0.086
     - 3.051
     - 0.173
   * - 6.0e16
     - 0.123
     - 3.646
     - 0.208
   * - 7.0e16
     - 0.167
     - 4.240
     - 0.243
   * - 8.0e16
     - 0.217
     - 4.815
     - 0.278
   * - 9.0e16
     - 0.273
     - 5.397
     - 0.313
   * - 1.0e17
     - 0.314
     - 5.578
     - 0.325
   * - 1.5e17
     - 0.690
     - 8.227
     - 0.491
   * - 2.0e17
     - 1.201
     - 10.79
     - 0.660
   * - 3.0e17
     - 2.602
     - 15.75
     - 1.005
   * - 4.0e17
     - 4.474
     - 20.52
     - 1.360
   * - 5.0e17
     - 6.790
     - 25.16
     - 1.726
   * - 6.0e17
     - 9.510
     - 29.65
     - 2.100
   * - 7.0e17
     - 12.64
     - 34.11
     - 2.488
   * - 8.0e17
     - 16.13
     - 38.44
     - 2.879
   * - 9.0e17
     - 20.00
     - 42.75
     - 3.285
   * - 1.0e18
     - 24.22
     - 47.01
     - 3.699
   * - 1.5e18
     - 50.28
     - 67.80
     - 5.912
   * - 2.0e18
     - 93.91
     - 88.02
     - 8.354
   * - 3.0e18
     - 170.3
     - 127.0
     - 13.87
   * - 4.0e18
     - 276.7
     - 164.1
     - 20.12
   * - 5.0e18
     - 396.1
     - 199.6
     - 26.98
   * - 6.0e18
     - 522.8
     - 233.4
     - 34.34

As for the change in refractive index, we consider the simple Lorentz model:

.. math::
    \Delta n(E) = -\frac{1}{2}  \frac{N e^2}{m_e \epsilon_0 (E/\hbar)^2 n_0}

Intervalence absorption
~~~~~~~~~~~~~~~~~~~~~~~

The bandstructure of InGaAsP has two valence bands (light and heavy holes) and a third band due to the spin-orbit coupling. When a photon hits such a material which is hevily P-doped, it can cause an electron in the spin-orbit valence band to jump to one of the hole bands. This effect can be quite intense in III-V at 1550nm because the :math:`\Delta_{so} < 0.8eV`. In IngaAsP materials, this effect has been thorougly characterized, and here we employ the model from :cite:t:`weber_optimization_1994`:

.. math::
    \Delta \alpha(E,N) = 4.252\times 10^{-20} \exp \left(-3.657 E\right) P

.. math::
    \Delta n(E,N) = -\frac{\hbar c \alpha_0}{\pi} \frac{1}{2E}\left(e^{-bE}E_i(bE) + e^{bE}E_1(bE)\right) P 

where :math:`\alpha_0 = 4.252\times 10^-20 m^2` and :math:`b = 3.657 eV^{-1}`. :math:`E_i` and :math:`E_1` are the exponential integrals.


Field effects
-------------

In InGaAsP materials you have a linear and a quadratic field effect.

Pockels effect
~~~~~~~~~~~~~~

The pockels effect follows the works of :cite:t:`adachi_internal_1983` and :cite:t:`adachi_linear_1984`. The electro-optic coefficient :math:`r_{41}` is given by the sum of the free and piezoelectric contributions. Namely:

.. note::
    For simplicity, here we will not write the :math:`\Delta E_{BGN}` contribution, but consider it implied.

.. math::
    r_{41, free} = -\frac{1}{\epsilon_r^2} \left(E_0 g\left(\frac{E}{E_g}\right) + F_0\right)

.. math::
    \begin{aligned}
    r_{41, \text{piezo}} &= 
    - \frac{1}{\epsilon_r^2} 
    \Bigg( 
        C \bigg( 
            -g\left(\frac{E}{E_g}\right) 
            + 4 \frac{E_g}{\Delta_{\text{so}}} 
            \bigg[ 
                f\left(\frac{E}{E_g}\right) \\
            &\quad - \left(\frac{E_g}{E_g + \Delta_{\text{so}}}\right)^{1.5} 
                f\left(\frac{E}{E_g + \Delta_{\text{so}}}\right) 
            \bigg] 
        \bigg) 
        + D 
    \Bigg) e_{14}
    \end{aligned}

where 

.. math::
    g(\chi) = \frac{1}{\chi^2} \left( 2 - (1 + \chi)^{-0.5} - (1 - \chi)^{-0.5} \right)

.. math::
    f(\chi) = \frac{1}{\chi^2} \left( 2 - (1 + \chi)^{0.5} - (1 - \chi)^{0.5} \right)

The constants :math:`E_0`, :math:`F_0`, :math:`C`, :math:`D` are found via interpolation as:

.. math::
    Q(x,y) = xyQ_{GaAs} + x(1-y)Q_{GaP} + y(1-x)Q_{InAs} + (1-x)(1-y)Q_{InP}

and the constants for the binaries are found in the table below:

.. list-table:: Constants for materials
    :align: center
    :header-rows: 1

   * - Constant
     - InP
     - GaP
     - GaAs
     - InAs
   * - E0 (m V^-1)
     - -42.06e-12
     - -83.31e-12
     - -71.48e-12
     - -30.23e-12
   * - F0 (m V^-1)
     - 91.32e-12
     - 16.60e-12
     - 123.16e-12
     - 197.88e-12
   * - C (m^2 N^-1)
     - -0.36e-10
     - -0.06e-10
     - -0.21e-10
     - -1.48e-10
   * - D (m^2 N^-1)
     - 2.60e-10
     - 1.92e-10
     - 2.12e-10
     - 2.32e-10

Kerr effect
~~~~~~~~~~~

The kerr effect is due to the Franz-Keldysh effect, which causes a ripple in the absorption band-edge, which translates to a quadratic change in the refractive index. Here, we follow the work of :cite:t:`adachi_quadratic_1984` and :cite:t:`maat_inp-based_2001`. Contrary to other materials, like ferroelectric materials like lithium niobate, the quadratic response to an electric field is complex. This is a consequence of the Franz Keldysh effect. The absorption for TE and TM polarization is given by:

.. math::
    \Delta \alpha (E)_{TE,TM} = A_{TE,TM}\lambda \frac{|E_ext|}{E_g - E} 10^{-B_{TE,TM}\frac{(E_g-E)^{3/2}}{|E_{ext}|}}

As for the :math:`S_{11}` and :math:`S_{12}`, they are given by:

.. math::
    \begin{aligned}
        S_{11} &= C_{TE} \frac{E^2}{\epsilon_r^2(E_g^2 - E^2)^2} \\
        S_{12} &= C_{TM} \frac{E^2}{\epsilon_r^2(E_g^2 - E^2)^2}
    \end{aligned}

and the constants :math:`C_{TE}` and :math:`C_{TM}`:

.. math::
    \begin{aligned}
        C_{TE} &= -3.10\times 10^{-18} eV^2 m^2 V^{-2} \\
        C_{TM} &= -5.60\times 10^{-18} eV^2 m^2 V^{-2}
    \end{aligned}


.. autoclass:: imodulator.ElectroOpticalModel.InGaAsPElectroOpticalModel
   :members:
   :special-members: __init__


References
----------

.. bibliography::
   :filter: docname in docnames
   :style: plain

