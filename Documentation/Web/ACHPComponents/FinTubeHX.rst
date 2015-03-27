.. _Fin-Tube-HX:

Fin-Tube Heat Exchangers
************************

Overview
--------

In ACHP, all of the heat exchangers that exchange heat with an air stream are of the fin-tube type.  In practice, other types of heat exchangers (shell-tube, microchannel, etc.) are also possible, but they are not included in ACHP as of the development of this documentation.

Mathematical Description
------------------------

.. plot:: MPLPlots/HXSchematics.py

There is some disagreement in literature as to the best way to describe the tube layouts.  The nomenclature used here is that there are a number of banks of tubes, where the banks are the vertical columns of tubes when viewed end-on with the air flow passing from left to right.  The above heat exchanger has a total of 6 banks of tubes, with 4 tubes per bank; there are 4 refrigerant circuits.

There are then a certain number of tubes per bank which form the core of the heat exchanger (not considering the circuiting).  The following figure defines the terms which describe the heat exchanger core, including the longitudinal spacing (also sometimes called bank-to-bank spacing), 

.. plot:: MPLPlots/HXTubesTerms.py

The empirical correlations for air-side heat transfer and pressure drop of fin-tube heat exchangers can be found in section :ref:`Air-Side-Correlations`, and the fluid-side correlations can be found in :ref:`Fluid-Side-Correlations`.

In order to use the heat exchanger in ACHP, complexities of circuiting are neglected.  In practice, it is uncommon for the number of tubes per bank to be divisible by the number of circuits.  As a result, some circuits can be longer than others.  For the purposes of ACHP, average circuit lengths are employed.  The average circuit length is given by

.. math::
    :label: eqAC1

    L_{tubes,total}=N_{tubes/bank}N_{bank}L_{tube}
    
    \overline{L_{circuit}}=L_{tubes,total}/N_{circuits}
    
This average circuit length is primarily required for the calculation of the fluid-side pressure drop.  Other parameters that are required are the total refrigerant-side volume :math:`V_{r,total}`, which can be given by

.. math::
    :label: eqAC2
    
    V_{r,total}=L_{tubes,total}\frac{\pi D_i^2}{4}
    
where :math:`\pi D_i^2/4` is the internal cross-sectional area of the tube.  The refrigerant side surface area is given by

.. math::
    :label: eqAC3
    
    A_{r,total}=\pi D_iL_{tubes,total}

.. _Air-Side-Correlations:

Air-Side Correlations and Other Calculations
--------------------------------------------

.. plot:: MPLPlots/FinTermDefinitions.py

Geometric Parameters
^^^^^^^^^^^^^^^^^^^^

Fins per meter [1/m]

.. math::
    :label: eqAC4
    
    FPM = FPI / 0.0254

Fin pitch (distance between centerlines of fins)

.. math::
    :label: eqAC5
    
    p_f = 1 / FPM

Spacing between fins

.. math::
    :label: eqAC6
    
    s = 1 / FPM - t

Height of heat exchanger [m]

.. math::
    :label: eqAC7

    H = P_t N_{tubes/bank}
    
:math:`A_{duct}` is the face area [m\ :sup:`2`\ ] equivalent to the duct cross-section

.. math::
    :label: eqAC8
    
    A_{duct} = H L_{tube}
    
Number of fins in the tube sheet [-]

.. math:: 
    :label: eqAC9
    
    N_{fin} = L_{tube} FPM
    
Secant of theta is the area enhancement factor [-].  It captures the increase in area due to the waviness of the fins 

.. math:: 
    :label: eqAC10
    
    \sec\theta = \frac{\sqrt{x_f^2 + p_d^2}}{x_f}

Duct cross-sectional area that is not fin or tube [m\ :sup:`2`\ ]

.. math:: 
    :label: eqAC11
    
    A_c = A_{duct} - t N_{fin} (H-DN_{tubes/bank}) - N_{tubes/bank}  D L_{tube}

Total outer area of the tubes [m\ :sup:`2`\ ]

.. math:: 
    :label: eqAC12
    
    A_{tube} = N_{tubes/bank} N_{bank} \pi D L_{tube}

Wetted Area of a single fin [m\ :sup:`2`\ ]

.. math:: 
    :label: eqAC13
    
    A_{1fin} = 2 (H P_l N_{bank} \sec\theta  - N_{tubes/bank}N_{bank} \pi D^2/4)

Total wetted area of the fins [m\ :sup:`2`\ ]

.. math:: 
    :label: eqAC14
    
    A_f = N_{fin} A_{1fin}

Total air-side area including tube and fins [m\ :sup:`2`\ ]

.. math:: 
    :label: eqAC15
    
    A_{a,total} = A_f + N_{tubes/bank} N_{bank} \pi D (L_{tube}-N_{fin}t)

Heat Transfer and Pressure Drop Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Evaluate the mass flow rate based on inlet conditions

.. math::
    :label: eqAC16
    
    \rho_{ha} = \frac{1 + W}{v_{ha}}
    
.. math::
    :label: eqAC17
    
    \dot m_{ha} = \dot V_{ha} \rho_{ha}

    u_{max} = \frac{\dot m_{ha}}{\rho_{ha} A_c}

Specific heat, thermal conductivity and viscosity based on humid air property correlations.

.. math::
    :label: eqAC18

    \mathrm{Pr} = \frac{c_{p,ha} \mu_{ha}}{k_{ha}}

Reynolds number based on the tube diameter:

.. math::
    :label: eqAC19
    
    \mathrm{Re}_D = \frac{\rho_{ha} u_{max} D}{\mu_{ha}}

Empirical Overall Correlations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Wavy-Louvered fins from Wang,Tsai,Lu [#WangTsaiLu]_:

    Colburn j-Factor:

    .. math::
        :label: eqAC20

        j = 16.06 \mathrm{Re}_D^{-1.02 (p_f / D) - 0.256} \left(\frac{A_{a,total}}{A_{tube}}\right)^{-0.601} N_{bank}^{-0.069}\left(\frac{p_f}{D}\right)^{0.84}
        
    Air-side mean heat transfer coefficient:

    .. math::
        :label: eqAC21

        \alpha_a = \frac{j \rho_{ha} u_{max} c_{p,a}}{\mathrm{Pr}^{2/3}}

    Air-side pressure drop friction factor:

    .. math::
        :label: eqAC22
            
        f_{a,total}=0.264(0.105+0.708\exp(-\mathrm{Re}_D/225.0))\mathrm{Re}_D^{-0.637}\left(\frac{A_{a,total}}{A_{tube}}\right)^{0.263}\left(\frac{p_f}{D}\right)^{-0.317}   (\mathrm{Re}_D>1000)

        f_{a,total}=0.768(0.0494+0.142\exp(-\mathrm{Re}_D/1180.0))\left(\frac{A_{a,total}}{A_{tube}}\right)^{0.0195}\left(\frac{p_f}{D}\right)^{-0.121}   (\mathrm{Re}_D\leq1000)
        
The air mass flux through the heat exchanger can be defined by

.. math::
    :label: eqAC23
    
    G_c=\frac{\dot m_{ha}}{A_c}

which yields the air-side pressure drop (neglecting entrance and exit pressure drops) of

.. math::
    :label: eqAC24
    
    \Delta p_{a}=\frac{A}{A_{tube}}\frac{G_c^2}{2\rho_{ha}}f_{a,total}
 
Surface Efficiency
^^^^^^^^^^^^^^^^^^

Thermal gradients in the fins result in a finned surface that does not perform as efficiently as if the entire finned surface was at the fluid temperature passing through the tubes.  Correlation (or explicit solution where possible) is required to be used to determine the temperature profile in the fins, and from that, the finned surface efficiency.

Circular fins
"""""""""""""

Schmidt [#Schmidt]_ developed a correction term for circular fins to calculate fin efficiency. Using the circular fin correlation of Schmidt yields the solution for the finned surface efficiency of 

.. math::
    :label: eqAC25
    
    \phi = \left(\frac{r_f}{r} - 1\right) \left[1 + 0.35 \ln\left(\frac{r_f}{r}\right)\right]
    
    \eta_f = \frac{\tanh(m r \phi)}{ m r \phi}
    
where :math:`r` is the outer diameter of the tube, and :math:`r_f` is the outer radius of the circular fin.  The parameter :math:`m` is given by

.. math::
    :label: eqAC28
    
    m = \sqrt{\frac{2 \alpha_a (c_s/c_p)}{k_{fin} t}} 

.. _Staggered-Fin-Efficiency:

Staggered tubes
"""""""""""""""

For the staggered tube bank configurations (like that in ACHP), hexagonal cells with adiabatic boundaries are formed around each tube, and for each cell the fin efficiency can be obtained.  Mathematical analysis is used to convert this problem into that like the circular fin.

.. plot:: MPLPlots/StaggeredSurfaceEfficiency.py

.. math::
    :label: eqAC26
    
    r = \frac{D}{2}
    
    X_D = \frac{\sqrt{P_l^2 + P_t^2 / 4}}{2}
    
    X_T = \frac{P_t}{2}
    
The effective radius ratio is given by

.. math::
    :label: eqAC27
    
    \frac{r_f}{r} = 1.27 \frac{X_T}{r} \sqrt{\frac{X_D}{X_T} - 0.3}
    
And the :math:`m` factor is given by

.. math::
    :label: eqAC27a
    
    m = \sqrt{\frac{2 \alpha_a (c_s/c_p)}{k_{fin} t}} 
    
:math:`c_s/c_p` is the correction for heat/mass transfer for a wetted surface.  If the fins are dry, this parameter is set to unity, which yields the standard definition for the parameter :math:`m` for a fin

Fin efficiency based on analysis by Hong and Webb [#HongWebb]_ for wet and dry fins with staggered fins

.. math::
    :label: eqAC29a
    
    \phi = \left(\frac{r_f}{r} - 1\right) \left[1 + \left(0.3+\left(\frac{mr(\frac{r_f}{r}-r)}{2.5}\right)^{1.5-\frac{1}{12}\frac{r_f}{r}}\left(0.26\left(\frac{r_f}{r}\right)^{0.3}-0.3\right)\right) \ln\left(\frac{r_f}{r}\right)\right]

.. math::
    :label: eqAC29b
    
    \eta_f = \frac{\tanh(m r \phi)}{m r \phi} \cos(0.1 m r \phi)

Plots
^^^^^

.. plot:: MPLPlots/FinCorrelationsPlts.py

Overall efficiency
""""""""""""""""""

Once the finned surface efficiency (:math:`\eta_f`) is known, the overall surface efficiency is given by

.. math::
    :label: eqAC30
    
    \eta_a = 1 - \frac{A_f}{A_{a,total}} (1 - \eta_f)

.. only:: html

    .. rubric:: References

.. [#HongWebb] Kwang Taek Hong and Ralph L. Webb, 1996, "Calculation of Fin Efficiency for Wet and Dry Fins" *HVAC&R Research*,v. 2 `Link to File <http://www.tandfonline.com/doi/abs/10.1080/10789669.1996.10391331>`_

.. [#WangTsaiLu] Chi-Chuan Wang and Yu-Min Tsai and Ding-Chong Lu, 1998, "Comprehensive Study of Convex-Louver and Wavy Fin-and-Tube Heat Exchangers", *Journal of Thermophysics and Heat Transfer*

.. [#Schmidt] Schmidt T.E. 1945-46. La Production Calorifique des Surfaces Munies D'ailettes. *Bulletin De L'Institut International Du Froid Annexe G-5*.

**Nomenclature**

===============================  ===================================================
Variable                         Description
===============================  ===================================================
:math:`A_{a,total}`              Total air-side area (fins+tubes) [m\ :superscript:`2`\ ]
:math:`A_f`                      Air-side area of fins [m\ :superscript:`2`\ ]
:math:`A_{1fin}`                 Air-side area of 1 fin [m\ :superscript:`2`\ ]
:math:`A_c`                      Cross-sectional area of duct not filled with fins and tubes [m\ :superscript:`2`\ ]
:math:`A_{r,total}`              Surface area on the refrigerant (or fluid) side [m\ :superscript:`2`\ ]
:math:`A_{duct}`                 Duct cross-sectional area [m\ :superscript:`2`\ ]
:math:`A_{tube}`                 Cross-sectional area of the tubes [m\ :superscript:`2`\ ]
:math:`c_{p,a}`                  Air specific heat [J/kg\ :sub:`ha`\ /K]
:math:`D_i`                      Interior diameter of tubes [m]
:math:`G_c`                      Air mass flux [kg/m\ :sup:`2`\ ]
:math:`f_{a,total}`              Air-side friction factor [-]
:math:`FPI`                      Fins per inch [1/in]
:math:`FPM`                      Fins per meter [1/m]
:math:`k_{fin}`                  Fin conductivity [W/m/K]
:math:`H`                        Height [m]
:math:`j`                        Colburn j-factor [-]
:math:`L_{tube}`                 Length of one tube [m]
:math:`L_{tubes,total}`          Sum total length of all tubes [m]
:math:`\overline{L_{circuit}}`   Average length of circuit [m]
:math:`m`                        Non-dimensional group in fin efficiency calculation [-]
:math:`\dot m_{ha}`               Mass flow rate of humid air [kg/s]
:math:`N_{tubes/bank}`           Number of tubes per bank [-]
:math:`N_{bank}`                 Number of banks of tubes [-]
:math:`N_{circuits}`             Number of circuits [-]
:math:`N_{fin}`                  Number of fins [-]
:math:`p_f`                      Twice the amplitude of fins corrugation [m]
:math:`p_l`                      Longitudinal bank-bank pitch (in the flow direction)[m]
:math:`p_t`                      Transverse pitch (perpindicular in the flow direction)[m]
:math:`\Delta p_a`               Air side pressure drop [Pa]
:math:`\mathrm{Pr}`              Prandtl number [-]
:math:`r`                        Outer radius of tube [m]
:math:`r_f`                      Circular fin radius [m]
:math:`\mathrm{Re}_D`            Reynolds number based on tube OD [-]
:math:`s`                        Spacing between fins [m]
:math:`\sec \theta`              Area increase factor [-]
:math:`t`                        Fin thickness [m]
:math:`u_{max}`                  Maximum air velocity [m/s]
:math:`v_{ha}`                   Density of humid air [m^3/kg]
:math:`\dot V_{ha}`              Volumetric flow rate of humid air [m\ :sup:`3`\ /s]
:math:`V_{r,total}`              Open volume on the refrigerant (or fluid) side [m\ :sup:`3`\ ]
:math:`W`                        Humidity ratio [-]
:math:`x_f`                      Half-wavelength of fin wave [m]
:math:`X_D`                      Diagonal distance in fin efficiency equation [m]
:math:`X_T`                      Transverse distance in fin efficiency equation [m]
:math:`\alpha_a`                 Air-side mean heat transfer coefficient [W/m\ :sup:`2`\ /K
:math:`\rho_{ha}`                Density of air [kg/m\ :sup:`3`\ ]
:math:`\phi`                     Surface efficiency parameter [-]
:math:`\eta_a`                   Overall finned surface efficiency [-]
:math:`\eta_f`                   Fin efficiency [-]
===============================  ===================================================