
.. |eNtu| replace:: :math:`\varepsilon-\mathrm{Ntu}`
.. |Ntu| replace:: :math:`\mathrm{Ntu}`
.. |UA| replace:: :math:`\mathrm{UA}`

.. _Heat-Mass-HX:

**************************************************
Heat Exchangers with and without Dehumidification
**************************************************

In the heat exchanger models for ACHP (for evaporator, condenser, cooling coil, etc.), air and a working fluid exchange heat.  In each part of the heat exchanger model, the working fluid phase (single-phase, two-phase) is assumed to be known.  Therefore, the given segment of the heat exchanger can be analyzed using the analysis which follows here.  In addition, depending on the surface temperatures in the heat exchangers, there exists the possibility that water moisture in the air might condense on the finned surface of the heat exchanger, which must be properly addressed.

Essentially the heat exchanger models can be broken down into these four possibilities:

* Working fluid is two-phase, there is no possibility of condensation of air on coil (always in condenser, possibly in evaporator and cooling coil)
* Working fluid is single-phase, there is no possibility of condensation of air on coil (always in condenser, possibly in evaporator and cooling coil)
* Working fluid is two-phase, there is a possibility of condensation of air on coil (evaporator and cooling coil)
* Working fluid is single-phase, there is a possibility of condensation of air on coil (evaporator and cooling coil)

Heat Exchanger Analysis Fundamentals
====================================

The analysis that follows is heavily based on the :math:`\varepsilon-\mathrm{Ntu}` method of analysis of heat exchangers.  This method assumes that the specific heat of both fluids are constant (or a reasonable average value can be found).  The beauty of the |eNtu| method is that if the heat transfer coefficients and areas can be obtained explicitly, and the inlet temperatures of both streams are known, the heat transfer rate can be obtained explictly.  This proves highly beneficial to the development of computationally-efficient code.  In many of the individual heat exchanger solvers, the equations are employed in a slightly different configuration, but the governing equations presented here describe the fundamental models.

Broadly, the way the |eNtu| method works is to first find the heat transfer conductance :math:`\mathrm{UA}` from the heat transfer network analysis, which yields the number of thermal units (|Ntu|), which gives you the effectiveness :math:`\varepsilon`, which allows for the calculation of everything else.

Analysis of fully-dry HX without Dehumidification
=================================================

Two-phase on working fluid side
-------------------------------
The combination of a two-phase working fluid and no possibility of moisture removal from the air yields the simplest analysis.  The capacitance rate of air by definition is the limiting capacitance rate (specific heat of two-phase mixture is infinite), and therefore, the |UA| are given by

.. math::
    :label: eqWDHE1

    \mathrm{UA}_a=\eta_a \alpha_a A_a 
    
    \mathrm{UA}_r=\alpha_r A_r
    
the overall heat transfrer conductance by

.. math::
    :label: eqWDHE2

    \frac{1}{\mathrm{UA}}=\frac{1}{1/\mathrm{UA}_a+1/\mathrm{UA}_r}

and the |Ntu| from

.. math::
    :label: eqWDHE3

    \mathrm{Ntu}=\frac{\mathrm{UA}}{\dot m_a c_{p,a}}
    
The effectiveness is obtained from 

.. math::
    :label: eqWDHE4
    
    \varepsilon=1-\exp(-\mathrm{Ntu})
    
because the ratio of capacitance rates is zero (:math:`C_r=0`).  Since the air-side is by definition the limiting capacitance rate, the heat transfer rate is given by

.. math::
    :label: eqWDHE5
    
    \dot Q = \varepsilon \dot m_a c_{p,a} (T_{h,i}-T_{c,i})
    
where :math:`T_{h,i}` is the inlet temperature of the hot stream, and :math:`T_{c,i}` is the inlet temperature of the cold stream.  By this definition, the sign of :math:`\dot Q` is defined to be positive.

Single-phase on working fluid side
----------------------------------

When there does not exist the possibility that the humid air stream can condense, the analysis of pure-heat transfer is relatively simple.  In the heat exchanger, the air will approach the working fluid intlet temperature, and the working fluid will exchange heat with the air, potentially, though not necessarily, bringing the temperature of the working fluid closer to that of the air.

The heat transfer conductances on air- and working fluid-sides are given by

.. math::
    :label: eqWDHE6
    
    \mathrm{UA}_a=\eta_a \alpha_a A_a 
    
    \mathrm{UA}_w=\alpha_w A_w
    
where the heat transfer coefficients :math:`\alpha` and areas :math:`A` are governed by the analysis in sections :ref:`Air-Side-Correlations`, :ref:`Fin-Tube-HX`, and :ref:`Fluid-Side-Correlations`.  The overall heat transfer conductance :math:`\mathrm{UA}` can be obtained from

.. math::
    :label: eqWDHE7
    
    \frac{1}{\mathrm{UA}}=\frac{1}{1/\mathrm{UA}_a+1/\mathrm{UA}_w}

Thus the number of thermal units can be determined by 

.. math::
    :label: eqWDHE8
    
    \mathrm{Ntu}=\frac{\mathrm{UA}}{C_{min}}
    
where :math:`C_{min}` is the minimum capacitance rate, given by

.. math::
    :label: eqWDHE9
    
    C_{min}=\min[\dot m_a c_{p,a}, \dot m_w c_{p,w}]
    
Thus with :math:`\mathrm{Ntu}` known, the effectiveness can be obtained based on the geometry of the heat exchanger, shown here

* Counterflow: 
    
.. math::
    :label: eqWDHE10
    
    \varepsilon=\frac{1-\exp[-\mathrm{Ntu}(1-C_r)]}{1-C_r\exp[-\mathrm{Ntu}(1-C_r)]}

* Crossflow (:math:`C_{a}<C_w` [air-side capacitance rate is minimum and air is treated as being unmixed because the fins block mixing of the air and the fluid is treated as being mixed]): 

.. math::
    :label: eqWDHE11
    
    \varepsilon= \dfrac{1}{C_r} (1 - \exp(-C_r (1 - \exp(-\mathrm{Ntu}))))\\
    
* Crossflow (:math:`C_{w}<C_a` [working-fluid-side capacitance rate is minimum and air is treated as being unmixed because the fins block mixing of the air and the fluid is treated as being mixed]): 

.. math::
    :label: eqWDHE12
    
    \varepsilon = 1-\exp\left(-\frac{1}{C_r}[1-\exp(-C_r\mathrm{Ntu})]\right)

Which yields the heat transfer rate of 

.. math::
    :label: eqWDHE13
    
    \dot Q = \varepsilon C_{min} (T_{h,i}-T_{c,i})

Analysis of Partial-wet/Partial-dry Heat Exchangers with the possibility of Dehumidification
============================================================================================

.. _PWPD-Two-Phase:

Two-phase on working fluid side
-------------------------------
In the analysis here, averaged heat transfer coefficients on the air-side (from sections :ref:`Air-Side-Correlations` and :ref:`Fin-Tube-HX`) and averaged heat transfer coefficients on the refrigerant side (from two-phase evaporation part of section :ref:`Fluid-Side-Correlations`) are used.  The slight wrinkle on the refrigerant side is that the average heat transfer coefficient is a function of refrigerant quality, but the analysis here assumes that the inlet and outlet qualities of the refrigerant are known (or are being iterated for).

Another complication is that for refrigerants that experience glide during phase change (azeotropic blends like R404a, R410a, R507c, R407a, etc.), the saturation temperature is not constant during evaporation.  If the pressure during phase change is assumed to be constant, there are different temperatures for the saturated liquid (bubble point), and the saturated vapor (dew point) for a given saturation pressure.  The analysis which follows for the two-phase section requires an average saturation temperature, which is simply defined to be :math:`T_{sat,r}=(T_{bubble,r}+T_{dew,r})/2` if the fluid has glide.  Otherwise, the conventional definition of the saturation temperature :math:`T_{sat,r}` is used.

Fully-dry
"""""""""
The heat transfer conductances on inner- (working fluid), and outer-sides (air) are given by

.. math::
    :label: eqWDHE13a
    
    \mathrm{UA}_o=\alpha_a\eta_aA_{a,two-phase} 
    
    \mathrm{UA}_i=\alpha_{r,two-phase}A_{r,two-phase}
    
and the number of thermal units outside is given by

.. math::
    :label: eqWDHE13b
    
    \mathrm{Ntu}_o=\frac{\mathrm{UA}_o}{\dot m_{a,two-phase}c_{p,a}}
    
In the fully dry analysis for the two-phase section, the :math:`\mathrm{UA}` value is obtained from 

.. math::
    :label: eqWDHE14
    
    \mathrm{UA}_{dry}=\frac{1}{\mathrm{UA}_i^{-1}+\mathrm{UA}_o^{-1}}
    
which, since :math:`C_{min}` is on the air-side (because the refrigerant is changing phase), yields 

.. math::
    :label: eqWDHE15
    
    \mathrm{Ntu}_{dry}=\frac{\mathrm{UA}_{dry}}{\dot m_{a,two-phase}c_{p,a}}
   
and the effectiveness

.. math::
    :label: eqWDHE16
    
    \varepsilon_{dry}=1-\exp(-\mathrm{Ntu}_{dry})

which yields the fully-dry heat transfer rate of

.. math::
    :label: eqWDHE17
    
    \dot Q_{dry}=\varepsilon_{dry}c_{p,a}(T_{a,i}-T_{sat,r})
    
which allows to calculate the inlet and outlet surface temperatures from

.. math::
    :label: eqWDHE18
    
    T_{s,i}=\frac{\mathrm{UA}_oT_{a,i}+\mathrm{UA}_iT_{sat,r}}{\mathrm{UA}_o+\mathrm{UA}_i}
    
    T_{s,o}=\frac{\mathrm{UA}_oT_{a,o}+\mathrm{UA}_iT_{sat,r}}{\mathrm{UA}_o+\mathrm{UA}_i}
    
For a given inlet air dewpoint temperature of :math:`T_{dp}`, if

* :math:`T_{s,o}` is above :math:`T_{dp}`, the coil is fully dry

* :math:`T_{s,i}` is below :math:`T_{dp}`, the coil is fully wet, and proceed to the analysis in the section :ref:`fully-wet`

* :math:`T_{dp}` is between :math:`T_{s,i}` and :math:`T_{s,o}`, the coil is partially wet, and proceed to the analysis in the section :ref:`partially-wet-dry-two-phase`.

If the coil is fully dry, then the factor :math:`f_{dry}` is equal to unity, and the outlet air enthalpy can be calculated from 

.. math::
    :label: eqWDHE19
    
    h_{a,o}=h_{a,i}-\dot Q_{dry}/\dot m_{a,two-phase}
    
And the heat transfer rate for the two-phase section of the evaporator for the given set of :math:`x_{o,two-phase}` and :math:`w_{two-phase}` yields

.. math::
    :label: eqWDHE20
    
    \dot Q_{two-phase}(w_{two-phase},x_{o,two-phase})=\dot Q_{dry}

.. _fully-wet:

Fully-wet
"""""""""

In the fully-wet analysis, the fully-dry section of the two-phase portion of the evaporator doesn't exist.  Therefore, the analysis which follows for the partially-wet/ partially-dry case can be employed directly with the inlet to the wet section given by the evaporator air inlet condition.  This can be expressed as 

.. math::
    :label: eqWDHE21
    
    h_{a,x}=h_{a,i}
    
    T_{a,x}=T_{a,i}
    
    f_{dry}=0
    
.. _partially-wet-dry-two-phase:

Partially-dry
"""""""""""""

In the partially-wet/ partially-dry analysis for a heat exchanger with two-phase working fluid, the configuration for the coil is described by this figure:

.. plot:: MPLPlots/EvaporatorWetDry.py

To begin the analysis for the dry part (this analysis is skipped if the coil is fully wetted), the temperature :math:`T_{a,x}` can be obtained directly from an energy balance at the point where the dewpoint is reached (though its location is not known *a priori*).  
    
.. math::
    :label: eqWDHE22
    
    T_{a,x} = T_{dp} + \frac{\mathrm{UA}_i}{\mathrm{UA}_o}(T_{dp} - T_{sat,r})
    
where the ratio of :math:`\mathrm{UA}` factors can be obtained from

.. math::
    :label: eqWDHE23
    
    \frac{\mathrm{UA}_i}{\mathrm{UA}_o}=\frac{\alpha_rA_{r,two-phase}}{\eta_{a}\alpha_{a}A_{a,two-phase}}
    
which yields an explicit solution for the effectiveness (because the inlet and outlet states are known on the air-side) :math:`\varepsilon_{dry}` of

.. math::
    :label: eqWDHE24
    
    \varepsilon_{dry}=\frac{T_{a,i}-T_{a,x}}{T_{a,i}-T_{sat,r}}
    
and the :math:`f_{dry}` factor is obtained from

.. math::
    :label: eqWDHE25
    
    f_{dry}=-\frac{1}{\mathrm{Ntu}_{dry}}\ln(1-\varepsilon_{dry})
    
by solving the equation

.. math::
    :label: eqWDHE26
    
    \varepsilon_{dry}=1-\exp(-f_{dry}\mathrm{Ntu}_{dry})
    
The enthalpy of the air at the temperature :math:`T_{a,x}` can be obtained from

.. math::
    :label: eqWDHE27
    
    h_{a,x}=h_{air}(T=T_{a,x},\omega=\omega_{a,i})
    
which yields the heat transfer rate in the dry part of 

.. math::
    :label: eq-PWPDQdry

    \dot Q_{dry}=\dot m_{a,two-phase}c_{p,a}(T_{a,i}-T_{a,x})
    
This completes the analysis for the dry part of the partially-wet/ partially-dry analysis.  The wet portion begins with the calculation of the factor :math:`c_s`, which is governed by the equation

.. math::
    :label: eqWDHE28
    
    c_s=\left. \frac{\partial h_{sat}}{\partial T} \right|_{T=T_{sat,r}}
    
Using this value for :math:`c_s`, a new wetted-fin-efficiency (:math:`\eta_a^*`) is calculated using the analysis in section :ref:`Staggered fin surface efficiency <Staggered-Fin-Efficiency>`.  This then yields the :math:`\mathrm{UA}` values 

.. math::
    :label: eqWDHE29
    
    \mathrm{UA}_i=\alpha_rA_{r,two-phase}
    
    \mathrm{UA}_o^*=\eta_{a}^*\alpha_{a}A_{a,two-phase}
    
and the overall wet UA value of

.. math::
    :label: eqWDHE30
    
    \mathrm{UA}_{wet}=\frac{1}{c_s/\mathrm{UA}_i+c_{p,a}/\mathrm{UA}_o^*}
    
which yields the :math:`\mathrm{Ntu}` value of

.. math::
    :label: eqWDHE31
    
    \mathrm{Ntu}_{wet}=\frac{\mathrm{UA}_{wet}}{\dot m_{a,two-phase}}
    
and the wet surface effectiveness of 

.. math::
    :label: eqWDHE32
    
    \varepsilon_{wet}=1-\exp[-(1-f_{dry})\mathrm{Ntu}_{wet}]

With the air saturation enthalpy :math:`h_{a,sat,r}` at the refrigerant saturation temperature :math:`T_{sat,r}` of 

.. math::
    :label: eqWDHE33
    
    h_{a,sat,r}=h_{sat}(T=T_{sat,r},\phi=1)
    
the heat transfer rate in the wet section is given by

.. math::
    :label: eq-PWPDQwet

    \dot Q_{wet}=\varepsilon_{wet}\dot m_{a,two-phase}(h_{a,i}-h_{a,sat,r})
    
and the outlet enthalpy can be given by

.. math::
    :label: eqWDHE34
    
    h_{a,o}=h_{a,x}-\dot Q_{wet}/\dot m_{a,two-phase}
    
The air in the wet section interacts with a surface that has an effective surface enthalpy of

.. math::
    :label: eqWDHE35
    
    h_{a,s,s,e}=h_{a,x}-\frac{h_{a,x}-h_{a,o}}{1-\exp[-(1-f_{dry})\mathrm{Ntu}_o]}
    
which is obtained by considering just the effectiveness of humid-air mass transfer with a fixed-enthalpy :math:`h_{a,s,s,e}` saturated surface.  Which yields an effective surface temperature :math:`T_{s,e}` obtained iteratively from 

.. math::
    :label: eqWDHE36
    
    h_{a,s,s,e}=h_{air}(T=T_{s,e},\phi=1)
    
This yields the air outlet temperature of

.. math::
    :label: eqWDHE37

    T_{a,o}=T_{s,e} + (T_{a,x}-T_{s,e})\exp[-(1-f_{dry})\mathrm{Ntu_o}]
    
where the air is assumed to transfer heat with an isothermal surface at :math:`T_{s,e}`.  The total heat transfer rate in the partially-wet/ partially-dry analysis is then given from

.. math::
    :label: eqWDHE38
    
    \dot Q_{two-phase}=\dot Q_{dry}+\dot Q_{wet}
    
where :math:`\dot Q_{dry}` comes from Equation :eq:`eq-PWPDQdry`, and :math:`\dot Q_{wet}` comes from Equation :eq:`eq-PWPDQwet`.

.. _PWPD-Single-Phase:

Single-phase on working fluid side
----------------------------------

Fully Dry Analysis
""""""""""""""""""

For the fully dry coil with a single-phase working fluid, (taken to be refrigerant here for the purposes of nomenclature) the internal and external :math:`\mathrm{Ntu}` for the dry analysis are given by

.. math::
    :label: eqWDHE39
    
    \mathrm{UA}_i=\alpha_rA_r

    \mathrm{Ntu}_i=\frac{\mathrm{UA}_i}{\dot m_{r} c_{p,r}}

    \mathrm{Ntu}_o=\frac{\eta_a \alpha_a A_a}{\dot m_ac_{p,a}}

where the :math:`\mathrm{Ntu}` are in general defined by :math:`\mathrm{Ntu}=\mathrm{UA}/C`.

The capacitance rate ratio (assuming the minimum capacitance to be on the air side) is

.. math::
    :label: eqWDHE40
    
    C^*=\frac{\dot m_a c_{p,a}}{\dot m_r c_{p,r}}

The overall number of thermal units for the entire heat exchanger when dry is

.. math::
    :label: eqWDHE41
    
    \mathrm{Ntu}_{dry}=\frac{\mathrm{Ntu}_o}{1+C^* \frac{\mathrm{Ntu}_o}{\mathrm{Ntu}_i}}=\frac{\frac{\mathrm{UA}_o}{C_a}}{1+\frac{\mathrm{UA}_o}{\mathrm{UA}_i}}=\frac{\frac{\mathrm{UA}_o\mathrm{UA}_i}{C_a}}{\mathrm{UA}_i+\mathrm{UA}_o}
 
Counterflow effectiveness relations are used for a non-zero :math:`C^*` since the refrigerant is not changing phase (this is assumed):

.. math::
    :label: eqWDHE42
    
    \varepsilon_{dry} = \frac{1 - \exp(-\mathrm{Ntu}_{dry} (1 - C^*))} {1 - C^* \exp(-\mathrm{Ntu}_{dry} (1 - C^*))}
 
The air outlet temperature is given by

.. math::
    :label: eqWDHE43
    
    T_{a,o,dry}=T_{a,i}-\varepsilon_{dry}(T_{a,i}-T_{r,i})
 
and the dry analysis glycol outlet temperature is given by

.. math::
    :label: eqWDHE44
    
    T_{r,o}=T_{r,i}-C^*(T_{a,i}-T_{a,o,dry})
 
The dry air outlet enthalpy is equal to (from overall energy balance over the heat exchanger and assuming constant specific heat on the refrigerant side)

.. math::
    :label: eqWDHE45
    
    h_{a,o}=h_{a,i}-\frac{\dot m_rc_{p,r}(T_{r,i}-T_{r,o})}{\dot m_a}
 
The surface temperature at the air outlet is determined from an energy balance at glycol inlet (air outlet) and solving for surface temperature

.. math::
    :label: eqWDHE46
    
    T_{s,o}=T_{r,i}+C^*\frac{\mathrm{Ntu}_{dry}}{\mathrm{Ntu}_i}(T_{a,o,dry}-T_{r,i})
 
Heat transfer for dry analysis

.. math::
    :label: eqWDHE47
    
    \dot Q_{dry}=\varepsilon_{dry}\dot m_ac_{p,a}(T_{a,i}-T_{r,i})
 
And dry fraction

.. math::	
    :label: eqWDHE48
    
    f_{dry}=1.0

Fully Wet Analysis
""""""""""""""""""

If :math:`T_{s,o} < T_{dp}` , then there is at least some portion of the coil which is wetted (potentially all the coil) because some part of the heat exchanger surface is below the dew point of the entering air. The heat exchanger is then assumed to be fully wet in the next step.  The possibility still remains that the coil may not be fully wet, but by trying the fully-wet analysis it can be determined whether the coil is fully or partially wetted. Initially, the outlet refrigerant temperature :math:`T_{r,o}` is not known, but it is determined iteratively, with a first guess given by the fully-dry outlet refrigerant temperature from the above analysis.  

The bounding outlet state for the air (the lowest enthalpy it can achieve) is to be saturated at the refrigerant inlet temperature, for which the enthalpy of the saturated air is given by

.. math::
    :label: eqWDHE49
    
    h_{a,sat,r,i}=h(T=T_{r,i},\phi=1.0)
 
The air saturation specific heat at mean refrigerant temperature is

.. math::
    :label: eqWDHE50
    
    c_s=\left. \frac{\partial h_{a,sat}}{\partial T_{r}} \right| _{T_r=\frac{T_{r,i}+T_{r,o}}{2}}
 
and is polynomial fit based on EES saturated water vapor enthalpy data.  The outlet glycol temperature is initially unknown (but is determined iteratively).  

The effective humid air mass flow ratio

.. math::
    :label: eqWDHE51
    
    m^*=\frac{\dot m_a}{\dot m_r\frac{c_{p,r}}{c_s}}
 
and the NTU for wet coil

.. math::
    :label: eqWDHE52
    
    \mathrm{Ntu}_{wet}=\frac{\mathrm{Ntu}_o}{1+m^* \frac{\mathrm{Ntu}_o}{\mathrm{Ntu}_i}}
 
counterflow effectiveness

.. math::
    :label: eqWDHE53
    
    \varepsilon_{wet} = \frac{1 - \exp(-\mathrm{Ntu}_{wet} (1 - m^*))} {1 - m^* \exp(-\mathrm{Ntu}_{wet} (1 - m^*))}
 
Air outlet enthalpy

.. math::
    :label: eqWDHE54
    
    h_{a,o}=h_{a,i}-\varepsilon_{wet}(h_{a,i}-h_{a,sat,r,i})
 
The amount of heat transfer is

.. math::
    :label: eqWDHE55
    
    \dot Q_{wet}=\varepsilon_{wet}\dot m_a(h_{a,i}-h_{a,sat,r,i})
 
and the outlet glycol temperature is

.. math::
    :label: eqWDHE56
    
    T_{r,o}=T_{r,i}+\frac{\dot m_a}{\dot m_r c_{p,r}}(h_{a,i}-h_{a,o})
 
the saturated air enthalpy at the outlet glycol temperature is

.. math::
    :label: eqWDHE57
    
    h_{a,sat,r,o}=h(T=T_{r,o},\phi=1.0)
 
and the surface temperature at the air inlet is (from local energy balance)

.. math::
    :label: eqWDHE58
    
    T_{s,i} = T_{r,o}+\frac{\dot m_a}{\dot m_rc_{p,r}}\frac{\mathrm{Ntu}_{wet}}{\mathrm{Ntu}_i}(h_{a,i}-h_{a,sat,r,o})
 
the effective surface enthalpy is given by

.. math::
    :label: eqWDHE59
    
    h_{a,s,sat,e}=h_{a,i}+\frac{h_{a,o}-h_{a,i}}{1-\exp(-\mathrm{Ntu}_o)}
 
from which it it possible to solve

.. math::
    :label: eqWDHE60
    
    h_{a,s,sat,e}=h(T=T_{s,e},\phi=1.0)
 
for :math:`T_{s,e}` which makes it possible to solve for the outlet air temperature

.. math::
    :label: eqWDHE61
    
    T_{a,o}=T_{s,e} + (T_{a,i}-T_{s,e})\exp(-\mathrm{Ntu}_o)

At this point, if :math:`T_{s,i}` is also below the dewpoint, the coil is truly full wetted, and the analysis is finished.  If on the other hand, :math:`T_{s,i}>T_{dp}`, then the coil is partially-wetted and the analysis in the next section is required.

Partially-wet/Partially-dry
"""""""""""""""""""""""""""

If the cold fluid is water, glycol, or other similar single-phase fluid, the partially-wet/partially-dry analysis is quite a bit more complicated than for two-phase refrigerant.  The full derivation is presented in section :ref:`PWPD-Algorithm`, and the code for this analysis is in :meth:`~DryWetSegment.DryWetSegment`.
    
.. |m3| replace:: m\ :sup:`3`\ 
.. |m2| replace:: m\ :sup:`2`\ 

**Nomenclature**

===============================  ===================================================
Variable                         Description
===============================  ===================================================
:math:`A_a`                      Air-side area [|m2|]
:math:`A_{a,two-phase}`          Air-side area in two-phase section [|m2|]
:math:`A_{r,two-phase}`          Refrigerant-side area in two-phase section [|m2|]
:math:`A_r`                      Refrigerant-side area [|m2|]
:math:`A_w`                      Water-side area [|m2|]
:math:`C_{min}`                  Minimum capacitance rate [W/K]
:math:`c_{p,a}`                  Specific heat of humid air per kg dry air [J/kg\ :sub:`da`\ ]
:math:`c_{p,w}`                  Specific heat of water [J/kg/K]
:math:`c_s`                      Saturation specific heat of humid air [J//kg\ :sub:`da`\ /K]
:math:`C_r`                      Ratio of capacitance rates [-]
:math:`C^*`                      Ratio of capacitance rates [-]
:math:`f_{dry}`                  Fraction of coil that is dry [-]
:math:`f_{wet}`                  Fraction of coil that is wet [-]
:math:`h_{a,i}`                  Enthalpy of humid air at inlet to coil [J/kg\ :sub:`da`]
:math:`h_{a,o}`                  Enthalpy of humid air at outlet from coil [J/kg\ :sub:`da`]
:math:`h_{a,x}`                  Enthalpy of humid air at wet-dry interface [J/kg\ :sub:`da`]
:math:`h_{a,s,s,e}`              Effective surface saturated enthalpy [J/kg\ :sub:`da`]
:math:`h_{a,sat,r,i}`            Enthalpy of saturated air at refrigerant inlet temperature [J/kg\ :sub:`da`]
:math:`h_{a,sat,r}`              Enthalpy of saturated air at refrigerant saturation temperature [J/kg\ :sub:`da`]
:math:`m^*`                      Effective wet ratio of mass flow rates [-]
:math:`\dot m_{a}`               Mass flow rate of dry air [kg/s]
:math:`\dot m_{a,two-phase}`     Mass flow rate of dry air in two-phase section [kg/s]
:math:`\mathrm{Ntu}`             Number of thermal units [-]
:math:`\mathrm{Ntu_o}`           Outside number of thermal units [-]
:math:`\mathrm{Ntu}_{dry}`       Dry number of thermal units [-]
:math:`\mathrm{Ntu}_{wet}`       Wet number of thermal units [-]
:math:`\dot Q`                   Heat transfer rate [W]
:math:`\dot Q_{dry}`             Heat transfer rate from dry analysis [W]
:math:`\dot Q_{two-phase}`       Heat transfer rate in two-phase section [W]
:math:`\dot Q_{wet}`             Heat transfer rate from wet analysis [W]
:math:`T_{a,i}`                  Air inlet drybulb temperature [K]
:math:`T_{a,o}`                  Air outlet drybulb temperature [K]
:math:`T_{a,o,dry}`              Air outlet drybulb temperature from dry analysis [K]
:math:`T_{a,x}`                  Air drybulb temperature at wet-dry interface [K]
:math:`T_{bubble,r}`             Bubble-point temperature of refrigerant [K]
:math:`T_{c,i}`                  Cold stream inlet temperature [K]
:math:`T_{dp}`                   Dewpoint temperature of air [K]
:math:`T_{dew,r}`                Dew-point temperature of refrigerant [K]
:math:`T_{h,i}`                  Hot stream inlet temperature [K]
:math:`T_{r,i}`                  Refrigerant inlet temperature [K]
:math:`T_{r,o}`                  Refrigerant outlet temperature [K]
:math:`T_{sat,r}`                Saturation temperature of refrigerant [K]
:math:`T_{s,e}`                  Effective surface temperature[K]
:math:`T_{s,i}`                  Surface temperature at refrigerant inlet [K]
:math:`T_{s,o}`                  Surface temperature at refrigerant outlet [K]
:math:`\mathrm{UA}`              Overall surface conductance [W/K]
:math:`\mathrm{UA}_{dry}`        Overall surface conductance for dry analysis [W/K]
:math:`\mathrm{UA}_{wet}`        Overall surface conductance for wet analysis [W/K]
:math:`\mathrm{UA}_i`            Inside thermal conductance [W/K]
:math:`\mathrm{UA}_o`            Outside thermal conductance [W/K]
:math:`\mathrm{UA}_o^*`          Effective outside thermal conductance for wetted coil [W/K]
:math:`\mathrm{UA}_a`            Air-side thermal conductance [W/K]
:math:`\mathrm{UA}_r`            Refrigerant-side thermal conductance [W/K]
:math:`\mathrm{UA}_w`            Water-side thermal conductance [W/K]
:math:`\alpha_a`                 Mean air-side heat transfer coefficient [W/|m2|/K]
:math:`\alpha_w`                 Mean water-side heat transfer coefficient [W/|m2|/K]
:math:`\alpha_r`                 Mean refrigerant-side heat transfer coefficient [W/|m2|/K]
:math:`\alpha_{r,two-phase}`     Mean refrigerant-side heat transfer coefficient in two-phase section [W/|m2|/K]
:math:`\eta_a`                   Air-side surface effectivess (dry analysis) [-]
:math:`\eta_a^*`                 Air-side surface effectivess (wet analysis) [-]
:math:`\varepsilon`              Effectiveness [-]
:math:`\varepsilon_{dry}`        Effectiveness from dry analysis [-]
:math:`\varepsilon_{wet}`        Effectiveness from wet analysis [-]
:math:`\omega`                   Humidity ratio [-]
:math:`\omega_i`                 Humidity ratio at air inlet [-]
===============================  ===================================================

.. automodule:: DryWetSegment 
    :members:
    :undoc-members: