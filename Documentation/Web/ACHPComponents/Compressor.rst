.. _Compressor:

============
Compressor
============

Overview
--------
A compressor is the heart of the air-conditioning or refrigeration system, compressing the refrigerant from the low-pressure side of the system up to the high-pressure side of the system.  Because it plays such a a significant role in the overall system efficiency, the accuracy of the compressor model is quite important.  


Mathematical Description
------------------------

The compressor is modeled based on a 10-coefficient ARI compressor map which is very commonly used to characterize the performance of compressors.  The map is based on a given amount of superheat along with input saturated suction and discharge pressures.  Though most everything else in the program is based on metric units, the standard in America is to generate the map based on superheat and saturated temperatures in degrees Fahrenheit.  Refer to the following figure to see the definitions of the temperatures and pressures:

.. plot::

    from CoolProp.CoolProp import Props
    from CoolProp.Plots import Ph
    import pylab
    Ph('R410A')
    ps=Props('P','T',260,'Q',1,'R410A')
    pd=Props('P','T',300,'Q',1,'R410A')
    hs=Props('H','T',260,'Q',1,'R410A')
    hd=Props('H','T',300,'Q',1,'R410A')
    pylab.plot(hs,ps,'bo',mfc='b')
    pylab.plot(hd,pd,'bo',mfc='b')
    pylab.gca().axhline(ps)
    pylab.gca().axhline(pd)
    pylab.gca().text(450,ps+30,'$p_s$',va='bottom',size=12)
    pylab.gca().text(450,pd+30,'$p_d$',va='bottom',size=12)

    pylab.gca().text(hs+1,ps-180,'$T_s$')
    pylab.gca().text(hd+1,pd+40,'$T_d$')
    pylab.gca().set_xlim(300,500)
    pylab.gca().set_ylim(0,6000)
    
Thus the map-based mass flow rate (in lbm/hr) and electrical power (in W), can be given by

.. math::
    :label: eqComp1

    \dot {m}_{map} = M_{1}+M_{2}T_s+M_{3}T_{d}+M_{4}T_{s}^2+M_{5}T_{s}T_d+M_{6}T_{d}^2+M_{7}T_{s}^3+M_{8}T_{d}T_{s}^2+M_{9}T_{d}^{2}T_{s}+M_{10}T_{d}^3

.. math::
    :label: eqComp2
    
    \dot W_{map} = P_{1}+P_{2}T_s+P_{3}T_{d}+P_{4}T_{s}^2+P_{5}T_{s}T_d+P_{6}T_{d}^2+P_{7}T_{s}^3+P_{8}T_{d}T_{s}^2+P_{9}T_{d}^{2}T_{s}+P_{10}T_{d}^3

where the saturated suction dewpoint temperature :math:`T_s` and saturated discharge dewpoint temperatures :math:`T_d` are in degrees Fahrenheit.  To be specific, they are the dew temperatures, which are the same as the saturated vapor temperatures for pure fluids.  The coefficients :math:`M_1, M_2, ...` are the mass flow map coefficients and :math:`P_1, P_2, ....` are the electrical power map coefficients.  In practice, the compressor is unlikely to operate at exactly the map superheat.  As a result, the map predictions must be corrected to better match the actual operating conditions.  The map correction is based on the method of Rice. et al. [#Rice]_ , which yields
   
.. math::
    :label: eqComp3
    
    \dot m_{actual} = \left[1 + 0.75 \left(\frac{v_{map}}{ v_{actual}} - 1\right ) \right ]  \dot m_{map}

where the subscripts *actual* refer to the properties evaluated at the actual superheat at the suction flange, and the *map* subscripts refer to the properties evaluated at the given map superheat.  Similarly, the electrical power correction is given by

.. math::
    :label: eqComp4
    
    \dot W_{actual} = \dot W_{map}  \frac{\dot m_{actual}}{ \dot m_{map}} \frac{ h_{2s,actual} - h_{1,actual}}{h_{2s,map} - h_{1,map}}

An energy balance over the compressor yields

.. math::
    :label: eqComp5
    
    \dot W_{actual}+\dot Q_{amb}+\dot m_{actual}(h_{1,actual}-h_{2,actual})=0


Since the electrical power is known from the corrected map, and the heat transfer can be expressed as a fraction of the electrical power by

.. math ::
    :label: eqComp6
    
    \dot Q_{amb}=-f_p\dot W_{actual}

the outlet enthalpy of the compressor can be therefore given by

.. math::
    :label: eqComp7
    
    h_{2,actual}=\frac{\dot W_{actual}(1-f_p)}{\dot m_{actual}}+h_{1,actual}
    
.. only:: html

    .. rubric:: References

.. [#Rice] Rice, C. K. and A. E. Dabiri, 1981. "A Compressor Simulation Model with Corrections for the Level of Suction Gas Superheat," ASHRAE Transactions, Vol. 87, Part 2, pp.771-782.

Nomenclature

.. |degF| replace:: :math:`^{\circ}\mathrm{F}`

===============================  ===================================================
Variable                         Description
===============================  ===================================================
:math:`f_p`                      Fraction of electrical power lost at heat transfer [-]
:math:`h_{1,actual}`             Enthalpy of refrigerant at actual superheat [J/kg]
:math:`h_{2s,actual}`            Isentropic enthalpy of refrigerant at discharge pressure using actual superheat [J/kg]
:math:`h_{1,map}`                Enthalpy of refrigerant at map superheat [J/kg]
:math:`h_{2s,map}`               Isentropic enthalpy of refrigerant at discharge pressure using map superheat [J/kg]
:math:`\dot m_{actual}`          Actual refrigerant mass flow rate [kg/s]
:math:`\dot m_{map}`             Refrigerant mass flow rate from map [lb\ :sub:`m`\ /hr]
:math:`M_1,M_2,...`              Mass flow map coefficients [varied]
:math:`P_1,P_2,...`              Electrical power map coefficients [varied]
:math:`p_d`                      Discharge dew pressure [Pa (absolute)]
:math:`p_s`                      Suction dew pressure [Pa (absolute)]
:math:`\dot Q_{amb}`             Ambient heat loss [W]
:math:`T_d`                      Discharge dew temperature [|degF|]
:math:`T_s`                      Suction dew temperature [|degF|]
:math:`v_{actual}`               Specific volume of refrigerant at actual superheat [m\ :superscript:`3`\ /kg]
:math:`v_{map}`                  Specific volume of refrigerant at map superheat [m\ :superscript:`3`\ /kg]
:math:`\dot W_{actual}`          Actual compressor electrical power [W]
:math:`\dot W_{map}`             Compressor electrical power from map [W]
===============================  ===================================================

Compressor Sample Code
----------------------
Minimal Component Test:

.. literalinclude:: ComponentTests/CompressorTest.py
    :language: python
    
If you open an IPython(x,y) shell in the root of the documentation (folder Documentation/Web relative to the main trunk), and run the commands below, you should get

.. ipython::

    In [1]: execfile('ACHPComponents/ComponentTests/CompressorTest.py')
    
If not, first stop should be the :ref:`FAQS`

Component Class Documentation
-----------------------------

.. py:module:: Compressor    
.. autoclass:: CompressorClass
    :members:
    :undoc-members:
    