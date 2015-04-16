
.. _Line-Set:

============
Line Set
============

Overview
--------

A line set is essentially just a set of tubes of a reasonably long length that allows for fluid to be moved from one place to another.  It is commonly insulated in order to decrease heat transfer between the tube and ambient.

Mathematical Description
------------------------

.. plot:: MPLPlots/LineSetCV.py

An energy balance over the control volume in the tube gives

.. math::
    :label: eqLS1
    
    \dot m_r c_{p,r}dT=UdA(T_{\infty}-T)
    
where if :math:`T_{\infty}` is greater than :math:`T`, the differential :math:`dT` is positive.  Therefore separation of variables yields

.. math::
    :label: eqLS2
    
    \frac{dT}{T_{\infty}-T}=\frac{UdA}{\dot m_r c_{p,r}}

and upon integration (with a :math:`u`-substitution)

.. math::
    :label: eqLS3

    -\left. \ln (T_{\infty}-T)\right|_{T_i}^{T_o}=\frac{UA}{\dot m_r c_{p,r}}
    
finally yields

.. math::
    :label: eqLS4

    \frac{T_o-T_{\infty}}{T_i-T_{\infty}}=\exp\left(-\frac{UA}{\dot m_r c_{p,r}}\right)
    
and the outlet temperature from

.. math::
    :label: eqLS5
    
    T_o=T_{\infty}+(T_i-T_{\infty})\exp\left(-\frac{UA}{\dot m_r c_{p,r}}\right)
    
heat transfer 

.. math::
    :label: eqLS6
    
    Q=\dot m_r c_{p,r} (T_o-T_i)
        
.. plot:: MPLPlots/LineSetCrossSection.py

Based on the concentric geometry, the heat transfer network terms are given by

.. math::
    :label: eqLS7
    
    R_{tube}=\frac{\ln(D_o/D_i)}{2\pi L k_{tube}}
    
.. math::
    :label: eqLS8
    
    R_{insul}=\frac{\ln[(D_o+2t_{insul})/D_o]}{2\pi L k_{insul}}
    
.. math::
    :label: eqLS9
    
    \mathrm{UA}_i=h_i\pi D_i L
    
    \mathrm{UA}_o=h_o\pi (D_o+2t_{insul}) L

and the overall heat conductance is given by

.. math::
    :label: eqLS10
    
    \mathrm{UA}=\frac{1}{\mathrm{UA}_i^{-1}+R_{tube}+R_{insul}+\mathrm{UA}_o^{-1}}

The pressure drop and refrigerant charge are evaluated from section :ref:`Single-Phase-Fluid-Correlations`.
    
.. |m3| replace:: m\ :sup:`3`\ 
.. |m2| replace:: m\ :sup:`2`\ 

**Nomenclature**

===============================  ===================================================
Variable                         Description
===============================  ===================================================
:math:`c_{p,r}`                  Refrigerant specific heat [J/kg/K]
:math:`D_i`                      Inner diameter of tube [m]
:math:`D_o`                      Outer diameter of tube [m]
:math:`k_{insul}`                Thermal conductivity of insulation [W/m/K]
:math:`L`                        Length of tube [m]
:math:`\dot m_r`                 Refrigerant mass flow rate [kg/s]
:math:`R_{tube}`                 Thermal resistance from tube [K/W]
:math:`R_{insul}`                Thermal resistance from insulation [K/W]
:math:`t_{insul}`                Thickness of insulation [m]
:math:`T_i`                      Inlet temperature [K]
:math:`T_o`                      Outlet temperature [K]
:math:`T_{\infty}`               Ambient temperature [K]
:math:`\mathrm{UA}`              Overall heat transfer conductance [W/K]
:math:`\mathrm{UA_i}`            Inner heat transfer conductance [W/K]
:math:`\mathrm{UA_o}`            Outer heat transfer conductance [W/K]
:math:`\alpha_i`                 Inner heat transfer coefficient [W/|m2|/K]
:math:`\alpha_o`                 Outer heat transfer coefficient [W/|m2|/K]
===============================  ===================================================

Component Class Documentation
-----------------------------

.. py:module:: LineSet    
.. autoclass:: LineSetClass
    :members:
    :undoc-members:
    