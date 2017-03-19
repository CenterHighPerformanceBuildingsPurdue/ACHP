.. _Evaporator:

Evaporator
==========

Overview
------------

In the evaporator, refrigerant enters at some vapor quality between 0% and 100%, and typically exits as a superheated vapor. The analysis for the evaporator shares many features with the analysis employed for the :ref:`Condenser` and the :ref:`Cooling-Coil`, but the evaporator's analysis is significantly more complicated due to the fact that there are moving boundaries on both the refrigerant and air sides.  On the refrigerant-side, the moving boundary is between two-phase refrigerant and superheated refrigerant, and on the air-side, there is a moving boundary between wet and dry parts of the coil.

Mathematical Description
------------------------

The evaporator is treated as being cross-counter-flow.  What this means practically is that the mass flow of air to and total air side surface area for each of the superheated and two-phase portions (assuming each exists), are given by 

.. math::
    :label: eqEvap1

    A_{a,superheat}=w_{superheat}A_{a,total}
    
    A_{a,two-phase}=w_{two-phase}A_{a,total}
    
    \dot m_{a,superheat}=w_{superheat}\dot m_{a,total}
    
    \dot m_{a,two-phase}=w_{two-phase}\dot m_{a,total}

which is analogous to the analysis employed for the condenser.  On the refrigerant-side, the surface area for superheated and two-phase sections are given by

.. math::
    :label: eqEvap1
    
    A_{r,superheat}=w_{superheat}A_{r,total}
    
    A_{r,two-phase}=w_{two-phase}A_{r,total}
    
Two-Phase Section
-----------------
In the two-phase section, the target heat transfer rate is given by

.. math::
    :label: eqEvap1
    
    \dot Q_{target}=\dot m_r(x_{r,o}-x_{r,i})h_{fg}
    
and the heat transfer rate from the partially-wet/partially-dry analysis in section :ref:`PWPD-Two-Phase` is given the name :math:`\dot Q_{PWPD}`.  To begin with, the first guess is that all the heat exchanger is in the two-phase region.  If using the outlet quality of saturated vapor (:math:`x_{r,o}=1`) with all the heat exchanger in the two-phase region, and :math:`\dot Q_{PWPD}` is greater than :math:`\dot Q_{target}`, too much of the heat exchanger area was given to the two-phase section, and there must be a superheated section as well.  Thus, the length fraction :math:`w_{two-phase}` is iteratively altered using a bounded solver (:math:`w_{two-phase}` is bounded between 0 and 1).  The remaining area is given to the superheated section, so the superheated section circuit length fraction :math:`w_{superheat}` is given by

.. math::
    :label: eqEvap1
    
    w_{superheat}=1-w_{two-phase}

Superheated Section
-------------------

In the superheated section, for a given :math:`w_{superheat}`, the analysis is exactly the same as the :ref:`Cooling-Coil`, and given by the section :ref:`PWPD-Single-Phase`.  As in the cooling coil, the working fluid is single-phase, and the air-side surface may be fully wet, partially wet, or fully dry, and all the inlet conditions for the superheated section are known.  The air inlet state to the superheated portion is assumed to be the same as for the two-phase portion because the superheated portion is typically rather small, and both superheated and two-phase portions should see approximately the same inlet air state.  The refrigerant properties are calculated using the same single-phase correlations as for the cooling coil.

Minimal Example Code
--------------------

.. literalinclude:: ComponentTests/EvaporatorTest.py

If you open an IPython shell in the root of the documentation (folder ``Documentation/Web`` relative to the main trunk), and run the commands below, you should get output something like

.. ipython::

    In [1]: %run 'ACHPComponents/ComponentTests/EvaporatorTest.py'
    
If not, first stop should be the :ref:`FAQS`

**Nomenclature**

===============================  ===================================================
Variable                         Description
===============================  ===================================================
:math:`A_{a,total}`              Total air side area (fins+tubes) [m\ :sup:`2`\ ]
:math:`A_{a,superheat}`          Air-side area for superheated portion of circuit [m\ :sup:`2`\ ]
:math:`A_{a,two-phase}`          Air-side area for two-phase portion of circuit [m\ :sup:`2`\ ]
:math:`A_{r,superheat}`          Refrigerant-side area for superheated portion of circuit [m\ :sup:`2`\ ]
:math:`A_{r,two-phase}`          Refrigerant-side area for two-phase portion of circuit [m\ :sup:`2`\ ]
:math:`h_{fg}`                   Latent heat of refrigerant [J/kg]
:math:`\dot m_{a,total}`         Mass flow rate of dry air [kg/s]
:math:`\dot m_{a,superheat}`     Mass flow rate of dry air in superheated section [kg/s]
:math:`\dot m_{a,two-phase}`     Mass flow rate of dry air in two-phase section [kg/s]
:math:`\dot m_r`                 Mass flow rate of refrigerant [kg/s]
:math:`\dot Q_{target}`          Target heat transfer rate [W]
:math:`\dot Q_{PWPD}`            Heat transfer rate from the partially wet/partially-dry analysis [W]
:math:`x_{r,i}`                  Refrigerant inlet quality [-]
:math:`x_{r,o}`                  Refrigerant outlet quality [-]
===============================  ===================================================

Component Class Documentation
-----------------------------

.. py:module:: ACHP.Evaporator    
.. autoclass:: EvaporatorClass
    :members:
    :undoc-members:
