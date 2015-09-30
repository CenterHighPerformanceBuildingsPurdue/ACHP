.. |Tevapcond| replace:: :math:`T_{evap}` and :math:`T_{cond}`
 
***************************
Secondary Loop Cycle Solver
***************************

Secondary Loop Cooling Mode
===========================

.. plot:: MPLPlots/CycleSolver/SLACSchematic.py

For secondary loops in cooling mode, there is an internal heat exchanger which physically separates the secondary loop and the refrigerant loop.

In this case, there are now three inputs, and three residuals (to be defined later).  The three inputs are :math:`T_{g,i,cc}`, |Tevapcond|.

The two loops can be solved separately, where for the secondary loop, the inlet temperature to the cooling coil is known, and the secondary working fluid's properties are independent of pressure.

To begin, the :ref:`Cooling-Coil` model is employed to calculate the heat transfer rate in the cooling coil, which gives the process from state point 5 to state point 6.

The :ref:`Line-Set` model is used to determine the heat transfer and pressure drop in the line going from state point 6 to state point 7.

The pump model is run to determine how much electrical power is consumed in the pump, which gives the state point 8, the glycol inlet to the internal heat exchanger.

The refrigerant loop is then solved.  The refrigerant superheat at the outlet of the IHX is imposed as an input for the cycle.  Thus state point 1 is known, and the :ref:`Compressor` model is used to calculate state point 2, the compressor mass flow rate, electrical power, etc..  

The :ref:`Condenser` model is then solved using the state point 2 as the inlet, and yielding the (hopefully) subcooled refrigerant outlet state point 3.  

Finally, the isenthalpic throttling process is used to determine the state point 4 at the refrigerant inlet to the IHX.

Lastly, the :ref:`Plate-Heat-Exchanger` model is run, using the inputs at state points 8 and 4, and yielding the refrigerant outlet state of state point 1'.

The residuals are given by an energy balance between state points 1 and 1', an energy balance on the secondary loop, and either matching the mass or the subcooling on the refrigerant side.  As in the DX system analysis, the residual vector is given by

.. math::
    :label: eqCS11
    
    \vec \Delta=[\vec \Delta_1, \vec \Delta_2, \vec \Delta_3]

with

.. math::
    :label: eqCS12

    \vec \Delta_1=\dot m_r(h_1-h_{1'})

    \vec \Delta_2=\left \lbrace \begin{array}{cc} m_r-m_{r,target} & (\mbox{Charge imposed}) \\ \Delta T_{sc,cond}-\Delta T_{sc,cond,target} & (\mbox{Subcooling imposed})\end{array} \right .
    
    \vec \Delta_3=\dot W_{pump}+\dot Q_{cc}-\dot Q_{IHX}+\dot Q_{\mathrm{suppply line}}+\dot Q_{\mathrm{return line}}

The set of :math:`T_{g,i,cc}`, |Tevapcond| which solve the residual equations are obtained by a multi-dimensional solver (see :ref:`Numerical-Methods-NDsolve`).

The common metrics of system efficiency are 

.. math::
    :label: eqCS7

    COP=\frac{\dot Q_{cc}}{\dot W_{comp}}
    
    COSP=\frac{\dot Q_{cc}-\dot W_{fan,cc}}{\dot W_{comp}+\dot W_{fan,cc}+\dot W_{fan,cond}+\dot W_{pump}}
    
.. literalinclude:: SampleCycles/SampleSLAC.py
    
which should yield the output, when run, of

.. ipython::
    :suppress:
    
    #This line will take the output from the script and write it to file
    In [1]: import sys; old=sys.stdout; f=open('ACHPModel/SampleCycles/SLACOutput.txt','w'); sys.stdout=f;execfile('ACHPModel/SampleCycles/SampleSLAC.py'); sys.stdout=old; f.close()
    
.. literalinclude:: SampleCycles/SLACOutput.txt

Secondary Loop Heating Mode
===========================
This section is left intentionally empty

Cycle Solver Class Documentation
================================

.. py:module:: ACHP.Cycle    
.. autoclass:: SecondaryCycleClass
    :members:
    :undoc-members:
    
    
.. |m3| replace:: m\ :sup:`3`\ 
.. |m2| replace:: m\ :sup:`2`\ 
.. |DTsctarg| replace:: :math:`\Delta T_{sc,cond,target}`

**Nomenclature**

===============================  ===================================================
Variable                         Description
===============================  ===================================================
:math:`COP`                      Coefficient of Performance [-]
:math:`COSP`                     Coefficient of System Performance [-]
:math:`h_1`                      Enthalpy at outlet of evaporator [J/kg]
:math:`h_{1'}`                   Enthalpy after going around the cycle [J/kg]
:math:`\dot m_r`                 Refrigerant mass flow rate [kg/s]
:math:`m_r`                      Model-predicted charge [kg]
:math:`m_{r,target}`             Target refrigerant charge in system [kg]
:math:`T_{cond}`                 Condensing (dewpoint) temperature [K]
:math:`T_{evap}`                 Evaporating (dewpoint) temperature [K]
:math:`T_{g,i,cc}`               Glycol inlet temperature to cooling coil [K]
:math:`\Delta T_{sh}`            IHX outlet superheat [K]
:math:`\Delta T_{sc,cond}`       Condenser outlet subcooling [K]
|DTsctarg|                       Condenser outlet subcooling target [K]
:math:`\dot Q_{cc}`              Cooling coil heat transfer rate [W]
:math:`\dot W_{fan,cc}`          Cooling coil fan power [W]
:math:`\dot W_{fan,cond}`        Condenser fan power [W]
:math:`\dot W_{comp}`            Compressor power input [W]
:math:`\vec \Delta`              Residual vector [varied]
===============================  ===================================================
