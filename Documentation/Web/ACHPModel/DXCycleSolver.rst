.. |Tevapcond| replace:: :math:`T_{evap}` and :math:`T_{cond}`

.. _DX-Cycle-Solver:

*****************************
Direct Expansion Cycle Solver
*****************************

In the cycle solver, the goal is to use all the physical component models together in order to obtain the performance of the combined cycle.  The process is slightly different for each of the 4 standard configurations of ACHP, but there are many common themes among the configurations.  ACHP has been designed to investigate direct expansion (DX) systems as well as systems employing a secondary working loop, in both cooling and heating modes.

A number of simplifying assumptions are employed at the cycle level:

* Pressure drops in each component are calculated, but not used directly(see :ref:`Pressure-Drop-Correction` for explanation)
* There is no charge inventory in the compressor shell
* The evaporator outlet refrigerant superheat is imposed (implies an expansion device with perfect superheat control)
* No pressure drops except for in the components considered

.. _Direct-Expansion-AC:

Direct Expansion Cooling Mode
=============================
The most common and straightforward system is the direct expansion cooling mode system.  The schematic of this cycle is shown in this figure:

.. plot:: MPLPlots/CycleSolver/DXACSchematic.py

In addition to the geometry of each component, required inputs are the superheat at the outlet of the evaporator :math:`\Delta T_{sh,evap}`, and either the total refrigerant charge :math:`m_r` or the refrigerant subcooling at the outlet of the condenser :math:`\Delta T_{sc,cond}`.

The primary independent variables in the cycle solver are the saturation temperatures of the refrigerant at evaporation and condensation, given by  :math:`T_{evap}` and :math:`T_{cond}` respectively.  For pure refrigerants the saturation temperature for a given pressure is the same for saturated liquid and saturated vapor.  For pseudo-pure fluids and blends, the refrigerant dew temperature is used which corresponds to saturated vapor.  

Before the full cycle-solver is run, the preconditioner in section :ref:`Cycle-Solver-Preconditioner` is used to get a good guess for the temperatures :math:`T_{evap}` and :math:`T_{cond}` using extremely simple cycle models.

Once the preconditioner has been run, preliminary values for :math:`T_{evap}` and :math:`T_{cond}` are known, and the first iteration of the cycle model may begin.  Execution of the cycle model follows in the same direction as the refrigerant flow.  The pressure drops through low-pressure and high-pressure parts of the system are assumed to be zero initially.

The cycle analysis begins with the vapor line which returns refrigerant vapor from the evaporator that is typically in the air duct back to the condensing unit outdoors.  In the first iteration of the cycle solver, the mass flow rate of refrigerant in the vapor line is not known, but the compressor map is used to provide a reasonable guess value for the mass flow rate of refrigerant.  The :ref:`Line-Set` model is used to calculate the process from point 1 to point 2.  The model is run with known flow rates and an inlet temperature of :math:`T_{evap}+\Delta T_{sh}` in order to calculate the state at the inlet to the compressor.

The compressor compresses refrigerant from state point 2 to state point 3.  The :ref:`Compressor` model is used which is based on an empirical correlation with superheat correction and it yields the outlet state 3 as well as the compressor mass flow rate, compressor electrical power, etc..  The compressor model requires the temperatures :math:`T_{evap}` and :math:`T_{cond}` as well as the compressor inlet superheat.

The condenser takes the superheated refrigerant at state point 3 and condenses it to a subcooled refrigerant at state point 4.  The :ref:`Condenser` model is used to calculate the process, and this condenser model is based on a moving-boundary model.  The condenser model requires :math:`T_{cond}` as an input, among others.

After the condenser, the subcooled refrigerant passes through the liquid-line that takes refrigerant from the condenser to the indoor coil that is inside the ductwork in the home.  The :ref:`Line-Set` model is used to model the flow that passes from state point 4 to state point 5.

The expansion device then expands the refrigerant from the high-pressure side of the system to the low-pressure side of the system.  In the current iteration of ACHP, no expansion device is included.  Thus the expansion device is just a constant-enthalpy throttling device that takes refrigerant from the high-side pressure at state point 5 to the low-side pressure at state point 6.

At the outlet of the expansion device, the refrigerant passes into the evaporator at some two-phase quality.  The :ref:`Evaporator` model is used to model the performance of the evaporator.  In the evaporator, refrigerant is heated by the air-stream (which cools the air stream) from state point 6 back to state point 1', which should be superheated (note: not necessarily so at intermediate iterations of cycle model) close to state point 1.  If the values of the cycle independent variables :math:`T_{evap}` and :math:`T_{cond}` have been well selected, the state points 1 and 1' will be coincident.

Once the cycle model has been run around the loop from state point 1 to state point 1', the residuals can be calculated.  The residuals are terms that when the cycle model has converged should all be equal to zero.

One of the residuals that is always active is an energy balance over the cycle.  You left from state point 1, and you should hopefully arrive back there if energy is conserved in the cycle.  Thus, the first residual is given by

.. math::
    :label: eqCS1

    \vec \Delta_1=\dot m_r(h_1-h_{1'})

where :math:`h_1` and :math:`h_{1'}` are the enthalpies of the refrigerant at state points 1 and 1' respectively.

Since there are two independent variables, there must be a second constraint, which in this case is a charge-level constraint.  Either the charge is constrained directly with an imposed charge level, or indirectly with an imposed refrigerant subcooling.  Thus the second residual can be given by

.. math::
    :label: eqCS2

    \vec \Delta_2=\left \lbrace \begin{array}{cc} m_r-m_{r,target} & (\mbox{Charge imposed}) \\ \Delta T_{sc,cond}-\Delta T_{sc,cond,target} & (\mbox{Subcooling imposed})\end{array} \right .

and a numerical solver is used to drive the residual vector :math:`\vec \Delta` to sufficiently close to zero by altering :math:`T_{evap}` and :math:`T_{cond}` (see :ref:`Numerical-Methods-NDsolve`).

.. _Pressure-Drop-Correction:

Pressure-Drop Correction
------------------------
After the cycle model has iterated to convergence, the pressure drops are then considered.  In reality, the pressure drop in each component results in a lower pressure at the inlet of the next component in the refrigerant loop.  In terms of modeling, coupling pressure drop and the component models causes great numerical difficulties.  The compromise that is employed instead is to run all the component models without pressure drop, but calculate the high-side pressure drop as the pressure drop of the condenser and liquid-line:

.. math::
    :label: eqCS3

    \Delta p_{high}=\Delta p_{cond}+\Delta p_{liquid-line}
    
and similarly, the low-side pressure drop is defined by

.. math::
    :label: eqCS4

    \Delta p_{low}=\Delta p_{evap}+\Delta p_{vapor-line}
    
These pressure drops are then employed to shift the saturation temperatures used in the compressor map in order to yield less refrigerant mass flow rate, and a higher compressor power.

The new effective compressor suction and discharge pressures are

.. math::
    :label: eqCS5

    p_{evap}^*=p(T_{evap})-\Delta p_{low}

    p_{cond}^*=p(T_{cond})+\Delta p_{high}
    
and the new, effective compressor dew temperatures are

.. math::
    :label: eqCS6

    T_{evap}^*=T(p=p_{evap},x=1)
    
    T_{cond}^*=T(p=p_{cond},x=1)
    
This calculated pressure drop is used until the model reaches convergence again, at which point the pressure drop tems are updated, and the model is run again.  This process continues until the imposed low- and high-side pressure drops are equal to the pressure drop terms calculated from the converged cycle model.

Cycle Performance Parameters
----------------------------

The common metrics of system efficiency are 

.. math::
    :label: eqCS7

    COP=\frac{\dot Q_{evap}}{\dot W_{comp}}
    
    COSP=\frac{\dot Q_{evap}-\dot W_{fan,evap}}{\dot W_{comp}+\dot W_{fan,evap}+\dot W_{fan,cond}}


.. literalinclude:: SampleCycles/SampleDXAC.py
    
which should yield the output, when run, of

.. ipython::
    :suppress: 
    
    #This line will take the output from the script and write it to file
    In [1]: import sys; old=sys.stdout; f=open('ACHPModel/SampleCycles/DXACOutput.txt','w'); sys.stdout=f;execfile('ACHPModel/SampleCycles/SampleDXAC.py'); sys.stdout=old; f.close()
    
.. literalinclude:: SampleCycles/DXACOutput.txt
     
Direct Expansion Heating Mode
=============================
The heat pump configuration of the system is as shown here:

.. plot:: MPLPlots/CycleSolver/DXHPSchematic.py

Physically, reversing valves are used to switch the mode of the system and the directions of the flows.  What was the condenser of the air conditioning system becomes the evaporator of the heat pump and *vice versa*, and the line sets are configured in a slightly different way.  Other than that, the analysis of the heat pump is directly analogous to that of the :ref:`Direct-Expansion-AC` system.

A :ref:`Preconditioner <Cycle-Solver-Preconditioner>` is used to get approximate values for |Tevapcond|, and using these values (which are iteratively modified using numerical methods), the solution for the cycle performance is found.

As with the cooling mode, the cycle analysis follows the refrigerant flow path around the loop.  

Beginning at the outlet of the evaporator, state point 1 is known because :math:`T_{evap}` and :math:`\Delta T_{sh}` are known.  Thus the compressor model is used directly to calculate the electrical power, refrigerant mass flow rate and state point 2.

The :ref:`Line-Set` model is then applied to the flow from the outlet of the compressor at state point 2 to the inlet of the condenser at state point 3.

The :ref:`Condenser` model is used to model the condensing process from state point 3 to a subcooled state at state point 4.  

The :ref:`Line-Set` model is used to model the flow of subcooled refrigerant at state point 4 back to the expansion device at state point 5.

As in cooling mode, the expansion device is assumed to be an ideal expansion device, which means that the working process is a constant-enthalpy expansion from state point 5 to state point 6.

The :ref:`Evaporator` model is then used to calculate the evaporation process of refrigerant from state point 6 to state point 1'.

As in cooling mode, the residual vector is given by

.. math::
    :label: eqCS8
    
    \vec \Delta=[\vec \Delta_1, \vec \Delta_2]

with

.. math::
    :label: eqCS9

    \vec \Delta_1=\dot m_r(h_1-h_{1'})

    \vec \Delta_2=\left \lbrace \begin{array}{cc} m_r-m_{r,target} & (\mbox{Charge imposed}) \\ \Delta T_{sc}-\Delta T_{sc,target} & (\mbox{Subcooling imposed})\end{array} \right .
    
The set of |Tevapcond| which solve the residual equations are obtained by a multi-dimensional solver (see :ref:`Numerical-Methods-NDsolve`).

Cycle Performance Parameters
----------------------------

The common metrics of system efficiency are 

.. math::
    :label: eq CS10

    COP=\frac{\dot Q_{cond}}{\dot W_{comp}}
    
    COSP=\frac{\dot Q_{cond}+\dot W_{fan,cond}}{\dot W_{comp}+\dot W_{fan,evap}+\dot W_{fan,cond}}
 
Minimal working example
-----------------------

.. literalinclude:: SampleCycles/SampleDXHP.py
    
which should yield the output, when run, of

.. ipython::
    :suppress:
    
    #This line will take the output from the script and write it to file
    In [1]: import sys; old=sys.stdout; f=open('ACHPModel/SampleCycles/DXHPOutput.txt','w'); sys.stdout=f;execfile('ACHPModel/SampleCycles/SampleDXHP.py'); sys.stdout=old; f.close()
    
.. literalinclude:: SampleCycles/DXHPOutput.txt

Cycle Solver Code Documentation
-------------------------------

.. py:module:: Cycle    
.. autoclass:: DXCycleClass
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
:math:`T_{evap}`                 Evaporating (dewpoint) temperature [K]
:math:`T_{cond}`                 Condensing (dewpoint) temperature [K]
:math:`T_{evap}^*`               Effective evaporating temperature [K]
:math:`T_{cond}^*`               Effective condensing temperature [K]
:math:`\Delta T_{sh}`            Evaporator outlet superheat [K]
:math:`\Delta T_{sc,cond}`       Condenser outlet subcooling [K]
|DTsctarg|                       Condenser outlet subcooling target [K]
:math:`\Delta p_{cond}`          Pressure drop in condenser [kPa]
:math:`\Delta p_{evap}`          Pressure drop in evaporator [kPa]
:math:`\Delta p_{high}`          Pressure drop on high-pressure side of system [kPa]
:math:`\Delta p_{liquid-line}`   Pressure drop in liquid line [kPa]
:math:`\Delta p_{low}`           Pressure drop on low-pressure side of system [kPa]
:math:`\Delta p_{vapor-line}`    Pressure drop in vapor line [kPa]
:math:`p_{evap}^*`               Effective evaporation saturation pressure [kPa]
:math:`p_{cond}^*`               Effective condensing saturation pressure [kPa]
:math:`\dot Q_{evap}`            Evaporator heat transfer rate [W]
:math:`\dot W_{fan,evap}`        Evaporator fan power [W]
:math:`\dot W_{fan,cond}`        Condenser fan power [W]
:math:`\dot W_{comp}`            Compressor power input [W]
:math:`\vec \Delta_1`            Residual vector [varied]
===============================  ===================================================