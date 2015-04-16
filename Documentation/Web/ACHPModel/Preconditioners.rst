.. _Cycle-Solver-Preconditioner:

Cycle Solver Preconditioner
***************************

The basic idea of a preconditioner is to use an extremely simple model in order to obtain good initial values for the solver used for the cycle in ACHP.  

Conventional system
===================

.. plot:: MPLPlots/Preconditioner/PreconditionerFlowChart.py

The following assumptions are employed in the preconditioner employed for the conventional heat pump and air conditioning models:

* Condenser air stream is the minimum capacitance rate in the condenser, and the limiting outlet state for the air in the condenser is :math:`T_{cond}`.  [This further assumes that there is no subcooled section, but this is ok, because there is also a pinch point at the saturated liquid point, so the basic analysis still works.  In the end, the preconditioner is only used to get an *approximate* solution anyway, so inaccuracy is acceptable, as long as it is reasonable.]
* Effectiveness of all heat exchangers is known, and fixed
* Compressor inlet superheat is fixed and known
* Compressor is adiabatic
* Evaporator is fully wet, fully dry, or a simple weighted mix of the two
* Line sets are not considered

Essentially the preconditioner operates with the same independent variables as ACHP, which are the dew temperatures of the refrigerant in the evaporator and condenser, given by the variables :math:`T_{evap}` and :math:`T_{cond}` respectively.

Compressor
----------

For a given set of :math:`T_{evap}`, :math:`T_{cond}`, and known superheat :math:`\Delta T_{sh}`, the inlet state and outlet pressure for the compressor is known.  Therefore, the compressor map can be used to predict the refrigerant mass flow rate as well as the compressor power.  This yields the values

.. math::

    \dot m_r=f_{map}(T_{evap},T_{cond},\Delta T_{sh})
    
    \dot W_{comp}=f_{map}(T_{evap},T_{cond},\Delta T_{sh})
    
Condenser 
---------

The heat transfer rate in the condenser can therefore be given by the value

.. math::

    \dot Q_{cond}=\varepsilon_{HX} \rho_{ha} \dot V_{ha,cond} c_{p,a}(T_{i,a,cond}-T_{cond})
    
The condenser heat transfer rate based on the imposed subcooling can also be given by

.. math::

    \dot Q_{cond,\Delta h}=\dot m_r (h_{r,o,comp}-h(T_{cond}-\Delta T_{sc},p_{cond}))
    
Theoretically these two heat transfer rate terms should match, and consistency is imposed by the numerical solver that is employed.

Evaporator
----------

As usual, the evaporator is the most complicated component due to the possibility of the evaporator coil being fully wet, fully dry, or partially wet and partially dry.  As in the full evaporator model, the evaporator is first considered to be fully dry, yielding the heat transfer rate of

.. math::

    \dot Q_{evap,dry}=\varepsilon_{HX} \rho_{ha} \dot V_{ha,evap} c_{p,a}(T_{i,a,evap}-T_{evap})
    
Then using the dry evaporator heat transfer analysis it is possible to determine the surface temperature.  The :math:`\mathrm{UA}` values can be obtained from

.. math::

    \mathrm{UA}_r = \alpha_r A_{r,total}

    \mathrm{UA}_a = \eta_a \alpha_a A_{a,total}
    
The outlet temperature of the air can be given from

.. math::

    T_{o,a,evap}=T_{i,a,evap}-\frac{\dot Q_{evap,dry}}{\dot m_{a,total}c_{p,a}}
    
which yields the air inlet surface temperature of

.. math::

    T_{s,a,i}=\frac{\mathrm{UA}_aT_{i,a,evap}+\mathrm{UA}_rT_{evap}}{\mathrm{UA}_a+\mathrm{UA}_r}
    
and the air outlet surface temperature of

.. math::

    T_{s,a,o}=\frac{\mathrm{UA}_aT_{o,a,evap}+\mathrm{UA}_rT_{evap}}{\mathrm{UA}_a+\mathrm{UA}_r}
    
If both :math:`T_{s,a,o}` and :math:`T_{s,a,i}` are above the dewpoint temperature of the entering air (:math:`T_{dp}`), the rate of heat transfer in the evaporator is equal to the dry-analysis heat transfer rate.  If both :math:`T_{s,a,o}` and :math:`T_{s,a,i}` are below the dewpoint of the entering air, the coil is entirely wet, for which the heat transfer rate can be obtained from

.. math::

    \dot Q_{evap,wet}=\varepsilon_{HX} \rho_{ha} \dot V_{ha,evap}(h_{a,i}-h_{a,s,evap})

where :math:`h_{a,s,evap}` is the saturated air enthalpy at :math:`T_{evap}` and :math:`h_{a,i}` is the enthalpy of the inlet air to the evaporator.

If the dewpoint of the inlet air is somewhere between :math:`T_{s,a,o}` and :math:`T_{s,a,i}`, the heat transfer rate in the evaporator is given by a simple weighting.  This yields the following solution for the evaporator:

.. math::

    \dot Q_{evap}=\left\lbrace \begin{array}{cc}\dot Q_{evap,dry} & T_{s,a,o} > T_{dp} \\ \dot Q_{evap,wet} & T_{s,a,i} < T_{dp} \\ \frac{T_{s,a,o}-T_{dp}}{T_{s,a,o}-T_{s,a,i}}\dot Q_{evap,wet}+ \left(1-\frac{T_{s,a,o}-T_{dp}}{T_{s,a,o}-T_{s,a,i}}\right)\dot Q_{evap,dry}& T_{s,a,i}> T_{dp} > T_{s,a,o}\end{array}\right.
    
Solution Method
---------------

The residuals to driven to zero are therefore an overall energy balance over the system, as well as matching :math:`\dot Q_{cond}` and :math:`\dot Q_{cond,\Delta h}`.  So the residual vector as a function of  :math:`T_{evap}` and :math:`T_{cond}` can be expressed as

.. math::

    \vec{\Delta}(T_{evap},T_{cond})=\left[ \begin{array}{c} \dot Q_{evap}+\dot W_{comp}+\dot Q_{cond} \\ \dot Q_{cond}+\dot Q_{cond,\Delta h} \end{array}  \right]
    
and a two-dimensional solver can be used to drive the norm of :math:`\vec{\Delta}` to sufficiently close to zero by altering :math:`T_{evap}` and :math:`T_{cond}`.

Heating Mode
------------
In heating mode, the system schematic remains exactly the same, and the same analysis is used, but the physical geometry of the evaporator and condenser are swapped.

Secondary Loop Systems
======================

.. plot:: MPLPlots/Preconditioner/PreconditionerFlowChartSL.py

The same basic structure is employed for the preconditioner for the secondary loop systems, except that one more variable must be determined by the preconditioner.  The preconditioner for the secondary loop system is used to determine the saturation temperatures :math:`T_{evap}` and :math:`T_{cond}`, as well as the cooling coil inlet temperature :math:`T_{g,i,cc}`.

The same exact analysis as for the DX preconditioner is employed for the compressor and condenser, and a very similar analysis is used for the cooling coil.  The cooling coil analysis mirrors that of the evaporator, as described here.

Cooling Coil
------------

The cooling coil is first considered to be fully dry, yielding the heat transfer rate of

.. math::

    \dot Q_{cc,dry}=\varepsilon_{HX} \rho_{ha} \dot V_{ha,cc} c_{p,a}(T_{a,i,cc}-T_{g,i,cc})
    
Then using the dry cooling coil heat transfer analysis it is possible to determine the surface temperature.  The :math:`\mathrm{UA}` values can be obtained from

.. math::

    \mathrm{UA}_g = \alpha_g A_g

    \mathrm{UA}_a = \eta_a \alpha_a A_a
    
The outlet temperature of the air can be given from

.. math::

    T_{a,o,cc}=T_{a,i,cc}-\frac{\dot Q_{cc,dry}}{\dot m_ac_{p,a}}
    
    T_{g,o,cc}=T_{g,i,cc}+\frac{\dot Q_{cc,dry}}{\dot m_g c_{p,g}}
    
which yields the air inlet surface temperature of

.. math::

    T_{s,a,i}=\frac{\mathrm{UA}_aT_{a,i,cc}+\mathrm{UA}_gT_{g,i,cc}}{\mathrm{UA}_a+\mathrm{UA}_g}
    
and the air outlet surface temperature of

.. math::

    T_{s,a,o}=\frac{\mathrm{UA}_aT_{a,o,cc}+\mathrm{UA}_gT_{g,o,cc}}{\mathrm{UA}_a+\mathrm{UA}_g}
    
If both :math:`T_{s,a,o}` and :math:`T_{s,a,i}` are above the dewpoint temperature of the entering air (:math:`T_{dp}`), the rate of heat transfer in the evaporator is equal to the dry-analysis heat transfer rate.  If both :math:`T_{s,a,o}` and :math:`T_{s,a,i}` are below the dewpoint of the entering air, the coil is entirely wet, for which the heat transfer rate can be obtained from

.. math::

    \dot Q_{cc,wet}=\varepsilon_{HX} \rho_{ha} \dot V_{ha,cc}(h_{a,i}-h_{a,s,cc})

where :math:`h_{a,s,cc}` is the saturated air enthalpy at :math:`T_{g,i}` and :math:`h_{a,i}` is the enthalpy of the inlet air to the cooling coil.

If the dewpoint of the inlet air is somewhere between :math:`T_{s,a,o}` and :math:`T_{s,a,i}`, the heat transfer rate in the cooling coil is given by a simple weighting.  This yields the following solution for the cooling coil:

.. math::

    \dot Q_{cc}=\left\lbrace \begin{array}{cc}\dot Q_{cc,dry} & T_{s,a,o} > T_{dp} \\ \dot Q_{cc,wet} & T_{s,a,i} < T_{dp} \\ \frac{T_{s,a,o}-T_{dp}}{T_{s,a,o}-T_{s,a,i}}\dot Q_{cc,wet}+ \left(1-\frac{T_{s,a,o}-T_{dp}}{T_{s,a,o}-T_{s,a,i}}\right)\dot Q_{cc,dry}& T_{s,a,i}> T_{dp} > T_{s,a,o}\end{array}\right.
    
Internal Heat Exchanger
-----------------------
Once the cooling coil code has been run, the glycol outlet temperature of the cooling coil can be obtained from

.. math::

    T_{g,o,cc}=T_{g,i,cc}+\dot Q_{cc}/\dot m_g
    
The heat transfer rate in the internal heat exchanger is then given by

.. math::

    \dot Q_{IHX}=\varepsilon_{HX}\dot m_g c_{p,g}(T_{g,o,cc}-T_{evap})
    
because the glycol is the limiting capacitance rate in the two-phase portion of the IHX.

Solution Method
---------------

The residuals to driven to zero are therefore an overall energy balance over the refrigerant loop, matching :math:`\dot Q_{cond}` and :math:`\dot Q_{cond,\Delta h}`, and an energy balance over the secondary loop.  So the residual vector as a function of  :math:`T_{evap}`, :math:`T_{cond}`, and :math:`T_{g,i,cc}` can be expressed as

.. math::

    \vec{\Delta}(T_{evap},T_{cond},T_{g,i,cc})=\left[ \begin{array}{c} \dot Q_{IHX}+\dot W_{comp}+\dot Q_{cond} \\ \dot Q_{cond}+\dot Q_{cond,\Delta h} \\ \dot Q_{cc}- \dot Q_{IHX} \end{array}  \right]
    
and a three-dimensional solver can be used to drive the norm of :math:`\vec{\Delta}` to sufficiently close to zero by altering :math:`T_{evap}`, :math:`T_{cond}`, and :math:`T_{g,i,cc}`.

The code for the preconditioners can be found in :download:`Preconditioners.py <../../../PyACHP/Preconditioners.py>`

Nomenclature

===============================  ===================================================
Variable                         Description
===============================  ===================================================
:math:`\alpha_g`                 Mean glycol heat transfer coefficient [W/:math:`\mathrm{m}^2`/K]
:math:`\alpha_r`                 Mean refrigerant heat transfer coefficient [W/:math:`\mathrm{m}^2`/K]
:math:`\alpha_a`                 Mean air heat transfer coefficient [W/:math:`\mathrm{m}^2`/K]
:math:`\vec{\Delta}`             Residual vector [W]
:math:`\eta_a`                   Overall air-side surface efficiency [-]
:math:`\varepsilon_{HX}`         Effectiveness of heat exchangers [-]
:math:`\rho_{ha}`                Density of humid air [kg\ :subscript:`da`\ /m\ :sup:`3`\ ]
:math:`A_{a,total}`              Total air-side surface area of evaporator (fins+tubes) [:math:`\mathrm{m}^2`]
:math:`A_{r,total}`              Total refrigerant-side surface area of evaporator [:math:`\mathrm{m}^2`]
:math:`c_{p,a}`                  Specific heat of humid air [J/kg\ :subscript:`da`\ /K]
:math:`h_{a,s,cc}`               Enthalpy of air saturated at :math:`T_{g,i,cc}` [J/kg\ :subscript:`da`\ ]
:math:`h_{a,s,sat}`              Enthalpy of air saturated at :math:`T_{evap}` [J/kg\ :subscript:`da`\ ]
:math:`h_{a,i}`                  Enthalpy of air inlet to evaporator [J/kg\ :subscript:`da`\ ]
:math:`h_{r,o,comp}`             Compressor outlet enthalpy [J/kg]
:math:`T_{a,i,cond}`             Condenser air inlet dry-bulb temperature [K]
:math:`T_{a,i,cc}`               Cooling coil air inlet dry-bulb temperature [K]
:math:`T_{a,o,cc}`               Cooling coil outlet dry-bulb temperature [K]
:math:`T_{a,i,evap}`             Evaporator air inlet dry-bulb temperature [K]
:math:`T_{a,o,evap}`             Evaporator air outlet dry-bulb temperature [K]
:math:`T_{g,i,cc}`               Cooling coil glycol inlet temperature [K]
:math:`T_{g,o,cc}`               Cooling coil glycol outlet temperature [K]
:math:`T_{dp}`                   Dewpoint temperature of humid air [K]
:math:`T_{evap}`                 Evaporator dew temperature [K]
:math:`T_{cond}`                 Condenser dew temperature [K]
:math:`T_{s,a,i}`                Surface temperature of air at air inlet [K]
:math:`T_{s,a,o}`                Surface temperature of air at air outlet [K]
:math:`\Delta T_{sh}`            Compressor suction superheat [K]
:math:`\Delta T_{sc}`            Condenser outlet subcooling [K]
:math:`p_{cond}`                 Condenser pressure [Pa (abs)]
:math:`\dot m_g`                 Mass flow rate of glycol [kg/s]
:math:`\dot m_r`                 Mass flow rate of refrigerant [kg/s]
:math:`\dot m_{a,total}`         Mass flow rate of dry air through evaporator [kg\ :subscript:`da`\ /s]
:math:`\dot Q_{evap}`            Evaporator heat transfer rate [W]
:math:`\dot Q_{evap,dry}`        Evaporator fully-dry heat transfer rate [W]
:math:`\dot Q_{evap,wet}`        Evaporator fully-wet heat transfer rate [W]
:math:`\dot Q_{cc}`              Cooling Coil heat transfer rate [W]
:math:`\dot Q_{cc,dry}`          Cooling Coil fully-dry heat transfer rate [W]
:math:`\dot Q_{cc,wet}`          Cooling coil fully-wet heat transfer rate [W]
:math:`\dot Q_{cond}`            Condenser heat transfer rate [W]
:math:`\dot Q_{cond,\Delta h}`   Condenser heat transfer rate from change in enthalpy [W]
:math:`\dot Q_{IHX}`             Internal Heat Exchanger heat transfer rate [W]
:math:`\mathrm{UA}_a`            Air-side :math:`\mathrm{UA}` value [W/K]
:math:`\mathrm{UA}_g`            Glycol-side :math:`\mathrm{UA}` value [W/K]
:math:`\mathrm{UA}_r`            Refrigerant-side :math:`\mathrm{UA}` value [W/K]
:math:`\dot V_{ha,cond}`         Volumetric flow rate of humid air in condenser [m\ :sup:`3`\ /s ]
:math:`\dot V_{ha,evap}`         Volumetric flow rate of humid air in evaporator [m\ :sup:`3`\ /s ]
:math:`\dot V_{ha,cc}`           Volumetric flow rate of humid air in cooling coil [m\ :sup:`3`\ /s ]
:math:`\dot W_{comp}`            Electrical power of compressor [W]
===============================  ===================================================