.. _Condenser:

Condenser
============

Overview
--------

The goal of a condenser is to take a refrigerant stream at a superheated state and cool it at nearly constant pressure so that it condenses to a subcooled refrigerant that can then be throttled to the low-pressure side of the system, exiting the throttling valve at low temperature.

Mathematical Description
------------------------

On the refrigerant side, the refrigerant path can be divided into as many as three sections - a superheated section, a two-phase section, and a subcooled section.  In practice, operating conditions (particularly low charge), may result in a condenser that does not have a subcooled section, and at extreme conditions, does not have a two-phase section either.  The following plot shows the possible configurations for the first two cases described:

.. plot:: MPLPlots/CondenserSections.py

In both cases, the schematic represents an averaged circuit on the refrigerant side.  Overall, the goal of solver for the condenser is to determine how much of the length of an averaged circuit is in the superheated, two-phase, and subcooled phases.  The fraction of the circuit length in each of these segments are given by :math:`w_{superheat}`, :math:`w_{two-phase}`, and :math:`w_{subcool}` respectively. The sum of these factors must be equal to unity, and if any of the sections do not exist, its respective circuit fraction parameter :math:`w` is set to zero.

Thus if the average length of one circuit is given by :math:`\overline{L_{circuit}}`, then the lengths of each segment are given by


.. math::
    :label: eqCond1
    
    L_{superheat}=w_{superheat}\overline{L_{circuit}}
    
    L_{two-phase}=w_{two-phase}\overline{L_{circuit}}
    
    L_{subcool}=w_{subcool}\overline{L_{circuit}}
    
and for consistency

.. math::
    :label: eqCond2
    
    w_{superheat}+w_{two-phase}+w_{subcool}=1
    
On the air side, the condenser is treated as being pure cross-flow.  This is a particularly good assumption for the condition where the condenser is composed of a single bank of tubes, which is a common configuration for condensers in the USA.  Thus, since the flow is pure crossflow, the air-side area, as well as the air flow rate for each of the superheated, two-phase, and subcooled sections can be assumed to be proportional to the area for that section (See below).

The analysis here makes a number of assumptions:

#. There is no condensation of moist air on the outside of the tubes since the temperature of refrigerant is above the dewpoint of the air stream
#. The flow is evenly balanced between all the circuits on the refrigerant side
#. The flow is evenly distributed on the air side

Air-side Considerations
^^^^^^^^^^^^^^^^^^^^^^^

On the air side, the air flow rate flow :math:`\dot m_{a,total}` is assumed to be evenly distributed across the face of the coil which has a total air-side area of :math:`A_{a,total}`. Therefore the areas and mass flow rates for each section are given by

.. math::
    :label: eqCond3
    
    A_{a,superheat}=w_{superheat}A_{a,total}
    
    A_{a,two-phase}=w_{two-phase}A_{a,total}
    
    A_{a,subcool}=w_{subcool}A_{a,total}
    
    \dot m_{a,superheat}=w_{superheat}\dot m_{a,total}
    
    \dot m_{a,two-phase}=w_{two-phase}\dot m_{a,total}
    
    \dot m_{a,subcool}=w_{subcool}\dot m_{a,total}

.. plot:: MPLPlots/CondenserFace.py

Condenser Algorithm
^^^^^^^^^^^^^^^^^^^

#. Explicitly solve the superheated section to determine :math:`w_{superheat}`

#. First assume that the refrigerant exits the two-phase portion at a quality of 0, and calculate the required :math:`w_{two-phase}` from the heat transfer analysis for the two-phase section.  Then either

    * If the value of :math:`w_{superheat} + w_{two-phase}` is less than 1, there is a subcooled section
    
        * The outlet quality of the two-phase section is therefore known to be 0, which yields the value for :math:`w_{two-phase}`
    
        * Use the remainder of the condenser area for the subcooled section (:math:`w_{subcool}=1-w_{two-phase}-w_{superheat}`), solve subcooled section for heat transfer rate, pressure drop, etc.
    
    * If the value of :math:`w_{superheat} + w_{two-phase}` is greater than 1, there is no subcooled section
    
        * Area fraction for two phase section is known to be :math:`w_{two-phase}=1-w_{superheat}`, iteratively calculate the outlet quality of the condenser
        
The analysis for each of these calculations is described in the sections that follow.

Superheated Section
^^^^^^^^^^^^^^^^^^^

In the superheated region, the inlet refrigerant temperature and the outlet refrigerant temperature (dew temperature) are known.  Therefore, assuming pure crossflow, the fraction of the area required for the superheated region can be explicitly obtained from

.. math::
    :label: eqCond4
    
    w_{superheat}=-\frac{\ln(1 - \Psi)}{ [1 - \exp(-\mathrm{UA}_{overall}/(c_{p,a}\dot m_{a,total}))]}\frac{\dot m_r c_{p,r}}{\dot m_{a,total} c_{p,a}}

where

.. math::
    :label: eqCond7a
    
    \Psi=\frac{(T_{r,i}-T_{dew,r})}{(T_{r,i}-T_{a,i})}
    
and where the derivation of this term can be obtained from :ref:`CondenserAreaDerivation`.  The superheated portion is assumed to exist.

Two-Phase Section
^^^^^^^^^^^^^^^^^

In the two-phase section of the heat exchanger, there are two basic possibilities.  Either the outlet of the two-phase section is at some two-phase quality (because there is no subcooled section) or the outlet of the two-phase section is at a quality of 0 (saturated liquid) because there is a subcooled section.  The first step is to assume that the outlet of the two-phase section is at a quality of 0 and calculate the required fraction of the circuit length.  Before the length fraction can be calculated, the average refrigerant side heat transfer coefficient is required, for which the analysis can be found in section :ref:`Shah Condensation <Shah-Condensation>`.  The average refrigerant-side heat transfer coefficient is a function of refrigerant outlet quality.

In either configuration of the two-phase section (solving for outlet quality or solving for :math:`w`), the effectiveness is known, since the parameter :math:`w` cancels out of the solution for :math:`\mathrm{Ntu}`, which means that the Ntu are independent of the length fraction :math:`w`.  Thus the effectiveness can be obtained directly from 

.. math::
    :label: eqCond5
    
    \varepsilon_{two-phase}=1-\exp\left(\frac{\mathrm{UA}_{overall}}{\dot m_{a,total}c_{p,a}}\right)
    
where the value of :math:`\mathrm{UA}_{overall}` is given by

.. math::
    :label: eqCond6
    
    \mathrm{UA}_{overall}=\dfrac{1}{  (\eta_a \alpha_a A_{a,total})^{-1} + (\alpha_{r,two-phase} A_{r,total})^{-1}}
    
For a given outlet quality of the two-phase section (:math:`x_{out,r,two-phase}`), the length fraction :math:`w_{two-phase}` can be obtained from

.. math::
    :label: w_2phase_condenser
    
    w_{two-phase}=-\frac{\dot m_r h_{fg}(1-x_{out,r,two-phase})}{\dot m_{a,total} c_{p,a}(T_{a,i}-T_{sat,r}) \varepsilon_{two-phase}}
    
Otherwise, if the quality is being iterated for, the residual to be driven to zero by altering the outlet quality of the two phase section  (:math:`x_{out,r,two-phase}`) is

.. math::
    :label: eqCond7
    
    \Delta=(1-w_{superheat})-w_{two-phase}(x_{out,r,two-phase})
    
where :math:`w_{two-phase}` is evaluated from equation :eq:`w_2phase_condenser`, and :math:`w_{superheat}` is known from its solution above.  Unfortunately, iterative methods are required since the value of the average refrigerant two-phase heat transfer coefficient is dependent on the two-phase section outlet quality.

Subcooled Section
^^^^^^^^^^^^^^^^^

In the subcooled section, if it exists, the available area is known (since :math:`w_{subcool}` is known), as are the inlet air and refrigerant temperatures.  The inlet refrigerant temperature to the subcooled section is the bubble temperature of the refrigerant, and the inlet air temperature is the same inlet air temperature as for all the other sections. Thus the :math:`\mathrm{UA}` value can be given by

.. math::
    :label: eqCond8
    
    \mathrm{UA}=\dfrac{w_{subcool}}{  (\eta_a \alpha_a A_{a,total})^{-1} + (\alpha_{r} A_{r,total})^{-1}}
    
and the minimum and maximum capacitance rates of air and subcooled refrigerant can be obtained from

.. math::
    :label: eqCond9
    
    C_{min}=\min[{\dot m_r c_{p,r},\dot m_{a,total}c_{p,a}w_{subcool}}]
    
    C_{max}=\max[{\dot m_r c_{p,r},\dot m_{a,total}c_{p,a}w_{subcool}}]
    
which give the :math:`\mathrm{Ntu}` from 

.. math::
    :label: eqCond10
    
    C_r=\frac{C_{min}}{C_{max}}
    
    \mathrm{Ntu}=\frac{\mathrm{UA}}{C_{min}}
    
which yields the effectiveness, if the minimum capacitance is on the air side, of

.. math::
    :label: eqCond11
    
    \varepsilon_{subcool} = \frac{1}{C_r} \left(1 - \exp(-C_r (1 - \exp(-\mathrm{Ntu})))\right)
    
or, if the minimum capacitance rate is on the refrigerant side, the effectiveness is given by

.. math::
    :label: eqCond12
    
    \varepsilon_{subcool} = 1 - \exp\left(-\frac{1}{C_r} (1 - \exp(-C_r \mathrm{Ntu}))\right)
    
The total heat transfer rate in the subcooled section is given by

.. math::
    :label: eqCond13
    
    \dot Q_{subcool}=-\varepsilon_{subcool}C_{min}(T_{bubble,r}-T_{i,a})
    
which is negative since heat is removed from the refrigerant.

Terminal Calculations
^^^^^^^^^^^^^^^^^^^^^

After all the sections of the condenser have been solved and all the requisite terms determined, the results from each section are grouped together, yielding the overall terms

.. math::
    :label: eqCond14
    
    \dot Q=\dot Q_{superheat}+\dot Q_{two-phase}+\dot Q_{subcool}
    
    \Delta p_{r}=\Delta p_{r,superheat}+\Delta p_{r,two-phase}+\Delta p_{r,subcool}
    
    m_{r}=m_{r,superheat}+m_{r,two-phase}+m_{r,subcool}
    
If there is a subcooled section, the condenser outlet subcooling is evaluated from

.. math::
    :label: eqCond15
    
    \Delta T_{sc}=T_{bubble,r}-T_{r,o}

or if there is no subcooled section, an effective subcooling amount is calculated from

.. math::
    :label: eqCond16
    
    \Delta T_{sc}=\frac{h_{fg}x_{out,r,two-phase}}{c_{p,dew}}
    
where :math:`c_{p,dew}` is the specific heat of saturated liquid.  This effective subcooling parameter is primarily needed to allow for continuous behavior during the solving process.  Hopefully the subcooled section exists at cycle model convergence.

Condenser Sample Code
---------------------

Minimal Component Test:

.. literalinclude:: ComponentTests/CondenserTest.py


If you open an IPython(x,y) shell in the root of the documentation (folder Documentation/Web relative to the main trunk), and run the commands below, you should get

.. ipython::

    In [1]: execfile('ACHPComponents/ComponentTests/CondenserTest.py')
    
If not, first stop should be the :ref:`FAQS`

**Nomenclature**

===============================  ===================================================
Variable                         Description
===============================  ===================================================
:math:`A_{a,total}`              Total air side area (fins+tubes) [m\ :sup:`2`\ ]
:math:`A_{a,superheat}`          Air-side area for superheated portion of circuit [m\ :sup:`2`\ ]
:math:`A_{a,two-phase}`          Air-side area for two-phase portion of circuit [m\ :sup:`2`\ ]
:math:`A_{a,subcool}`            Air-side area for subcooled portion of circuit [m\ :sup:`2`\ ]
:math:`A_{r,total}`              Total refrigerant side area of all tubes [m\ :sup:`2`\ ]
:math:`c_{p,r}`                  Specific heat of refrigerant [J/kg/K]
:math:`c_{p,dew}`                  Specific heat of refrigerant at dewpoint temperature [J/kg/K]
:math:`c_{p,a}`                  Specific heat of humid air per kg of dry air [J/kg\ :sub:`da`\ /K]
:math:`C_{min}`                  Minimum capacitance rate [W/K]
:math:`C_{max}`                  Maximum capacitance rate [W/K]
:math:`C_{r}`                    Ratio of capacitance rates [-]
:math:`h_{fg}`                   Refrigerant latent heat [J/kg]
:math:`\overline{L_{circuit}}`   Average circuit length [m]
:math:`\dot m_{a,total}`         Mass flow rate of dry air [kg/s]
:math:`\dot m_{a,superheat}`     Mass flow rate of dry air in superheated section [kg/s]
:math:`\dot m_{a,two-phase}`     Mass flow rate of dry air in two-phase section [kg/s]
:math:`\dot m_{a,subcool}`       Mass flow rate of dry air in subcooled section [kg/s]
:math:`\dot m_r`                 Mass flow rate of refrigerant [kg/s]
:math:`m_r`                      Mass of refrigerant charge (total) [kg]
:math:`m_{r,superheat}`          Mass of refrigerant charge in superheated section [kg]
:math:`m_{r,two-phase}`          Mass of refrigerant charge in two-phase section [kg]
:math:`m_{r,subcool}`            Mass of refrigerant charge in subcooled section [kg]
:math:`\mathrm{Ntu}`             Number of thermal units [-]
:math:`L_{superheat}`            Length of superheated portion of circuit [m]
:math:`L_{two-phase}`            Length of two-phase portion of circuit [m]
:math:`L_{subcool}`              Length of subcooled portion of circuit [m]
:math:`\Delta p_r`               Refrigerant-side pressure drop (total) [Pa]
:math:`\Delta p_{r,superheat}`   Refrigerant-side pressure drop in superheated section [Pa]
:math:`\Delta p_{r,two-phase}`   Refrigerant-side pressure drop in two-phase section [Pa]
:math:`\Delta p_{r,subcool}`     Refrigerant-side pressure drop in subcooled section [Pa]
:math:`\dot Q`                   Total heat transfer rate [W]
:math:`\dot Q_{superheat}`       Heat transfer rate in superheated section [W]
:math:`\dot Q_{two-phase}`       Heat transfer rate in two-phase section [W]
:math:`\dot Q_{subcool}`         Heat transfer rate in subcooled section [W]
:math:`\Delta T_{sc}`            Subcooling amount [K]
:math:`T_{a,i}`                  Air inlet temperature [K]
:math:`T_{bubble,r}`             Bubble-point temperature of refrigerant [K]
:math:`T_{dew,r}`                Dew-point temperature of refrigerant [K]
:math:`T_{r,i}`                  Refrigerant inlet temperature [K]
:math:`T_{sat,r}`                Refrigerant saturation temperature [K]
:math:`\mathrm{UA}_{overall}`    Overall heat transfer conductance [W/K]
:math:`w_{superheat}`            Fraction of length of superheated portion of circuit [-]
:math:`w_{two-phase}`            Fraction of length of two-phase portion of circuit [-]
:math:`w_{subcool}`              Fraction of length of subcooled portion of circuit [-]
:math:`x_{out,r,two-phase}`      Quality of refrigerant at outlet of two-phase section [-]
:math:`\eta_a`                   Overall surface effectiveness on air side [-]
:math:`\alpha_{a}`               Air side heat transfer coefficient [W/m\ :sup:`2`\ /K]
:math:`\alpha_{r,two-phase}`     Refrigerant side heat transfer coefficient [W/m\ :sup:`2`\ /K]
:math:`\Delta`                   Residual term [-]
:math:`\varepsilon_{two-phase}`  Effectiveness of heat transfer in the two-phase section [-]
:math:`\varepsilon_{subcool}`    Effectiveness of heat transfer in the subcooled section [-]
:math:`\Psi`                     Effectiveness term from derivation [-]
===============================  ===================================================

Component Class Documentation
-----------------------------

.. py:module:: Condenser    
.. autoclass:: CondenserClass
    :members:
    :undoc-members:
    