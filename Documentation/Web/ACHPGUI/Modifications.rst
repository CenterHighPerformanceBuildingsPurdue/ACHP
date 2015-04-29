Making Modifications
--------------------

While it is nice to be able to see the performance of the default system, you probably want to make some changes for your own system.  To do that, it is necessary to understand how ACHP is configured.

ACHP was built to be able to handle both conventional four-component air-conditioning and heat pump systems as well as secondary loop systems that use a secondary working fluid between the primary working fluid and the indoor coil.  This type of system has been proposed, and is used, in large industrial scale systems.

Comments
^^^^^^^^

In each of the input screens, if you hover your mouse over the text which describes a field, the tooltip that pops up will have further description of the field.

The user interface of ACHP is broken up into three parts - Input, Solvers, and Output.  The tabs at the top of the screen allow you to move between the different parts.

Main Inputs
^^^^^^^^^^^

A number of different types of refrigerants can be used.  The types of refrigerants that are available in ACHP are 

* Working fluids that are implemented in `CoolProp <http://coolprop.sourceforge.net>`_.  A list of working fluids implmented in CoolProp can be found at the CoolProp website.  These refrigerants are listed at the top of the list, before the REFPROP refrigerants
* Pure refrigerants, and defined blends from REFPROP.  These fluids are listed as starting with "REFPROP-"
* Arbitrary blends of fluids from REFPROP.  This functionality is achieved by creating a string which describes the blend.  For instance, the fluid name "REFPROP-MIX:R32[0.697615]&R125[0.302385]" is the refrigerant R410A which is composed of a R32 and R125 mole fractions of 0.697615 and 0.302385 respectively.  The syntax of the string must be followed, but can be extended if more fluids comprise the mixture.  For instance, the fluid name "REFPROP-MIX:R125[0.35782]&R134a[0.038264]&R143a[0.60392]" would be the fluid name string for R404A

Since no expansion device model is employed in ACHP, the superheat is imposed as a model input.  The superheat is that at the outlet of the evaporator, or in the case of the secondary loop system in cooling mode, the outlet of the internal heat exchanger.

Either the subcooling or the charge can be imposed by selecting one of the variables and providing its value, and the other one is then solved for.  The convergence characteristics for subcooling imposed are slightly better than for charge imposed due to the formulation of the preconditioner.  At the very least, if charge is desired to be imposed, a reasonable value for the charge should be obtained by fixing the subcooling for one point.

On the main screen of ACHP you can select whether the system is operating in cooling mode, or in heating mode, and whether a secondary loop is employed, or whether the system is a direct expansion system.  As the mode and secondary loop/DX options are changed, the cycle schematic will also change.  

If the cycle is a secondary loop cycle, the secondary working fluid can be selected from the dropdown box.

Heat Exchanger Inputs
^^^^^^^^^^^^^^^^^^^^^

In both secondary loop systems and direct expansion systems, there are two heat exchangers that are used to transfer heat with the indoor and outdoor air streams.  The input screens for both of the heat exchangers are identical.  

Currently, there are three type of fins coded -  plain fin, herringbone fin and wavy lanced fin.  If other fin types are needed, they must be coded into ACHP.  

The inputs are divided into three subgroups - fins, tubes and air.  In the fins group, you can select the fin pitch, waviness parameters of the fin, as well as the fin thickness, and fin material thermal conductivity.

In the tubes section, the geometry of the tubes and the circuits is defined.  At any time, clicking on the **Show Circuits...** button will pop up a screen that will show an end-on view of the heat exchanger in order to demonstrate the construction of the currently defined heat exchanger.  Clicking on the **Select** button will open up a small window that will allow you to select from standard tube dimensions.

Pump and Internal Heat Exchanger Inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the pump and internal heat exchanger inputs pane, it is possible to set inputs for the pump and the internal heat exchanger.  The pump inputs are straightforward, but the internal heat exchanger inputs need some description.

In cooling mode, either the coaxial heat exchanger or the plate heat exchanger can be used.  In heating mode, the coaxial heat exchanger model has not been modified to handle condensation on the refrigerant side, so therefore it is not possible to use the coaxial heat exchanger in heating mode.

Compressor Inputs
^^^^^^^^^^^^^^^^^
The compressor model is based on a 10-coefficient compressor model as in ANSI/AHRI standard 540-2004

.. math::

	\dot {m} = M_{1}+M_{2}T_s+M_{3}T_{d}+M_{4}T_{s}^2+M_{5}T_{s}T_d+M_{6}T_{d}^2+M_{7}T_{s}^3+M_{8}T_{d}T_{s}^2+M_{9}T_{d}^{2}T_{s}+M_{10}T_{d}^3
	
.. math::

	Power = P_{1}+P_{2}T_s+P_{3}T_{d}+P_{4}T_{s}^2+P_{5}T_{s}T_d+P_{6}T_{d}^2+P_{7}T_{s}^3+P_{8}T_{d}T_{s}^2+P_{9}T_{d}^{2}T_{s}+P_{10}T_{d}^3

where :math:`T_s` and :math:`T_d` are refrigerant dew point temperatures in degrees Fahrenheit, the mass flow rate :math:`\dot m` is in lbm/hr and the power is in Watts.  The coefficients can be read in from a comma-delimited file by clicking the **Load Coefficients** button where the first column in the file are the coefficients :math:`M_1,M_2,...` and the second column are the coefficients :math:`P_1,P_2,...`.  If the compressor coefficients are loaded from a file, it will over-write the coefficients for the compressor present in the user interface.  A sample compressor coefficient file would be::

	117.316,-461.3
	5.094,-18.6
	-1.593,46.9
	4.48E-02,-0.21
	-2.14E-02,0.43
	1.04E-02,-0.44
	7.90E-05,2.25E-04
	-5.73E-05,2.37E-03
	1.79E-04,-3.32E-03
	-8.08E-05,2.50E-03

Information about the compressor map can be obtained by clicking on the **Info** button

If the primary refrigerant is changed, an appropriate compressor map must be employed for the given primary refrigerant.

Line Set Inputs
^^^^^^^^^^^^^^^
The line set allows the condensing unit and the indoor coil to be physically separated.  

In the line set inputs pane you can set inputs for the line set.  The supply line set goes from the outdoor unit to the indoor coil, and the return line goes from the indoor coil to the outdoor unit.
