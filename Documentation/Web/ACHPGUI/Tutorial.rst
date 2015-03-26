Tutorial
========

Getting Started
---------------

In order to get started using ACHP, start the ACHP.exe executable which after showing a splash screen, should boot up into a screen that looks like

.. image:: /ACHPGUI/images/ACHPMain.png
	:width: 100%
	
This is where a number of the important parameters can be changed.

The default configuration file that ships with ACHP is loaded from the file *Default.cfg* and is a direct expansion air conditioning system that is derived from the experimental and modeling work of Bo Shen [#fBoShen]_.  At any time, the file *Default.cfg* can be overwritten which changes the parameters that will be loaded into ACHP when it starts.

To run the model, press the F5 key, or go into the menu to Solve->Solve (F5)

Once the model has run (you can watch the console for intermediate output while it is running), all the output fields are then populated.

In the output are things like

Main Output Screen

.. image:: /ACHPGUI/images/MainOutput.png
	:width: 100%
	
Temperature-Entropy plots

.. image:: /ACHPGUI/images/TsPlot.png
	:width: 100%
	
Pressure-Enthalpy plots

.. image:: /ACHPGUI/images/PhPlot.png
	:width: 100%
	
Output screens for each component

.. image:: /ACHPGUI/images/ComponentOutput.png
	:width: 100%

To save the data from a run to a file, click the **Write to File...** button on the main output screen.  It will prompt you for a location to save the file.  The file comprises the all the data that is output on the output screens of ACHP.  It is a comma-separated file that you can open up in Excel.

.. rubric:: Footnotes

.. [#fBoShen] Bo Shen, 2006, "Improvement and Validation of Unitary Air Conditioner and Heat Pump Simulation Models at Off-Design Conditions" Final Report 1173-RP `Link to file <http://rp.ashrae.biz/page/rp-1173.pdf>`_