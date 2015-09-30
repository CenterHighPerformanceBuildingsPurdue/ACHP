.. _Cooling-Coil:

Cooling Coil
============

Overview
--------

The goal of a cooling coil is to cool and/or dehumidify an air stream by passing the air stream over a finned set of tubes that have cooler glycol (or similar liquid) flowing through them.  In principle the exact same physical heat exchanger can be used to heat an air stream.

Mathematical Description
------------------------
The cooling coil is a particular case of the Heat Exchanger modeling presented in section :ref:`Heat-Mass-HX`.  The water (or similar liquid) passing through the heat exchanger is single phase, and the humid air is assumed to be moist enough that condensation on the heat exchanger surfaces is possible.  Therefore, the partially-wet/partially-dry analysis of section :ref:`PWPD-Single-Phase` must be used.  The use of this analysis still allows for the possibility that the coil is fully wet or fully dry. The entire area of the heat exchanger is employed in order to calculate :math:`A_{a,total}` and :math:`\dot m_{a,total}`.

Cooling Coil Sample Code
------------------------
Minimal Component Test:

.. literalinclude:: ComponentTests/CoolingCoilTest.py
    :language: python
    
If you open an IPython(x,y) shell in the root of the documentation (folder Documentation/Web relative to the main trunk), and run the commands below, you should get

.. ipython::

    In [1]: execfile('ACHPComponents/ComponentTests/CoolingCoilTest.py')
    
If not, first stop should be the :ref:`FAQS`

Component Class Documentation
-----------------------------

.. py:module:: ACHP.CoolingCoil    
.. autoclass:: CoolingCoilClass
    :members:
    :undoc-members:
