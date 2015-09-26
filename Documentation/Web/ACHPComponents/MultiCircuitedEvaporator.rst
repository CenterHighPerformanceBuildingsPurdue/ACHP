==========================
Multi-Circuited Evaporator
==========================

In an evaporator, there is always the possibility that the distribution of refrigerant or air may be unequal between the circuits, which usually results in lower evaporator performance.  In the case of equal distribution between circuits, the analysis of the :ref:`Evaporator` can be used, otherwise the multi-circuited evaporator is needed.  The multi-circuited evaporator (MCE) analyzed as being a set of evaporators, each comprising one circuit, that are then each fed some (not necessarily even) distribution of refrigerant and air.
The tubes per banks are distributed among the different circuits. If the tubes per bank is divisible by number of circuits, all the circuits will have the same number of tubes. Otherwise, the circuits are ordered from fewer to more if not evenly distributed.

Mathematical Description
------------------------

The MCE model is capable of handling the following types of maldistribution:

* Volumetric air flow mal-distribution

* Air-side inlet condition mal-distribution

* Refrigerant mass-flow-per-circuit mal-distribution
 
* Refrigerant quality-per-circuit mal-distribution

For the MCE, there are :math:`N_{circuits}` circuits, and each circuit is treated as an individual evaporator, since then the analysis for each circuit can be provided by the analysis for the conventional evaporator with one circuit.

For the MCE, the total refrigerant mass flow rate :math:`\dot m_r` is known as an input (which arises from the compressor map).  The total refrigerant mass flow rate per circuit can be given by

.. math::

    \dot m_{r,i} = \gamma _i \dot m_r

where :math:`\gamma_i` is the mass flow distribution factor for the i-th circuit of the evaporator.  The sum of the indices is given by

.. math::

    \sum_{i=0}^{N_{circuits}-1}[\gamma_i]=1

and if the flow is equally distributed, all the terms :math:`\gamma_i` are equal.

If the inlet refrigerant quality is not balanced between circuits, the refrigerant vapor mal-distribution can be given by a set of weighting parameters that distribute refrigerant vapor among the circuits.  The total amount of vapor entering the evaporator can be given by

.. math::
    
    \dot m_v=x\dot m_r
    
and the amount of vapor entering the i-th circuit can be given by 

.. math::

    \dot m_{v,i}=\xi_i \dot m_v
    
where the factors :math:`\xi_i` also must sum to unity.  The inlet quality for each circuit is then equal to 

.. math::

    x_i=\frac{\dot m_{v,i}}{\dot m_{r,i}}
    
where all the :math:`x_i` values must be greater than zero.  The evaporator component model takes enthalpy as the inlet, which can be calculated from

.. math::

    h_i=h(p_{evap},x_i)
    
using the CoolProp property routines.

MCE Sample Code
---------------------

Minimal Component Test:

.. literalinclude:: ComponentTests/MCETest.py

If you open an IPython(x,y) shell in the root of the documentation (folder Documentation/Web relative to the main trunk), and run the commands below, you should get

.. ipython::

    In [1]: execfile('ACHPComponents/ComponentTests/MCETest.py')
    
If not, first stop should be the :ref:`FAQS`

Component Class Documentation
-----------------------------

.. py:module:: ACHP.MultiCircuitEvaporator    
.. autoclass:: MultiCircuitEvaporatorClass
    :members:
    :undoc-members:
