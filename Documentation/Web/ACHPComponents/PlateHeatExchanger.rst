.. |Ntu| replace:: :math:`\mathrm{Ntu}`

.. _Plate-Heat-Exchanger:

Plate-Heat-Exchangers
*********************

Overview and Geometry
=====================
The motivating factor that drives the use of brazed-plate heat exchangers is that they are a highly-compact heat exchanger that allows for excellent heat transfer between two fluids with very well controlled pressure drop.  They tend to be slightly more expensive than equivalent coaxial type heat exchangers due to their more exacting manufacturing requirements.  But they can be easily altered to add more plates to give more surface area for increased heat transfer area and lower pressure drop.  The trade-off as usual is that adding plates to decrease the pressure drop also results in a decrease in heat transfer coefficient, which means that each m\ :sup:`2`\  of surface area in the HX becomes less useful.

In the most basic configuration of a plate heat exchanger, hot and cold streams in pure counterflow alternate through the stack of plates.  From the side view, a simplified schematic of the BPHE is something like this:

.. plot:: MPLPlots/PHE/PHESideView.py

In practice, it is sometimes useful to have one of the stream do multiple passes for one pass of the other stream, but this capability is not included in the BPHE model as of this time.  This is commonly used when the capacitance rates are very different, a sufficient heat transfer rate cannot be achieved for one fluid, or when one of the fluids is particularly sensitive to pressure drop. All of these issues are particularly strongly felt for the flow of gases.

In practice, what you get when you put it all together is a heat exchanger that looks something like from face-on:

.. plot:: MPLPlots/PHE/PHEOverallGeometry.py

Robust gasketing of the plates is required to ensure that the fluid phases do not mix at the inlets and outlets of the plates.

Each of the plates that form the internal surface of the BPHE are formed of plates with a wavy shape, and the edge of the plate looks something like this:

.. plot:: MPLPlots/PHE/PHEWavyPlate.py

Based on this geometry, the hydraulic diameter :math:`d_h` is defined by

.. math::
    :label: eqPHE1

    d_h=\frac{4\hat a}{\Phi}

where the parameter :math:`\Phi` is the ratio of the actual area to the planar area enclosed by the edges of the plate (bounded by :math:`L` and :math:`B` in the figure above).  Typically the plates do not have exactly a sinusoidal profile, but, to a decent approximation, their profile is sinusoidal, which results in the value of :math:`\Phi` of

.. math::
    :label: eqPHE2    
    
    \Phi=\frac{1}{6}\left(1+\sqrt{1+X^2}+4\sqrt{1+X^2/2} \right) 
    
where the wavenumber :math:`X` is given by

.. math::
    :label: eqPHE3
    
    X=2\pi\hat a/\Lambda

When the plates are put together to form a stack, the plates are alternated, and as a result, chevron-shaped flow paths are formed, which have the effect of yielding highly mixed flow, resulting in good heat transfer coefficients

A stack of :math:`N_{plates}` plates form the heat exchanger.  There are a total of :math:`N_{channels}` channels formed between the plates.  If :math:`N_{channels}` is evenly divisible by 2, both fluids have the same number of channels.  If :math:`N_{channels}` is not evenly divisible by 2, one stram must have one extra channel.  The outer two plates don't provide any heat transfer as they are just used to maintain the channel structure for the outermost channels.  There are therefore a total of :math:`N_{plates}-2` active plates, for which the active area of one side of one plate is equal to 

.. math::
    :label: eqPHE4
    
    A_{1p}=B_pL_p\Phi
    
Each channel of a fluid gets two sides of this area, which yields the cold- and hot-side wetted areas of

.. math::
    :label: eqPHE5
    
    A_h= 2N_{channels,h}A_{1p}
    
    A_c= 2N_{channels,c}A_{1p}
    
.. math::
    :label: eqPHE6
    
    V_h = N_{channels,h}B_p2\hat a
    
    V_c = N_{channels,c}B_p2\hat a
    
and the mass flow rate of the hot and cold fluids per circuit (needed for correlations), is equal to

.. math::
    :label: eqPHE7
    
    \dot m_{h,ch}=\dot m_{h}/N_{channels,h}
    
    \dot m_{c,ch}=\dot m_{c}/N_{channels,c}

Heat Transfer and Pressure Drop Correlations
============================================

Single-Phase Flow
-----------------

When there is single-phase flow on one side of the heat exchanger, the analysis of Martin [#Martin]_ from the VDI Heat Atlas is employed.  The nomenclature used here mirrors that of Martin.  Using this analysis, the pressure drop and heat transfer coefficients for the fluid flowing between the plates can be calculated.

The Reynolds number for the flow through the channel between two plates is given by

.. math::
    :label: eqPHE8
    
    \mathrm{Re}=\frac{\rho w d_h}{\eta}
    
where the velocity per channel is given by

.. math::
    :label: eqPHE9
    
    w=\frac{\dot m_{ch}}{2\hat{a}\rho B_p}
    
The pressure drop and heat transfer coefficients are usual a function of the Reynolds number, and if the flow is laminar (:math:`\mathrm{Re}<2000`), the 
factors :math:`\zeta_0` and :math:`\zeta_{1,0}` are given by

.. math::
    :label: eqPHE10
    
    \zeta_0=\frac{64}{\mathrm{Re}}
    
    \zeta_{1,0}=\frac{597}{\mathrm{Re}}+3.85
        
and if the flow is turbulent (:math:`\mathrm{Re}\geq 2000`), the factors :math:`\zeta_0` and :math:`\zeta_{1,0}` are given by

.. math::
    :label: eqPHE11
    
    \zeta_0=\frac{1}{(1.8\ln(\mathrm{Re})-1.5)^2}
    
    \zeta_{1,0}=\frac{39}{\mathrm{Re}^{0.289}}

The friction factor :math:`\zeta` is obtained from

.. math::
    :label: eqPHE12
    
    \frac{1}{\sqrt{\zeta}}=\frac{\cos\varphi}{\sqrt{b\tan\varphi+c\sin\varphi+\zeta_0/\cos\varphi}}+\frac{1-\cos\varphi}{\sqrt{\zeta_1}}
    
where the factor :math:`\zeta_1` is given by

.. math::
    :label: eqPHE13
    
    \zeta_1=a\zeta_{1,0}
    
and the factors :math:`a` and :math:`b` and :math:`c` given by Martin are

.. math::
    :label: eqPHE14
    
    a=3.8

    b=0.18
        
    c=0.36
        
The Hagen number is defined by

.. math::
    :label: eqPHE15
    
    \mathrm{Hg}=\frac{\zeta\mathrm{Re}^2}{2}=\frac{\rho\Delta pd_h^3}{\eta^2L_p}

which gives the value for the pressure drop.

The Nusselt number is obtained from

.. math::
    :label: eqPHE16
    
    \mathrm{Nu}=c_q\mathrm{Pr}^{1/3}(\eta/\eta_w)^{1/6} [2 \mathrm{Hg} \sin(2\phi)]^{q}
    
where the recommended values of the constants :math:`c_q` and :math:`q` from Martin are 0.122 and 0.374 respectively.  Finally the overall heat transfer coefficient is obtained from

.. math::
    :label: eqPHE17
    
    \alpha=\frac{k\mathrm{Nu}}{d_h}

Two-Phase Evaporating Flow
--------------------------

When the fluid flow is evaporating, it is quite a bit more difficult to determine the best model to use.  There are contradictory conclusions drawn in literature as to what type of heat transfer is occurring.

It seems like the most accepted view, though is open to debate, is that the flow is governed by nucleate boiling within the channels, and as a result, nucleate pool boiling relations are employed in order to calculate the heat transfer coefficient.  This model has some features which are not well-suited to implementation into the PHE model.  For one, there is no quality dependence on heat transfer coefficient, which yields un-physically high values of heat transfer coefficient at high quality (should go to the saturated vapor gas heat transfer coefficient at pure vapor).

In spite of these shortcomings, the pool boiling correlation of Cooper [#Cooper]_ was used.  This yields a simple form of the solution for the heat transfer coefficient.  The heat transfer coefficient is obtained from

.. math::
    :label: eqPHE18
    
    \alpha=55 (p^*)^{0.12-0.2 \log_{10}(R_p)}(-\log_{10}(p^*))^{-0.55} (q")^{0.67} M^{-0.5}

where

.. math::
    :label: eqPHE19
    
    p^*=p/p_{crit}
    
    q"=\dot Q/A
    
and :math:`M` is the molar mass (kg/kmol) of the fluid and :math:`R_p` is the relative roughness of the surface.

In order to calculate the pressure drop in evaporating flow in the PHE channel, the frictional pressure drop is calculated using the Lockhart-Martinelli two-phase pressure drop correlation from section :ref:`Lockhart-Martinelli<Lockhart-Martinelli>` with the value of the parameter C of 4.67 as recommended by Claesson [#Claesson]_.  The accelerational pressure change is given from the same section.

Two-Phase Condensing Flow
-------------------------
The available models for condensing flow in PHE share many of the shortcomings of the evaporating flow models.  There is a paucity of good data available, since most of the know-how is controlled by the major PHE manufacturers.  That said, many researchers have studied this topic, but the parameter space (geometrically and thermodynamically) is quite vast.  This is still a topic that could do with further study.

Longo has conducted studies [#Longo2004]_ [#Longo2010]_ [#Longo2010b]_ that looked at condensation in PHE, and from these studies it can be seen that at low equivalent Reynolds number (:math:`\mathrm{Re}_{eq}<1750`), the j-factor is nominally constant at a value of 60, and above that, it is linear with equivalent Reynolds number, so the j-factor can be given by

.. math::
    :label: eqPHE20
    
    j=\left\lbrace \begin{array}{cc} 60 & \mathrm{Re}_{eq}<1750 \\ \frac{75-60}{3000-1750}(\mathrm{Re}_{eq}-1750)+60 & \mathrm{Re}_{eq}\geq 1750 \end{array} \right.
    
where the equivalent Reynolds number is defined by

.. math::
    :label: eqPHE21
    
    \mathrm{Re}_{eq}=\frac{G\left[(1-\overline x)+\overline x \sqrt{\frac{\rho_L}{\rho_V}}\right] d_h}{\eta_L}

which finally yields the heat transfer coefficient

.. math::
    :label: eqPHE22
    
    \alpha=\frac{jk\mathrm{Pr}^{1/3}}{d_h}
    
Mathematical Description
========================

With the set of required correlations defined, it is now possible to analyze the plate heat exchanger for a range of different configurations.  The PHE model is constructed to be general enough that it can handle any phase of fluids entering into the heat exchanger.  The basic idea behind the PHE model is a two-step process:

#. Determine the bounding heat transfer rate (100% effectiveness) limited by taking each fluid to the inlet temperature of the other fluid.  This is the most amount of heat transfer possible.  Also watch out for internal pinch points
#. Since the heat transfer rate is now bounded between zero and the maximum, iterate to find the actual heat transfer rate in the heat exchanger in order to yield the "right-size" heat exchanger (see description below)

Bounds on Heat Transfer Rate
----------------------------

Since the PHE is pure counterflow, the coldest possible temperature that the hot stream can achieve is the inlet temperature of the cold stream, and similarly, the hottest temperature that the cold stream can achieve is the inlet temperature of the cold stream.  The inlet enthalpies of the hot stream :math:`h_{h,i}` and the cold stream :math:`h_{c,i}` allow to calculate a preliminary value for the upper bound on the heat transfer rate:

.. math::
    :label: eqPHE23
    
    \dot Q_{max,h}=\dot m_h [h_{h,i}-h(T=T_{c,i},p=p_{h,i},Ref_h)]
    
    \dot Q_{max,c}=\dot m_c [h(T=T_{h,i},p=p_{c,i},Ref_c)-h_{c,i}]
    
    \dot Q_{max,\varepsilon=1}=\max[\dot Q_{max},\dot Q_{min}]
    
Using this preliminary bound on the heat transfer rate, it is then possible to determine the enthalpies and temperatures of both fluids at each of their phase transitions(if they exist).

In the case of an evaporator that cools a water stream, there is no possibility of temperature inversion within the heat exchanger because the refrigerant enters at some quality greater than zero, and there is no possibility that the isobar of the refrigerant could intersect the isobar of the hot water.  The following figure shows this configuration:

.. _PHEQmaxEvapcells:

.. plot:: MPLPlots/PHE/PHEQmaxEvapcells.py

Since the heat transfer rate is known, it can be readily determined whether any of the phase transitions can be physically reached.  In the configuration shown here, there are two regions; in cell 1 the hot fluid (water) is single phase and the cold fluid (refrigerant) is evaporating, and in cell 2, the hot fluid is still single-phase, and the cold fluid is now single phase as well.

On the other hand, if the refrigerant were condensing, and entering at some subcooling amount greater than zero, for instance 10 K, the analysis is slightly different.  In this case, it is entirely possible that there could be temperature inversion at the heat transfer rate given by :math:`\dot Q_{max,\varepsilon=1}`, as shown in the folllowing physically impossible figure:

.. plot:: MPLPlots/PHE/PHEQmaxCondcells.py

Thus, a new maximum heat transfer rate :math:`\dot Q_{max}` can be determined that is less than :math:`\dot Q_{max,\varepsilon=1}` whereby the temperatures of the two streams are equated at the possible pinch point, which resembles something like the following figure:

.. plot:: MPLPlots/PHE/PHEQmaxCondPinchedcells.py 

In this case, the water is limiting the heat transfer rate, and the maximum heat transfer rate can be given by taking the water all the way to the dew temperature of the refrigerant, and using the known heat transfer rate in cell 3.  The cold-stream pinch enthalpy is given by

.. math::
    :label: eqPHE24
    
    h_{pinch}=h(T=T_{dew,h},p=p_c,Ref_c)

Since the inlet enthalpy and outlet enthalpy (saturated vapor) of the hot refrigerant are known in cell 3, the heat transfer rate in cell 3 is known from

.. math::
    :label: eqPHE25
    
    \dot Q_{cell 3}= \dot m_{h}(h_{h,i}-h(T=T_{dew,h},x=1,Ref_h))
    
and the new limiting heat transfer rate can be given by

.. math::
    :label: eqPHE26
    
    \dot Q_{max}= \dot m_c(h_{pinch}-h_{c,i})+\dot Q_{cell 3}

where the contribution :math:`\dot m_c(h_{pinch}-h_{c,i})` is from heating up the cold fluid to the pinch point temperature. 

Calculation of Heat Transfer Rate
---------------------------------
Now that the physical bounds on the heat transfer rate in the PHE have been determined, it is now possible to finish analyzing the PHE performance.  For a given :math:`\dot Q<\dot Q_{max}`, there are a number of different cells, and in each one, at least one of the fluids has a phase transition.  In the degenerate case that both fluids are single-phase throughout the PHE, there is only one cell, and no phase transitions anywhere in the heat exchanger.

The discussion that follows here assumes that the heat transfer rate :math:`\dot Q` is known, but in practice, it is iteratively obtained by a bounded 1-D solver because :math:`\dot Q` is known to be between 0 and :math:`\dot Q_{max}`.

For a given :math:`\dot Q`, the outlet enthalpies are known, which begins the process of buiding enthalpy vectors for both streams.  The outlet enthalpies for each stream are given by

.. math::
    :label: eqPHE27
    
    h_{h,o}=h_{h,i}-\dot Q/\dot m_h
    
    h_{c,o}=h_{c,i}+\dot Q/\dot m_c
    
which yields the initial enthalpy vectors (ordered from low to high enthalpy) of 

.. math::
    :label: eqPHE28
    
    \vec h_h = [h_{h,o} , h_{h,i}]

    \vec h_c = [h_{c,i} , h_{c,o}]
    
To these enthalpy vectors are now added any phase transitions that exist; a phase transition exists if its corresponding saturation enthalpy is between the inlet and outlet enthalpies of the fluid.  With each phase transition enthalpy comes a partner enthalpy of the other stream.  This set of enthalpy vectors then define the enthalpies of both streams at each cell edge.  For instance, in the case shown in Figure :ref:`Evaporator<PHEQmaxEvapcells>`, there is one phase transition where the refrigerant transitions between two-phase and superheated vapor.  The enthalpy of the cold stream at the transition point is given by

.. math::
    :label: eqPHE29
    
    h_{PT}=h(T=T_{dew,c},x=1,\mathrm{Ref}_c)
    
and the enthalpy of the hot stream at the phase transition :math:`h_{PT}^*` can be obtained by an energy balance over cell 2, which yields

.. math::
    :label: eqPHE30
    
    \dot m _h (h_{h,i}-h_{PT})=\dot m_c (h_{c,o}-h_{PT}^*)
    
or

.. math::
    :label: eqPHE31
    
    h_{PT}^*=h_{c,o}-\frac{\dot m _h (h_{h,i}-h_{PT})}{\dot m_c}
    
and now the enthalpy vectors are given by the values

.. math::
    :label: eqPHE32
    
    \vec h_h = [h_{h,o} , h_{PT}, h_{h,i}]

    \vec h_c = [h_{c,i} , h_{PT}^*,h_{c,o}]

If there are multiple phase transitions on each side, the same method is applied, where the phase transition enthalpies and their partner enthalpies are obtained by an energy balance on the new cell that is formed, working from the outer edges of the enthalpy vectors towards the inside since the outlet enthalpies of both streams are known and can be used in the energy balances to back out partner enthalpies.

For a given value of :math:`\dot Q`, each of the enthalpy vectors has the same length of :math:`N_{cell}+1`, which then form the enthalpy boundaries for the :math:`N_{cell}` cells.

In each cell, first the phase of each fluid must be determined.  Each fluid will have the same phase throughout the entire cell (that was the whole point in the first place!).  The average enthalpy of each fluid in the cell can be used to determine the phase of the each fluid in the cell.  Our goal now is to determine how much of the physical length of the heat exchanger is required to obtain the given duty in each cell.  The required physical heat exchanger length of the cell :math:`w` can be given by

.. math::
    :label: eqPHE33
    
    L_i=w_iL

where all the :math:`w_i` parameters must sum to unity (1.0).  

In a given cell, the heat transfer rate is known because this is how the enthalpy vectors have been constructed.  The heat transfer rate in the cell can be given by

.. math::
    :label: eqPHE34
    
    \dot Q_i=\dot m_r (\vec h_{h,i+1}- \vec h_{h,i})

So long as at least one of the fluids in the cell is single-phase, the effectiveness in the cell can be defined by

.. math::
    :label: eqPHE35
    
    \varepsilon=\frac{\dot Q_{i}}{C_{min}(T_{h,i,cell}-T_{c,i,cell})}
    
where :math:`T_{h,i,cell}` and :math:`T_{c,i,cell}` are the hot fluid and cold fluid inlet temperatures to the cell.  The minimum capacitance rate :math:`C_{min}` is by definition on the single-phase-fluid side.  In the single-phase/two-phase cell case, the minimum capacitance rate is given by

.. math::
    :label: eqPHE36
    
    C_{min}=\dot m_{\mathrm{single-phase}}c_{p,\mathrm{single-phase}}
    
and since the flow is pure counter-flow, the |Ntu| can be obtained directly from 

.. math::
    :label: eqPHE37
    
    \mathrm{Ntu}=\frac{\varepsilon}{1-\varepsilon} \qquad (C_r=0)
    
If both fluids are single phase, the minimum capacitance rate can be obtained from

.. math::
    :label: eqPHE38
    
    C_{min}=\min[\dot m_h c_{p,h},\dot m_c c_{p,c}]
    
    C_{max}=\max[\dot m_h c_{p,h},\dot m_c c_{p,c}]
    
    C_{r}=C_{min}/C_{max}
    
which yields the |Ntu| for the single-phase/single-phase cell with pure counterflow of 

.. math::
    :label: eqPHE39
    
    \mathrm{Ntu}=\frac{1}{C_r-1}\ln\left(\frac{\varepsilon-1}{\varepsilon C_{r}-1} \right) \qquad (C_r>0)

and the required heat conductance can be obtained from

.. math::
    :label: eqPHE40
    
    \mathrm{UA}_{req}=\mathrm{Ntu} C_{min}
    
The actual heat transfer conductance in the call can be given by

.. math::
    :label: eqPHE41
    
    \mathrm{UA}_{actual}=\frac{1}{\alpha cA_c}+\frac{t}{kA}+\frac{1}{\alpha_h A_h}
    
where the areas are based on the total wetted area of the heat exchanger and local heat transfer coefficients (:math:`\alpha_h,\alpha_c`) for the cell are employed.  The fraction of the heat exchanger that would be required for the given thermal duty in the cell can be obtained from 

.. math::
    :label: eqPHE42
    
    w_i=\frac{\mathrm{UA}_{req}}{\mathrm{UA}_{actual}}

Determination of Thermal Duty
-----------------------------
Finally, the heat transfer rate in the PHE is obtained through iterative methods.  The value of :math:`\dot Q` is known to be between zero and :math:`\dot Q_{max}`, and the residual to be driven to zero by a numerical solver is

.. math::
    :label: eqPHE43
    
    \Delta=1-\sum_i [w_i] 
    
which will be zero if :math:`\dot Q` has been appropriately found.

.. only:: html

    .. rubric References

.. [#Longo2004] Longo, G.; Gasparella, A. & Sartori, R. (2004), Experimental heat transfer coefficients during refrigerant vaporisation and condensation inside herringbone-type plate heat exchangers with enhanced surfaces, *International Journal of Heat and Mass Transfer* 47, 4125-4136.

.. [#Longo2010] Longo, G. (2010), Heat transfer and pressure drop during HFC refrigerant saturated vapour condensation inside a brazed plate heat exchanger, *International Journal of Heat and Mass Transfer* 53, 1079-1087.

.. [#Longo2010b] Longo, G. (2010), Heat transfer and pressure drop during hydrocarbon refrigerant condensation inside a brazed plate heat exchanger, *International Journal of Refrigeration* 33, 944-953.

.. [#Claesson] Claesson, J., 2004. Thermal and Hydraulic Performance of Compact Brazed Plate Heat Exchangers Operating as Evaporators in Domestic Heat Pumps. PhD Thesis. KTH.

.. [#Cooper] Cooper, M.G., 1984. Heat Flow Rates in Saturated Nucleate Pool Boiling-A Wide-Ranging Examination Using Reduced Properties. Advances in Heat Transfer, 16, pp.157-239.

.. [#Martin] Holger Martin, VDI Heat Atlas 2010, Chapter B6: Pressure Drop and Heat Transfer in Plate Heat Exchangers

.. A couple of replacements to save typing

.. |m3| replace:: m\ :sup:`3`\ 
.. |m2| replace:: m\ :sup:`2`\ 
.. |Qmaxe1| replace:: :math:`\dot Q_{max,\varepsilon=1}`

**Nomenclature**

===============================  ===================================================
Variable                         Description
===============================  ===================================================
:math:`\hat a`                   Amplitude of plate corrugation [m]
:math:`a`                        Constant for equation [-]
:math:`b`                        Constant for equation [-]
:math:`c`                        Constant for equation [-]
:math:`c_q`                      Constant for equation [-]
:math:`\alpha`                   Heat transfer coefficient [W/|m2|/K]
:math:`A_{1p}`                   Area of one plate [|m2|]
:math:`A`                        Area on given side [|m2|]
:math:`A_h`                      Area on hot side [|m2|]
:math:`A_c`                      Area on cold side [|m2|]
:math:`B`                        Width of wetted section [m]
:math:`B_p`                      Port-port centerline distance [m]
:math:`c_{p,single-phase}`       Specific heat of single-phase fluid [J/kg/K]
:math:`c_{p,c}`                  Specific heat of cold fluid [J/kg/K]
:math:`c_{p,h}`                  Specific heat of hot fluid [J/kg/K]
:math:`C_{min}`                  Minimum capacitance rate [W/K]
:math:`C_{max}`                  Maximum capacitance rate [W/K]
:math:`C_{r}`                    Ratio of capacitance rates [W/K]
:math:`G`                        Refrigerant mass flux [kg/|m2|/s]
:math:`\mathrm{Hg}`              Hagen number [-]
:math:`d_h`                      Hydraulic diameter [m]
:math:`h_{c,i}`                  Enthalpy of cold stream at inlet [J/kg/K]
:math:`h_{h,i}`                  Enthalpy of hot stream at inlet [J/kg/K]
:math:`h_{c,o}`                  Enthalpy of cold stream at outlet [J/kg/K]
:math:`h_{h,o}`                  Enthalpy of hot stream at outlet [J/kg/K]
:math:`h_{pinch}`                Enthalpy of stream at pinch [J/kg/K]
:math:`\vec h_c`                 Vector of cold stream enthalpies [J/kg/K]
:math:`\vec h_h`                 Vector of hot stream enthalpies [J/kg/K]
:math:`h_{PT}`                   Enthalpy at phase transition [J/kg/K]
:math:`h_{PT}^*`                 Complementary enthalpy at phase transition [J/kg/K]
:math:`j`                        Colburn j-factor [-]
:math:`k`                        Thermal conductivity of fluid [W/m/K]
:math:`L_i`                      Length of a given cell [m]
:math:`L`                        Length of wetted section [|m2|]
:math:`L_p`                      Port-port centerline distance [m]
:math:`\dot m_c`                 Total mass flow rate of cold fluid [kg/s]
:math:`\dot m_h`                 Total mass flow rate of hot fluid [kg/s]
:math:`\dot m_{c,ch}`            Mass flow rate of cold fluid per channel [kg/s]
:math:`\dot m_{h,ch}`            Mass flow rate of hot fluid per channel [kg/s]
:math:`M`                        Molar mass [kg/kmol]
:math:`N_{channels}`             Number of channels [-]
:math:`N_{channels,c}`           Number of channels on the cold side [-]
:math:`N_{channels,h}`           Number of channels on the hot side [-]
:math:`N_{plates}`               Number of plates [-]
:math:`N_{cell}`                 Number of cells [-]
:math:`\mathrm{Nu}`              Nusselt number [-]
:math:`p^*`                      Reduced pressure [-]
:math:`p_{crit}`                 Critical pressure [kPa]
:math:`p`                        Saturation pressure [kPa]
:math:`p_c`                      Pressure of cold stream [kPa]
:math:`p_{c,i}`                  Inlet pressure of cold stream [kPa]
:math:`p_{h,i}`                  Inlet pressure of hot stream [kPa]
:math:`\mathrm{Pr}`              Prandtl number [-]
:math:`q`                        Constant for equation [-]
:math:`q"`                       Heat transfer flux [W/|m2|]
:math:`\dot Q`                   Heat transfer rate [W]
:math:`\dot Q_{i}`               Heat transfer rate in the given cell [W]
:math:`\dot Q_{max,c}`           Cold stream max heat transfer rate [W]
:math:`\dot Q_{cell3}`           Heat transfer rate in the highest-enthalpy cell [W]
:math:`\dot Q_{max,h}`           Hot stream max heat transfer rate [W]
:math:`\dot Q_{max}`             Maximum heat transfer rate [W]
|qmaxe1|                         Max heat transfer rate taking each stream to inlet temp of opposite fluid [W]
:math:`\mathrm{Re}`              Reynolds number [-]
:math:`\mathrm{Re}_{eq}`         Equivalent Reynolds number [-]
:math:`T_{h,i}`                  Hot stream inlet temperature to PHE [K]
:math:`T_{c,i}`                  Cold stream inlet temperature to PHE [K]
:math:`T_{h,i,cell}`             Hot stream inlet temperature to cell [K]
:math:`T_{c,i,cell}`             Cold stream inlet temperature to cell [K]
:math:`T_{dew,h}`                Dewpoint temperature of refrigerant [K]
:math:`\mathrm{UA}_{reg}`        Required conductace [W/K]
:math:`\mathrm{UA}_{actual}`     Actual conductance available [W/K]
:math:`V_c`                      Total volume on cold side [|m3|]
:math:`V_h`                      Total volume on hot side [|m3|]
:math:`w`                        Velocity of fluid in channel [m/s]
:math:`w_i`                      Fraction of total length for given cell [-]
:math:`x`                        Refrigerant quality [-]
:math:`\bar x`                   Average refrigerant quality [-]
:math:`X`                        Wave number [-]
:math:`\eta`                     Viscosity of the fluid [Pa-s]
:math:`\eta_L`                   Viscosity of saturated liquid [Pa-s]
:math:`\eta_w`                   Viscosity at the wall temperature [Pa-s]
:math:`\varphi`                  Plate inclination angle [rad]
:math:`\Lambda`                  Plate corrugation wavelength [m]
:math:`\Phi`                     Area increase factor [-]
:math:`\rho`                     Fluid density [kg/|m3|]
:math:`\rho_L`                   Saturated liquid fluid density [kg/|m3|]
:math:`\rho_V`                   Saturated vapor fluid density [kg/|m3|]
:math:`\zeta`                    Friction factor [-]
:math:`\zeta_0`                  Friction factor [-]
:math:`\zeta_{1,0}`              Friction factor [-]
:math:`\zeta_1`                  Friction factor [-]
:math:`\varepsilon`              Effectiveness [-]
===============================  ===================================================