.. _PWPD-Algorithm:

Algorithm of partial-wet-partial-dry evaporator
=============================================================

This derivation was provided by Howard Cheung but was originally carried out by Jim Braun.  Editorial modifications have been made by Ian Bell.

Introduction
------------

Partial-wet-partial-dry analysis is the technique given in Braun [#Braun]_ used to estimate the performance of an air-to-refrigerant evaporating coil with a higher accuracy for estimating the latent heat transfer than completely wet or dry analysis. The methodology is to divide the heat exchanger into two sections by locating the point on the heat exchanger surface where the dewpoint is reached. For the part with surface temperature below dewpoint, they are lumped together for a completely wet analysis. The remaining one is lumped together as a dry section.

To illustrate the idea, the document gives an example on a steady counterflow air-to-refrigerant heat exchanger in which air is cooled and moisture is condensed on the heat exchanger surface. The governing equations are listed based on the assumptions, and an algorithm is derived to solve the implicit mathematical model for the analysis.

Background
----------

In the following analysis, the following assumptions are employed:

* Steady state
* Counterflow heat exchanger
* Single-phase fluid flow on both side
* Constant specific heat throughout the entire heat exchanger
* Constant heat transfer coefficient for air-to-surface and surface-to-refrigerant heat transfer
* Coil completely covered with condensate for wet section
* Unity Lewis number

The analysis can be divided into three parts:

#. Completely dry analysis
#. Completely wet analysis
#. Partial-wet-partial-dry analysis

Initially, the dry analysis is implemented. Should the surface temperature at the air outlet be higher than the dewpoint of the air, the dry coil assumption is accepted. Otherwise the wet analysis is carried out on the entire heat exchanger and the surface temperature at the air inlet is examined. If the temperature is lower than dewpoint, the completely wet coil assumption is accepted. If both assumptions are rejected, there must exist a point on the heat exchanger where the temperature is at the dewpoint and the partial-wet-partial-dry analysis is implemented.

Completely Dry Analysis
^^^^^^^^^^^^^^^^^^^^^^^

To conduct the competely dry analysis, a simple :math:`\varepsilon-\mathrm{Ntu}` method on a counterflow heat exchanger is used. The governing equations of the :math:`\varepsilon-\mathrm{Ntu}` method are listed as Equations :eq:`eq-Q_rdry` to :eq:`eq-Q_adry`.

.. math::
    :label: eq-Q_rdry
    
    \dot{Q}_{dry} = \dot{m}_r c_{p,r} (T_{r,out} - T_{r,in})

.. math::
    :label: eq-Q_dry

    \dot{Q}_{dry} = \varepsilon C_{min} (T_{a,in}-T_{r,in}) 

.. math::
    :label: eq-epsilon_d

    \varepsilon = \frac{1 - \exp(-\mathrm{Ntu}_{dry}(1-C_{ratio}))}{1 - C_{ratio} \exp(-\mathrm{Ntu}_{dry}(1-C_{ratio}))}

.. math::
    :label: eq-Ntu

    \mathrm{Ntu}_{dry} = \frac{\mathrm{UA}_o}{C_{min}}

.. math::
    :label: eq-UA

    \frac{1}{\mathrm{UA}_o} = \frac{1}{U_a A_a}+\frac{1}{U_r A_r}

.. math::
    :label: eq-C_min

    C_{min} = min(\dot{m}_a c_{p,a}, \dot{m}_r c_{p,r})

.. math::
    :label: eq-C_ratio

    C_{ratio} = \frac{C_{min}}{max(\dot{m}_a c_{p,a},\dot{m}_r c_{p,r})} 

.. math::
    :label: eq-Q_adry

    \dot{Q} = \dot{m}_a c_{p,a} (T_{a,in} - T_{a,out}) 

These equations can be solved analytically for the heat exchanger performance. After solving the heat exchanger, one may find the temperature on the surface of the heat exchanger at the air outlet by Equation :eq:`Ts_dry`.

.. math::
    :label: Ts_dry
    
    U_a A_a (T_{a,out} - T_{s,a,out}) = U_{r} A_r (T_{s,a,out} - T_{r,in}) 

If the temperature :math:`T_{s,a,out}` is higher than dewpoint of inlet air, the coil is said to be dry and the heat exchanger performance analysis is completed. Otherwise completely wet analysis is conducted.

Completely Wet Analysis
^^^^^^^^^^^^^^^^^^^^^^^

In the case of wet analysis, the unity Lewis number is used such that the temperatures in the :math:`\varepsilon-\mathrm{Ntu}` method are converted to the corresponding air-water mixture enthalpies to account for the condensation of moisture from air on the heat exchanger surface. The :math:`\varepsilon-\mathrm{Ntu}` method is modified to form governing equations listed from Equation :eq:`eq-Q_w` to :eq:`eq-Q_wa`.

.. math::
    :label: eq-Q_w

    \dot{Q}_{wet} = \varepsilon^* \dot{m}_{min} (h_{a,in}-h_{sat,r,in})

.. math::
    :label: eq-epsilon_wet

    \varepsilon^* = \frac{1 - \exp(-Ntu^*(1-\dot{m}_{ratio}))}{1 - \dot{m}_{ratio} \exp(-Ntu^*(1-\dot{m}_{ratio}))}

.. math::
    :label: eq-Ntu_wet
    
    Ntu^* = \frac{UA^*_o}{\dot{m}_{min}}

.. math::
    :label: eq-m_min

    \dot{m}_{min} = min(\dot{m}_a,\dot{m}_r \frac{c_{p,r}}{c_{s}} )

.. math::
    :label: eq-mstar

    \dot{m}_{ratio} = \frac{\dot{m}_{min}}{max(\dot{m}_a,\dot{m}_r \frac{c_{p,r}}{c_{s}})}

.. math::
    :label: eq-c_s

    c_s = \frac{dh_{sat}}{dT}|_{T = \frac{T_{r,in}+T_{r,out}}{2}}

.. math::
    :label: eq-UAstar

    \frac{1}{UA^*_o} = \frac{c_{p,a}}{U_a^* A_a}+\frac{c_s}{U_{r} A_r}

.. math::
    :label: eq-h_sse

    h_{s,sat,eff} = h_{a,in} + \frac{h_{a,out}-h_{a,in}}{1 - \exp(-\frac{U_a^* A_a}{\dot{m}_a})}

.. math::
    :label: eq-T_aout

    T_{a,out} = T_{s,eff} + (T_{a,in} - T_{s,eff})\exp(-\frac{U_a A_a }{\dot{m}_a c_{p,a}})

.. math::
    :label: eq-Q_wa

    \dot{Q}_{wet} = \dot{m}_a (h_{a,in} - h_{a,out})

These equations are formed explicitly and can be solved analytically for the performance of the heat exchanger. The wet coil assumption is verified by comparing the heat exchanger surface temperature with the dewpoint of inlet air. The temperature can be calculated by Equation :eq:`eq-T_swet` which the definition of :math:`c_{s,local}` is given in Equation :eq:`eq-c_slocal`.

.. math::
    :label: eq-T_swet
    
    U_r A_r (T_{s,in} - T_{r,out}) = UA^*_o(c_{s,local})(h_{a,in} - h_{sat,r,out})

.. math::
    :label: eq-c_slocal
    
    c_{s,local} = \frac{dh_{sat}}{dT}|_{T = \frac{T_{a,in}+T_{r,out}}{2}}

If the temperature is lower than the dewpoint, the completely wet coil assumption is accepted. Otherwise the calculation will proceed to the next part: partial-wet-partial-dry analysis.

Partial-Wet-Partial-Dry Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The partial-wet-partial-dry analysis divides the heat exchanger into two regions: a wet region and a dry region. The heat transfer rate can be divided into two parts as shown in Equation :eq:`eq-Q_pdpw`:

.. math::
    :label: eq-Q_pdpw

    \dot{Q} = \dot{Q}_{f,dry} + \dot{Q}_{f,wet}

Both sections can be addressed based on :math:`\varepsilon-\mathrm{Ntu}` method except that a separate set of governing equations are used for each section. The dry section can be described as lists of equations from Equation :eq:`eq-Qep_dry` to :eq:`eq-ep_fdry`.

.. math::
    :label: eq-Qep_dry
    
    \dot{Q}_{f,dry} = \varepsilon_{f,dry} C_{min} (T_{a,in} - T_{r,x})

.. math::
    :label: eq-Qaf_dry

    \dot{Q}_{f,dry} = \dot{m}_a c_{p,a} (T_{a,in}-T_{a,x})

.. math::
    :label: eq-Qrf_dry

    \dot{Q}_{f,dry} = \dot{m}_r c_{p,r} (T_{r,out} - T_{r,x})

.. math::
    :label: eq-ep_fdry

    \varepsilon_{f,dry}=\frac{1-\exp(-f_{dry}\mathrm{Ntu}_{dry}(1-C_{ratio}))}{1-C_{ratio}\exp(-f_{dry}\mathrm{Ntu}_{dry}(1-C_{ratio}))}

The wet region is governed by a similar set of equation with :math:`\varepsilon-\mathrm{Ntu}` method as listed from Equation :eq:`eq-Qep_wet` to :eq:`eq-ep_fwet`.

.. math::
    :label: eq-Qep_wet

    \dot{Q}_{f,wet} = \varepsilon_{f,wet} \dot{m}_{min} (h_{a,x} - h_{sat,r,in})

.. math::
    :label: eq-Q_afwet

    \dot{Q}_{f,wet} = \dot{m}_a (h_{a,x}-h_{a,out})

.. math::
    :label: eq-Q_rfwet

    \dot{Q}_{f,wet} = \dot{m}_r c_{pr} (T_{r,x} - T_{r,in})

.. math::
    :label: eq-h_fsse

    h_{f,s,sat,eff} = h_{a,x} + \frac{h_{a,out}-h_{a,x}}{1 - \exp(-\frac{(1-f_{dry})U^*_a A_a}{\dot{m}_a})}

.. math::
    :label: eq-T_faout

    T_{a,out} = T_{f,s,eff} + (T_{a,x} - T_{f,s,eff})\exp(-\frac{U_a A_a }{\dot{m}_a c_{p,a}})

.. math::
    :label: eq-ep_fwet

    \varepsilon_{f,wet}=\frac{1-\exp(-(1-f_{dry})\mathrm{Ntu}_{wet}(1-\dot{m}^*))}{1-\dot{m}^*\exp(-(1-f_{dry})\mathrm{Ntu}_{wet}(1-\dot{m}^*))}

At the intersection between the dry and wet region, the heat transfer is governed as Equation :eq:`eq-HT_pdpw`.

.. math::
    :label: eq-HT_pdpw

    \mathrm{UA}_o(T_{a,x} - T_{r,x}) = U_a A_a (T_{a,x}-T_{dp})

Unlike the previous analyses, these equations cannot be solved analytically because only the inlet conditions of refrigerant and air are known. To solve equations iteratively, a bounded solver on :math:`f_{dry}` can be used because :math:`f_{dry}` is proved to be between 0 and 1 from the previous analysis. One way is to calculate the refrigerant outlet temperature based on the dry region only and on both wet and dry region and compare the error between the two methods. Should the error close to zero, the heat exchanger is solved. In this case, one can derive a function which when the solution is reached, the function equals to zero as Equation :eq:`eq-solution`.

.. math::
    :label: eq-solution

    g(x_{true}) = 0


The function :math:`g` is defined as Equation :eq:`eq-sol_T` as a function of :math:`f_{dry}` and a solution of :math:`f_{dry}` is said to be found if :math:`g` equals some value very close to zero.

.. math::
    :label: eq-sol_T

    g(f_{dry}) = T_{r,out}(f_{dry})|_{\mbox{from dry region only}} - T_{r,out}(f_{dry})|_{\mbox{from both regions}}

The following subsections describe how to find :math:`T_{r,out}(f_{dry})|_{\mbox{from dry region only}}` and :math:`T_{r,out}(f_{dry})|_{\mbox{from both regions}}`.

Solving Refrigerant Outlet Temperature with Dry Region Only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the dry region only, a simplification is first conducted on Equation :eq:`eq-ep_fdry`.

.. math::
    :label: eq-step01

    B = \exp(-f_{dry}\mathrm{Ntu}_{dry}(1-C_{ratio}))
    
.. math::
    :label: eq-step02

    B = \frac{1-\varepsilon_{f,dry}}{1-C_{ratio}\varepsilon_{f,dry}}

A combination of of Equations :eq:`eq-Qep_dry` and  :eq:`eq-Qrf_dry` is done.

.. math::
    :label: eq-step2b

    \dot{m}_r c_{p,r} (T_{r,out}-T_{r,x}) = \varepsilon_{f,dry} C_{min} (T_{a,in} - T_{r,x}) \\

.. math::
    :label: eq-step03
    
    T_{r,x} = \frac{\dot{m}_r c_{p,r}T_{r,out}-\varepsilon_{f,dry} C_{min}T_{a,in}}{\dot{m}_r c_{p,r}-\varepsilon_{f,dry} C_{min}}

The terms in Equation :eq:`eq-HT_pdpw` can be rearranged to Equation :eq:`eq-step04`.

.. math::
    :label: eq-step04
    
    T_{a,x} = \frac{U_a A_a T_{r,x}-\mathrm{UA}_o T_{dp}}{U_aA_a-\mathrm{UA}_o}


The Equations :eq:`eq-Qep_dry` and :eq:`eq-Qaf_dry` can also be combined together to form Equation :eq:`eq-step05b` through Equation :eq:`eq-step05`.

.. math::
    :label: eq-step05

    \dot{m}_a c_{p,a} (T_{a,in} - T_{a,x}) = \varepsilon_{f,dry} C_{min} (T_{a,in} - T_{r,x}) 

.. math::
    :label: eq-step05b

    \varepsilon_{f,dry} = \frac{\dot{m}_a c_{p,a}}{C_{min}} \frac{T_{a,in} - T_{a,x}}{T_{a,in} - T_{r,x}}

Equation in :eq:`eq-step03` is re-arranged to form Equation :eq:`eq-step06` through the step in Equation :eq:`eq-step06b`.

.. math::
    :label: eq-step06b

    T_{a,in} - T_{r,x} = T_{a,in} - \frac{\dot{m}_r c_{p,r}T_{r,out}-\varepsilon_{f,dry} C_{min}T_{a,in}}{\dot{m}_r c_{p,r}-\varepsilon_{f,dry} C_{min}} 
    
.. math::
    :label: eq-step06
    
    T_{a,in} - T_{r,x} = \frac{\dot{m}_r c_{p,r}}{\dot{m}_r c_{p,r}-\varepsilon_{f,dry}C_{min}}(T_{a,in} - T_{r,out}) 
    

Equation :eq:`eq-step04` can also be arranged in a similar form as Equation :eq:`eq-step08`.

.. math::
    :label: eq-step08a
    
    T_{a,in} - T_{a,x} = T_{a,in} - \frac{U_aA_a T_{r,out}-\mathrm{UA}_o T_{dp}}{U_aA_a-\mathrm{UA}_o} 

.. math::
    :label: eq-step08
    
    T_{a,in} - T_{a,x} = \frac{U_aA_a(T_{a,in} - T_{r,x}) - \mathrm{UA}_o(T_{a,in} - T_{dp})}{U_aA_a - \mathrm{UA}_o}

Dividing Equation :eq:`eq-step08` by Equation :eq:`eq-step06` can construct Equation :eq:`eq-step09`.

.. math::
    :label: eq-step09
    
    \frac{T_{a,in} - T_{a,x}}{T_{a,in} - T_{r,x}} = \frac{U_aA_a-\mathrm{UA}_o\frac{T_{a,in} - T_{dp}}{T_{a,in} - T_{r,x}}}{U_aA_a-\mathrm{UA}_o}

The terms in Equation :eq:`eq-step06` can be arranged again as Equation :eq:`eq-step10` so the left-hand side of the equation is the same as Equation :eq:`eq-step09`.

.. math::
    :label: eq-step10
    
    \frac{T_{a,in} - T_{a,x}}{T_{a,in} - T_{r,x}} = \frac{\dot{m}_r c_{p,r} - \varepsilon_{f,dry} C_{min}}{\dot{m}_r c_{p,r}} \frac{T_{a,in}-T_{dp}}{T_{a,in} - T_{r,out}}

An :math:`\mathrm{Ntu}_o` can be defined with :math:`\mathrm{UA}_o` as Equation :eq:`eq-step11`.

.. math::
    :label: eq-step11
    
    \mathrm{Ntu}_o = \frac{\mathrm{UA}_o}{\dot{m}_a c_{p,a}}


The :math:`\mathrm{Ntu}_o` from Equation :eq:`eq-step11` is substituted into Equation :eq:`eq-step09` together with :eq:`eq-Ntu` to form :eq:`eq-step12` so that the equation can be defined in terms of :math:`\mathrm{Ntu}`\ s rather than :math:`\mathrm{UA}`\ s.

.. math::
    :label: eq-step12
    
    \frac{T_{a,in} - T_{a,x}}{T_{a,in} - T_{r,x}} = \frac{C_{min}\mathrm{Ntu}_{dry}-\dot{m}_a c_{p,a} \mathrm{Ntu}_o\frac{T_{a,in} - T_{dp}}{T_{a,in} - T_{r,x}}}{C_{min}\mathrm{Ntu}_{dry}-\dot{m}_a c_{p,a} \mathrm{Ntu}_o}

For convenience, another dimensionless variable :math:`\Gamma` is defined in Equation :eq:`eq-step13`.

.. math::
    :label: eq-step13

    \Gamma = \frac{T_{a,in} - T_{a,x}}{T_{a,in} - T_{r,x}}

By combining Equations :eq:`eq-step10` and :eq:`eq-step12` together, one can express :math:`\Gamma` in Equation :eq:`eq-step13` as Equation :eq:`eq-step14`.

.. math::
    :label: eq-step14

    \Gamma = \frac{C_{min}\mathrm{Ntu}_{dry} - \mathrm{Ntu}_o \dot{m}_a c_{p,a} \frac{\dot{m}_r c_{p,r} - \varepsilon_{f,dry} C_{min}}{\dot{m}_r c_{p,r}}\frac{T_{a,in}-T_{dp}}{T_{a,in} - T_{r,out}}}{C_{min}\mathrm{Ntu}_{dry} - \mathrm{Ntu}_o \dot{m}_a c_{p,a}} 

Equation :eq:`eq-step13` can also be formulated as Equation :eq:`eq-step15` with Equation :eq:`eq-step06`.

.. math::
    :label: eq-step15
    
    \Gamma = \frac{C_{min}\varepsilon_{f,dry}}{\dot{m}_a c_{p,a}}

Equations :eq:`eq-step14` and :eq:`eq-step15` can be equated together and by rearranging the subject as :math:`\varepsilon_{f,dry}`, one can form Equation :eq:`eq-step16`.

.. math::
    :label: eq-step16
    
    \varepsilon_{f,dry} = \frac{\dot{m}_r c_{p,r} \dot{m}_a c_{p,a} [C_{min}\mathrm{Ntu}_{dry}(T_{a,in} - T_{r,out}) - \dot{m}_a c_{p,a} \mathrm{Ntu}_o(T_{a,in} - T_{dp})]}{\dot{m}_r c_{p,r} C^2_{min} \mathrm{Ntu}_{dry}(T_{a,in} - T_{r,out}) - C_{min} \mathrm{Ntu}_o \dot{m}_a c_{p,a}[\dot{m}_r c_{p,r}(T_{a,in} - T_{r,out})+\dot{m}_a c_{p,a} (T_{a,in} - T_{dp})]} 

Further calculation can change the subject of Equation :eq:`eq-step16` as the form of the right hand side of Equation :eq:`eq-step02` to establish Equation :eq:`eq-step17`.

.. math::
    :label: eq-step17

    \frac{1-\varepsilon_{f,dry}}{1-C_{ratio}\varepsilon_{f,dry}} = \frac{\Xi_1+ \Xi_{02}}{\Xi_{03} + \Xi_{04}}

    \Xi_{01} = (\dot{m}_a c_{p,a})^2 \mathrm{Ntu}_o (\dot{m}_r c_{p,r} - C_{min})(T_{a,in} - T_{dp})

    \Xi_{02} = C_{min}\dot{m}_r c_{p,r} [\mathrm{Ntu}_{dry}(C_{min}-\dot{m}_a c_{p,a})+ \mathrm{Ntu}_o \dot{m}_a c_{p,a}](T_{a,in} - T_{r,out})

    \Xi_{03} = C_{min} \dot{m}_r c_{p,r}[\mathrm{Ntu}_{dry}(C_{min} - C_{ratio} \dot{m}_a c_{p,a}) - \mathrm{Ntu}_o \dot{m}_a c_{p,a}](T_{a,in} - T_{r,out})

    \Xi_{04} = (\dot{m}_a c_{p,a})^2 \mathrm{Ntu}_o [C_{ratio} \dot{m}_r c_{p,r} - C_{min}] (T_{a,in} - T_{dp})

Equation :eq:`eq-step02` can also be rearranged as Equation :eq:`eq-step18`.

.. math::
    :label: eq-step18
    
    f_{dry} = -\frac{1}{(1-C_{ratio})\mathrm{Ntu}_{dry}}\ln B

Another term can be defined as :math:`K` in Equation :eq:`eq-step19`.

.. math::
    :label: eq-step19
    
    K = (1-C_{ratio})\mathrm{Ntu}_{dry}

The subject in Equation :eq:`eq-step17` can be written as :math:`B` from Equation :eq:`eq-step02` to Equation :eq:`eq-step20`.

.. math::
    :label: eq-step20
    
    B  = \frac{\Xi_1+ \Xi_{02}}{\Xi_{03} + \Xi_{04}}

using the definitions from Equation :eq:`eq-step17`. Equation :eq:`eq-step20` can be simplified should :math:`C_{ratio}` and :math:`C_{min}` be known. If :math:`C_{min} = \dot{m}_a c_{p,a}`, Equation :eq:`eq-step20` can be written as Equation :eq:`eq-step21`.

.. math::
    :label: eq-step21
    
    B = \frac{(T_{dp} - T_{r,out})+C_{ratio}(T_{a,in}-T_{dp})}{[1-\frac{\mathrm{Ntu}_{dry}}{\mathrm{Ntu}_o}(1-C_{ratio})](T_{a,in}-T_{r,out})}

With the definitions of :math:`K` and :math:`f_{dry}` in Equations :eq:`eq-step18` and :eq:`eq-step19`, one can write Equation :eq:`eq-step21` into Equation :eq:`eq-step22`.

.. math::
    :label: eq-step22

    f_{dry} = -\frac{1}{K}\ln\frac{(T_{dp} - T_{r,out})+C_{ratio}(T_{a,in}-T_{dp})}{[1-\frac{K}{\mathrm{Ntu}_o}](T_{a,in}-T_{r,out})}

The refrigerant outlet temperature in this case can be calculated as Equation :eq:`eq-step24` from Equation :eq:`eq-step22`.

.. math::
    :label: eq-step24
    
    T_{r,out} = \frac{T_{dp} - \exp(-Kf_{dry})(1-\frac{K}{\mathrm{Ntu}_o})T_{a,in} + C_{ratio} (T_{a,in} - T_{dp})}{1-\exp(-Kf_{dry})(1-\frac{K}{\mathrm{Ntu}_o})}

Similarly, when :math:`C_{min} = \dot{m}_{r} c_{p,r}`, Equation :eq:`eq-step20` will be written as Equation :eq:`eq-step25`.

.. math::
    :label: eq-step25
    
    B = \frac{C_{ratio}[1+ \frac{\mathrm{Ntu}_{dry}}{\mathrm{Ntu}_o}(1-C_{ratio})](T_{ain}-T_{r,out})}{C_{ratio}(T_{dp} - T_{r,out})+(T_{a,in}-T_{dp})}

Similar derivation can be made on :math:`f_{dry}` and :math:`T_{r,out}` in this case to form Equations :eq:`eq-step26` and :eq:`eq-step27`.

.. math::
    :label: eq-step26

    f_{dry} = -\frac{1}{K}\ln\frac{C^* [1 + \frac{K}{\mathrm{Ntu}_o}](T_{a,in}-T_{r,out})}{C_{ratio}(T_{dp} - T_{r,out})+(T_{a,in}-T_{dp})}
    
.. math::
    :label: eq-step27
    
    T_{r,out} = \frac{\exp(-Kf_{dry})[T_{a,in} + (C_{ratio}-1)T_{dp}] - C_{ratio}(1+\frac{K}{\mathrm{Ntu}_o})T_{a,in}}{C_{ratio}\exp(-Kf_{dry}) - C_{ratio}(1+\frac{K}{\mathrm{Ntu}_o})}


:math:`T_{r,out}(f_{dry})|_{\mbox{from dry region only}}` in Equation :eq:`eq-sol_T` can then be solved by Equations :eq:`eq-step19`, :eq:`eq-step24` and :eq:`eq-step27`, depending on the value of :math:`C_{min}`.


Solving Refrigerant Outlet Temperature with Both Regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While the method to find :math:`T_{r,out}(f_{dry})|_{\mbox{from dry region only}}` is depicted in the previous section, the solution of :math:`T_{r,out}(f_{dry})|_{\mbox{from both regions}}` is described below. Equations :eq:`eq-Qep_dry` and :eq:`eq-Qaf_dry` can be combined to yield

.. math::
    :label: eq-step001

    T_{a,x} = T_{a,i} - \varepsilon_{f,dry} \frac{C_{min}}{\dot{m}_a c_{p,a}} (T_{a,in} - T_{r,x})

Because the specific heat of air is taken to be constant, :math:`\Delta h=c_p\Delta T`.  Thus, the enthalpy of the air at the wet-dry interface can be given by

.. math::
    :label: eq-step001a

    h_{a,x} = h_{a,i} - \varepsilon_{f,dry} \frac{C_{min}}{\dot{m}_a} (T_{a,in} - T_{r,x})
    
Through the use of Equation :eq:`eq-step001a`, the :math:`\varepsilon-\mathrm{Ntu}` equation on the wet region Equation :eq:`eq-Qep_wet` can be expressed as

.. math::
    :label: eq-step002
    
    \dot{Q}_{f,wet} = \varepsilon_{f,wet} \dot{m}_{min} (h_{a,in} - \varepsilon_{f,dry} \frac{C_{min}}{\dot{m}_a} (T_{a,in} - T_{r,x})  - h_{sat,r,in})

The value of :math:`\dot{Q}_{f,wet}` can be substituted from Equation :eq:`eq-Q_rfwet` which yields the value for :math:`T_{r,x}` of

.. math::
    :label: eq-step003

    T_{r,x} = \frac{T_{r,in} + \varepsilon_{f,wet}\frac{\dot{m}_{min}}{\dot{m}_r c_{p,r}}(h_{a,in} - h_{s,r,in} -  \varepsilon_{f,dry}\frac{C_{min}}{\dot{m}_a} T_{a,in})}{1-\varepsilon_{f,wet} \varepsilon_{f,dry}\frac{C_{min}\dot{m}_{min}}{\dot{m}_r c_{p,r} \dot{m}_a} }


:math:`T_{r,out}` from this method can be written (by combining Equations :eq:`eq-Qep_dry` and :eq:`eq-Qrf_dry`) as 

.. math::
    :label: eq-step004

    T_{r,out} = \varepsilon_{f,dry} \frac{C_{min}}{\dot{m}_r c_{p,r}} T_{a,in} + (1-\varepsilon_{f,dry} \frac{C_{min}}{\dot{m}_r c_{p,r}}) T_{r,x}

By solving Equations :eq:`eq-ep_fdry`, :eq:`eq-ep_fwet`, :eq:`eq-step003` and :eq:`eq-step004`, one can find :math:`T_{r,out}(f_{dry})|_{\mbox{from both regions}}` in Equation :eq:`eq-sol_T`. The :math:`g(f_{dry})` in Equation :eq:`eq-sol_T` can be found for diffferent :math:`f_{dry}` and the one which gives a function value close to zero is the numerical solution of :math:`f_{dry}`. With the value of :math:`f_{dry}`, all other variables in the partial-dry-partial-wet analysis can be computed and the heat exchanger performance can be solved.

**Nomenclature**

=============================   ===================================================================
Variable                        Description
=============================   ===================================================================
:math:`A`                       Surface Area [m\ :sup:`2`\ ]
:math:`B`                       Dimensionless variable [--]
:math:`C`                       Capacity Rate [W/K]
:math:`c_{p}`                   Specific Heat Capacity [J/kg/K]
:math:`c_s`                     Analogous specific heat capacity for air-water enthalpy [J/kg/K]
:math:`f`                       Proportion of dry section [--]
:math:`g`                       Function to be solved [any]
:math:`h`                       Air-water mixture enthalpy [J/kg\ :sub:`ha`\ ]
:math:`K`                       Dimensionless variable [--]
:math:`\dot{m}`                 Mass Flow Rate [kg/s]
:math:`\mathrm{Ntu}`            Number of transfer unit [--]
:math:`\mathrm{Ntu_o}`          Overall Number of transfer units [--]
:math:`\dot{Q}`                 Heat Transfer Rate [W]
:math:`T`                       Temperature [K]
:math:`U`                       Heat Transfer Coefficient [W/m\ :sup:`2`\ /K]
:math:`\mathrm{UA}_o`           Overall Heat Conductance [W/K]
:math:`\mathrm{UA}^*_o`         Overall Heat and Mass Transfer Conductance [W/K]
:math:`x_{true}`                Solution of :math:`g(x) = 0` [any]
:math:`\varepsilon`             Heat Exchanger Effectiveness [--]
:math:`\Gamma`                  Dimensionless variable [--]
=============================   ===================================================================

=============================   ===================================================================
Subscript/Superscript           Description
=============================   ===================================================================
:math:`a`                       Of air side
:math:`dp`                      Of dewpoint
:math:`dry`                     Of dry region
:math:`eff`                     Effective
:math:`f`                       Of partial condition
:math:`in`                      At inlet
:math:`local`                   Localized
:math:`min`                     Minimum
:math:`out`                     At outlet
:math:`r`                       Of refrigerant side; in the case of :math:`h`, it means the air-water enthalpy at the temperature of the refrigerant
:math:`ratio`                   Of ratio
:math:`s`                       Of surface
:math:`sat`                     At saturation
:math:`wet`                     Of wet region
:math:`x`                       At intersection
:math:`*`                       Of/Adjusted for mass transfer
=============================   ===================================================================

.. [#Braun] Braun, J. E., 1988. Methodologies for the Design and Control of Central Cooling Plants. Ph.D. thesis, University of Wisconsin - Madison.
