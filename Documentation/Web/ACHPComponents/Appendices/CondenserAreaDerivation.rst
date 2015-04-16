
.. _CondenserAreaDerivation:

Condenser Area Derivation
*************************

For a configuration like the cross-flow Condenser, the area on both refrigerant and air sides is proportional to a parameter :math:`w`, the length fraction for a given circuit.  In addition, the air mass flow rate is proportional to the parameter :math:`w`. As a result, independent of whether the minimum capacitance rate is on the air- or refrigerant-side, the same result for the area fraction :math:`w` is obtained, as shown from the derivation below.

For this derivation, the inlet temperatures of both streams are known, and the outlet stream of the refrigerant is known.  In addition, the mass flow rates of refrigerant and air are known.  Therefore, the actual amount of heat transfer is also known.

.. math::
    :label: eqCAD1

    \mathrm{UA} =w\dfrac{1}{  (\eta_a \alpha_a A_{a,total})^{-1} + (\alpha_{r} A_{r,total})^{-1}}=w\mathrm{UA}_{overall}\\
    
    C_{min} = \min[\dot{m}_r c_{p,r}, w \dot{m}_{a,total} c_{p,a}]

Minimum capacitance rate on air side
====================================

If :math:`C_{min}` on air side:

.. math::
    :label: eqCAD2
    
    \mathrm{Ntu} = \frac{\mathrm{UA}}{wc_{p,a}  \dot m_{a,total}}=\frac{\mathrm{UA}_{overall}}{c_{p,a}\dot m_{a,total}}
 
:math:`w` cancels out, leaving :math:`\mathrm{Ntu}` independent of :math:`w`.  Energy balance yields

.. math::
    :label: eqCAD3
    
    \dot m_rc_{p,r}(T_{r,i}-T_{r,o})=\varepsilon C_{min} (T_{r,i}-T_{a,i})
 
.. math::
    :label: eqCAD4
    
    \dot m_rc_{p,r}(T_{r,i}-T_{r,o})=\varepsilon w c_{p,a}  \dot m_{a,total} (T_{r,i}-T_{a,i})
 
.. math::
    :label: eqCAD5
    
    \frac{(T_{r,i}-T_{r,o})}{(T_{r,i}-T_{a,i})}=\varepsilon C_r
 
LHS is constant, call it :math:`\Psi`.  The minimum capacitance rate is on the air side (:math:`C_{min}=  w\dot m_{a,total} c_{p,a}` ), cross-flow with :math:`C_{max}` mixed (ref.) and :math:`C_{min}` unmixed (air) yields

.. math::
    :label: eqCAD6
    
    \varepsilon = \displaystyle\frac{1}{C_r} (1 - \exp(-C_r (1 - \exp(-\mathrm{Ntu}))))\\

Now solve for :math:`C_r`

.. math::
    :label: eqCAD7
    
    \Psi = (1 - \exp(-C_r (1 - \exp(-\mathrm{Ntu}))))\\
 
.. math::
    :label: eqCAD8
    
    \exp(-C_r (1 - \exp(-\mathrm{Ntu}))) = 1 - \Psi\\
 
.. math::
    :label: eqCAD9
    
    -C_r (1 - \exp(-\mathrm{Ntu})) = \ln(1 - \Psi)\\
 
.. math::
    :label: eqCAD10
    
    C_r  = -\frac{\ln(1 - \Psi)}{(1 - \exp(-\mathrm{Ntu}))}
 
Coming back to the definition of :math:`C_r` as the ratio of capacitance rates, you can get :math:`w` from

.. math::
    :label: eqCAD11
    
    w=C_r\frac{\dot m_rc_{p,r}}{\dot m_{a,total} c_{p,a}}

and since :math:`C_r` is already known, you obtain

.. math::
    :label: eqCAD12
    
    w=-\frac{\ln(1 - \Psi)}{(1 - \exp(-\mathrm{Ntu}))} \frac{\dot m_rc_{p,r}}{\dot m_{a,total} c_{p,a}}
    
.. math::
    :label: eqCAD13
    
    w=-\frac{\ln(1 - \Psi)}{(1 - \exp(-\mathrm{UA}_{overall}/(c_{p,a}\dot m_{a,total})))} \frac{\dot m_rc_{p,r}}{\dot m_{a,total} c_{p,a}}
 

Minimum capacitance rate on refrigerant side
============================================

If :math:`C_{min}` on refrigerant side:

.. math::
    :label: eqCAD14
    
    \mathrm{Ntu} = \frac{\mathrm{UA}}{\dot m_r c_{p,r}}=\frac{w\mathrm{UA}_{overall}}{\dot m_r c_{p,r}}
 
.. math::
    :label: eqCAD15
    
    C_r\mathrm{Ntu} = \frac{\dot m_r c_{p,r}}{w\dot m_{a,total} c_{p,a}}\frac{w\mathrm{UA}_{overall}}{\dot m_r c_{p,r}}
 
.. math::
    :label: eqCAD16
    
    C_r\mathrm{Ntu} = \frac{\mathrm{UA}_{overall}}{\dot m_{a,total} c_{p,a}}
 

Energy balance yields

.. math::
    :label: eqCAD17
    
    \dot m_rc_{p,r}(T_{r,i}-T_{r,o})=\varepsilon C_{min} (T_{r,i}-T_{a,i})
 
.. math::
    :label: eqCAD18
    
    \dot m_rc_{p,r}(T_{r,i}-T_{r,o})=\varepsilon\dot m_rc_{p,r} (T_{r,i}-T_{a,i})
 
.. math::
    :label: eqCAD19
    
    \varepsilon=\frac{(T_{r,i}-T_{r,o})}{(T_{r,i}-T_{a,i})}
 
Right-hand-side is also equal to :math:`\Psi` from above.  Effectiveness with :math:`C_{min}` mixed (ref.) and :math:`C_{max}` unmixed (air) yields

.. math::
    :label: eqCAD20
    
    \varepsilon = 1 - \exp(-\displaystyle\frac{1}{C_r}  (1 - \exp(-C_r \mathrm{Ntu})))
 
.. math::
    :label: eqCAD21
    
    \exp(-\displaystyle\frac{1}{C_r}  (1 - \exp(-C_r \mathrm{Ntu})))  = 1 - \varepsilon
 
.. math::
    :label: eqCAD22
    
    -\dfrac{1}{C_r}  (1 - \exp(-C_r \mathrm{Ntu}))  = \ln(1 - \epsilon)
 
.. math::
    :label: eqCAD23
    
    {C_r}=-\frac{ (1 - \exp(-C_r \mathrm{Ntu}))}{\ln(1 - \varepsilon)}=\frac{\dot m_r c_{p,r}}{w\dot m_{a,total} c_{p,a}}
 
.. math::
    :label: eqCAD24
    
    w=-\frac{\ln(1 - \varepsilon)\dot m_r c_{p,r}}{ (1 - \exp(-C_r \mathrm{Ntu}))\dot m_{a,total} c_{p,a}}
 
.. math::
    :label: eqCAD25
    
    w=-\frac{\ln(1 - \Psi)\dot m_r c_{p,r}}{ [1 - \exp(-\mathrm{UA}_{overall}/(c_{p,a}\dot m_{a,total}))]\dot m_{a,total} c_{p,a}}
    
Thus both assuming that the minimum capacitance rate is on the air- or refrigerant-sides yields exactly the same solution, which conveniently allows for an explicit solution independent of whether the air-side is the limiting capacitance rate or not