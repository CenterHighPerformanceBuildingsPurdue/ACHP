from __future__ import division, print_function, absolute_import
from CoolProp.CoolProp import HAPropsSI 
from scipy.optimize import fsolve, minimize
from . import Correlations
from math import pi
from .Solvers import MultiDimNewtRaph
import numpy as np
import CoolProp as CP

def DXPreconditioner(Cycle,epsilon=0.96):
    
    #Assume the heat exchangers are highly effective
    #Condensing heat transfer rate from enthalpies
    rho_air=1.1
    Cp_air =1005 #[J/kg/K]
        
    #AbstractState
    AS = Cycle.AS
    
    def OBJECTIVE(x):
        Tevap=x[0]
        Tcond=x[1]
        DT_sh=abs(x[2])
        
        #Use fixed effectiveness to get a guess for the condenser capacity
        Qcond=epsilon*Cycle.Condenser.Fins.Air.Vdot_ha*rho_air*(Cycle.Condenser.Fins.Air.Tdb-Tcond)*Cp_air
        
        AS.update(CP.QT_INPUTS,1.0,Tevap)
        pevap=AS.p() #[pa]
        AS.update(CP.QT_INPUTS,1.0,Tcond)
        pcond=AS.p() #[pa]
        Cycle.Compressor.pin_r=pevap
        Cycle.Compressor.pout_r=pcond
        Cycle.Compressor.Tin_r=Tevap+DT_sh
        Cycle.Compressor.Ref=Cycle.Ref
        Cycle.Compressor.Oil=Cycle.Oil
        Cycle.Compressor.AS=Cycle.AS
        Cycle.Compressor.Calculate()
        W=Cycle.Compressor.W
        
        # Evaporator fully-dry analysis
        Qevap_dry=epsilon*Cycle.Evaporator.Fins.Air.Vdot_ha*rho_air*(Cycle.Evaporator.Fins.Air.Tdb-Tevap)*Cp_air
        
        #Air-side heat transfer UA
        Evap=Cycle.Evaporator
        Evap.mdot_r=Cycle.Compressor.mdot_r
        Evap.psat_r=pevap
        Evap.AS=Cycle.AS
        Evap.Initialize()
        UA_a=Evap.Fins.h_a*Evap.Fins.A_a*Evap.Fins.eta_a
        Tin_a=Evap.Fins.Air.Tdb
        Tout_a=Tin_a+Qevap_dry/(Evap.Fins.mdot_da*Evap.Fins.cp_da)
        #Refrigerant-side heat transfer UA
        UA_r=Evap.A_r_wetted*Correlations.ShahEvaporation_Average(0.5,0.5,AS,Evap.G_r,Evap.ID,Evap.psat_r,Qevap_dry/Evap.A_r_wetted,Evap.Tbubble_r,Evap.Tdew_r)
        #Get wall temperatures at inlet and outlet from energy balance
        T_so_a=(UA_a*Evap.Tin_a+UA_r*Tevap)/(UA_a+UA_r)
        T_so_b=(UA_a*Tout_a+UA_r*Tevap)/(UA_a+UA_r)
        
        Tdewpoint=HAPropsSI('D','T',Cycle.Evaporator.Fins.Air.Tdb, 'P',101325, 'R',Evap.Fins.Air.RH)
        
        #Now calculate the fully-wet analysis
        #Evaporator is bounded by saturated air at the refrigerant temperature.
        h_ai=HAPropsSI('H','T',Cycle.Evaporator.Fins.Air.Tdb, 'P',101325, 'R', Cycle.Evaporator.Fins.Air.RH) #[J/kg_da]
        h_s_w_o=HAPropsSI('H','T',Tevap, 'P',101325, 'R', 1.0) #[J/kg_da]
        Qevap_wet=epsilon*Cycle.Evaporator.Fins.Air.Vdot_ha*rho_air*(h_ai-h_s_w_o)
        
        #Coil is either fully-wet, fully-dry or partially wet, partially dry
        if T_so_a>Tdewpoint and T_so_b>Tdewpoint:
            #Fully dry, use dry Q
            f_dry=1.0
        elif T_so_a<Tdewpoint and T_so_b<Tdewpoint:
            #Fully wet, use wet Q
            f_dry=0.0
        else:
            f_dry=1-(Tdewpoint-T_so_a)/(T_so_b-T_so_a)
        Qevap=f_dry*Qevap_dry+(1-f_dry)*Qevap_wet
        
        if Cycle.ImposedVariable == 'Subcooling': #if Subcooling impose
            AS.update(CP.PT_INPUTS,pcond,Tcond-Cycle.DT_sc_target)
            h_target = AS.hmass() #[J/kg]
            Qcond_enthalpy=Cycle.Compressor.mdot_r*(Cycle.Compressor.hout_r-h_target)
        else: #otherwise, if Charge impose
            AS.update(CP.PT_INPUTS,pcond,Tcond-5)
            h_target = AS.hmass() #[J/kg]
            Qcond_enthalpy=Cycle.Compressor.mdot_r*(Cycle.Compressor.hout_r-h_target)
        
        Qevap_enthalpy=Cycle.Compressor.mdot_r*(Cycle.Compressor.hin_r-h_target)
        
        resids=[Qevap+W+Qcond,Qcond+Qcond_enthalpy,Qevap-Qevap_enthalpy]#,f_dry]
        return resids
    
    Tevap_init=Cycle.Evaporator.Fins.Air.Tdb-15
    Tcond_init=Cycle.Condenser.Fins.Air.Tdb+8
    DT_sh_init=Cycle.Evaporator.DT_sh
    x=fsolve(OBJECTIVE,[Tevap_init,Tcond_init,DT_sh_init])
    DT_evap=Cycle.Evaporator.Fins.Air.Tdb-x[0]
    DT_cond=x[1]-Cycle.Condenser.Fins.Air.Tdb
    DT_sh=abs(x[2])
    
    return DT_evap-1, DT_cond+1, DT_sh

def SecondaryLoopPreconditioner(Cycle,epsilon=0.9):
    rho_air= 1.1 #[kg/m^2]
    Cp_air = 1005 #[J/kg-K]

    #AbstractState
    AS = Cycle.AS
    #AbstractState for SecLoopFluid
    AS_SLF = Cycle.AS_SLF
        
    def OBJECTIVE(x):
        Tevap=x[0]
        Tcond=x[1]
        Tin_CC=x[2]
        if Cycle.Mode=='AC':
            #Condenser heat transfer rate
            Qcond=epsilon*Cycle.Condenser.Fins.Air.Vdot_ha*rho_air*(Cycle.Condenser.Fins.Air.Tdb-Tcond)*Cp_air
            
            #Compressor power
            AS.update(CP.QT_INPUTS,1.0,Tevap)
            pevap=AS.p() #[Pa]
            AS.update(CP.QT_INPUTS,1.0,Tcond)
            pcond=AS.p() #[Pa]
            Cycle.Compressor.pin_r=pevap
            Cycle.Compressor.pout_r=pcond
            Cycle.Compressor.Tin_r=Tevap+Cycle.Compressor.DT_sh
            Cycle.Compressor.Ref=Cycle.Ref
            Cycle.Compressor.AS=Cycle.AS
            Cycle.Compressor.Calculate()
            W=Cycle.Compressor.W
            
            Qcoolingcoil_dry=epsilon*Cycle.CoolingCoil.Fins.Air.Vdot_ha*rho_air*(Cycle.CoolingCoil.Fins.Air.Tdb-Tin_CC)*Cp_air
            
            # Air-side heat transfer UA
            CC=Cycle.CoolingCoil
            #CC.AS_g=AS_SLF
            CC.Initialize()
            UA_a=CC.Fins.h_a*CC.Fins.A_a*CC.Fins.eta_a
            #Air outlet temp from dry analysis
            Tout_a=CC.Tin_a-Qcoolingcoil_dry/(CC.Fins.mdot_a*CC.Fins.cp_a)
            
            # Refrigerant side UA
            f,h,Re=Correlations.f_h_1phase_Tube(Cycle.Pump.mdot_g/CC.Ncircuits, CC.ID, Tin_CC, CC.pin_g, CC.AS_g)
            UA_r=CC.A_g_wetted*h
            #Glycol specific heat
            AS_SLF.update(CP.PT_INPUTS,Cycle.Pump.pin_g,Tin_CC)
            cp_g=AS_SLF.cpmass() #[J/kg-K]
            #Refrigerant outlet temp
            Tout_CC=Tin_CC+Qcoolingcoil_dry/(Cycle.Pump.mdot_g*cp_g)
            
            #Get wall temperatures at inlet and outlet from energy balance
            T_so_a=(UA_a*CC.Tin_a+UA_r*Tout_CC)/(UA_a+UA_r)
            T_so_b=(UA_a*Tout_a+UA_r*Tin_CC)/(UA_a+UA_r)
            
            Tdewpoint=HAPropsSI('D','T',CC.Fins.Air.Tdb,'P',101325,'R',CC.Fins.Air.RH)
            #Now calculate the fully-wet analysis
            #Evaporator is bounded by saturated air at the refrigerant temperature.
            h_ai= HAPropsSI('H','T',CC.Fins.Air.Tdb,'P',101325,'R',CC.Fins.Air.RH)
            h_s_w_o=HAPropsSI('H','T',Tin_CC,'P',101325,'R',1.0)
            Qcoolingcoil_wet=epsilon*CC.Fins.Air.Vdot_ha*rho_air*(h_ai-h_s_w_o)
            
            #Coil is either fully-wet, fully-dry or partially wet, partially dry
            if T_so_a>Tdewpoint and T_so_b>Tdewpoint:
                #Fully dry, use dry Q
                f_dry=1.0
            elif T_so_a<Tdewpoint and T_so_b<Tdewpoint:
                #Fully wet, use wet Q
                f_dry=0.0
            else:
                f_dry=1-(Tdewpoint-T_so_a)/(T_so_b-T_so_a)
            Qcoolingcoil=f_dry*Qcoolingcoil_dry+(1-f_dry)*Qcoolingcoil_wet
        
            Tin_IHX=Tin_CC+Qcoolingcoil/(Cycle.Pump.mdot_g*cp_g)
            QIHX=epsilon*Cycle.Pump.mdot_g*cp_g*(Tin_IHX-Tevap)
            
            AS.update(CP.PT_INPUTS,pcond,Tcond-Cycle.DT_sc_target)
            h_target = AS.hmass() #[J/kg]
            Qcond_enthalpy=Cycle.Compressor.mdot_r*(Cycle.Compressor.hout_r-h_target)
            resids=[QIHX+W+Qcond,Qcond+Qcond_enthalpy,Qcoolingcoil-QIHX]
            return resids
        
        elif Cycle.Mode=='HP':
            #Evaporator heat transfer rate
            Qevap=epsilon*Cycle.Evaporator.Fins.Air.Vdot_ha*rho_air*(Cycle.Evaporator.Fins.Air.Tdb-Tevap)*Cp_air
            
            #Compressor power
            AS.update(CP.QT_INPUTS,1.0,Tevap)
            pevap=AS.p() #[pa]
            AS.update(CP.QT_INPUTS,1.0,Tcond)
            pcond=AS.p() #[pa]
            Cycle.Compressor.pin_r=pevap
            Cycle.Compressor.pout_r=pcond
            Cycle.Compressor.Tin_r=Tevap+Cycle.Evaporator.DT_sh
            Cycle.Compressor.Ref=Cycle.Ref
            Cycle.Compressor.AS=Cycle.AS
            Cycle.Compressor.Calculate()
            W=Cycle.Compressor.W
            
            #Evaporator will be dry
            Qcoolingcoil=epsilon*Cycle.CoolingCoil.Fins.Air.Vdot_ha*rho_air*(Tin_CC-Cycle.CoolingCoil.Fins.Air.Tdb)*Cp_air
            
            #Glycol specifi heat
            AS_SLF.update(CP.PT_INPUTS,Cycle.Pump.pin_g,Tin_CC)
            cp_g=AS_SLF.cpmass() #[J/kg/K]
            
            Tin_IHX=Tin_CC-Qcoolingcoil/(Cycle.Pump.mdot_g*cp_g)
            QIHX=epsilon*Cycle.Pump.mdot_g*cp_g*(Tin_IHX-Tcond)
            
            AS.update(CP.PT_INPUTS,pcond,Tcond-Cycle.DT_sc_target)
            h_target = AS.hmass() #[J/kg]
            QIHX_enthalpy=Cycle.Compressor.mdot_r*(Cycle.Compressor.hout_r-h_target)
            
            resids=[QIHX+W+Qevap,QIHX+QIHX_enthalpy,Qcoolingcoil+QIHX]
            return resids
          
    solverFunc=fsolve
    if Cycle.Mode=='AC':
        Tevap_init=Cycle.CoolingCoil.Fins.Air.Tdb-15
        Tcond_init=Cycle.Condenser.Fins.Air.Tdb+8
        Tin_CC=Tevap_init+1
        #First try using the fsolve algorithm
        try:
            x=fsolve(OBJECTIVE,[Tevap_init,Tcond_init,Tin_CC])
        except:
            #If that doesnt work, try the Mult-Dimensional Newton-raphson solver
            try:
                x=MultiDimNewtRaph(OBJECTIVE,[Tevap_init,Tcond_init,Tin_CC])
            except:
                x=[Tevap_init,Tcond_init,284]
        DT_evap=Cycle.CoolingCoil.Fins.Air.Tdb-x[0]
        DT_cond=x[1]-Cycle.Condenser.Fins.Air.Tdb
        Tin_CC=x[2]
    elif Cycle.Mode=='HP':
        Tevap_init=Cycle.Evaporator.Fins.Air.Tdb-8
        Tcond_init=Cycle.CoolingCoil.Fins.Air.Tdb+15
        Tin_CC=Tcond_init-1
        x=solverFunc(OBJECTIVE,[Tevap_init,Tcond_init,Tin_CC])
        DT_evap=Cycle.Evaporator.Fins.Air.Tdb-x[0]
        Tin_CC=x[2]
        DT_cond=x[1]-Tin_CC
    else:
        raise ValueError()
        
    return DT_evap-2,DT_cond+2,Tin_CC
    
