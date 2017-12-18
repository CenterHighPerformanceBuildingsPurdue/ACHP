from __future__ import division, print_function, absolute_import
from math import log,exp
from CoolProp.CoolProp import HAPropsSI, cair_sat
from ACHP.FinCorrelations import WavyLouveredFins, HerringboneFins, PlainFins
from ACHP.MicroFinCorrelations import MultiLouveredMicroFins

class DWSVals():
    """ 
    Empty Class for passing data with DryWetSegment
    """
    def __init__(self):
        #Don't do anything
        pass

def DryWetSegment(DWS):
    """
    Generic solver function for dry-wet mixed surface conditions for a given element.
    Can handle superheated, subcooled and two-phase regions.
    Does not handle the pressure drops, only HT required to get the dry/wet interface
    """
    
    #List of required parameters
    RequiredParameters=['Tin_a','h_a','cp_da','eta_a','A_a','pin_a','RHin_a','Tin_r','pin_r','h_r','cp_r','A_r','Rw','mdot_r','Fins','FinsType']
    
    #Check that all the parameters are included, raise exception otherwise
    for param in RequiredParameters:
        if not hasattr(DWS,param):
            raise AttributeError("Parameter "+param+" is required for DWS class in DryWetSegment")
    
    #Retrieve values from structures defined above
    Tin_a=DWS.Tin_a
    if DWS.h_a<0.000000001:
        print ("Warning: Dws.h_a was constrained to 0.001, original value: ", DWS.h_a)
        h_a=0.000000001
    else:
        h_a=DWS.h_a
    cp_da=DWS.cp_da
    eta_a=DWS.eta_a  #from fin correlations, overall airside surface effectiveness
    A_a=DWS.A_a
    pin_a=DWS.pin_a
    RHin_a=DWS.RHin_a
    mdot_da=DWS.mdot_da

    Tin_r=DWS.Tin_r
    pin_r=DWS.pin_r
    if DWS.h_r<0.000000001:
        print ("Warning: Dws.h_r was constrained to 0.001, original value: ", DWS.h_r)
        h_r=0.000000001
    else:
        h_r=DWS.h_r
    cp_r=DWS.cp_r
    A_r=DWS.A_r
    mdot_r=DWS.mdot_r
    Rw=DWS.Rw
    
    #Calculate the dewpoint (amongst others)
    omega_in=HAPropsSI('W','T',Tin_a,'P',pin_a,'R',RHin_a)
    Tdp=HAPropsSI('D','T',Tin_a,'P',pin_a,'W',omega_in)
    hin_a=HAPropsSI('H','T',Tin_a,'P',pin_a,'W',omega_in) #[J/kg_da]
    
    # Internal UA between fluid flow and outside surface (neglecting tube conduction)
    UA_i=h_r*A_r #[W/K], from Shah or f_h_1phase_Tube-fct -> Correlations.py
    # External UA between wall and free stream
    UA_o=eta_a*h_a*A_a #[W/K], from fin correlations
    # wall UA
    UA_w=1/Rw
    # Internal Ntu
    Ntu_i=UA_i/(mdot_r*cp_r)   #[-]
    # External Ntu (multiplied by eta_a since surface is finned and has lower effectiveness)
    Ntu_o=eta_a*h_a*A_a/(mdot_da*cp_da) #[-]
    

    if DWS.IsTwoPhase: #(Two-Phase analysis)
        UA=1/(1/(h_a*A_a*eta_a)+1/(h_r*A_r)+1/UA_w); #overall heat transfer coefficient
        Ntu_dry=UA/(mdot_da*cp_da); #Number of transfer units
        epsilon_dry=1-exp(-Ntu_dry);  #since Cr=0, e.g. see Incropera - Fundamentals of Heat and Mass Transfer, 2007, p. 690
        Q_dry=epsilon_dry*mdot_da*cp_da*(Tin_a-Tin_r);
        Tout_a=Tin_a-Q_dry/(mdot_da*cp_da); #outlet temperature, dry fin

        T_so_a=(UA_o*Tin_a+UA_i*Tin_r)/(UA_o+UA_i); #inlet surface temperature (neglect wall thermal conductance)
        T_so_b=(UA_o*Tout_a+UA_i*Tin_r)/(UA_o+UA_i);  #outlet surface temperature (neglect wall thermal conductance)

        if T_so_b>Tdp:
            #All dry, since surface at outlet dry
            f_dry=1.0
            Q=Q_dry #[W]
            Q_sensible=Q #[W]
            hout_a=hin_a-Q/mdot_da #[J/kg_da]
            # Air outlet humidity ratio
            DWS.omega_out = omega_in #[kg/kg]
        else:
            if T_so_a<Tdp:
                #All wet, since surface at inlet wet
                f_dry=0.0
                Q_dry=0.0
                T_ac=Tin_a #temp at onset of wetted wall
                h_ac=hin_a #enthalpy at onset of wetted surface
            else:
                # Partially wet and dry (i.e T_so_b<Tdp<T_so_a)

                # Air temperature at the interface between wet and dry surface
                # Based on equating heat fluxes at the wall which is at dew point UA_i*(Tw-Ti)=UA_o*(To-Tw)
                T_ac = Tdp + UA_i/UA_o*(Tdp - Tin_r)
                # Dry effectiveness (minimum capacitance on the air side by definition)
                epsilon_dry=(Tin_a-T_ac)/(Tin_a-Tin_r)
                # Dry fraction found by solving epsilon=1-exp(-f_dry*Ntu) for known epsilon from above equation
                f_dry=-1.0/Ntu_dry*log(1.0-epsilon_dry)
                # Enthalpy, using air humidity at the interface between wet and dry surfaces, which is same humidity ratio as inlet
                h_ac=HAPropsSI('H','T',T_ac,'P',pin_a,'W',omega_in) #[J/kg_da]
                # Dry heat transfer
                Q_dry=mdot_da*cp_da*(Tin_a-T_ac)

            # Saturation specific heat at mean water temp (c_s : partial derivative dh_sat/dT @ Tsat_r)
            c_s=cair_sat(Tin_r)*1000  #[J/kg-K]
            # Find new, effective fin efficiency since cs/cp is changed from wetting
            # Ratio of specific heats [-]
            DWS.Fins.Air.cs_cp=c_s/cp_da
            DWS.Fins.WetDry='Wet'
            
            #Compute the fin efficiency based on the user choice of FinsType
            if DWS.FinsType == 'WavyLouveredFins':
                WavyLouveredFins(DWS.Fins)
            elif DWS.FinsType == 'HerringboneFins':
                HerringboneFins(DWS.Fins)
            elif DWS.FinsType == 'PlainFins':
                PlainFins(DWS.Fins)
            elif DWS.FinsType == 'MultiLouveredMicroFins':
                MultiLouveredMicroFins(DWS.Fins)
            
            eta_a_wet=DWS.Fins.eta_a_wet
            UA_o=eta_a_wet*h_a*A_a
            Ntu_o=eta_a_wet*h_a*A_a/(mdot_da*cp_da)
                
            # Wet analysis overall Ntu for two-phase refrigerant
            # Minimum capacitance rate is by definition on the air side
            # Ntu_wet is the NTU if the entire two-phase region were to be wetted
            UA_wet=1/(c_s/UA_i+cp_da/UA_o+c_s/UA_w)
            Ntu_wet=UA_wet/(mdot_da)
            # Wet effectiveness [-]
            epsilon_wet=1-exp(-(1-f_dry)*Ntu_wet)
            # Air saturated at refrigerant saturation temp [J/kg]
            h_s_s_o=HAPropsSI('H','T',Tin_r, 'P',pin_a, 'R', 1.0) #[kJ/kg_da]
            
            # Wet heat transfer [W]
            Q_wet=epsilon_wet*mdot_da*(h_ac-h_s_s_o)
            # Total heat transfer [W]
            Q=Q_wet+Q_dry
            # Air exit enthalpy [J/kg]
            hout_a=h_ac-Q_wet/mdot_da
            # Saturated air temp at effective surface temp [J/kg_da]
            h_s_s_e=h_ac-(h_ac-hout_a)/(1-exp(-(1-f_dry)*Ntu_o))
            # Effective surface temperature [K]
            T_s_e = HAPropsSI('T','H',h_s_s_e,'P',pin_a,'R',1.0)
            # Outlet dry-bulb temp [K]
            Tout_a = T_s_e+(T_ac-T_s_e)*exp(-(1-f_dry)*Ntu_o)
            #Sensible heat transfer rate [kW]
            Q_sensible=mdot_da*cp_da*(Tin_a-Tout_a)
        #Outlet is saturated vapor
        Tout_r=DWS.Tdew_r 
            
    else: #(Single-Phase analysis)
        #Overall UA
        UA = 1 / (1 / (UA_i) + 1 / (UA_o) + 1 / (UA_w));
        # Min and max capacitance rates [W/K]
        Cmin = min([cp_r * mdot_r, cp_da * mdot_da])
        Cmax = max([cp_r * mdot_r, cp_da * mdot_da])
        # Capacitance rate ratio [-]
        C_star = Cmin / Cmax
        # Ntu overall [-]
        Ntu_dry = UA / Cmin
        
        if Ntu_dry<0.0000001:
            print("warning:  NTU_dry in dry wet segment was negative. forced it to positive value of 0.001!")
            Ntu_dry=0.0000001

        # Counterflow effectiveness [-]
        #epsilon_dry = ((1 - exp(-Ntu_dry * (1 - C_star))) / 
        #   (1 - C_star * exp(-Ntu_dry * (1 - C_star))))
        
        #Crossflow effectiveness (e.g. see Incropera - Fundamentals of Heat and Mass Transfer, 2007, p. 662)
        if (cp_r * mdot_r)<(cp_da * mdot_da):
            epsilon_dry= 1-exp(-C_star**(-1)*(1-exp(-C_star*(Ntu_dry))))
            #Cross flow, single phase, cmax is airside, which is unmixed
        else:
            epsilon_dry=(1/C_star)*(1-exp(-C_star*(1-exp(-Ntu_dry))))
            #Cross flow, single phase, cmax is refrigerant side, which is mixed

        # Dry heat transfer [W]
        Q_dry = epsilon_dry*Cmin*(Tin_a-Tin_r)
        # Dry-analysis air outlet temp [K]
        Tout_a_dry=Tin_a-Q_dry/(mdot_da*cp_da)
        # Dry-analysis outlet temp [K]
        Tout_r=Tin_r+Q_dry/(mdot_r*cp_r)
        # Dry-analysis air outlet enthalpy from energy balance [J/kg]
        hout_a=hin_a-Q_dry/mdot_da
        # Dry-analysis surface outlet temp [K] (neglect wall thermal conductance)
        Tout_s=(UA_o*Tout_a_dry+UA_i*Tin_r)/(UA_o+UA_i)
        # Dry-analysis surface inlet temp [K] (neglect wall thermal conductance)
        Tin_s=(UA_o*Tin_a+UA_i*Tout_r)/(UA_o+UA_i)
        # Dry-analysis outlet refrigerant temp [K]
        Tout_r_dry=Tout_r
        # Dry fraction [-]
        f_dry=1.0
        # Air outlet humidity ratio [-]
        DWS.omega_out = omega_in

        # If inlet surface temp below dewpoint, whole surface is wetted 
        if Tin_s<Tdp:
            isFullyWet=True
        else:
            isFullyWet=False

        if Tout_s<Tdp or isFullyWet:
            # There is some wetting, either the coil is fully wetted or partially wetted 

            # Loop to get the correct c_s 
            # Start with the inlet temp as the outlet temp
            x1=Tin_r+1 #Lowest possible outlet temperature
            x2=Tin_a-1 #Highest possible outlet temperature
            eps=1e-8
            iter=1
            change=999
            while ((iter<=3 or change>eps) and iter<100):
                if (iter==1):
                    Tout_r=x1;
                if (iter>1):
                    Tout_r=x2;

                Tout_r_start=Tout_r;
                # Saturated air enthalpy at the inlet water temperature [J/kg]
                h_s_w_i=HAPropsSI('H','T',Tin_r,'P', pin_a, 'R', 1.0) #[J/kg_da]
                # Saturation specific heat at mean water temp [J/kg]
                c_s=cair_sat((Tin_r+Tout_r)/2.0)*1000
                # Ratio of specific heats [-]
                DWS.Fins.Air.cs_cp=c_s/cp_da
                # Find new, effective fin efficiency since cs/cp is changed from wetting
                # Based on the user choice of FinsType
                if DWS.FinsType == 'WavyLouveredFins':
                    WavyLouveredFins(DWS.Fins)
                elif DWS.FinsType == 'HerringboneFins':
                    HerringboneFins(DWS.Fins)
                elif DWS.FinsType == 'PlainFins':
                    PlainFins(DWS.Fins)
                elif DWS.FinsType == 'MultiLouveredMicroFins':
                    MultiLouveredMicroFins(DWS.Fins)
                # Effective humid air mass flow ratio
                m_star=mdot_da/(mdot_r*(cp_r/c_s))
                #compute the new Ntu_owet
                Ntu_owet = eta_a*h_a*A_a/(mdot_da*cp_da)
                m_star = min([cp_r * mdot_r/c_s, mdot_da])/max([cp_r * mdot_r/c_s, mdot_da])
                mdot_min = min([cp_r * mdot_r/c_s, mdot_da])
                # Wet-analysis overall Ntu [-] (neglect wall thermal conductance)
                Ntu_wet=Ntu_o/(1+m_star*(Ntu_owet/Ntu_i))
                if(cp_r * mdot_r> c_s * mdot_da):
                    Ntu_wet=Ntu_o/(1+m_star*(Ntu_owet/Ntu_i))
                else:
                    Ntu_wet=Ntu_i/(1+m_star*(Ntu_i/Ntu_owet))
                    
                # Counterflow effectiveness for wet analysis
                epsilon_wet = ((1 - exp(-Ntu_wet * (1 - m_star))) / 
                    (1 - m_star * exp(-Ntu_wet * (1 - m_star))))
                # Wet-analysis heat transfer rate
                Q_wet = epsilon_wet*mdot_min*(hin_a-h_s_w_i)
                # Air outlet enthalpy [J/kg_da]
                hout_a=hin_a-Q_wet/mdot_da
                # Water outlet temp [K]
                Tout_r = Tin_r+mdot_da/(mdot_r*cp_r)*(hin_a-hout_a)
                # Water outlet saturated surface enthalpy [J/kg_da]
                h_s_w_o=HAPropsSI('H','T',Tout_r, 'P',pin_a, 'R', 1.0) #[J/kg_da]
                #Local UA* and c_s
                UA_star = 1/(cp_da/(eta_a*h_a*A_a)+cair_sat((Tin_a+Tout_r)/2.0)*1000*(1/(h_r*A_r)+1/UA_w))
                # Wet-analysis surface temperature [K]
                Tin_s = Tout_r + UA_star/(h_r*A_r)*(hin_a-h_s_w_o)
                # Wet-analysis saturation enthalpy [J/kg_da]
                h_s_s_e=hin_a+(hout_a-hin_a)/(1-exp(-Ntu_owet))
                # Surface effective temperature [K]
                T_s_e=HAPropsSI('T','H',h_s_s_e,'P',pin_a,'R',1.0)
                # Air outlet temp based on effective temp [K]
                Tout_a=T_s_e + (Tin_a-T_s_e)*exp(-Ntu_o)
                #Sensible heat transfer rate [W]
                Q_sensible = mdot_da*cp_da*(Tin_a-Tout_a)
                # Error between guess and recalculated value [K]
                errorToutr=Tout_r-Tout_r_start;
                    
                if(iter>500):
                    print("Superheated region wet analysis T_outr convergence failed")
                    DWS.Q=Q_dry;
                    return
                if iter==1:
                    y1=errorToutr
                if iter>1:
                    y2=errorToutr
                    x3=x2-y2/(y2-y1)*(x2-x1)
                    change=abs(y2/(y2-y1)*(x2-x1))
                    y1=y2; x1=x2; x2=x3
                if hasattr(DWS,'Verbosity') and DWS.Verbosity>7:
                    print("Fullwet iter %d Toutr %0.5f dT %g" %(iter,Tout_r,errorToutr))
                #Update loop counter
                iter+=1
                
            # Fully wetted outlet temperature [K]
            Tout_r_wet=Tout_r
            # Dry fraction
            f_dry=0.0
            
            if (Tin_s>Tdp and not isFullyWet):

                #Partially wet and partially dry with single-phase on refrigerant side
                
                """
                -----------------------------------------------------------
                                            |
                * Tout_a   <----            * T_a,x                 <---- * Tin_a
                                            |
                 ____________Wet____________|              Dry 
                ----------------------------o T_dp ------------------------
                                            |
                * Tin_r    ---->            * T_w,x                 ----> * Tout_r               
                                            |
                -----------------------------------------------------------
                """

                iter = 1
                
                # Now do an iterative solver to find the fraction of the coil that is wetted
                x1=0.0001
                x2=0.9999
                eps=1e-8
                while ((iter<=3 or error>eps) and iter<100):
                    if iter==1:
                        f_dry=x1
                    if iter>1:
                        f_dry=x2
                    
                    K=Ntu_dry*(1.0-C_star)
                    expk = exp(-K*f_dry)
                    if cp_da*mdot_da < cp_r*mdot_r:
                        Tout_r_guess = (Tdp + C_star*(Tin_a - Tdp)-expk*(1-K/Ntu_o)*Tin_a)/(1-expk*(1-K/Ntu_o))
                    else:
                        Tout_r_guess = (expk*(Tin_a+(C_star-1)*Tdp)-C_star*(1+K/Ntu_o)*Tin_a)/(expk*C_star-C_star*(1+K/Ntu_o))

                    # Wet and dry effective effectivenesses
                    epsilon_dry = ((1 - exp(-f_dry*Ntu_dry * (1 - C_star))) / 
                        (1 - C_star * exp(-f_dry*Ntu_dry * (1 - C_star))))
                    epsilon_wet = ((1 - exp(-(1-f_dry)*Ntu_wet * (1 - m_star))) / 
                        (1 - m_star * exp(-(1-f_dry)*Ntu_wet * (1 - m_star))))
                    
                    # Temperature of "water" where condensation begins
                    T_w_x=(Tin_r+(mdot_min)/(cp_r * mdot_r)*epsilon_wet*(hin_a-h_s_w_i-epsilon_dry*Cmin/mdot_da*Tin_a))/(1-(Cmin*mdot_min)/(cp_r * mdot_r * mdot_da)*epsilon_wet*epsilon_dry)
                    # Temperature of air where condensation begins [K]
                    # Obtained from energy balance on air side
                    T_a_x = Tin_a - epsilon_dry*Cmin*(Tin_a - T_w_x)/(mdot_da*cp_da)
                    # Enthalpy of air where condensation begins
                    h_a_x = hin_a - cp_da*(Tin_a - T_a_x)
                    # New "water" temperature (stored temporarily to be able to build change
                    Tout_r=(Cmin)/(cp_r * mdot_r)*epsilon_dry*Tin_a+(1-(Cmin)/(cp_r * mdot_r)*epsilon_dry)*T_w_x
                    # Difference between initial guess and outlet 
                    error=Tout_r-Tout_r_guess
                    
                    if(iter>500):
                        print("Superheated region wet analysis f_dry convergence failed")
                        DWS.Q=Q_dry
                        return
                    if iter==1:
                        y1=error
                    if iter>1:
                        y2=error
                        x3=x2-y2/(y2-y1)*(x2-x1)
                        change=abs(y2/(y2-y1)*(x2-x1))
                        y1=y2; x1=x2; x2=x3;
                    if hasattr(DWS,'Verbosity') and DWS.Verbosity>7:
                        print("Partwet iter %d Toutr_guess %0.5f diff %g f_dry: %g"%(iter,Tout_r_guess,error,f_dry))
                    #Update loop counter
                    iter+=1
                
                # Wet-analysis saturation enthalpy [J/kg]
                h_s_s_e=h_a_x+(hout_a-h_a_x)/(1-exp(-(1-f_dry)*Ntu_owet))
                # Surface effective temperature [K]
                T_s_e=HAPropsSI('T','H',h_s_s_e,'P',pin_a,'R',1.0)
                # Air outlet temp based on effective surface temp [K]
                Tout_a=T_s_e + (T_a_x-T_s_e)*exp(-(1-f_dry)*Ntu_o)
                # Heat transferred [W]
                Q=mdot_r*cp_r*(Tout_r-Tin_r)
                # Dry-analysis air outlet enthalpy from energy balance [J/kg]
                hout_a=hin_a-Q/mdot_da
                #Sensible heat transfer rate [kW]
                Q_sensible = mdot_da*cp_da*(Tin_a-Tout_a)
            else:
                Q=Q_wet
        else:
            # Coil is fully dry
            Tout_a=Tout_a_dry
            Q=Q_dry
            Q_sensible=Q_dry
          

    DWS.f_dry=f_dry
    DWS.omega_out=HAPropsSI('W','T',Tout_a,'P',101325,'H',hout_a)
    DWS.RHout_a=HAPropsSI('R','T',Tout_a,'P',101325,'W',DWS.omega_out)
    DWS.Tout_a=Tout_a
    DWS.Q=Q
    DWS.Q_sensible=Q_sensible
    
    DWS.hout_a=hout_a
    DWS.hin_a=hin_a
    DWS.Tout_r=Tout_r
    DWS.Twall_s=Tout_r - Q/UA_i #inner wall temperature for gas cooler model
