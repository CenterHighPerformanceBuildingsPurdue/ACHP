from __future__ import division, print_function, absolute_import
from CoolProp.CoolProp import PropsSI
import CoolProp as CP
from ACHP.convert_units import *

class VICompressorTelloClass():
    """
    Vapor injection Compressor Model based on Tello-Oquendo Model  <http://www.sciencedirect.com/science/article/pii/S0140700716304042>
    
    Required Parameters:
        
    ===========   ==========  ========================================================================
    Variable      Units       Description
    ===========   ==========  ========================================================================
    M             g/s         A numpy-like list of compressor map coefficients for suction mass flow
    P             Watts       A numpy-like list of compressor map coefficients for electrical power
    K             --          A numpy-like list of compressor map coefficients for injection mass flow
    Ref           N/A         A string representing the refrigerant
    Tin_r         K           Refrigerant inlet temperature
    Tinj_r        K           Refrigerant injection temperature
    pin_r         Pa          Refrigerant suction pressure (absolute)
    pout_r        Pa          Refrigerant discharge pressure (absolute)
    pinj_r        Pa          Refrigerant injection pressure (absolute)
    fp            --          Fraction of electrical power lost as heat to ambient
    Vdot_ratio    --          Displacement Scale factor
    ===========   ==========  ========================================================================
    
    All variables are of double-type unless otherwise specified
        
    """
    def __init__(self,**kwargs):
        #Load up the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    def Update(self,**kwargs):
        #Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items with indices:
                [0] Description of value
                
                [1] Units of value
                
                [2] The value itself
        """
        
        return [ 
            ('M1','-',self.M[0]),
            ('M2','-',self.M[1]),
            ('M3','-',self.M[2]),
            ('M4','-',self.M[3]),
            ('M5','-',self.M[4]),
            ('M6','-',self.M[5]),
            ('M7','-',self.M[6]),
            ('M8','-',self.M[7]),
            ('M9','-',self.M[8]),
            ('M10','-',self.M[9]),
            ('P1','-',self.P[0]),
            ('P2','-',self.P[1]),
            ('P3','-',self.P[2]),
            ('P4','-',self.P[3]),
            ('P5','-',self.P[4]),
            ('P6','-',self.P[5]),
            ('P7','-',self.P[6]),
            ('P8','-',self.P[7]),
            ('P9','-',self.P[8]),
            ('P10','-',self.P[9]),
            ('P11','-',self.P[10]),
            ('K0','-',self.K[0]),
            ('K1','-',self.K[1]),
            ('Heat Loss Fraction','-',self.fp),
            ('Displacement scale factor','-',self.Vdot_ratio),
            ('Power','W',self.W),
            ('Suction mass flow rate','kg/s',self.mdot_r),
            ('Injection mass flow rate','kg/s',self.mdot_inj),
            ('Total mass flow rate','kg/s',self.mdot_tot),
            ('Inlet Temperature','K',self.Tin_r),
            ('Injection Temperature','K',self.Tinj_r),
            ('Outlet Temperature','K',self.Tout_r),
            ('Inlet Enthalpy','J/kg',self.hin_r),
            ('Injection Enthalpy','J/kg',self.hinj_r),
            ('Outlet Enthalpy','J/kg',self.hout_r),
            ('Inlet Pressure','Pa',self.pin_r),
            ('Injection Pressure','Pa',self.pinj_r),
            ('Outlet Pressure','Pa',self.pout_r),
            ('Inlet Entropy','J/kg-K',self.sin_r),
            ('Injection Entropy','J/kg-K',self.sinj_r),
            ('Outlet Entropy','J/kg-K',self.sout_r),
            ('Overall isentropic efficiency','-',self.eta_oi),
            ('Pumped flow rate','m**3/s',self.Vdot_pumped),
            ('Ambient heat loss','W',self.Q_amb)
         ]
        
    def Calculate(self):
        #AbstractState
        AS = self.AS
        
        #Local copies of coefficients
        P=self.P
        M=self.M
        K=self.K
        
        #Calculate suction superheat and dew temperatures
        AS.update(CP.PQ_INPUTS, self.pin_r, 1.0)
        self.Tsat_s_K=AS.T() #[K]
        
        AS.update(CP.PQ_INPUTS, self.pout_r, 1.0)
        self.Tsat_d_K=AS.T() #[K]
        #self.DT_sh_K=self.Tin_r-self.Tsat_s_K
        
        AS.update(CP.PQ_INPUTS, self.pinj_r, 1.0)
        self.Tsat_inj_K=AS.T() #[K]
        #self.DT_sh_inj_K=self.Tinj_r-self.Tsat_inj_K
        
        AS.update(CP.PT_INPUTS, self.pinj_r, self.Tinj_r)
        self.hinj_r=AS.hmass() #[J/kg]
        self.sinj_r=AS.smass() #[J/kg-K]
        
        AS.update(CP.PT_INPUTS, self.pin_r, self.Tin_r)
        self.hin_r=AS.hmass() #[J/kg]
        self.sin_r=AS.smass() #[J/kg-K]
        self.vin_r = 1 / AS.rhomass() #[m**3/kg]
        
        #Convert saturation temperatures in K to F
        Tsat_s = self.Tsat_s_K * 9/5 - 459.67
        Tsat_d = self.Tsat_d_K * 9/5 - 459.67
        Tsat_inj = self.Tsat_inj_K * 9/5 - 459.67
    
        #Apply the modified 11 coefficient AHRI map to saturation temps in F
        power_map = P[0] + P[1] * Tsat_s + P[2] * Tsat_d + P[3] * Tsat_s**2 + P[4] * Tsat_s * Tsat_d + P[5] * Tsat_d**2 + P[6] * Tsat_s**3 + P[7] * Tsat_d * Tsat_s**2 + P[8] * Tsat_d**2*Tsat_s + P[9] * Tsat_d**3 + P[10] * Tsat_inj
        #power in Watts
        
        mdot_map = M[0] + M[1] * Tsat_s + M[2] * Tsat_d + M[3] * Tsat_s**2 + M[4] * Tsat_s * Tsat_d + M[5] * Tsat_d**2 + M[6] * Tsat_s**3 + M[7] * Tsat_d * Tsat_s**2 + M[8] * Tsat_d**2*Tsat_s + M[9] * Tsat_d**3
        #mass in g/s convert go kg/s
        mdot_map /= 1000.0 
        
        # Add more mass flow rate to scale
        mdot_map*=self.Vdot_ratio
        power_map*=self.Vdot_ratio
         
        #injection mass flow rate
        mdot_inj = K[0] * mdot_map + K[1] * mdot_map * (self.pinj_r/self.pin_r)
        
        #P1 = self.pin_r
        #P2 = self.pout_r

        #T1_actual = self.Tsat_s_K + self.DT_sh_K
        #T1_map = self.Tsat_s_K + 20 * 5 / 9
    
        #AS.update(CP.PT_INPUTS, P1, T1_map)
        #v_map = 1 / AS.rhomass() #[m^3/kg]
        #s1_map = AS.smass() #[J/kg-K]
        #h1_map = AS.hmass() #[J/kg]
        
        #AS.update(CP.PT_INPUTS, P1, T1_actual)
        #s1_actual = AS.smass() #[J/kg-K]
        #h1_actual = AS.hmass() #[J/kg]
        #v_actual = 1 / AS.rhomass() #[m^3/kg]
        #F = 0.75
        mdot = mdot_map # for now use the mass directly from the map ###(1 + F * (v_map / v_actual - 1)) * mdot_map
        power = power_map #for now use the power directly from  the map without any corrections
        
        #The following enthalpies rae defined based on ASHRAE 23.1 (2015) --NOT yet published  
        #h_31s = specific enthalpy of refrigerant at injection pressure following an isentropic compression of the refrigerant from inlet of the compressor
        AS.update(CP.PSmass_INPUTS, self.pinj_r, self.sin_r)
        h_31s = AS.hmass() #[J/kg]
        
        #h_22s = specific enthalpy of refrigerant vapor after mixing of the intermediate pressure flow at injection state (i.e., state number 5 in ASHRAE 23.1) with partially compressed flow at 31s
        h_22s = (mdot*h_31s+ mdot_inj*self.hinj_r)/(mdot+mdot_inj)
        
        AS.update(CP.HmassP_INPUTS, h_22s, self.pinj_r)
        s_22s = AS.smass() #[J/kg-k]
        
        #h_32s = specific enthalpy of refrigerant vapor at compressor discharge pressure following an isentropic compression of the refrigerant vapor from state point h_22s
        AS.update(CP.PSmass_INPUTS, self.pout_r, s_22s)
        h_32s = AS.hmass() #[J/kg]

        #isentropic discharge enthalpy
        AS.update(CP.PSmass_INPUTS, self.pout_r, self.sin_r)
        h_41s = AS.hmass() #[J/kg] 
        
        #isentropic efficiency (assume same heat loss fraction for each stage, assume adiabatic mixing at state 22)
        self.eta_oi = (mdot*(h_31s-self.hin_r)+ (mdot+mdot_inj)*(h_32s-h_22s))/power
        
        #enthalpy after first stage compression
        h_31 = self.hin_r + (h_31s-self.hin_r)/(self.eta_oi*(1+self.fp)) #including the heat loss
        
        #enthalpy after mixing h_31 and h_inj
        h_22 = (mdot*h_31 + mdot_inj*self.hinj_r)/(mdot+mdot_inj)
        
        #discharge enthalpy (i.e., h_32 in ASHRAE 23.1)
        self.hout_r = h_22 + (h_32s-h_22s)/(self.eta_oi*(1+self.fp)) #including the heat loss
        
        
        #AS.update(CP.PSmass_INPUTS, P4, self.sin_r)
        #h2s_map = AS.hmass() #[J/kg]        
        
        #AS.update(CP.PSmass_INPUTS, P2, s1_actual)
        #h2s_actual = AS.hmass() #[J/kg]
    
        #Shaft power based on 20F superheat calculation from fit overall isentropic efficiency
        #power = power_map * (mdot / mdot_map) * (h2s_actual - h1_actual) / (h2s_map - h1_map)
    
        #h2 = power * (1 - self.fp) / mdot + h1_actual
        #isentropic efficiency
        #self.eta_oi=mdot*(h2s_actual-h1_actual)/(power)
        
        AS.update(CP.HmassP_INPUTS, self.hout_r, self.pout_r)
        self.Tout_r = AS.T() #[K]
        self.sout_r = AS.smass() #[J/kg-K]        
        
        self.mdot_r = mdot
        self.mdot_inj = mdot_inj
        self.mdot_tot = mdot + mdot_inj
        self.W=power
        self.CycleEnergyIn=power*(1-self.fp)
        self.Vdot_pumped= mdot*self.vin_r
        self.Q_amb=-self.fp*power

        #For Plotting:
        self.h_21 = self.hin_r
        self.h_22 = h_22
        self.h_31 = h_31
        self.h_32 = self.hout_r
        self.h_5 = self.hinj_r
        self.h_22s = h_22s
        self.h_31s = h_31s
        self.h_32s = h_32s     
        self.h_41s = h_41s #this should not be used for plotting VI Tello compressor  

        
if __name__=='__main__':        
    import numpy as np
    
    #Abstract State        
    Ref = 'R407C'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref)
    
    Tin_dew = 273.15 #assume 5K superheat
    pin_r = PropsSI('P','T',Tin_dew,'Q',1,Ref)
    pout_r = PropsSI('P','T',333.15,'Q',1,Ref)
    pinj_r = np.sqrt(pin_r * pout_r)
    Tinj_dew = PropsSI('T','P',pinj_r,'Q',1,Ref)

    kwds={
        'K':[-0.3845,0.3296],
        'M':[133.3171898,0.508718380832,-2.15889885692,0.00847246179835,0.009495493298,0.0170835511659,3.65431994895E-05,6.660136064E-06,-4.719716435E-05,-4.61719969253E-05],
        'P':[1.003,0.9998 ,1.134 ,0.9032 ,-0.5003 ,0.5136 ,-0.001366 ,-0.009186 ,0.005756 ,-0.00156 ,1.053],
        'AS':AS,
        'Tin_r':Tin_dew+5,
        'pin_r':pin_r,
        'pout_r':pout_r,
        'pinj_r':pinj_r,
        'Tinj_r':Tinj_dew+5,
        'fp':0.12, #Fraction of electrical power lost as heat to ambient
        'Vdot_ratio': 1.0, #Displacement Scale factor
    }
    Comp=VICompressorTelloClass(**kwds)
    Comp.Calculate()
    print (Comp.W,'W')
    print (Comp.Tout_r,'K')
    print (Comp.mdot_r,'kg/s')
    print (Comp.mdot_inj,'kg/s')
    print (Comp.Q_amb, 'W')
    print (Comp.eta_oi*100, '%')
