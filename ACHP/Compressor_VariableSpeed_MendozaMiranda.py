from __future__ import division, print_function, absolute_import
from math import fabs
from CoolProp.CoolProp import PropsSI
from ACHP.OilPropLib import *
import CoolProp as CP

class VariableSpeedCompressorClass():
    """
    Compressor Model based on J.M. Mendoza-Miranda et al. "Comparative evaluation of
    R1234yf, R1234ze(E) and R450A as alternatives to R134a in a variable-speed recip comp",
    Energy 114(2016), 753-766
    
    Applicable to reciprocating, rotary, and scroll compressors.

    Required Parameters:
        
    ===========    ==========  ========================================================================
    Variable       Units       Description
    ===========    ==========  ========================================================================
    a_etav          --         A numpy-like list of compressor map coefficients for volumetric eff.
    a_etais         --         A numpy-like list of compressor map coefficients for isentropic eff.
    a_etaoi          --         A numpy-like list of compressor map coefficients for overall eff.
    Ref            N/A         A string representing the refrigerant
    Oil            N/A         A string representing the lubricant oil  
    Tin_r          K           Refrigerant inlet temperature
    pin_r          Pa          Refrigerant suction pressure (absolute)
    pout_r         Pa          Refrigerant discharge pressure (absolute)
    Vdisp          m^3         Displacement volume
    Vdot_ratio     --          Displacement Scale factor
    V_oil_sump     m^3         Total volume of oil sump inside the compressor shell
    shell_pressure N/A         A string defining the shell pressure of the compressor
    N              rpm         Compressor rotational speed
    ===========   ==========  ========================================================================

        
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
            ('Heat Loss Fraction','-',self.fp),
            ('Displacement scale factor','-',self.Vdot_ratio),
            ('Power','W',self.W),
            ('Mass flow rate','kg/s',self.mdot_r),
            ('Inlet Temperature','K',self.Tin_r),
            ('Outlet Temperature','K',self.Tout_r),
            ('Inlet Enthalpy','J/kg',self.hin_r),
            ('Outlet Enthalpy','J/kg',self.hout_r),
            ('Overall isentropic efficiency','-',self.eta_oi),
            ('Pumped flow rate','m^3/s',self.Vdot_pumped),
            ('Ambient heat loss','W',self.Q_amb),
            ('Refrigerant change in oil sump','kg',self.Charge)
         ]
        
    def Calculate(self):
        #AbstractState
        AS = self.AS
        
        # Reference molar mass of R134a
        MM_R134a = 102.032 #[kg/kmol]
        
        # Molar mass of working fluid
        MM = AS.molar_mass() #[kg/mol]
        MM *= 1000 #[kg/kmol]
        
        # Ambient conditions
        Tamb = self.Tamb
        
        # Compressor suction conditions
        pin = self.pin_r
        Tin = self.Tin_r
        
        AS.update(CP.PT_INPUTS, pin, Tin)
        rhoin = AS.rhomass() #[kg/m^3]
        hin = AS.hmass() #[J/kg]
        sin = AS.smass() #[J/kg-K]
        
        # Compressor discharge conditions
        pout = self.pout_r
    
        AS.update(CP.PSmass_INPUTS, pout, sin)
        hout_is = AS.hmass() #[J/kg]
        Tout_is = AS.T() #[K]
    
        # Reference rotational speed
        N_r = 9.334 #[rps]
        
        # Rotational speed
        N = self.N/60 #[rps]

        # Compressor displacement with scale factor
        Vdisp = self.Vdisp*self.Vdot_ratio
        
        # Isentropic enthalpy difference
        deltah_is = hout_is - hin
        
        #Local copies of coefficients
        a_etav=self.a_etav
        a_etais=self.a_etais
        a_etaoi=self.a_etaoi
        
        # pi-groups isentropic efficiency
        pi2_etais = pout/pin
        pi3_etais = N_r/N
        pi4_etais = (N**3*Vdisp)/deltah_is
        pi5_etais = deltah_is*rhoin/pin 
        pi6_etais = 1/fabs((Tin + Tout_is)/2 - Tamb)
        pi7_etais = MM_R134a/MM
        
        # pi-groups volumetric efficiency
        pi2_etav = pout/pin
        pi3_etav = (rhoin/pin)**1.5*N**3*Vdisp
        pi4_etav = N_r/N
        pi5_etav = MM_R134a/MM
    
        # Expressions for volumetric efficiency, isentropic efficiency and 
        # overall isentropic efficiency
        eta_v = a_etav[0]*pi2_etav**a_etav[1]*pi4_etav**a_etav[2]*pi5_etav**a_etav[3]
        eta_is = a_etais[0]*pi2_etais**a_etais[1]*pi3_etais**a_etais[2]*pi4_etais**a_etais[3]*pi6_etais**a_etais[4]
        eta_oi = a_etaoi[0]*pi2_etais**a_etaoi[1]*pi3_etais**a_etaoi[2]*pi6_etais**a_etaoi[3]*pi7_etais**a_etaoi[4]
    
        # Refrigerant mass flow rate
        mdot = eta_v*rhoin*Vdisp*N #[kg/s]
        
        # Discharge enthalpy
        hout = hin + (hout_is - hin)/eta_is
        
        # Discharge temperature and specific entropy
        AS.update(CP.HmassP_INPUTS, hout, pout)
        Tout = AS.T() #[K]
        sout = AS.smass() #[J/kg-K]         
        
        # Power input
        Wdot_el = mdot*(hout_is - hin)/eta_oi
        
        # Heat loss to ambient
        self.Q_amb = mdot*(hout - hin) - Wdot_el
        
        # Outputs
        self.hin_r = hin
        self.sin_r = sin
        self.hout_r = hout
        self.Tout_r = Tout
        self.sout_r = sout
        self.mdot_r = mdot
        self.W = Wdot_el
        self.eta_v = eta_v
        self.eta_is = eta_is
        self.eta_oi = eta_oi
        self.fp = -self.Q_amb/Wdot_el
        self.CycleEnergyIn=Wdot_el*(1-self.fp)
        self.Vdot_pumped= mdot/rhoin
        
        # Estimate refrigerant dissolved in the oil sump
        T_ave = (Tin + self.Tout_r)/2
        if self.shell_pressure == 'high-pressure':
            p_shell = pout
        elif self.shell_pressure == 'low-pressure':
            p_shell = pin
        
        # Solubility fraction
        self.x_Ref,error = Solubility_Ref_in_Liq(self.Ref,self.Oil,T_ave,p_shell/1000)
        
        AS.update(CP.PT_INPUTS, p_shell, T_ave)
        rho_shell = AS.rhomass() #[kg/m^3]

        rhomass_oil = rho_oil(self.Oil,T_ave-273.15)
        self.m_oil = self.V_oil_sump*rhomass_oil
        
        # Amount of refrigerant dissolved in the oil sump
        self.Charge = self.m_oil*self.x_Ref/(1-self.x_Ref)
        
        
if __name__=='__main__':
    
    """
    Example from: J.M. Mendoza-Miranda et al. "Comparative evaluation of
    R1234yf, R1234ze(E) and R450A as alternatives to R134a in a variable-speed recip comp",
    Energy 114(2016), 753-766
    
    ---------------------------------------------------------------------------
    'a_etav':[1,-0.2678,-0.0106,0.7195], #Volumetric eff. coeff.
    'a_etais':[1,0.0753,0.2183,0.0015,0.0972], #Isentropic eff. coeff.
    'a_etaoi':[1,-0.1642,0.2050,0.0659,0.7669], #Overall eff. coeff.
    
    Compressor geometry:
    Displacement (cm^3): 68.1
    Number cylinders: 2
    Piston diameter (m): 0.085
    Stroke (m): 0.060

    Applicative range fro compressor efficiencies proposed:
    Evaporating temperature (K): 260-280
    Condensing temperature (K): 310-330
    Suction temperature (K): 270-300
    Compressor type: Reciprocating
    Compressor rotation speed (rpm): 400-600
    Refrigerants: R134a,R1234yf,R1234ze(E),R450A
    Ambient temperature (K): 288-300
    ---------------------------------------------------------------------------
    
    
    Additional compressors used for validation of the pi-correlations:
    
    Dataset 1
    Compressor type: Reciprocating
    Displacement (cm^3): 170
    Compressor rotation speed (rpm): 1000-3000
    Nominal compressor power consumption: N/A
    Evaporating temperature (K): 259-279
    Condensing temperature (K): 315-355
    Suction temperature (K): 262-285
    Discharge temperature (K): 323-382
    Refrigerant type: R134a
    'a_etav':[1,-0.2678,-0.0106,0.7195], #Volumetric eff. coeff.
    'a_etais':[1,0.0753,0.2183,0.0015,0.0972], #Isentropic eff. coeff.
    'a_etaoi':[1,-0.1642,0.2050,0.0659,0.7669], #Overall eff. coeff.
    
    Dataset 2
    Compressor type: Rotary
    Displacement (cm^3): 15.4
    Compressor rotation speed (rpm): 2850
    Nominal compressor power consumption: 0.55 W
    Evaporating temperature (K): 266-288
    Condensing temperature (K): 313-318
    Suction temperature (K): 270-293
    Discharge temperature (K): 333-350
    Refrigerant type: R134a/R1234yf
    'a_etav':[1.2,-0.2678,-0.0106,0.7195], #Volumetric eff. coeff.
    'a_etais':[0.85,0.0753,0.2183,0.0015,0.0972], #Isentropic eff. coeff.
    'a_etaoi':[1,-0.1642,0.2050,0.0659,0.7669], #Overall eff. coeff.    

    Dataset 3
    Compressor type: Scroll
    Displacement (cm^3): 100.6
    Compressor rotation speed (rpm): 2900
    Nominal compressor power consumption: 3.4 kW
    Evaporating temperature (K): 247-274
    Condensing temperature (K): 316-342
    Suction temperature (K): 290-292
    Discharge temperature (K): 364-415
    Refrigerant type: R22
    'a_etav':[1.35,-0.2678,-0.0106,0.7195], #Volumetric eff. coeff.
    'a_etais':[1,0.0753,0.2183,0.0015,0.0972], #Isentropic eff. coeff.
    'a_etaoi':[1,-0.1642,0.2050,0.0659,0.7669], #Overall eff. coeff.    
    
    """
    
    #Abstract State 
    Ref = 'R134a'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref)
    
    Tsat_ev = 270 #[K]
    DT_sh = 7 + 8 #[K]
    
    AS.update(CP.QT_INPUTS, 1, Tsat_ev)
    pin_r = AS.p() #[Pa]
    
    Tsat_cd = 320 #[K]
    DT_sc = 8 #[K]

    AS.update(CP.QT_INPUTS, 1, Tsat_cd)
    pout_r = AS.p() #[Pa]
    
    # Compressor inlet temperature
    Tin_r = Tsat_ev + DT_sh
    
    for i in range(1):
        kwds={
              'AS': AS, #Abstract state
              'Ref': Ref,
              'Tamb': 290,
              'Tin_r':Tin_r,
              'pin_r':pin_r,
              'pout_r':pout_r,
              'Vdisp': 0.000681, #Displacement volume, m^3
              'Vdot_ratio': 1.0, #Displacement Scale factor
              'N': 500, #Rotational speed, rpm
              'a_etav':[1,-0.2678,-0.0106,0.7195], #Volumetric eff. coeff.
              'a_etais':[1,0.0753,0.2183,0.0015,0.0972], #Isentropic eff. coeff.
              'a_etaoi':[1,-0.1642,0.2050,0.0659,0.7669], #Overall eff. coeff.
              'shell_pressure': 'low-pressure',
              'Oil': 'POE32',
              'V_oil_sump': 0.0,
              }
        Comp=VariableSpeedCompressorClass(**kwds)
        Comp.Calculate()
        print ('Power:', Comp.W,'W')
        print ('Mass flow rate:',Comp.mdot_r,'kg/s')
        print ('Volumetric efficiency:', Comp.eta_v, '-')
        print ('Isentropic efficiency:', Comp.eta_is, '-')
        print ('Overall efficiency:', Comp.eta_oi, '-')
        print ('Refrigerant dissolved in oil sump:', Comp.Charge,'kg')