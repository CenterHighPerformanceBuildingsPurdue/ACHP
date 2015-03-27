from __future__ import division #Make integer 3/2 give 1.5 in python 2.x
from CoolProp.CoolProp import PropsSI #,T_hp, h_sp

class CompressorClass():
    """
    Compressor Model based on 10-coefficient Model from `ANSI/AHRI standard 540 <http://www.ahrinet.org/App_Content/ahri/files/standards%20pdfs/ANSI%20standards%20pdfs/ANSI-ARI-540-2004%20latest.pdf>`_
    
    Required Parameters:
        
    ===========   ==========  ========================================================================
    Variable      Units       Description
    ===========   ==========  ========================================================================
    M             Ibm/hr      A numpy-like list of compressor map coefficients for mass flow
    P             Watts       A numpy-like list of compressor map coefficients for electrical power
    Ref           N/A         A string representing the refrigerant
    Tin_r         K           Refrigerant inlet temperature
    pin_r         Pa          Refrigerant suction pressure (absolute)
    pout_r        Pa          Refrigerant discharge pressure (absolute)
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
            ('Heat Loss Fraction','-',self.fp),
            ('Displacement scale factor','-',self.Vdot_ratio),
            ('Power','W',self.W),
            ('Mass flow rate','kg/s',self.mdot_r),
            ('Inlet Temperature','K',self.Tin_r),
            ('Outlet Temperature','K',self.Tout_r),
            ('Inlet Enthalpy','J/kg',self.hin_r),
            ('Outlet Enthalpy','J/kg',self.hout_r),
            ('Overall isentropic efficiency','-',self.eta_oi),
            ('Pumped flow rate','m^3/s',self.Vdot_pumped)
            ('Ambient heat loss','W',self.Q_amb)
         ]
        
    def Calculate(self):
        #Local copies of coefficients
        P=self.P
        M=self.M
        
        #Calculate suction superheat and dew temperatures
        self.Tsat_s_K=PropsSI('T','P',self.pin_r,'Q',1.0,self.Ref)
        self.Tsat_d_K=PropsSI('T','P',self.pout_r,'Q',1.0,self.Ref)
        self.DT_sh_K=self.Tin_r-self.Tsat_s_K
        
        #Convert saturation temperatures in K to F
        Tsat_s = self.Tsat_s_K * 9/5 - 459.67
        Tsat_d = self.Tsat_d_K * 9/5 - 459.67
    
        #Apply the 10 coefficient ARI map to saturation temps in F
        power_map = P[0] + P[1] * Tsat_s + P[2] * Tsat_d + P[3] * Tsat_s**2 + P[4] * Tsat_s * Tsat_d + P[5] * Tsat_d**2 + P[6] * Tsat_s**3 + P[7] * Tsat_d * Tsat_s**2 + P[8] * Tsat_d**2*Tsat_s + P[9] * Tsat_d**3
        mdot_map = M[0] + M[1] * Tsat_s + M[2] * Tsat_d + M[3] * Tsat_s**2 + M[4] * Tsat_s * Tsat_d + M[5] * Tsat_d**2 + M[6] * Tsat_s**3 + M[7] * Tsat_d * Tsat_s**2 + M[8] * Tsat_d**2*Tsat_s + M[9] * Tsat_d**3
    
        # Convert mass flow rate to kg/s from lbm/h
        mdot_map *= 0.000125998 
    
        # Add more mass flow rate to scale
        mdot_map*=self.Vdot_ratio
        power_map*=self.Vdot_ratio
    
        P1 = self.pin_r
        P2 = self.pout_r
        T1_actual = self.Tsat_s_K + self.DT_sh_K
    
        v_map = 1 / PropsSI('D', 'T', self.Tsat_s_K + 20.0/9.0*5.0, 'P', P1, self.Ref)
        v_actual = 1 / PropsSI('D', 'T', self.Tsat_s_K + self.DT_sh_K, 'P', P1, self.Ref)
        F = 0.75
        mdot = (1 + F * (v_map / v_actual - 1)) * mdot_map
    
        T1_map = self.Tsat_s_K + 20 * 5 / 9
        s1_map = PropsSI('S', 'T', T1_map, 'P', P1, self.Ref)
        h1_map = PropsSI('H', 'T', T1_map, 'P', P1, self.Ref)
        h2s_map = PropsSI('H','P',P2,'S',s1_map,self.Ref)
        #h2s_map = h_sp(self.Ref, s1_map, P2, T1_map + 20) #+20 for guess value
    
        s1_actual = PropsSI('S', 'T', T1_actual, 'P', P1, self.Ref)
        h1_actual = PropsSI('H', 'T', T1_actual, 'P', P1, self.Ref)
        h2s_actual = PropsSI('H','P',P2,'S',s1_actual,self.Ref)
        #h2s_actual = h_sp(self.Ref, s1_actual, P2, T1_actual + 20) #+20 for guess value
    
        #Shaft power based on 20F superheat calculation from fit overall isentropic efficiency
        power = power_map * (mdot / mdot_map) * (h2s_actual - h1_actual) / (h2s_map - h1_map)
    
        h2 = power * (1 - self.fp) / mdot + h1_actual #/1000
        self.eta_oi=mdot*(h2s_actual-h1_actual)/(power) #/1000
        self.Tout_r = PropsSI('T','H',h2,'P',P2,self.Ref)
        #self.Tout_r = T_hp(self.Ref, h2, P2, T1_map + 20) #Plus 20 for guess value for discharge temp
        self.sout_r = PropsSI('S','T',self.Tout_r,'P',P2,self.Ref) #* 1000
        self.sin_r = PropsSI('S','T',self.Tin_r,'P',P1,self.Ref) #* 1000
        self.hout_r = h2 #* 1000
        self.hin_r = h1_actual #* 1000
        self.mdot_r=mdot
        self.W=power
        self.CycleEnergyIn=power*(1-self.fp)
        self.Vdot_pumped=mdot/PropsSI('D','T',self.Tin_r,'P',P1,self.Ref)
        self.Q_amb=-self.fp*power
        
if __name__=='__main__':        
    for i in range(1):
        kwds={
              'M':[217.3163128,5.094492028,-0.593170311,4.38E-02,-2.14E-02,1.04E-02,7.90E-05,-5.73E-05,1.79E-04,-8.08E-05],
              'P':[-561.3615705,-15.62601841,46.92506685,-0.217949552,0.435062616,-0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03],
              'Ref':'R134a',
              'Tin_r':280,
              'pin_r':PropsSI('P','T',279,'Q',1.0,'R134a'),
              'pout_r':PropsSI('P','T',315,'Q',1.0,'R134a'),
              'fp':0.15, #Fraction of electrical power lost as heat to ambient
              'Vdot_ratio': 1.0 #Displacement Scale factor
              }
        Comp=CompressorClass(**kwds)
        Comp.Calculate()
        print Comp.W,'W'
        print Comp.Vdot_pumped,'m^3/s'
        print Comp.Q_amb, 'W'