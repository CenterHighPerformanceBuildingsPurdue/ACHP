from __future__ import division
from CoolProp.CoolProp import PropsSI, FluidsList #,IsFluidType
from Correlations import f_h_1phase_Tube,TrhoPhase_ph
from math import log,pi,exp

class LineSetClass():
    def __init__(self,**kwargs):
        #Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    def Update(self,**kwargs):
        #Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        return [
            ('Length of tube','m',self.L),
            ('Supply line OD','m',self.OD),
            ('Supply line ID','m',self.ID),
            ('Tube Conductivity','W/m-K',self.k_tube),
            ('Insulation thickness','m',self.t_insul),
            ('Insulation conductivity','W/m-K',self.k_insul),
            ('Air overall HTC','W/m^2-K',self.h_air),
            ('Air Temperature','K',self.T_air),
            ('Q Total','W',self.Q),
            ('Pressure drop ','Pa',self.DP),
            ('Reynolds # Fluid','-',self.Re_fluid),
            ('Mean HTC Fluid','W/m^2-K',self.h_fluid),
            ('Charge','kg',self.Charge),
            ('Inlet Temperature','K',self.Tin),
            ('Outlet Temperature','K',self.Tout)
         ]
    
    def Calculate(self):
        if not 'INCOMP' in self.Ref: #IsFluidType(self.Ref,'Brine'):
            #Figure out the inlet state
            self.Tbubble=PropsSI('T','P',self.pin,'Q',0.0,self.Ref)
            self.Tdew=PropsSI('T','P',self.pin,'Q',1.0,self.Ref)
        else:
            #It is a brine
            self.Tbubble = None
            self.Tdew = None
        self.Tin,self.rhoin,self.Phasein=TrhoPhase_ph(self.Ref,self.pin,self.hin,self.Tbubble,self.Tdew)
        self.f_fluid, self.h_fluid, self.Re_fluid=f_h_1phase_Tube(self.mdot, self.ID, self.Tin, self.pin, self.Ref)
        # Specific heat capacity [J/kg-K]                        
        cp=PropsSI('C','T',self.Tin,'P',self.pin,self.Ref) #*1000
        # Density [kg/m^3]
        rho=PropsSI('D','T',self.Tin, 'P', self.pin, self.Ref)
    
        #Thermal resistance of tube
        R_tube=log(self.OD/self.ID)/(2*pi*self.L*self.k_tube)
        #Thermal resistance of insulation
        R_insul=log((self.OD+2.0*self.t_insul)/self.OD)/(2*pi*self.L*self.k_insul);
        #Convective UA for inside the tube 
        UA_i=pi*self.ID*self.L*self.h_fluid;
        #Convective UA for the air-side
        UA_o=pi*(self.OD+2*self.t_insul)*self.L*self.h_air;
        
        #Avoid the possibility of division by zero if h_air is zero
        if UA_o<1e-12:
            UA_o=1e-12
    
        #Overall UA value
        UA=1/(1/UA_i+R_tube+R_insul+1/UA_o)
        
        #Outlet fluid temperature [K]
        self.Tout=self.T_air-exp(-UA/(self.mdot*cp))*(self.T_air-self.Tin)
        #Overall heat transfer rate [W] 
        self.Q=self.mdot*cp*(self.Tout-self.Tin)
        self.hout=self.hin+self.Q/self.mdot
        
        #Pressure drop calculations for superheated refrigerant
        v=1./rho
        G=self.mdot/(pi*self.ID**2/4.0)
        #Pressure gradient using Darcy friction factor
        dpdz=-self.f_fluid*v*G**2/(2.*self.ID) #Pressure gradient
        self.DP=dpdz*self.L
        
        #Charge in Line set [kg]
        self.Charge=pi*self.ID**2/4.0*self.L*rho