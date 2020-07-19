from __future__ import division, print_function, absolute_import
from CoolProp.CoolProp import PropsSI
from ACHP.Correlations import f_h_1phase_Tube,TrhoPhase_ph
from math import log,pi,exp
import CoolProp as CP

class LineSetClass():
    def __init__(self,**kwargs):
        # Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    def Update(self,**kwargs):
        # Load the parameters passed in
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
        # AbstractState
        if hasattr(self,'Backend'): #check if backend is given
            AS = CP.AbstractState(self.Backend, self.Ref)
            if hasattr(self,'MassFrac'):
                AS.set_mass_fractions([self.MassFrac])
        else: #otherwise, use the defualt backend
            AS = CP.AbstractState('HEOS', self.Ref)
        self.AS = AS
        
        if not 'IncompressibleBackend' in AS.backend_name():
            # Figure out the inlet state
            AS.update(CP.PQ_INPUTS, self.pin, 0.0)
            self.Tbubble=AS.T() #[K]
            AS.update(CP.PQ_INPUTS, self.pin, 1.0)
            self.Tdew=AS.T() #[K]
        else:
            # It is a brine
            self.Tbubble = None
            self.Tdew = None
        
        self.Tin,self.rhoin,self.Phasein=TrhoPhase_ph(self.AS,self.pin,self.hin,self.Tbubble,self.Tdew)
        ###Solver shows TwoPhase in the first iteration, the following if statement just to avoid ValueError with CoolProp for pseudo-pure refrigerants
        if self.Phasein =='TwoPhase':
            self.f_fluid, self.h_fluid, self.Re_fluid=f_h_1phase_Tube(self.mdot, self.ID, self.Tin-1, self.pin, self.AS)
            AS.update(CP.PT_INPUTS, self.pin, self.Tin-1)
            # Specific heat capacity [J/kg-K]                        
            cp=AS.cpmass()
            # Density [kg/m^3]
            rho=AS.rhomass()
        else: #Single phase
            self.f_fluid, self.h_fluid, self.Re_fluid=f_h_1phase_Tube(self.mdot, self.ID, self.Tin, self.pin, self.AS)
            AS.update(CP.PT_INPUTS, self.pin, self.Tin)
            # Specific heat capacity [J/kg-K]                        
            cp=AS.cpmass()
            # Density [kg/m^3]
            rho=AS.rhomass()

        # Thermal resistance of tube
        R_tube=log(self.OD/self.ID)/(2*pi*self.L*self.k_tube)
        # Thermal resistance of insulation
        R_insul=log((self.OD+2.0*self.t_insul)/self.OD)/(2*pi*self.L*self.k_insul);
        # Convective UA for inside the tube 
        UA_i=pi*self.ID*self.L*self.h_fluid;
        # Convective UA for the air-side
        UA_o=pi*(self.OD+2*self.t_insul)*self.L*self.h_air;
        
        # Avoid the possibility of division by zero if h_air is zero
        if UA_o<1e-12:
            UA_o=1e-12
    
        # Overall UA value
        UA=1/(1/UA_i+R_tube+R_insul+1/UA_o)
        
        # Outlet fluid temperature [K]
        self.Tout=self.T_air-exp(-UA/(self.mdot*cp))*(self.T_air-self.Tin)
        # Overall heat transfer rate [W] 
        self.Q=self.mdot*cp*(self.Tout-self.Tin)
        self.hout=self.hin+self.Q/self.mdot
        
        # Pressure drop calculations for single phase refrigerant
        v=1./rho
        G=self.mdot/(pi*self.ID**2/4.0)
        # Pressure gradient using Darcy friction factor
        dpdz=-self.f_fluid*v*G**2/(2.*self.ID) #Pressure gradient
        self.DP=dpdz*self.L
        
        # Charge in Line set [kg]
        self.Charge=pi*self.ID**2/4.0*self.L*rho


class LineSetOptionClass():
    def __init__(self,**kwargs):
        # Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    def Update(self,**kwargs):
        # Load the parameters passed in
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
                
        # AbstractState
        AS = self.AS
        if hasattr(self,'MassFrac'):
            AS.set_mass_fractions([self.MassFrac])
        elif hasattr(self, 'VoluFrac'):
            AS.set_volu_fractions([self.VoluFrac])
        

        if not 'IncompressibleBackend' in AS.backend_name():
            # Figure out the inlet state
            AS.update(CP.PQ_INPUTS, self.pin, 0.0)
            self.Tbubble=AS.T() #[K]
            AS.update(CP.PQ_INPUTS, self.pin, 1.0)
            self.Tdew=AS.T() #[K]
            
            # Check for supercritical state
            if self.pin > AS.p_critical():
                # Supercritical
                self.Tbubble= None
                self.Tdew= None
            
        else:
            # It is a brine or incompressible
            self.Tbubble = None
            self.Tdew = None
        
        if hasattr(self,'LineSetOption'): #check if LineSetOption is given
        
            if self.LineSetOption == 'On':
                
                self.Tin,self.rhoin,self.Phasein=TrhoPhase_ph(self.AS,self.pin,self.hin,self.Tbubble,self.Tdew)
                ### Solver shows TwoPhase in the first iteration, the following if statement just to avoid ValueError with CoolProp for pseudo-pure refrigerants
                if self.Phasein =='Supercritical': #TO DO: Need to be UPDATED with Petterson et al. (2000) correlation for transcritical CO2
                    self.f_fluid, self.h_fluid, self.Re_fluid=f_h_1phase_Tube(self.mdot, self.ID, self.Tin, self.pin, self.AS)
                    AS.update(CP.PT_INPUTS, self.pin, self.Tin)
                    # Specific heat capacity [J/kg-K]                        
                    cp=AS.cpmass()
                    # Density [kg/m^3]
                    rho=AS.rhomass()
                    # Specific entropy [J/kg-K]
                    sin=AS.smass()
                elif self.Phasein =='TwoPhase':
                    self.f_fluid, self.h_fluid, self.Re_fluid=f_h_1phase_Tube(self.mdot, self.ID, self.Tin-1, self.pin, self.AS)
                    AS.update(CP.PT_INPUTS, self.pin, self.Tin-1)
                    # Specific heat capacity [J/kg-K]                        
                    cp=AS.cpmass()
                    # Density [kg/m^3]
                    rho=AS.rhomass()
                    # Specific entropy [J/kg-K]
                    sin=AS.smass()
                else: #Single phase
                    self.f_fluid, self.h_fluid, self.Re_fluid=f_h_1phase_Tube(self.mdot, self.ID, self.Tin, self.pin, self.AS)
                    AS.update(CP.PT_INPUTS, self.pin, self.Tin)
                    # Specific heat capacity [J/kg-K]                        
                    cp=AS.cpmass()
                    # Density [kg/m^3]
                    rho=AS.rhomass()
                    # Specific entropy [J/kg-K]
                    sin=AS.smass()
        
                # Inlet specific entropy
                self.sin=sin
                # Thermal resistance of tube
                R_tube=log(self.OD/self.ID)/(2*pi*self.L*self.k_tube)
                # Thermal resistance of insulation
                R_insul=log((self.OD+2.0*self.t_insul)/self.OD)/(2*pi*self.L*self.k_insul);
                # Convective UA for inside the tube 
                UA_i=pi*self.ID*self.L*self.h_fluid;
                # Convective UA for the air-side
                UA_o=pi*(self.OD+2*self.t_insul)*self.L*self.h_air;
                
                # Avoid the possibility of division by zero if h_air is zero
                if UA_o<1e-12:
                    UA_o=1e-12
            
                # Overall UA value
                UA=1/(1/UA_i+R_tube+R_insul+1/UA_o)

                # Outlet fluid temperature [K]
                self.Tout=self.T_air-exp(-UA/(self.mdot*cp))*(self.T_air-self.Tin)
                # Overall heat transfer rate [W] 
                self.Q=self.mdot*cp*(self.Tout-self.Tin)
                self.hout=self.hin+self.Q/self.mdot
                
                # Pressure drop calculations for single phase refrigerant
                v=1./rho
                G=self.mdot/(pi*self.ID**2/4.0)
                # Pressure gradient using Darcy friction factor
                dpdz=-self.f_fluid*v*G**2/(2.*self.ID) #Pressure gradient
                self.DP=dpdz*self.L

                # Outlet specific entropy
                AS.update(CP.HmassP_INPUTS , self.hout, self.pin+self.DP)
                self.sout=AS.smass()
                
                # Charge in Line set [kg]
                self.Charge=pi*self.ID**2/4.0*self.L*rho

            else: # LineSetOption == 'Off'
                print('lineSetOption is off for '+str(self.Name))
                self.Tout=self.Tin
                
                # Inlet specific entropy
                AS.update(CP.HmassP_INPUTS , self.hin, self.pin)
                self.sin=AS.smass()
                
                # Overall heat transfer rate [W] 
                self.Q= 0.0 
                self.hout=self.hin

                self.Re_fluid = 0.0
                self.h_fluid = 0.0
                
                # Pressure gradient using Darcy friction factor
                self.DP= 0.0

                # Outlet specific entropy
                AS.update(CP.HmassP_INPUTS , self.hout, self.pin+self.DP)
                self.sout=AS.smass()
                
                # Charge in Line set [kg]
                self.Charge= 0.0       

        else:
            print(str(self.Name)+ ' is neglected')
            self.Tout=self.Tin
            
            # Inlet specific entropy
            AS.update(CP.HmassP_INPUTS , self.hin, self.pin)
            self.sin=AS.smass()
                
            # Overall heat transfer rate [W] 
            self.Q= 0.0 
            self.hout=self.hin

            self.Re_fluid = 0.0
            self.h_fluid = 0.0
                
            # Pressure gradient using Darcy friction factor
            self.DP= 0.0
 
            # Outlet specific entropy
            AS.update(CP.HmassP_INPUTS , self.hout, self.pin+self.DP)
            self.sout=AS.smass()
            
            # Charge in Line set [kg]
            self.Charge= 0.0
        
if __name__=='__main__':
    #Abstract State        
    Ref = 'R410A'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref)
    
    kwargs={
            'L':7.6,
            'k_tube':0.19,
            't_insul':0.02,
            'k_insul':0.036,
            'T_air':297,
            'AS': AS,
            'h_air':0.0000000001,
            'Tin': PropsSI('T','P',500000,'Q',0,Ref)-10,
            'pin': 500000,           #Pressure of the fluid at the inlet i.e 500kPa
            'hin': PropsSI('H','P',500000,'T',PropsSI('T','P',500000,'Q',0,Ref)-10,Ref),    #Enthalpy of the fluid at the inlet i.e subcooled 10K below bubble
            'mdot': 0.03,#fluid mass flow rate
            'LineSetOption': 'On',
            'Name': 'LineSet',
            }
    

    LineSetSupply = LineSetOptionClass(**kwargs)
    LineSetSupply.OD=0.009525
    LineSetSupply.ID=0.007986
    LineSetSupply.Calculate()
    

    LineSetReturn = LineSetOptionClass(**kwargs)
    LineSetReturn.OD=0.01905
    LineSetReturn.ID=0.017526
    LineSetReturn.Calculate()    

    print (LineSetSupply.OutputList())
    print (LineSetReturn.OutputList())