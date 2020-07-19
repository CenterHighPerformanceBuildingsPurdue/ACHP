from __future__ import division, print_function, absolute_import

from math import pi,exp,log,sqrt,tan,cos,sin
from scipy.constants import g
import numpy as np
from CoolProp.CoolProp import PropsSI
import CoolProp as CP
from ACHP.Correlations import NaturalConv_HTC

class LiquidReceiverClass():
        
    """
    The liquid receiver acts as a damper in the cycle, absorbing the the mass flow rate 
    fluctuations. More concretely, a different explanation can be given.
    When the liquid receiver gets subcooled or saturated liquid at its top, it can be assumed to be
    in thermodynamic equilibrium at each time, because liquid and vapor have the same pressure when
    they enter it (indeed, if the reservoir isn't full, the vapor contained in it must be saturated, as it is in
    presence of liquid). In the inferior part of the tank, the mix of saturated and subcooled liquid (already
    present) allows the working fluid to exit it in a subcooled liquid state. The saturation pressure and
    temperature then reign then in the superior part of the reservoir. Thus, with this component, the
    charge fluctuations are litteraly absorbed, put to an equilibrium value, and the subcooling becomes
    null (this fact can't be stated in the presence of non-condensable gases).
    
    level = (h_v_sat - h)/(h_v_sat - h_l_sat)*(rho/rho_l_sat)
    
    """

    def __init__(self,**kwargs):
        # Load up the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def Update(self,**kwargs):
        # Update the parameters passed in
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
            ('Liquid Receiver Total Volume','m3',self.V_tank),
            ('Liquid Receiver Total Charge','kg',self.Charge),
            ('Inlet Temperature','K',self.Tin_r),
            ('Outlet Temperature','K',self.Tout_r),
            ('Inlet Pressure','Pa',self.pin_r),
            ('Inlet Density', 'kg/m3',self.rhoin_r),
            ('Outlet Pressure','Pa',self.pout_r)
           ]
           
    def Calculate(self):

        # AbstractState
        AS = self.AS
        
        # Inlet conditions
        P1 = self.pin_r
        T1 = self.Tin_r
        AS.update(CP.PT_INPUTS, self.pin_r, self.Tin_r)
        self.rhoin_r = AS.rhomass() #[kg/m^3]
        self.sin_r = AS.smass() #[J/kg-K]
        self.hin_r = AS.hmass() #[J/kg]
        
        # Static pressure (rho*g*h) between inlet and outlet of the tank"
        self.pout_r = self.pin_r #+ (rho_in*g*self.h_ports)/1000

        # Adiabatic
        self.Q = 0.0 #[W]
        
        # No temperature gradient is observed in the reservoir.
        self.Tout_r = self.Tin_r 
        self.sout_r = self.sin_r 

        """
        "Additional notes on liquid receiver:"
        "x_ex_tank=0"	"due to the presence of non condensable gas  (air, due to leakage) in the working fluid, 
        "the liquid at the exit of the tank is not saturated..."

        #h_su_tank=h_ex_cd
        #V_ex_tank = m_dot/rho_ex_tank 
        """
    
        # Refrigerant charge in the accumulator 
        """
        The tank is characterized by an internal diameter and heigth (ID,h)
        and by the maximum level of refrigerant inside
        """
        self.V_tank = pi*self.ID**2/4.0*self.level_r
        self.Charge = self.V_tank*self.rhoin_r

class SuctionAccumulatorClass():
        
    """
    Accumulator at the compressor suction side
    The inlet state is typically superheated
    
    Q > 0 if added to the tank (Tamb-T_w) > 0
    """

    def __init__(self,**kwargs):
        # Load up the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def Update(self,**kwargs):
        # Update the parameters passed in
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
            ('Accumulator Total Volume','m3',self.V_tank),
            ('Accumulator Total Charge','kg',self.Charge),
            ('Inlet Temperature','K',self.Tin_r),
            ('Outlet Temperature','K',self.Tout_r),
            ('Inlet Pressure','Pa',self.pin_r),
            ('Inlet Density', 'kg/m3',self.rhoin_r),
            ('Outlet Pressure','Pa',self.pout_r)
           ]
           
    def Calculate(self):

        # AbstractState
        AS = self.AS
        
        # Inlet conditions
        P1 = self.pin_r
        T1 = self.Tin_r
        AS.update(CP.PT_INPUTS, self.pin_r, self.Tin_r)
        self.rhoin_r = AS.rhomass() #[kg/m^3]
        self.sin_r = AS.smass() #[J/kg-K]
        self.hin_r = AS.hmass() #[J/kg]
        
        # Static pressure (rho*g*h) between inlet and outlet of the tank"
        self.pout_r = self.pin_r #+ (rho_in*g*self.h_ports)/1000
        
        # Non-adiabatic tank
        Ref = 'Air'
        AS_film = CP.AbstractState('HEOS',Ref)
        T_w = self.Tin_r # [K]
        T = self.Tamb # [K]
        P_film = 101325 #[Pa]
        
        # Vertical side shell
        h_c_vertical = NaturalConv_HTC(AS_film,'vertical_plate',T_w,T,P_film,self.h_tank,D_pipe=None,PlateNum=None)
        A_shell_vertical = pi*self.OD*self.h_tank # m^2
        Q_vertical = h_c_vertical*A_shell_vertical*(T-T_w) 

        # Top Surface
        h_c_topshell = NaturalConv_HTC(AS_film,'horizontal_plate',T_w,T,P_film,self.OD,PlateNum='upper_surface')
        A_shell_top = pi*self.OD**2/4 # m2
        Q_topshell = h_c_topshell*A_shell_top*(T-T_w)                     

        # Bottom Surface
        h_c_bottomshell = NaturalConv_HTC(AS_film,'horizontal_plate',T_w,T,P_film,self.OD,PlateNum='lower_surface')
        A_shell_bottom = pi*self.OD**2/4 # m2
        Q_bottomshell = h_c_bottomshell*A_shell_bottom*(T-T_w)  
        
        # Total heat transfer
        self.Q = Q_vertical + Q_topshell + Q_bottomshell # [W]
        
        # Tank energy balance
        self.hout_r = self.hin_r + self.Q
        
        # No temperature gradient is observed in the reservoir.
        AS.update(CP.HmassP_INPUTS,self.hout_r,self.pout_r)
        self.Tout_r = AS.T() #[K]
        self.rhoout_r = AS.rhomass() #[kg/m3]
        self.sout_r = AS.smass() #[J/kg-K]

        # Use average density to calculate charge
        rho_tank = (self.rhoin_r + self.rhoout_r)/2
        
        # Refrigerant charge in the buffer tank 
        self.V_tank = pi*self.ID**2/4.0*self.h_tank
        self.Charge = self.V_tank*rho_tank
        
if __name__=='__main__':

    """
    Example liquid receiver
    """
    #Abstract State        
    Ref = 'R134a'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref)
    kwds={
            'AS': AS,
            'pin_r':527374.817,
            'Tin_r':15.48+273.15,
            'ID':0.3,
            'level_r': 0.5,
            #'h_ports':0.5
            }
    LiquidReceiver=LiquidReceiverClass(**kwds)
    LiquidReceiver.Calculate()
        
    print ('Receiver Volume [cm3]', LiquidReceiver.V_tank*1e6)
    print  ('Charge [kg]',LiquidReceiver.Charge)
    
    """
    Example buffer tank
    """
    #Abstract State        
    Ref = 'R410A'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref)
    kwds={
            'AS': AS,
            'pin_r':1054208.4,
            'Tin_r':19.07+273.15,
            'Tamb':25+273.15,
            'ID':0.115,
            'OD': 0.1236,
            'h_tank': 0.274,
            }
    Accumulator=SuctionAccumulatorClass(**kwds)
    Accumulator.Calculate()
        
    print ('Accumulator Volume [cm3]', Accumulator.V_tank*1e6)
    print  ('Charge [kg]',Accumulator.Charge)
    print ('Heat Loss [W]', Accumulator.Q)
    print ('Tin [K]', Accumulator.Tin_r)
    print ('Tout [K]',Accumulator.Tout_r)