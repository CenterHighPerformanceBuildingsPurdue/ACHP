from __future__ import division, absolute_import, print_function
from CoolProp.CoolProp import PropsSI
import CoolProp as CP

from math import pi,exp,log,sqrt,tan,cos,sin

class ExpDevClass():
    """
     Expansion device model
     based on the work developed by Bo Shen in ACMODEL
    """
    def __init__(self,**kwargs):
        #Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def Update(self,**kwargs):
        #Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    
#     def OutputList(self): #TODO: fix this list of outputs
#         """
#             Return a list of parameters for this component for further output
#             
#             It is a list of tuples, and each tuple is formed of items:
#                 [0] Description of value
#                 [1] Units of value
#                 [2] The value itself
#         """
#         
#         return [
#             ('Effective Length','m',self.Lp),
#             ('Wetted area','m^2',self.A_h_wetted),
#             ('Outlet Superheat','K',self.Tin_c-self.Tdew_c),
#             ('Q Total','W',self.Q),
#             ('Q Superheat Hot','W',self.Q_superheated_h),
#             ('Q Two-Phase Hot','W',self.Q_2phase_h),
#             ('Q Subcooled Hot','W',self.Q_subcooled_h),
#             ('Q Superheat Cold','W',self.Q_superheated_c),
#             ('Q Two-Phase Cold','W',self.Q_2phase_c),
#             ('Q Subcooled Cold','W',self.Q_subcooled_c),
#             ('Inlet hot stream temp','K',self.Tin_h),
#             ('Outlet hot stream temp','K',self.Tout_h),
#             ('Inlet cold stream temp','K',self.Tin_c),
#             ('Outlet cold stream temp','K',self.Tout_c),
#             ('Charge Total','kg',self.Charge_h),
#             ('Charge Superheat','kg',self.Charge_superheated_h),
#             ('Charge Two-Phase','kg',self.Charge_2phase_h),
#             ('Charge Subcool','kg',self.Charge_subcooled_h),
#             ('Charge Total','kg',self.Charge_c),
#             ('Charge Superheat','kg',self.Charge_superheated_c),
#             ('Charge Two-Phase','kg',self.Charge_2phase_c),
#             ('Charge Subcool','kg',self.Charge_subcooled_c),
#             ('Hot HTC Superheat','W/m^2-K',self.h_superheated_h),
#             ('Hot HTC Two-Phase','W/m^2-K',self.h_2phase_h),
#             ('Hot HTC Subcool','W/m^2-K',self.h_subcooled_h),
#             ('Cold Mean HTC Superheat','W/m^2-K',self.h_superheated_c),
#             ('Cold Mean HTC Ref. Two-Phase','W/m^2-K',self.h_2phase_c),
#             ('Cold Mean HTC Ref. Subcool','W/m^2-K',self.h_subcooled_c),
#             ('Pressure Drop Hot','Pa',self.DP_h),
#             ('Pressure Drop Hot superheated','Pa',self.DP_superheated_h),
#             ('Pressure Drop Hot 2 phase','Pa',self.DP_2phase_h),
#             ('Pressure Drop Hot subcooled','Pa',self.DP_subcooled_h),
#             ('Pressure Drop Cold','Pa',self.DP_c),
#             ('Pressure Drop Cold superheated','Pa',self.DP_superheated_c),
#             ('Pressure Drop Cold 2 phase','Pa',self.DP_2phase_c),
#             ('Pressure Drop Cold subcooled','Pa',self.DP_subcooled_c),
#             ('Area fraction Superheat Hot','-',self.w_superheated_h),
#             ('Area fraction Two-Phase Hot','-',self.w_2phase_h),
#             ('Area fraction Subcooled Hot','-',self.w_subcooled_h),
#             ('Area fraction Superheat Cold','-',self.w_superheated_c),
#             ('Area fraction Two-Phase Cold','-',self.w_2phase_c),
#             ('Area fraction Subcooled Cold','-',self.w_subcooled_c)
#          ]
        
    def Initialize(self):
        
        #AbstractState
        assert hasattr(self,'AS'), 'Please specify the Abstract State'
        
        #If the user doesn't include the ExpType, fail
        assert hasattr(self,'ExpType'), 'Please specify the type of the expansion device'
        
    def Calculate(self):
        
        #Initialize
        self.Initialize()
        #AbstractState
        AS = self.AS
        
        if self.ExpType == 'Ideal':
            #===================================================================
            # No information about expansion device is given
            #===================================================================
            #inlet state
            if self.pin_r > AS.p_critical(): #Supercritical
                AS.update(CP.HmassP_INPUTS, self.hin_r, self.pin_r)
                self.sin_r = AS.smass() #[J/kg-K]
                self.Tin_r = AS.T() #[K]
            else: #other refrigerants  
                AS.update(CP.PQ_INPUTS, self.pin_r, 0.0)
                Tbubble_in = AS.T() #[K]
                h_l_in = AS.hmass() #[J/kg]
                s_l_in = AS.smass() #[J/kg-K]
                AS.update(CP.PQ_INPUTS, self.pin_r, 1.0)
                Tdew_in = AS.T() #[K]
                h_v_in = AS.hmass() #[J/kg]
                s_v_in = AS.smass() #[J/kg-K]
                
                self.xin_r = (self.hin_r-h_l_in)/(h_v_in-h_l_in)
                if (self.xin_r>0.999):
                    print ("ExpDev :: Upstream state in the expansion device is superheated")
                    raise
                if (self.xin_r>0.0 and self.xin_r<1.0): #two-phase state at the inlet
                    self.sin_r = self.xin_r*s_v_in+(1-self.xin_r)*s_l_in #[J/kg-K]
                    self.Tin_r = self.xin_r*Tdew_in+(1-self.xin_r)*Tbubble_in #[K]
                else: #liquid state at the inlet
                    AS.update(CP.HmassP_INPUTS, self.hin_r, self.pin_r)
                    self.sin_r = AS.smass() #[J/kg-K]
                    self.Tin_r = AS.T() #[K]
                
            #outlet state (assume h = constant)
            self.hout_r = self.hin_r #[J/kg]
            
            AS.update(CP.PQ_INPUTS, self.pout_r, 0.0)
            Tbubble_out = AS.T() #[K]
            h_l_out = AS.hmass() #[J/kg]
            s_l_out = AS.smass() #[J/kg-K]
            AS.update(CP.PQ_INPUTS, self.pout_r, 1.0)
            Tdew_out = AS.T() #[K]
            h_v_out = AS.hmass() #[J/kg]
            s_v_out = AS.smass() #[J/kg-K]
            
            #outlet state (two-phase)
            self.xout_r = (self.hout_r-h_l_out)/(h_v_out-h_l_out) #[-]
            self.Tout_r = self.xout_r*Tdew_out+(1-self.xout_r)*Tbubble_out #[K]
            self.sout_r = self.xout_r*s_v_out+(1-self.xout_r)*s_l_out #[J/kg-K]
        
        else:
            raise

                  
if __name__=='__main__':
    #Abstract State
    Ref = 'R410A'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref)
    
    print('Example for Ideal expansion device')
    params={
        'AS':AS,
        'ExpType':'Ideal',     #expansion device type
        'pin_r': PropsSI('P','T',60+273.15,'Q',0,Ref), #upsteam pressure
        'hin_r': PropsSI('H','P',PropsSI('P','T',60+273.15,'Q',0,Ref),'Q',0,Ref), #upstream enthalpy
        'pout_r': PropsSI('P','T',10+273.15,'Q',0,Ref), #downstream pressure        
    }
        
    Exp=ExpDevClass(**params)
    Exp.Calculate()
    print('Tout =',Exp.Tout_r,'[K]')
    print('hout =',Exp.hout_r,'[J/kg]')
    print('xout =',Exp.xout_r,'[-]')
    print()
    
