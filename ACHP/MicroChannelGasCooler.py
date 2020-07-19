from __future__ import division, print_function, absolute_import
from math import pi,log,exp
from CoolProp.CoolProp import HAPropsSI
from ACHP.Correlations import f_h_1phase_MicroTube, Petterson_supercritical_average
from ACHP.MicroFinCorrelations import MultiLouveredMicroFins, MicroFinInputs, IsFinsClass
from ACHP.DryWetSegment import DWSVals, DryWetSegment
from ACHP.ACHPTools import ValidateFields

from scipy.optimize import brentq, fsolve
import CoolProp as CP

class FinVals():
    def __init__(self):
        pass
    
class MicroChannelGasCoolerClass():
    def __init__(self,**kwargs):
        #Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)      
    
    def Update(self,**kwargs):
        #Update the parameters passed in
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
            ('Volumetric flow rate of humid air','m^3/s',self.Fins.Air.Vdot_ha),
            ('Inlet Dry bulb temp','K',self.Tin_a),
            ('Inlet Air pressure','Pa',self.Fins.Air.p),
            ('Inlet Air Relative Humidity','-',self.Fins.Air.RH),
            ('Tubes per bank','-',self.Fins.Tubes.NTubes),
            ('Number of banks','-',self.Fins.Tubes.Nbank),
            ('Number of passes','-',self.Fins.Tubes.Npass),
            ('Number of ports','-',self.Fins.Tubes.Nports),
            ('Length of tube','m',self.Fins.Tubes.Ltube),
            ('Tube width','m',self.Td),
            ('Tube height','m',self.Ht),
            ('Tube spacing','m',self.Fins.Tubes.b),
            ('Tube thickness','m',self.Fins.Tubes.tw),
            ('Wall port thickness','m',self.Fins.Tubes.twp),
            ('Channel aspect ratio','-',self.Fins.Tubes.beta),
            ('Tube Conductivity','W/m-K',self.Fins.Tubes.kw),
            ('Fin length','m',self.Fins.Fins.Lf),
            ('Fin thickness','m',self.Fins.Fins.t),
            ('Fin Conductivity','W/m-K',self.Fins.Fins.k_fin),
            ('Louver angle','degree',self.Fins.Louvers.Lalpha),
            ('Louver pitch','m',self.Fins.Louvers.lp),
            ('Louver cut length','m',self.Fins.Llouv),
            ('Fins Type','-',self.FinsType),
            ('Q Total','W',self.Q),
            ('Q Supercritical','W',self.Q_supercritical),
            ('Q Supercritical_liquid','W',self.Q_supercrit_liq),
            ('Inlet Temp','K',self.Tin_r),
            ('Outlet Temp','K',self.Tout_r),
            ('Pressure Drop Total','Pa',self.DP_r),
            ('Pressure Drop Supercritical','Pa',self.DP_r_supercritical),
            ('Pressure Drop Supercritical_liquid','Pa',self.DP_r_supercrit_liq),
            ('Charge Total','kg',self.Charge),
            ('Charge Supercritical','kg',self.Charge_supercritical),
            ('Charge Supercritical_liquid','kg',self.Charge_supercrit_liq),
            ('Mean HTC Superheat','W/m^2-K',self.h_r_supercritical),
            ('Mean HTC Supercritical_liquid','W/m^2-K',self.h_r_supercrit_liq),
            ('Wetted Area Fraction Supercritical','-',self.w_supercritical),
            ('Wetted Area Fraction Supercritical_liquid','-',self.w_supercrit_liq),
            ('Mean Air HTC','W/m^2-K',self.Fins.h_a*self.h_a_tuning),
            ('Surface Effectiveness','-',self.Fins.eta_a),
            ('Air-side area (fin+tubes)','m^2',self.Fins.A_a),
            ('Mass Flow rate of Dry Air','kg/s',self.Fins.mdot_da),
            ('Mass Flow rate of Humid Air','kg/s',self.Fins.mdot_ha),
            ('Pressure Drop Air-side (core only)','Pa',self.Fins.dP_a),
            ('Pressure Drop Air-side (total)','Pa',self.dP_a),
            ('Number of Circuits','-',self.Ncircuits),
            ('Approach temperature degree','K',self.DT_app),
        ]
        

    def Initialize(self):
        #Only validate the first time
        if not hasattr(self,'IsValidated'):
            self.Fins.Validate()
            reqFields=[
               ('Fins',IsFinsClass,None,None),
               ('FinsType',str,None,None),
               ('mdot_r',float,0.00001,20),
               ('Tin_r',float,200,500),
               ('psat_r',float,0.01,200000000)
               ]
            optFields=['Verbosity','AS','h_a_tuning','h_r_tuning','DP_tuning']
            ValidateFields(self.__dict__,reqFields,optFields)
            self.IsValidated=True
        
        #set tuning factors to 1 in case not given by user
        if not hasattr(self,'h_a_tuning'):
            self.h_a_tuning = 1
        if not hasattr(self,'h_r_tuning'):
            self.h_r_tuning = 1
        if not hasattr(self,'DP_tuning'):
            self.DP_tuning = 1
            
        #AbstractState
        assert hasattr(self,'AS'), 'Please specify the Abstract State'
        
        # Retrieve some parameters from nested structures 
        # for code compactness
        self.Ltube=self.Fins.Tubes.Ltube    #tube length
        self.NTubes=self.Fins.Tubes.NTubes  #number of tube (per bank)
        self.Nbank=self.Fins.Tubes.Nbank    #number of banks
        self.Tin_a=self.Fins.Air.Tdb        #inlet air temperature
        self.Pin_a =self.Fins.Air.p         #inlet air pressure
        self.RHin_a=self.Fins.Air.RH        #inlet air relative humidity
        self.Td=self.Fins.Tubes.Td          #tube outside width (depth)
        self.Ht=self.Fins.Tubes.Ht          #tube outside height (major diameter)
        self.b=self.Fins.Tubes.b            #tube spacing
        self.tw=self.Fins.Tubes.tw          #tube thickness
        self.Npass=self.Fins.Tubes.Npass    #Number of passes on ref-side (per bank)
        self.kw=self.Fins.Tubes.kw          #thermal conductivity of tube wall
        self.Nports=self.Fins.Tubes.Nports  #Number of rectangular ports
        self.twp=self.Fins.Tubes.twp        #Port wall thickness
        self.beta=self.Fins.Tubes.beta      #channel (port) aspect ratio (=width/height)
        
        # Define Number of circuits (=number of tubes per pass)
        self.Ncircuits = self.NTubes/self.Npass 
        # Calculate an effective length of circuit if circuits are 
        # not all the same length
        TotalLength=self.Ltube*self.NTubes*self.Nbank
        self.Lcircuit=TotalLength/self.Ncircuits
                                                             
        # Volume of refrigerant = rectangle of tube + circular part at the ends - thickness between ports
        self.V_r = ((self.Td-self.Ht)*(self.Ht-2.0*self.tw) + (pi/4.0) * (self.Ht - 2.0*self.tw)**2 - (self.Ht-2.0*self.tw)*self.twp*(self.Nports-1)) * self.Lcircuit * self.Ncircuits
        # Tube wetted area = tube straight length + circular shape at the ends - horizontal port thickness  + vertical thickness between ports
        self.A_r_wetted = (2.0*(self.Td - self.Ht) + pi*(self.Ht-2.0*self.tw) - 2.0*self.twp*(self.Nports-1) + 2.0*(self.Ht-2.0*self.tw)*(self.Nports-1)) * self.Lcircuit * self.Ncircuits
        # Free-flow area on refrigerant-side = area of rectangle tube + circular parts at end - area thickness between ports
        self.A_c = ((self.Td-self.Ht)*(self.Ht-2.0*self.tw) + (pi/4.0) * (self.Ht - 2.0*self.tw)**2 - self.twp*(self.Ht-2.0*self.tw)*(self.Nports-1)) * self.Ncircuits
        # Hydraulic diameter on ref-side
        self.Dh = 4*self.A_c*self.Lcircuit/self.A_r_wetted
        # Mass flux ref-side
        self.G_r = self.mdot_r / self.A_c
        
        # Total conduction area (exclude port's thickness)
        self.Aw = 2 * (self.Td - self.twp*(self.Nports-1)) * self.Lcircuit * self.Ncircuits
        # Thermal resistance at the wall
        self.Rw = self.tw /(self.kw*self.Aw)
        
        # Define known parameters
        self.AS.update(CP.PT_INPUTS, self.psat_r, self.Tin_r)
        self.hin_r=self.AS.hmass() #[J/kg]
        self.sin_r=self.AS.smass() #[J/kg-K]
        
        # Define critical pressure and temperature
        self.Pcr=self.AS.p_critical() #[Pa]
        self.Tcr=self.AS.T_critical() #[K]
        
        #critical enthalpy at defined pressure
        self.AS.update(CP.PT_INPUTS, self.psat_r, self.Tcr)
        self.hcr = self.AS.hmass() #[J/kg]
        
        #triple temperature
        self.Ttriple = self.AS.Ttriple()
        
        self.Fins.Air.RHmean=self.Fins.Air.RH
        
        #Update with user FinType
        if self.FinsType == 'MultiLouveredMicroFins':
            MultiLouveredMicroFins(self.Fins)
        
        self.mdot_ha=self.Fins.mdot_ha #[kg_ha/s]
        self.mdot_da=self.Fins.mdot_da #[kg_da/s]
        
        
    def Calculate(self):
        
        #Initialize
        self.Initialize()
        AS = self.AS
        
        #assume we have all supercritical region
        self.Tout_r_cr=self.Tcr
        #give an intial guess for the inner wall temperature
        self.T_w = (self.Tout_r_cr+self.Tin_a)/2
        #If we have already used too much of the HX (max possible sum of w is 1.0)
        if self._Supercritical_Forward(1.0)<0:
            self.existsSubcooled=False
            self.w_supercritical=1.0
            def OBJECTIVE(Tout_r_cr):
                self.Tout_r_cr=Tout_r_cr
                AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r_cr)
                hout = AS.hmass()
                Q_target=self.mdot_r*(hout-self.hin_r)
                self._Supercritical_Forward(self.w_supercritical)
                return self.Q_supercritical-Q_target
            brentq(OBJECTIVE,self.Tin_r,self.Tcr)
            #Zero out all the supercritical_liquid parameters
            self.Q_supercrit_liq=0.0
            self.DP_r_supercrit_liq=0.0
            self.Charge_supercrit_liq=0.0
            self.w_supercrit_liq=0.0
            self.h_r_supercrit_liq=0.0
            self.Re_r_supercrit_liq=0.0
            self.Tout_a_supercrit_liq=0.0
            self.fdry_supercrit_liq=0.0    
        else:
            #By definition then we have a supercritical_liquid portion, solve for it
            self.existsSubcooled=True 
            self.w_supercritical=brentq(self._Supercritical_Forward,0.00000000001,0.9999999999)
            def OBJECTIVE(Tout_r_sc):
                self.Tout_r_sc=Tout_r_sc
                AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r_sc)
                hout = AS.hmass()
                Q_target=self.mdot_r*(hout-self.hcr)
                self._Supercrit_liq_Forward(1-self.w_supercritical)
                return self.Q_supercrit_liq-Q_target
            brentq(OBJECTIVE,self.Tcr,self.Ttriple)
            
        #Overall calculations
        self.Q=self.Q_supercritical+self.Q_supercrit_liq
        self.DP_r=(self.DP_r_supercritical+self.DP_r_supercrit_liq)*self.DP_tuning #correcting the pressure drop
        self.Charge=self.Charge_supercritical+self.Charge_supercrit_liq
        
        #Average air outlet temperature (area fraction weighted average) [K]
        self.Tout_a=self.w_supercritical*self.Tout_a_supercritical+self.w_supercrit_liq*self.Tout_a_supercrit_liq
        
        #Outlet enthalpy obtained from energy balance
        self.hout_r=self.hin_r+self.Q/self.mdot_r
        
        AS.update(CP.HmassP_INPUTS, self.hout_r, self.psat_r)
        self.Tout_r = AS.T() #[K]
        self.sout_r = AS.smass() #[J/kg-K]
        
        #Approach temperature
        self.DT_app = self.Tout_r - self.Tin_a
        
        self.hmean_r=self.w_supercritical*self.h_r_supercritical+self.w_supercrit_liq*self.h_r_supercrit_liq
        self.UA_r=self.hmean_r*self.A_r_wetted
        self.UA_a=(self.Fins.h_a*self.h_a_tuning)*self.Fins.A_a*self.Fins.eta_a
        self.UA_w=1/self.Rw
        
        #Upadte air-side pressure drop based on the outlet air temperature
        #the air-side pressure drop here include momentum, contraction and expansion effects
        #Objective function
        def OBJECTIVE_PD(x):
            Pair_o = x[0]
            W = x[1]
            if W < 0: #to ensure that the humidty ratio is 
                print ('Microchannel GasCooler -- Humidity ratio for air pressure drop is less than zero. Humidity ratio is set to 0.0')
                W = 0.0
            v_da=HAPropsSI('V','T',self.Tout_a,'P',Pair_o,'W',W)
            W_new = HAPropsSI('W','T',self.Tout_a,'P',Pair_o,'V',v_da)
            
            #outlet air density
            rho_o = 1 / v_da*(1+W_new) #[m^3/kg_ha]
            #mean air density
            rho_m = pow(0.5*(1/self.Fins.rho_i_air + 1/rho_o),-1)
            #air-side pressure drop including momentum, expansion and contraction effects
            DeltaP_air = self.Fins.G_air**2/2/self.Fins.rho_i_air * ((1 - self.Fins.sigma**2 + self.Fins.Kc_tri) + 2*(self.Fins.rho_i_air/rho_o - 1) + self.Fins.f_a*self.Fins.A_a/self.Fins.A_a_c *(self.Fins.rho_i_air/rho_m) - (1 - self.Fins.sigma**2 -self.Fins.Ke_tri)*(self.Fins.rho_i_air/rho_o))
            
            resids=[(self.Pin_a-Pair_o)-DeltaP_air, W-W_new]     
            return resids
        
        #Initial guesses
        P_init = self.Pin_a
        w_init = HAPropsSI('W','T',self.Tin_a,'P',self.Pin_a,'R',self.RHin_a)
        #solve for outlet air pressure and outlet humidity ratio
        x=fsolve(OBJECTIVE_PD,[P_init,w_init])
        #update the air-side pressure drop
        self.dP_a = self.Pin_a - x[0]
        
    def _Supercritical_Forward(self,w_supercritical):
        # **********************************************************************
        #                      SUPERCRITICAL PART 
        # **********************************************************************
        
        #AbstractState
        AS = self.AS
        
        DWS=DWSVals() #DryWetSegment structure (only dry-analysis, single phase is used)
    
        # Store temporary values to be passed to DryWetSegment
        DWS.Fins=self.Fins
        DWS.FinsType = self.FinsType                                            
        DWS.A_a=self.Fins.A_a*w_supercritical
        DWS.cp_da=self.Fins.cp_da
        DWS.eta_a=self.Fins.eta_a
        DWS.h_a=self.Fins.h_a*self.h_a_tuning  #Heat transfer coefficient, not enthalpy
        DWS.mdot_da=self.mdot_da*w_supercritical
        DWS.pin_a=self.Fins.Air.p

        DWS.Tin_a=self.Tin_a
        DWS.RHin_a=self.Fins.Air.RH
    
        DWS.Tin_r=self.Tin_r
        DWS.A_r=self.A_r_wetted*w_supercritical
        DWS.Rw=self.Rw/w_supercritical
        DWS.pin_r=self.psat_r
        DWS.mdot_r=self.mdot_r
        DWS.IsTwoPhase=False
        
        
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r_cr)
        hout = AS.hmass() #[J/kg]
        #Target heat transfer to go from inlet temperature to iterative outlet temperature
        Q_target=self.mdot_r*(hout-self.hin_r)
        
        if Q_target>0:
            raise ValueError('Q_target in Gas cooler must be negative')
        
        # This block calculates the average refrigerant heat transfer coefficient, average friction factor, average specific heat, and average density
        h_r, f_r_supercritical, DWS.cp_r, rho_supercritical = Petterson_supercritical_average(self.Tout_r_cr, self.Tin_r, self.T_w, self.AS, self.G_r, self.Dh, 0, self.Dh/self.Lcircuit, 0, self.psat_r, -Q_target/DWS.A_r);
        DWS.h_r=h_r*self.h_r_tuning #correct the supercritical HTC
        
        #Compute Fins Efficiency based on FinsType 
        DryWetSegment(DWS)
        
        self.T_w = DWS.Twall_s #inner surface wall temperature (refrigerant)
        self.Q_supercritical=DWS.Q
        self.h_r_supercritical=DWS.h_r
        self.fdry_supercritical=DWS.f_dry
        self.Tout_a_supercritical=DWS.Tout_a
        
        #Pressure drop calculations for supercritical refrigerant
        v_r=1./rho_supercritical
        #Pressure gradient using Darcy friction factor
        dpdz_r=-f_r_supercritical*v_r*self.G_r**2/(2*self.Dh) #Pressure gradient
        self.DP_r_supercritical=dpdz_r*self.Lcircuit*w_supercritical
        #charge for the supercritical portion
        self.Charge_supercritical = w_supercritical * self.V_r * rho_supercritical
        
        if self.Verbosity>7:
            print(w_supercritical,DWS.Q,Q_target,"w_supercritical,DWS.Q,Q_target")
        
        return Q_target-DWS.Q
    
    
    def _Supercrit_liq_Forward(self,w_supercrit_liq):
        # **********************************************************************
        #                      SUPERCRITICAL_LIQUID PART 
        # **********************************************************************
        self.w_supercrit_liq=w_supercrit_liq
        
        if self.w_supercrit_liq<0:
            raise ValueError('w_supercrit_liq in Gas cooler cannot be less than zero')
        
        #AbstractState
        AS = self.AS
        
        DWS=DWSVals() #DryWetSegment structure
    
        # Store temporary values to be passed to DryWetSegment
        DWS.A_a=self.Fins.A_a*w_supercrit_liq
        DWS.cp_da=self.Fins.cp_da
        DWS.eta_a=self.Fins.eta_a
        DWS.h_a=self.Fins.h_a*self.h_a_tuning  #Heat transfer coefficient
        DWS.mdot_da=self.mdot_da*w_supercrit_liq
        DWS.pin_a=self.Fins.Air.p
        DWS.Fins=self.Fins
        DWS.FinsType = self.FinsType           
    
        # Inputs on the air side to two phase region are inlet air again
        DWS.Tin_a=self.Tin_a
        DWS.RHin_a=self.Fins.Air.RH
    
        DWS.Tin_r=self.Tcr
        DWS.A_r=self.A_r_wetted*w_supercrit_liq
        DWS.Rw=self.Rw/w_supercrit_liq
        
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tcr-1)
        DWS.cp_r=AS.cpmass() #[J/kg-K] 
        
        DWS.pin_r=self.psat_r
        DWS.mdot_r=self.mdot_r
        DWS.IsTwoPhase=False
        
    
#         # Friction factor and HTC in the refrigerant portions.
#         # Average fluid temps are used for the calculation of properties 
#         # Average temp of refrigerant is average of sat. temp and outlet temp
#         # Secondary fluid is air over the fins
#         self.f_r_supercrit_liq, self.h_r_supercrit_liq, self.Re_r_supercrit_liq=f_h_1phase_MicroTube(
#           self.G_r, self.Dh, self.Tcr-1.0, self.psat_r, self.AS,
#           "Single")
#         
#         # Average Refrigerant heat transfer coefficient
#         DWS.h_r=self.h_r_supercrit_liq*self.h_r_tuning #correct the supercritical liquid HTC
#         
#         #Run DryWetSegment
#         DryWetSegment(DWS)
#         
#         AS.update(CP.PT_INPUTS, self.psat_r, (self.Tcr + DWS.Tout_r) / 2)
#         rho_supercrit_liq=AS.rhomass() #[kg/m^3]
#         self.Charge_supercrit_liq = self.w_supercrit_liq * self.V_r * rho_supercrit_liq
#     
#         #Pressure drop calculations for supercrit_liqed refrigerant
#         v_r=1/rho_supercrit_liq
#         #Pressure gradient using Darcy friction factor
#         dpdz_r=-self.f_r_supercrit_liq*v_r*self.G_r**2/(2*self.Dh)  #Pressure gradient
#         self.DP_r_supercrit_liq=dpdz_r*self.Lcircuit*self.w_supercrit_liq
#         
#         # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
#         # Heat is removed here from the refrigerant since it is condensing
#         self.Q_supercrit_liq=DWS.Q
#         self.fdry_supercrit_liq=DWS.f_dry
#         self.Tout_a_supercrit_liq=DWS.Tout_a
#         self.Tout_r=DWS.Tout_r

        AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r_sc)
        hout = AS.hmass() #[J/kg]
        #Target heat transfer to go from inlet temperature (critical) to iterative outlet temperature
        Q_target=self.mdot_r*(hout-self.hcr)
        
        if Q_target>0:
            raise ValueError('Q_target in Gas cooler must be negative')
        
        # Friction factor and HTC in the refrigerant portions.
        # Average fluid temps are used for the calculation of properties 
        # Average temp of refrigerant is average of sat. temp and outlet temp
        # Secondary fluid is air over the fins
        h_r, self.f_r_supercrit_liq, DWS.cp_r, rho_supercrit_liq = Petterson_supercritical_average(self.Tout_r_sc, self.Tcr, self.T_w, self.AS, self.G_r, self.Dh, 0, self.Dh/self.Lcircuit, 0, self.psat_r, -Q_target/DWS.A_r);
        #use the same correction factor of supercritical region 
        DWS.h_r=h_r*self.h_r_tuning #correct the supercritical liquid HTC
        
        #Run DryWetSegment
        DryWetSegment(DWS)
        
        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
        # Heat is removed here from the refrigerant since it is condensing
        self.T_w = DWS.Twall_s #inner surface wall temperature (refrigerant)
        self.Q_supercrit_liq=DWS.Q
        self.h_r_supercrit_liq=DWS.h_r
        self.fdry_supercrit_liq=DWS.f_dry
        self.Tout_a_supercrit_liq=DWS.Tout_a
        self.Tout_r=DWS.Tout_r
        
        self.Charge_supercrit_liq = self.w_supercrit_liq * self.V_r * rho_supercrit_liq
    
        #Pressure drop calculations for supercrit_liqed refrigerant
        v_r=1/rho_supercrit_liq
        #Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_supercrit_liq*v_r*self.G_r**2/(2*self.Dh)  #Pressure gradient
        self.DP_r_supercrit_liq=dpdz_r*self.Lcircuit*self.w_supercrit_liq
        
        return Q_target-DWS.Q
    
         
def SampleMicroChannelGasCooler():
    Fins=MicroFinInputs()
    Fins.Tubes.NTubes=61.354           #Number of tubes (per bank for now!)
    Fins.Tubes.Nbank=1                 #Number of banks (set to 1 for now!)
    Fins.Tubes.Npass=3                 #Number of passes (per bank-averaged)
    Fins.Tubes.Nports=1                #Number of rectangular ports
    Fins.Tubes.Ltube=0.30213           #length of a single tube
    Fins.Tubes.Td=0.0333               #Tube outside width (depth)
    Fins.Tubes.Ht= 0.002               #Tube outside height (major diameter)
    Fins.Tubes.b=0.00635               #Tube spacing     
    Fins.Tubes.tw=0.0003               #Tube wall thickness     
    Fins.Tubes.twp=0.0003              #Port (channel) wall thickness     
    Fins.Tubes.beta=1                  #Port (channel) aspect ratio (=width/height)
    Fins.Tubes.kw=117                  #wall thermal conductivity
    
    Fins.Fins.FPI=11.0998              #Fin per inch
    Fins.Fins.Lf=0.0333                #Fin length
    Fins.Fins.t=0.000152               #Fin thickness
    Fins.Fins.k_fin=117                #Fin thermal conductivity
    
    Fins.Air.Vdot_ha=0.281             #rated volumetric flowrate (m^3/s)
    Fins.Air.Tmean=29.4+273.15   
    Fins.Air.Tdb=29.4+273.15           #Dry Bulb Temperature
    Fins.Air.p=101325                  #Air pressure in Pa
    Fins.Air.RH=0.5                    #Relative Humidity
    Fins.Air.RHmean=0.5
    Fins.Air.FanPower=160    
    
    Fins.Louvers.Lalpha=20             #Louver angle, in degree
    Fins.Louvers.lp=0.001              #Louver pitch
    Fins.Louvers.Llouv=0.005737        #Louver cut length
    
    #Abstract State
    Ref = 'R744'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref)
    
    params={
        'AS': AS,
        'mdot_r': 0.076,
        'Tin_r': 110.6+273.15,
        'psat_r': 11000000,
        'Fins': Fins,
        'FinsType': 'MultiLouveredMicroFins',
        'Verbosity':0,
        'h_a_tuning':1,
        'h_r_tuning':1,
        'DP_tuning':1,
    }
    GasCool=MicroChannelGasCoolerClass(**params)
    GasCool.Calculate()
    return GasCool
    
if __name__=='__main__':
    #This runs if you run this file directly
    from time import time
    t1=time()
    GasCool=SampleMicroChannelGasCooler()

    #print(GasCool.OutputList())
    
    print('Heat transfer rate in gas cooler is', GasCool.Q,'W')
    print('Heat transfer rate in gas cooler (supercritical section) is',GasCool.Q_supercritical,'W')
    print('Heat transfer rate in gas cooler (supercritical_liquid section) is',GasCool.Q_supercrit_liq,'W')
    print('Fraction of circuit length in supercritical section is',GasCool.w_supercritical)
    print('Fraction of circuit length in supercritical_liquid section is',GasCool.w_supercrit_liq)
    print('Refrigerant outlet temperature is',GasCool.Tout_r-273.15, 'C')
    print('Air outlet temperature is',GasCool.Tout_a-273.15, 'C')
    print ('Took '+str(time()-t1)+' seconds to run Gas Cooler model')