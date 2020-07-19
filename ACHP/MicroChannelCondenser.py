from __future__ import division, print_function, absolute_import
from math import pi,log,exp
from CoolProp.CoolProp import HAPropsSI
from ACHP.Correlations import f_h_1phase_MicroTube,KM_Cond_Average,TwoPhaseDensity,AccelPressureDrop 
from ACHP.MicroFinCorrelations import MultiLouveredMicroFins, MicroFinInputs, IsFinsClass
from scipy.optimize import brentq, fsolve
from ACHP.ACHPTools import ValidateFields
import CoolProp as CP

class FinVals():
    def __init__(self):
        pass
    
class MicroCondenserClass():
    def __init__(self,**kwargs):
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
            ('Fins per inch','1/in',self.Fins.Fins.FPI),
            ('Fin length','m',self.Fins.Fins.Lf),
            ('Fin thickness','m',self.Fins.Fins.t),
            ('Fin Conductivity','W/m-K',self.Fins.Fins.k_fin),
            ('Louver angle','degree',self.Fins.Louvers.Lalpha),
            ('Louver pitch','m',self.Fins.Louvers.lp),
            ('Louver cut length','m',self.Fins.Llouv),
            ('Fins Type','-',self.FinsType),
            ('Q Total','W',self.Q),
            ('Q Superheat','W',self.Q_superheat),
            ('Q Two-Phase','W',self.Q_2phase),
            ('Q Subcool','W',self.Q_subcool),
            ('Inlet Temp','K',self.Tin_r),
            ('Outlet Temp','K',self.Tout_r),
            ('Pressure Drop Total','Pa',self.DP_r),
            ('Pressure Drop Superheat','Pa',self.DP_r_superheat),
            ('Pressure Drop Two-Phase','Pa',self.DP_r_2phase),
            ('Pressure Drop Subcool','Pa',self.DP_r_subcool),
            ('Charge Total','kg',self.Charge),
            ('Charge Superheat','kg',self.Charge_superheat),
            ('Charge Two-Phase','kg',self.Charge_2phase),
            ('Charge Subcool','kg',self.Charge_subcool),
            ('Mean HTC Superheat','W/m^2-K',self.h_r_superheat),
            ('Mean HTC Two-phase','W/m^2-K',self.h_r_2phase),
            ('Mean HTC Subcool','W/m^2-K',self.h_r_subcool),
            ('Wetted Area Fraction Superheat','-',self.w_superheat),
            ('Wetted Area Fraction Two-phase','-',self.w_2phase),
            ('Wetted Area Fraction Subcool','-',self.w_subcool),
            ('Mean Air HTC','W/m^2-K',self.Fins.h_a*self.h_a_tuning),
            ('Surface Effectiveness','-',self.Fins.eta_a),
            ('Air-side area (fin+tubes)','m^2',self.Fins.A_a),
            ('Mass Flow rate of Dry Air','kg/s',self.Fins.mdot_da),
            ('Mass Flow rate of Humid Air','kg/s',self.Fins.mdot_ha),
            ('Pressure Drop Air-side (core only)','Pa',self.Fins.dP_a),
            ('Pressure Drop Air-side (total)','Pa',self.dP_a),
            ('Subcooling','K',self.DT_sc),
            ('Number of Circuits','-',self.Ncircuits)
        ]
        
    def Update(self,**kwargs):
        #Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def Calculate(self):
        #Only validate the first time
        if not hasattr(self,'IsValidated'):
            self.Fins.Validate()
            reqFields=[
               ('Fins',IsFinsClass,None,None),
               ('FinsType',str,None,None),
               ('mdot_r',float,0.00001,20),
               ('Tin_r',float,200,500),
               ('psat_r',float,0.01,20000000)
               ]
            optFields=['Verbosity','AS','h_a_tuning','h_tp_tuning','DP_tuning']
            ValidateFields(self.__dict__,reqFields,optFields)
            self.IsValidated=True
        
        #set tuning factors to 1 in case not given by user
        if not hasattr(self,'h_a_tuning'):
            self.h_a_tuning = 1
        if not hasattr(self,'h_tp_tuning'):
            self.h_tp_tuning = 1
        if not hasattr(self,'DP_tuning'):
            self.DP_tuning = 1
            
        #AbstractState
        AS = self.AS
        
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

        ## Bubble and dew temperatures (same for fluids without glide)
        AS.update(CP.PQ_INPUTS, self.psat_r, 0.0)
        self.Tbubble=AS.T() #[K]
        self.h_l = AS.hmass() #[J/kg]
        self.cp_satL = AS.cpmass() #[J/kg-K]
        
        AS.update(CP.PQ_INPUTS, self.psat_r, 1.0)
        self.Tdew=AS.T() #[K]
        self.h_v = AS.hmass() #[J/kg]
        
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
        
        #Definitely have a superheated portion
        self._Superheat_Forward()
        #Maybe have a full two-phase section
        #First try to run with a full two-phase section from quality of 1 to quality of 0
        self._TwoPhase_Forward()
        #If we have already used too much of the HX (max possible sum of w is 1.0)
        if self.w_2phase+self.w_superheat>1:
            #There is no subcooled portion, solve for outlet quality
            brentq(self._TwoPhase_Forward,0.0000001,0.9999999)
            #Zero out all the subcooled parameters
            self.Q_subcool=0.0
            self.DP_r_subcool=0.0
            self.Charge_subcool=0.0
            self.w_subcool=0.0
            self.h_r_subcool=0.0
            self.existsSubcooled=False
        else:
            #By definition then we have a subcooled portion, solve for it
            self.existsSubcooled=True 
            self._Subcool_Forward()
        
        #Overall calculations
        self.Q=self.Q_superheat+self.Q_2phase+self.Q_subcool
        self.DP_r=self.DP_r_superheat+self.DP_r_2phase+self.DP_r_subcool
        self.DP_r=self.DP_r*self.DP_tuning #correcting the pressure drop
        self.Charge=self.Charge_2phase+self.Charge_subcool+self.Charge_superheat
        
        #define known parameters
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tin_r)
        self.hin_r=AS.hmass() #[J/kg]
        self.sin_r=AS.smass() #[J/kg-K]
        
        if self.existsSubcooled==True:
            AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r)
            self.hout_r=AS.hmass() #[J/kg]
            self.sout_r=AS.smass() #[J/kg-K]
        else:
            self.Tout_r=self.xout_2phase*self.Tdew+(1-self.xout_2phase)*self.Tbubble
            AS.update(CP.QT_INPUTS, 0.0, self.Tbubble)
            h_l = AS.hmass() #[J/kg]
            s_l = AS.smass() #[J/kg-K]
            AS.update(CP.QT_INPUTS, 1.0, self.Tdew)
            h_v = AS.hmass() #[J/kg]
            s_v = AS.smass() #[J/kg-K]
            self.hout_r=h_l+self.xout_2phase*(h_v-h_l)
            self.sout_r=s_l+self.xout_2phase*(s_v-s_l)
            #Use the effective subcooling
            self.DT_sc=self.DT_sc_2phase
        
        #Calculate the mean outlet air temperature [K]
        self.Tout_a=self.Tin_a-self.Q/(self.Fins.cp_da*self.Fins.mdot_da)
        self.hmean_r=self.w_2phase*self.h_r_2phase+self.w_superheat*self.h_r_superheat+self.w_subcool*self.h_r_subcool
        self.UA_r=self.hmean_r*self.A_r_wetted
        self.UA_a=(self.Fins.h_a*self.h_a_tuning)*self.Fins.A_a*self.Fins.eta_a
        self.UA_w=1/self.Rw
        
        #Upadte air-side pressure drop based on the outlet air temperature
        #the air-side pressure drop here include momentum, contraction and expansion effects
        #Objective function
        def OBJECTIVE(x):
            Pair_o = x[0]
            W = x[1]
            if W < 0: #to ensure that the humidty ratio is 
                print ('Microchannel Condensder -- Humidity ratio for air pressure drop is less than zero. Humidity ratio is set to 0.0')
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
        x=fsolve(OBJECTIVE,[P_init,w_init])
        #update the air-side pressure drop
        self.dP_a = self.Pin_a - x[0]

        
    def _Superheat_Forward(self):  
        # **********************************************************************
        #                      SUPERHEATED PART 
        # **********************************************************************
        #AbstractState
        AS=self.AS
        #Dew temperature for constant pressure cooling to saturation
        Tdew=self.Tdew
        Tbubble=self.Tbubble
        
        # Average fluid temps are used for the calculation of properties 
        # Average temp of refrigerant is average of sat. temp and outlet temp		
        # Secondary fluid is air over the fins
        self.f_r_superheat, self.h_r_superheat, self.Re_r_superheat=f_h_1phase_MicroTube(self.G_r, self.Dh, (Tdew+self.Tin_r)/2.0, self.psat_r, self.AS, "Single")
        
        AS.update(CP.PT_INPUTS, self.psat_r, (Tdew+self.Tin_r)/2)    
        cp_r = AS.cpmass() #[J/kg-K]

        #Compute Fins Efficiency based on FinsType
        if self.FinsType == 'MultiLouveredMicroFins': 
            MultiLouveredMicroFins(self.Fins)
            
        self.mdot_da=self.Fins.mdot_da
        
        # Cross-flow in the superheated region.  
        # Using effectiveness-Ntu relationships for cross flow with non-zero Cr.
        UA_overall = 1. / (1. / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning) + 1. / (self.h_r_superheat * self.A_r_wetted) + self.Rw)
        epsilon_superheat=(Tdew-self.Tin_r)/(self.Tin_a-self.Tin_r)
        Ntu=UA_overall/(self.mdot_da*self.Fins.cp_da)
        if epsilon_superheat>1.0:
            epsilon_superheat=1.0-1e-12
        self.w_superheat=-log(1-epsilon_superheat)*self.mdot_r*cp_r/((1-exp(-Ntu))*self.mdot_da*self.Fins.cp_da)
              
        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
        # Heat is removed here from the refrigerant since it is being cooled
        self.Q_superheat = self.mdot_r * cp_r * (Tdew-self.Tin_r)
        
        AS.update(CP.PT_INPUTS, self.psat_r, (self.Tin_r+Tdew)/2.0)
        rho_superheat=AS.rhomass() #[kg/m^3]
        #Pressure drop calculations for superheated refrigerant
        v_r=1./rho_superheat;
        #Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_superheat*v_r*self.G_r**2/(2.*self.Dh) #Pressure gradient
        self.DP_r_superheat=dpdz_r*self.Lcircuit*self.w_superheat
        self.Charge_superheat = self.w_superheat * self.V_r * rho_superheat

        #Latent heat needed for pseudo-quality calc
        AS.update(CP.QT_INPUTS, 0.0, Tbubble)
        h_l = AS.hmass() #[J/kg]
        AS.update(CP.QT_INPUTS, 1.0, Tdew)
        h_v = AS.hmass() #[J/kg]
        h_fg = h_v - h_l #[J/kg]
        self.xin_r=1.0+cp_r*(self.Tin_r-Tdew)/h_fg
        
    def _TwoPhase_Forward(self,xout_r_2phase=0.0):
        # **********************************************************************
        #                      TWO-PHASE PART 
        # **********************************************************************
        """
            xout_r_2phase: quality of refrigerant at end of two-phase portion
                default value is 0.0 (full two phase region)
        """
        ## Bubble and dew temperatures (same for fluids without glide)
        Tbubble=self.Tbubble
        Tdew=self.Tdew
        ## Mean temperature for use in HT relationships
        Tsat_r=(Tbubble+Tdew)/2
        
        h_l = self.h_l #[J/kg]
        h_v = self.h_v #[J/kg]
        h_fg = h_v - h_l #[J/kg]
        
        # This block calculates the average frictional pressure drop griendient
        # and average refrigerant heat transfer coefficient by
        # integrating the local heat transfer coefficient between 
        # a quality of 1.0 and the outlet quality
        DPDZ_frict_2phase, h_r_2phase =KM_Cond_Average(xout_r_2phase,1.0,self.AS,self.G_r,self.Dh,Tbubble,Tdew,self.psat_r,self.beta)
        
        self.h_r_2phase=h_r_2phase*self.h_tp_tuning

        UA_overall = 1 / (1 / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning) + 1 / (self.h_r_2phase * self.A_r_wetted) + self.Rw);
        self.epsilon_2phase=1-exp(-UA_overall/(self.mdot_da*self.Fins.cp_da));
        self.w_2phase=-self.mdot_r*h_fg*(1.0-xout_r_2phase)/(self.mdot_da*self.Fins.cp_da*(self.Tin_a-Tsat_r)*self.epsilon_2phase);

        #Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
        #Heat is removed here from the refrigerant since it is condensing
        self.Q_2phase = self.epsilon_2phase * self.Fins.cp_da* self.mdot_da * self.w_2phase * (self.Tin_a-Tsat_r);
        
        self.xout_2phase=xout_r_2phase
        
        # Frictional pressure drop component
        DP_frict=DPDZ_frict_2phase*self.Lcircuit*self.w_2phase
        #Accelerational pressure drop component    
        DP_accel=-AccelPressureDrop(self.xout_2phase,1.0,self.AS,self.G_r,Tbubble,Tdew,slipModel='Zivi')*self.Lcircuit*self.w_2phase
        # Total pressure drop is the sum of accelerational and frictional components (neglecting gravitational effects)
        self.DP_r_2phase=DP_frict+DP_accel
    
        rho_average=TwoPhaseDensity(self.AS,self.xout_2phase,1.0,self.Tdew,self.Tbubble,slipModel='Zivi')
        self.Charge_2phase = rho_average * self.w_2phase * self.V_r    
        
        if self.Verbosity>7:
            print ('2phase cond resid', self.w_2phase-(1-self.w_superheat))
            print ('h_r_2phase',self.h_r_2phase)
        
        #Calculate an effective pseudo-subcooling based on the equality
        cp_satL=self.cp_satL
        self.DT_sc_2phase=-self.xout_2phase*h_fg/(cp_satL)
            
        #If the quality is being solved for, the length of the two-phase and subcooled
        # sections should add to the length of the HX.  Return the residual
        return self.w_2phase-(1-self.w_superheat)
    
    
    def _Subcool_Forward(self):
        # **********************************************************************
        #                      SUBCOOLED PART 
        # **********************************************************************
        self.w_subcool=1-self.w_2phase-self.w_superheat
        
        if self.w_subcool<0:
            raise ValueError('w_subcool in Condenser cannot be less than zero')
        #AbstractState
        AS = self.AS
        # Bubble temperature
        Tbubble=self.Tbubble
        
        # Based on the the construction of the cycle model there is guaranteed to be a 
        # two-phase portion of the heat exchanger
        A_a_subcool = self.Fins.A_a * self.w_subcool
        mdot_da_subcool = self.mdot_da * self.w_subcool
        A_r_subcool =  self.A_r_wetted * self.w_subcool
    
        # Friction factor and HTC in the refrigerant portions.
        # Average fluid temps are used for the calculation of properties 
        # Average temp of refrigerant is average of sat. temp and outlet temp
        # Secondary fluid is air over the fins
        self.f_r_subcool, self.h_r_subcool, self.Re_r_subcool=f_h_1phase_MicroTube(
          self.G_r, self.Dh, Tbubble-1.0, self.psat_r, self.AS,
          "Single")
        
        AS.update(CP.PT_INPUTS, self.psat_r, Tbubble-1)
        cp_r = AS.cpmass() #[J/kg-K]
    
        # Cross-flow in the subcooled region.
        R_a=1. / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning)
        R_r=1. / (self.h_r_subcool * self.A_r_wetted)
        UA_subcool = self.w_subcool / (R_a + R_r + self.Rw)
        Cmin=min([self.mdot_da*self.Fins.cp_da*self.w_subcool,self.mdot_r*cp_r])
        Cmax=max([self.mdot_da*self.Fins.cp_da*self.w_subcool,self.mdot_r*cp_r])
        Cr=Cmin/Cmax
        NTU=UA_subcool/Cmin
        
        if(self.mdot_da*self.Fins.cp_da*self.w_subcool>self.mdot_r*cp_r):
            #Minimum capacitance rate on refrigerant side
            epsilon_subcool = 1. - exp(-1. / Cr * (1. - exp(-Cr * NTU)))
        else:
            #Minimum capacitance rate on air side
            epsilon_subcool = 1 / Cr * (1 - exp(-Cr * (1 - exp(-NTU))))
        
        #Effectiveness for both fluids unmixed:
        #epsilon_subcool = 1 - exp(1/Cr * NTU**0.22 * (exp(-Cr*NTU**0.78)-1))
        
           
        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
        # Heat is removed here from the refrigerant since it is condensing
        self.Q_subcool=-epsilon_subcool*Cmin*(Tbubble-self.Tin_a)
        self.DT_sc=-self.Q_subcool/(self.mdot_r*cp_r)
        self.Tout_r=Tbubble-self.DT_sc
        
        AS.update(CP.PT_INPUTS, self.psat_r, (Tbubble + self.Tout_r) / 2)
        rho_subcool=AS.rhomass() #[kg/m^3]
        self.Charge_subcool = self.w_subcool * self.V_r * rho_subcool
    
        #Pressure drop calculations for subcooled refrigerant
        v_r=1/rho_subcool
        #Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_subcool*v_r*self.G_r**2/(2*self.Dh)  #Pressure gradient
        self.DP_r_subcool=dpdz_r*self.Lcircuit*self.w_subcool
        
def SampleMicroCondenser(AS,T=95):
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
    
    Fins.Air.Vdot_ha=1.05              #Air volume flow rate in m^3/s
    Fins.Air.Tmean=298 
    Fins.Air.Tdb=298                   #Air inlet temperature, K
    Fins.Air.p=100000                  #Air pressure in Pa
    Fins.Air.RHmean=0.5
    Fins.Air.RH=0.5                    #Air inlet relative humidity
    Fins.Air.FanPower=854.9           #Fan power, Watts
    
    Fins.Louvers.Lalpha=20             #Louver angle, in degree
    Fins.Louvers.lp=0.001              #Louver pitch
    Fins.Louvers.Llouv=0.005737        #Louver cut length
    
    params={
        'AS': AS,
        'mdot_r': 0.0683,
        'Tin_r': T+273.15,
        'psat_r': 3500000, 
        'Fins': Fins,
        'FinsType': 'MultiLouveredMicroFins',
        'Verbosity':0,
        'h_a_tuning':1,
        'h_tp_tuning':1,
        'DP_tuning':1
    }
    MicroCond=MicroCondenserClass(**params)
    MicroCond.Calculate()
    return MicroCond
    
if __name__=='__main__':
    #This runs if you run this file directly
    Ref = 'R410A'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref) #Abstract State        
    MicroCond=SampleMicroCondenser(AS, 95)
    #print (MicroCond.OutputList())
    print ('Heat transfer rate in condenser is', MicroCond.Q,'W')
    print ('Heat transfer rate in condenser (superheat section) is',MicroCond.Q_superheat,'W')
    print ('Heat transfer rate in condenser (twophase section) is',MicroCond.Q_2phase,'W')
    print ('Heat transfer rate in condenser (subcooled section) is',MicroCond.Q_subcool,'W')
    print ('Fraction of circuit length in superheated section is',MicroCond.w_superheat)
    print ('Fraction of circuit length in twophase section is',MicroCond.w_2phase)
    print ('Fraction of circuit length in subcooled section is',MicroCond.w_subcool) 