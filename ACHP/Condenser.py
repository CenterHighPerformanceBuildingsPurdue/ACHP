from __future__ import division, print_function, absolute_import
from math import pi,log,exp
from ACHP.Correlations import f_h_1phase_Tube,ShahCondensation_Average,LMPressureGradientAvg,TwoPhaseDensity,AccelPressureDrop 
from ACHP.FinCorrelations import WavyLouveredFins,FinInputs,IsFinsClass, HerringboneFins, PlainFins
from scipy.optimize import brentq
from ACHP.ACHPTools import ValidateFields
import CoolProp as CP

class FinVals():
    def __init__(self):
        pass
    
class CondenserClass():
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
            ('Tubes per bank','-',self.Fins.Tubes.NTubes_per_bank),
            ('Number of banks','-',self.Fins.Tubes.Nbank),
            ('Number circuits','-',self.Fins.Tubes.Ncircuits),
            ('Length of tube','m',self.Fins.Tubes.Ltube),
            ('Tube OD','m',self.OD),
            ('Tube ID','m',self.ID),
            ('Tube Long. Pitch','m',self.Fins.Tubes.Pl),
            ('Tube Transverse Pitch','m',self.Fins.Tubes.Pt),
            ('Tube Conductivity','W/m-K',self.Fins.Tubes.kw),
            ('Fins per inch','1/in',self.Fins.Fins.FPI),
            ('Fin waviness pd','m',self.Fins.Fins.Pd),
            ('Fin waviness xf','m',self.Fins.Fins.xf),
            ('Fin thickness','m',self.Fins.Fins.t),
            ('Fin Conductivity','W/m-K',self.Fins.Fins.k_fin),
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
            ('Mean Air HTC','W/m^2-K',self.Fins.h_a),
            ('Surface Effectiveness','-',self.Fins.eta_a),
            ('Air-side area (fin+tubes)','m^2',self.Fins.A_a*self.h_a_tuning),
            ('Mass Flow rate of Dry Air','kg/s',self.Fins.mdot_da),
            ('Mass Flow rate of Humid Air','kg/s',self.Fins.mdot_ha),
            ('Pressure Drop Air-side','Pa',self.Fins.dP_a),
            ('Subcooling','K',self.DT_sc)
        ]
        
    def Update(self,**kwargs):
        # Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def Calculate(self):
        # Only validate the first time
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
        
        # Set tuning factors to 1 in case not given by user
        if not hasattr(self,'h_a_tuning'):
            self.h_a_tuning = 1
        if not hasattr(self,'h_tp_tuning'):
            self.h_tp_tuning = 1
        if not hasattr(self,'DP_tuning'):
            self.DP_tuning = 1
                  
        # AbstractState
        AS = self.AS
        
        # Retrieve some parameters from nested structures 
        # for code compactness
        self.ID=self.Fins.Tubes.ID
        self.OD=self.Fins.Tubes.OD
        self.Ltube=self.Fins.Tubes.Ltube
        self.NTubes_per_bank=self.Fins.Tubes.NTubes_per_bank
        self.Nbank=self.Fins.Tubes.Nbank
        self.Ncircuits=self.Fins.Tubes.Ncircuits
        self.Tin_a=self.Fins.Air.Tdb
        self.kw=self.Fins.Tubes.kw          #thermal conductivity of tube wall
        
        ## Bubble and dew temperatures (same for fluids without glide)
        AS.update(CP.PQ_INPUTS, self.psat_r, 0.0)
        self.Tbubble=AS.T() #[K]
        self.h_l = AS.hmass() #[J/kg]
        self.cp_satL = AS.cpmass() #[J/kg-K]
        
        AS.update(CP.PQ_INPUTS, self.psat_r, 1.0)
        self.Tdew=AS.T() #[K]
        self.h_v = AS.hmass() #[J/kg]
        
        # Calculate an effective length of circuit if circuits are 
        # not all the same length
        TotalLength=self.Ltube*self.NTubes_per_bank*self.Nbank
        self.Lcircuit=TotalLength/self.Ncircuits
        
        self.V_r = pi * self.ID**2 / 4.0 * self.Lcircuit * self.Ncircuits
        self.A_r_wetted = pi * self.ID * self.Ncircuits * self.Lcircuit
        self.G_r = self.mdot_r/(self.Ncircuits*pi*self.ID**2/4.0)    
        
        # Thermal resistance at the wall
        self.Rw = log(self.OD/self.ID)/(2*pi*self.kw*self.Lcircuit*self.Ncircuits)
        
        # Define known parameters
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tin_r)
        self.hin_r=AS.hmass() #[J/kg]
        self.sin_r=AS.smass() #[J/kg-K]
        
        # Definitely have a superheated portion
        # First try to run with a full superheated section from Tin_r to Tdew
        self._Superheat_Forward(T_sh_out_r = self.Tdew)
        # Maybe have a full two-phase section
        # First try to run with a full two-phase section from quality of 1 to quality of 0
        self._TwoPhase_Forward()
        # If we have already used too much of the HX (max possible sum of w is 1.0)
        if self.w_superheat>=1:
            # There is partial superheated region
            brentq(self._Superheat_Forward,self.Tdew,self.Tin_r)
            # Zero out all the two-phase parameters
            self.Q_2phase=0.0
            self.DP_r_2phase=0.0
            self.Charge_2phase=0.0
            self.w_2phase=0.0
            self.h_r_2phase=0.0
            # Zero out all the subcooled parameters
            self.Q_subcool=0.0
            self.DP_r_subcool=0.0
            self.Charge_subcool=0.0
            self.w_subcool=0.0
            self.h_r_subcool=0.0
            self.existsSubcooled=False                 
        elif self.w_2phase+self.w_superheat>1:
            # There is no subcooled portion, solve for outlet quality
            brentq(self._TwoPhase_Forward,0.0000001,0.9999999)
            # Zero out all the subcooled parameters
            self.Q_subcool=0.0
            self.DP_r_subcool=0.0
            self.Charge_subcool=0.0
            self.w_subcool=0.0
            self.h_r_subcool=0.0
            self.existsSubcooled=False
        else:
            # By definition then we have a subcooled portion, solve for it
            self.existsSubcooled=True 
            self._Subcool_Forward()
        
        # Overall calculations
        self.Q=self.Q_superheat+self.Q_2phase+self.Q_subcool
        self.DP_r=self.DP_r_superheat+self.DP_r_2phase+self.DP_r_subcool
        self.DP_r=self.DP_r*self.DP_tuning #correcting the pressure drop
        self.Charge=self.Charge_2phase+self.Charge_subcool+self.Charge_superheat
        
        if self.existsSubcooled==True:
            # Outlet condition is in the subcooled zone 
            AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r)
            self.hout_r=AS.hmass() #[J/kg]
            self.sout_r=AS.smass() #[J/kg-K]
        elif self.w_superheat>=1:
            # Outlet condition is in the superheated zone 
            AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r)
            self.hout_r=AS.hmass() #[J/kg]
            self.sout_r=AS.smass() #[J/kg-K]   
            self.DT_sc = 0.0 #[K]                 
        else:
            # Outlet conditions is in the two-phase zone 
            self.Tout_r=self.xout_2phase*self.Tdew+(1-self.xout_2phase)*self.Tbubble
            AS.update(CP.QT_INPUTS, 0.0, self.Tout_r)
            h_l=AS.hmass() #[J/kg]
            s_l=AS.smass() #[J/kg-K]
            AS.update(CP.QT_INPUTS, 1.0, self.Tout_r)
            h_v=AS.hmass() #[J/kg]
            s_v=AS.smass() #[J/kg-K]
            self.hout_r= h_l + self.xout_2phase * (h_v - h_l)
            self.sout_r= s_l + self.xout_2phase * (s_v - s_l)
            # Use the effective subcooling
            self.DT_sc=self.DT_sc_2phase
        
        # Calculate the mean outlet air temperature [K]
        self.Tout_a=self.Tin_a-self.Q/(self.Fins.cp_da*self.Fins.mdot_da)
        self.hmean_r=self.w_2phase*self.h_r_2phase+self.w_superheat*self.h_r_superheat+self.w_subcool*self.h_r_subcool
        self.UA_r=self.hmean_r*self.A_r_wetted
        self.UA_a=(self.Fins.h_a*self.h_a_tuning)*self.Fins.A_a*self.Fins.eta_a
        self.UA_w=1/self.Rw
        
    def _Superheat_Forward(self,T_sh_out_r):
        # **********************************************************************
        #                      SUPERHEATED PART 
        # **********************************************************************
        # AbstractState
        AS = self.AS
        # Dew temperature for constant pressure cooling to saturation
        Tdew= T_sh_out_r #self.Tdew
        Tbubble=self.Tbubble
        
        # Average fluid temps are used for the calculation of properties 
        # Average temp of refrigerant is average of sat. temp and outlet temp		
        # Secondary fluid is air over the fins
        self.f_r_superheat, self.h_r_superheat, self.Re_r_superheat=f_h_1phase_Tube(self.mdot_r / self.Ncircuits, self.ID, 
            (Tdew+self.Tin_r)/2.0, self.psat_r, self.AS, "Single");
        
        AS.update(CP.PT_INPUTS, self.psat_r, (Tdew+self.Tin_r)/2)    
        cp_r = AS.cpmass() #[J/kg-K]
        rho_superheat= AS.rhomass() #[kg/m^3]

        # Compute Fins Efficiency based on FinsType 
        if self.FinsType == 'WavyLouveredFins':
            WavyLouveredFins(self.Fins)
        elif self.FinsType == 'HerringboneFins':
            HerringboneFins(self.Fins)
        elif self.FinsType == 'PlainFins':
            PlainFins(self.Fins)
            
        self.mdot_da=self.Fins.mdot_da
        
        # Cross-flow in the superheated region.  
        # Using effectiveness-Ntu relationships for cross flow with non-zero Cr.
        UA_overall = 1 / (1 / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning) + 1 / (self.h_r_superheat * self.A_r_wetted) + self.Rw)
        epsilon_superheat=(Tdew-self.Tin_r)/(self.Tin_a-self.Tin_r)
        Ntu=UA_overall/(self.mdot_da*self.Fins.cp_da)
        if epsilon_superheat>1.0:
            epsilon_superheat=1.0-1e-12
        self.w_superheat=-log(1-epsilon_superheat)*self.mdot_r*cp_r/((1-exp(-Ntu))*self.mdot_da*self.Fins.cp_da)
              
        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
        # Heat is removed here from the refrigerant since it is being cooled
        self.Q_superheat = self.mdot_r * cp_r * (Tdew-self.Tin_r)
        self.DT_sh = -self.Q_superheat/(self.mdot_r*cp_r)
        self.Tout_r = self.Tin_r-self.DT_sh
        
        # Pressure drop calculations for superheated refrigerant
        v_r=1./rho_superheat
        # Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_superheat*v_r*self.G_r**2/(2*self.ID) #Pressure gradient
        self.DP_r_superheat=dpdz_r*self.Lcircuit*self.w_superheat
        self.Charge_superheat = self.w_superheat * self.V_r * rho_superheat

        # Latent heat needed for pseudo-quality calc
        AS.update(CP.QT_INPUTS, 0.0, Tbubble)
        h_l = AS.hmass() #[J/kg]
        AS.update(CP.QT_INPUTS, 1.0, Tdew)
        h_v = AS.hmass() #[J/kg]
        h_fg = h_v - h_l #[J/kg]
        self.xin_r=1.0+cp_r*(self.Tin_r-Tdew)/h_fg
        
        return self.w_superheat-1
        
    def _TwoPhase_Forward(self,xout_r_2phase=0.0):
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
        
        # This block calculates the average refrigerant heat transfer coefficient by
        # integrating the local heat transfer coefficient between 
        # a quality of 1.0 and the outlet quality
        h_r_2phase=ShahCondensation_Average(xout_r_2phase,1.0,self.AS,self.G_r,self.ID,self.psat_r,Tbubble,Tdew);
        self.h_r_2phase=h_r_2phase*self.h_tp_tuning
        
        UA_overall = 1 / (1 / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning) + 1 / (self.h_r_2phase * self.A_r_wetted) +self.Rw);
        self.epsilon_2phase=1-exp(-UA_overall/(self.mdot_da*self.Fins.cp_da));
        self.w_2phase=-self.mdot_r*h_fg*(1.0-xout_r_2phase)/(self.mdot_da*self.Fins.cp_da*(self.Tin_a-Tsat_r)*self.epsilon_2phase);

        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
        # Heat is removed here from the refrigerant since it is condensing
        self.Q_2phase = self.epsilon_2phase * self.Fins.cp_da* self.mdot_da * self.w_2phase * (self.Tin_a-Tsat_r);
        
        self.xout_2phase=xout_r_2phase
        
        # Frictional pressure drop component
        DP_frict=LMPressureGradientAvg(self.xout_2phase,1.0,self.AS,self.G_r,self.ID,Tbubble,Tdew)*self.Lcircuit*self.w_2phase
        # Accelerational pressure drop component    
        DP_accel=-AccelPressureDrop(self.xout_2phase,1.0,self.AS,self.G_r,Tbubble,Tdew,slipModel='Zivi')*self.Lcircuit*self.w_2phase
        # Total pressure drop is the sum of accelerational and frictional components (neglecting gravitational effects)
        self.DP_r_2phase=DP_frict+DP_accel
    
        rho_average=TwoPhaseDensity(self.AS,self.xout_2phase,1.0,self.Tdew,self.Tbubble,slipModel='Zivi')
        self.Charge_2phase = rho_average * self.w_2phase * self.V_r    
        
        if self.Verbosity>7:
            print ('2phase cond resid', self.w_2phase-(1-self.w_superheat))
            print ('h_r_2phase',self.h_r_2phase)
        
        # Calculate an effective pseudo-subcooling based on the equality
        #     cp*DT_sc=-dx*h_fg
        cp_satL = self.cp_satL #[J/kg-K]
        self.DT_sc_2phase=-self.xout_2phase*h_fg/(cp_satL)
            
        # If the quality is being solved for, the length of the two-phase and subcooled
        # sections should add to the length of the HX.  Return the residual
        return self.w_2phase-(1-self.w_superheat)
    
    
    def _Subcool_Forward(self):
        self.w_subcool=1-self.w_2phase-self.w_superheat
        
        if self.w_subcool<0:
            raise ValueError('w_subcool in Condenser cannot be less than zero')
        
        # AbstractState
        AS = self.AS
        ## Bubble and dew temperatures (same for fluids without glide)
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
    
        self.f_r_subcool, self.h_r_subcool, self.Re_r_subcool=f_h_1phase_Tube(
          self.mdot_r / self.Ncircuits, self.ID, Tbubble-1.0, self.psat_r, self.AS,
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
            # Minimum capacitance rate on refrigerant side
            epsilon_subcool = 1. - exp(-1. / Cr * (1. - exp(-Cr * NTU)))
        else:
            # Minimum capacitance rate on air side
            epsilon_subcool = 1 / Cr * (1 - exp(-Cr * (1 - exp(-NTU))))
            
        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
        # Heat is removed here from the refrigerant since it is condensing
        self.Q_subcool=-epsilon_subcool*Cmin*(Tbubble-self.Tin_a)
        self.DT_sc=-self.Q_subcool/(self.mdot_r*cp_r)
        self.Tout_r=Tbubble-self.DT_sc
        
        AS.update(CP.PT_INPUTS, self.psat_r, (Tbubble + self.Tout_r) / 2)
        rho_subcool=AS.rhomass() #[kg/m^3]
        self.Charge_subcool = self.w_subcool * self.V_r * rho_subcool
    
        # Pressure drop calculations for subcooled refrigerant
        v_r=1/rho_subcool
        # Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_subcool*v_r*self.G_r**2/(2*self.ID)  #Pressure gradient
        self.DP_r_subcool=dpdz_r*self.Lcircuit*self.w_subcool
        
def SampleCondenser(AS,T=41.37):
    Fins=FinInputs()
    Fins.Tubes.NTubes_per_bank=41       #number of tubes per bank or row
    Fins.Tubes.Nbank=1                  #number of banks or rows
    Fins.Tubes.Ncircuits=5              #number of circuits
    Fins.Tubes.Ltube=2.286              #one tube length
    Fins.Tubes.OD=0.007
    Fins.Tubes.ID=0.0063904
    Fins.Tubes.Pl=0.0191                #distance between center of tubes in flow direction                                                
    Fins.Tubes.Pt=0.0222                #distance between center of tubes orthogonal to flow direction
    Fins.Tubes.kw=237                   #Wall thermal conductivity
    
    Fins.Fins.FPI=25                    #Number of fins per inch
    Fins.Fins.Pd=0.001                  #2* amplitude of wavy fin
    Fins.Fins.xf=0.001                  #1/2 period of fin
    Fins.Fins.t=0.00011                 #Thickness of fin material
    Fins.Fins.k_fin=237                 #Thermal conductivity of fin material
    
    Fins.Air.Vdot_ha=1.7934             #rated volumetric flowrate
    Fins.Air.Tmean=308.15   
    Fins.Air.Tdb=308.15                 #Dry Bulb Temperature
    Fins.Air.p=101325                   #Air pressure in Pa
    Fins.Air.RH=0.51                    #Relative Humidity
    Fins.Air.RHmean=0.51
    Fins.Air.FanPower=160    
    
    params={
        'AS': AS, #Abstract State
        'mdot_r': 0.0708,
        'Tin_r': T+20+273.15,
        'psat_r': 2500076.19, #PropsSI('P','T',T+273.15,'Q',1.0,'R410A') 
        'Fins': Fins,
        'FinsType': 'HerringboneFins',  #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'Verbosity':0,
        'h_a_tuning':1,
        'h_tp_tuning':1,
        'DP_tuning':1,
    }
    Cond=CondenserClass(**params)
    Cond.Calculate()
    return Cond
    
def SampleIntercooler(AS,T=41.37):
    # Only superheated region
    Fins=FinInputs()
    Fins.Tubes.NTubes_per_bank=41       #number of tubes per bank or row
    Fins.Tubes.Nbank=1                  #number of banks or rows
    Fins.Tubes.Ncircuits=5              #number of circuits
    Fins.Tubes.Ltube=3                  #one tube length
    Fins.Tubes.OD=0.007
    Fins.Tubes.ID=0.0063904
    Fins.Tubes.Pl=0.0191                #distance between center of tubes in flow direction                                                
    Fins.Tubes.Pt=0.0222                #distance between center of tubes orthogonal to flow direction
    Fins.Tubes.kw=237                   #Wall thermal conductivity
    
    Fins.Fins.FPI=25                    #Number of fins per inch
    Fins.Fins.Pd=0.001                  #2* amplitude of wavy fin
    Fins.Fins.xf=0.001                  #1/2 period of fin
    Fins.Fins.t=0.00011                 #Thickness of fin material
    Fins.Fins.k_fin=237                 #Thermal conductivity of fin material
    
    Fins.Air.Vdot_ha=8.998             #rated volumetric flowrate
    Fins.Air.Tmean=48.9+273.15   
    Fins.Air.Tdb=48.9+273.15                 #Dry Bulb Temperature
    Fins.Air.p=101325                   #Air pressure in Pa
    Fins.Air.RH=0.51                    #Relative Humidity
    Fins.Air.RHmean=0.51
    Fins.Air.FanPower=160    
    
    params={
        'AS': AS, #Abstract State
        'mdot_r': 1.209,
        'Tin_r': T+273.15,
        'psat_r': 716700, 
        'Fins': Fins,
        'FinsType': 'HerringboneFins',  #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'Verbosity':0,
        'h_a_tuning':1,
        'h_tp_tuning':1,
        'DP_tuning':1,
    }
    Cond=CondenserClass(**params)
    Cond.Calculate()
    return Cond    

if __name__=='__main__':
    # This runs if you run this file directly
    Ref = 'R410A'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref) #Abstract State
    Cond=SampleCondenser(AS,43.3)
    print (Cond.OutputList())
    
    print ('Heat transfer rate in condenser is', Cond.Q,'W')
    print ('Heat transfer rate in condenser (superheat section) is',Cond.Q_superheat,'W')
    print ('Heat transfer rate in condenser (twophase section) is',Cond.Q_2phase,'W')
    print ('Heat transfer rate in condenser (subcooled section) is',Cond.Q_subcool,'W')
    print ('Fraction of circuit length in superheated section is',Cond.w_superheat)
    print ('Fraction of circuit length in twophase section is',Cond.w_2phase)
    print ('Fraction of circuit length in subcooled section is',Cond.w_subcool) 
    print ()
    
    
    #This runs the intercooler
    Ref = 'R134A'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref) #Abstract State
    Cond=SampleIntercooler(AS,85.72)
    print (Cond.OutputList())
    
    print ('Heat transfer rate in condenser is', Cond.Q,'W')
    print ('Heat transfer rate in condenser (superheat section) is',Cond.Q_superheat,'W')
    print ('Heat transfer rate in condenser (twophase section) is',Cond.Q_2phase,'W')
    print ('Heat transfer rate in condenser (subcooled section) is',Cond.Q_subcool,'W')
    print ('Fraction of circuit length in superheated section is',Cond.w_superheat)
    print ('Fraction of circuit length in twophase section is',Cond.w_2phase)
    print ('Fraction of circuit length in subcooled section is',Cond.w_subcool)
    print ('Refrigerant Inlet Temperature is',Cond.Tin_r - 273.15) 
    print ('Refrigerant Outlet Temperature is',Cond.Tout_r - 273.15)