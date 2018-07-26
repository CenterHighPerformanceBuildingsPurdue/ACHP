from __future__ import division, print_function, absolute_import
from math import pi,log,exp

from ACHP.Correlations import f_h_1phase_Tube, Petterson_supercritical_average
from ACHP.FinCorrelations import WavyLouveredFins,FinInputs,IsFinsClass, HerringboneFins, PlainFins
from ACHP.DryWetSegment import DWSVals, DryWetSegment
from ACHP.ACHPTools import ValidateFields

from scipy.optimize import brentq
import CoolProp as CP

class FinVals():
    def __init__(self):
        pass
    
class GasCoolerClass():
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
            ('Pressure Drop Air-side','Pa',self.Fins.dP_a),
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
               ('psat_r',float,0.01,20000000)
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
        
        # Define critical pressure and temperature
        self.Pcr=AS.p_critical() #[Pa]
        self.Tcr=AS.T_critical() #[K]
        
        #critical enthalpy at defined pressure
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tcr)
        self.hcr = AS.hmass() #[J/kg]
        
        #triple temperature
        self.Ttriple = AS.Ttriple()
        
        self.Fins.Air.RHmean=self.Fins.Air.RH
        
        #Update with user FinType
        if self.FinsType == 'WavyLouveredFins':
            WavyLouveredFins(self.Fins)
        elif self.FinsType == 'HerringboneFins':
            HerringboneFins(self.Fins)
        elif self.FinsType == 'PlainFins':
            PlainFins(self.Fins)
        
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
        self.DP_r=(self.DP_r_supercritical+self.DP_r_supercrit_liq)*self.DP_tuning
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
        h_r, f_r_supercritical, DWS.cp_r, rho_supercritical = Petterson_supercritical_average(self.Tout_r_cr, self.Tin_r, self.T_w, self.AS, self.G_r, self.ID, 0, self.ID/self.Lcircuit, self.mdot_r / self.Ncircuits, self.psat_r, -Q_target/DWS.A_r);
        DWS.h_r = h_r*self.h_r_tuning #Correct refrigerant HTC with tuning factor
        
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
        dpdz_r=-f_r_supercritical*v_r*self.G_r**2/(2*self.ID) #Pressure gradient
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
#         self.f_r_supercrit_liq, self.h_r_supercrit_liq, self.Re_r_supercrit_liq=f_h_1phase_Tube(
#           self.mdot_r / self.Ncircuits, self.ID, self.Tcr-1.0, self.psat_r, self.AS,
#           "Single")
#         
#         # Average Refrigerant heat transfer coefficient
#         DWS.h_r=self.h_r_supercrit_liq*self.h_r_tuning
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
#         dpdz_r=-self.f_r_supercrit_liq*v_r*self.G_r**2/(2*self.ID)  #Pressure gradient
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
        h_r, self.f_r_supercrit_liq, DWS.cp_r, rho_supercrit_liq = Petterson_supercritical_average(self.Tout_r_sc, self.Tcr, self.T_w, self.AS, self.G_r, self.ID, 0, self.ID/self.Lcircuit, self.mdot_r / self.Ncircuits, self.psat_r, -Q_target/DWS.A_r);
        DWS.h_r = h_r*self.h_r_tuning #Correct refrigerant HTC with tuning factor
        
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
        dpdz_r=-self.f_r_supercrit_liq*v_r*self.G_r**2/(2*self.ID)  #Pressure gradient
        self.DP_r_supercrit_liq=dpdz_r*self.Lcircuit*self.w_supercrit_liq
        
        return Q_target-DWS.Q
    
         
def SampleGasCooler(AS):
    Fins=FinInputs()
    Fins.Tubes.NTubes_per_bank=18       #number of tubes per bank or row
    Fins.Tubes.Nbank=3                  #number of banks or rows
    Fins.Tubes.Ncircuits=1              #number of circuits
    Fins.Tubes.Ltube=0.61           #one tube length
    Fins.Tubes.OD=7.9/1000
    Fins.Tubes.ID=7.5/1000
    Fins.Tubes.Pl=19/1000                #distance between center of tubes in flow direction                                                
    Fins.Tubes.Pt=25/1000                #distance between center of tubes orthogonal to flow direction
    Fins.Tubes.kw=237                  #wall thermal conductivity (i.e pipe material)
    
    Fins.Fins.FPI=1/(1.5/1000/0.0254)   #Number of fins per inch
    Fins.Fins.Pd=0.001                  #2* amplitude of wavy fin
    Fins.Fins.xf=0.001                  #1/2 period of fin
    Fins.Fins.t=0.13/1000                #Thickness of fin material
    Fins.Fins.k_fin=237                 #Thermal conductivity of fin material
    
    Fins.Air.Vdot_ha=0.281                #rated volumetric flowrate (m^3/s)
    Fins.Air.Tmean=29.4+273.15   
    Fins.Air.Tdb=29.4+273.15            #Dry Bulb Temperature
    Fins.Air.p=101325                   #Air pressure in Pa
    Fins.Air.RH=0.5                     #Relative Humidity
    Fins.Air.RHmean=0.5
    Fins.Air.FanPower=160    
    
    params={
        'AS': AS,
        'mdot_r': 0.076,
        'Tin_r': 110.6+273.15,
        'psat_r': 11000000,
        'Fins': Fins,
        'FinsType': 'WavyLouveredFins',  #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'Verbosity':0,
        'h_a_tuning':1,
        'h_r_tuning':1,
        'DP_tuning':1,
    }
    Cond=GasCoolerClass(**params)
    Cond.Calculate()
    return Cond
    
if __name__=='__main__':
    #This runs if you run this file directly
    from time import time
    t1=time()      
    Ref = 'R744'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref) #Abstract State  
    Cond=SampleGasCooler(AS)

    #print(Cond.OutputList())
    
    print('Heat transfer rate in gas cooler is', Cond.Q,'W')
    print('Heat transfer rate in gas cooler (supercritical section) is',Cond.Q_supercritical,'W')
    print('Heat transfer rate in gas cooler (supercritical_liquid section) is',Cond.Q_supercrit_liq,'W')
    print('Fraction of circuit length in supercritical section is',Cond.w_supercritical)
    print('Fraction of circuit length in supercritical_liquid section is',Cond.w_supercrit_liq)
    print('Refrigerant outlet temperature is',Cond.Tout_r-273.15, 'C')
    print('Air outlet temperature is',Cond.Tout_a-273.15, 'C')
    print ('Took '+str(time()-t1)+' seconds to run Gas Cooler model')