from __future__ import division, print_function, absolute_import
from math import pi,log,exp

from ACHP.Correlations import f_h_1phase_Tube, f_h_cp_supercritical
from ACHP.FinCorrelations import WavyLouveredFins,FinInputs,IsFinsClass, HerringboneFins, PlainFins
from ACHP.ACHPTools import ValidateFields

from scipy.optimize import brentq
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
            ('Fins per inch','1/in',self.Fins.Fins.FPI),
            ('Fin waviness pd','m',self.Fins.Fins.Pd),
            ('Fin waviness xf','m',self.Fins.Fins.xf),
            ('Fin thickness','m',self.Fins.Fins.t),
            ('Fin Conductivity','W/m-K',self.Fins.Fins.k_fin),
            ('Fins Type','-',self.FinsType),
            ('Q Total','W',self.Q),
            ('Q Supercritical','W',self.Q_supercritical),
            ('Q Supercritical_liquid','W',self.Q_subcool),
            ('Inlet Temp','K',self.Tin_r),
            ('Outlet Temp','K',self.Tout_r),
            ('Pressure Drop Total','Pa',self.DP_r),
            ('Pressure Drop Supercritical','Pa',self.DP_r_supercritical),
            ('Pressure Drop Supercritical_liquid','Pa',self.DP_r_subcool),
            ('Charge Total','kg',self.Charge),
            ('Charge Supercritical','kg',self.Charge_supercritical),
            ('Charge Supercritical_liquid','kg',self.Charge_subcool),
            ('Mean HTC Superheat','W/m^2-K',self.h_r_supercritical),
            ('Mean HTC Supercritical_liquid','W/m^2-K',self.h_r_subcool),
            ('Wetted Area Fraction Supercritical','-',self.w_supercritical),
            ('Wetted Area Fraction Supercritical_liquid','-',self.w_subcool),
            ('Mean Air HTC','W/m^2-K',self.Fins.h_a),
            ('Surface Effectiveness','-',self.Fins.eta_a),
            ('Air-side area (fin+tubes)','m^2',self.Fins.A_a),
            ('Mass Flow rate of Dry Air','kg/s',self.Fins.mdot_da),
            ('Mass Flow rate of Humid Air','kg/s',self.Fins.mdot_ha),
            ('Pressure Drop Air-side','Pa',self.Fins.dP_a),
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
               ('Ref',str,None,None),
               ('Fins',IsFinsClass,None,None),
               ('FinsType',str,None,None),
               ('mdot_r',float,0.00001,20),
               ('Tin_r',float,200,500),
               ('psat_r',float,0.01,20000000)
               ]
            optFields=['Verbosity','Backend']
            ValidateFields(self.__dict__,reqFields,optFields)
            self.IsValidated=True
        
        #AbstractState
        if hasattr(self,'Backend'): #check if backend is given
            AS = CP.AbstractState(self.Backend, self.Ref)
        else: #otherwise, use the defualt backend
            AS = CP.AbstractState('HEOS', self.Ref)
        self.AS = AS
        
        # Retrieve some parameters from nested structures 
        # for code compactness
        self.ID=self.Fins.Tubes.ID
        self.OD=self.Fins.Tubes.OD
        self.Ltube=self.Fins.Tubes.Ltube
        self.NTubes_per_bank=self.Fins.Tubes.NTubes_per_bank
        self.Nbank=self.Fins.Tubes.Nbank
        self.Ncircuits=self.Fins.Tubes.Ncircuits
        self.Tin_a=self.Fins.Air.Tdb

        # Calculate an effective length of circuit if circuits are 
        # not all the same length
        TotalLength=self.Ltube*self.NTubes_per_bank*self.Nbank
        self.Lcircuit=TotalLength/self.Ncircuits
        
        self.V_r = pi * self.ID**2 / 4.0 * self.Lcircuit * self.Ncircuits
        self.A_r_wetted = pi * self.ID * self.Ncircuits * self.Lcircuit
        self.G_r = self.mdot_r/(self.Ncircuits*pi*self.ID**2/4.0)    
        
        # Define known parameters
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tin_r)
        self.hin_r=AS.hmass() #[J/kg]
        self.sin_r=AS.smass() #[J/kg-K]
        
        # Define critical pressure and temperature
        self.Pcr=AS.p_critical() #[Pa]
        self.Tcr=AS.T_critical() #[K]
        
        
        #try to run with a full supercritical section from Temperature of Tin_r to Tcr
        self._Supercritical_Forward(self.Tcr)

        #If we have already used too much of the HX (max possible sum of w is 1.0)
        if self.w_supercritical>1:
            #There is no supercritical_liquid portion, solve for outlet Temperature
            brentq(self._Supercritical_Forward,self.Tcr,self.Tin_r)
            #Zero out all the supercritical_liquid parameters
            self.Q_subcool=0.0
            self.DP_r_subcool=0.0
            self.Charge_subcool=0.0
            self.w_subcool=0.0
            self.h_r_subcool=0.0
            self.existsSubcooled=False
        else:
            #By definition then we have a supercritical_liquid portion, solve for it
            self.existsSubcooled=True 
            self._Subcool_Forward()
        
        #Overall calculations
        self.Q=self.Q_supercritical+self.Q_subcool
        self.DP_r=self.DP_r_supercritical+self.DP_r_subcool
        self.Charge=self.Charge_supercritical+self.Charge_subcool
        
        if self.existsSubcooled==True: #exit at supercritical_liquid region
            AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r)
            self.hout_r=AS.hmass() #[J/kg]
            self.sout_r=AS.smass() #[J/kg-K]
        else: #still in supercritical region
            self.Tout_r=self.Tout_r_cr
            AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r)
            self.hout_r=AS.hmass() #[J/kg]
            self.sout_r=AS.smass() #[J/kg-K]

        #Calculate the mean outlet air temperature [K]
        self.Tout_a=self.Tin_a-self.Q/(self.Fins.cp_da*self.Fins.mdot_da)
        self.hmean_r=self.w_supercritical*self.h_r_supercritical+self.w_subcool*self.h_r_subcool
        self.UA_r=self.hmean_r*self.A_r_wetted
        self.UA_a=self.Fins.h_a*self.Fins.A_a*self.Fins.eta_a
        
        
    def _Supercritical_Forward(self,Tout_r_cr):
        """
            Tout_r_cr: Temperature of refrigerant at end of supercritical portion
                default value is Tcr (full supercritical region)
        """
        #AbstractState
        AS = self.AS
        
        # This block calculates the average refrigerant heat transfer coefficient, average friction factor, average specific heat, and average density
        self.h_r_supercritical, self.f_r_supercritical, cp_r, rho_supercritical = f_h_cp_supercritical(Tout_r_cr, self.Tin_r, self.AS, self.OD, self.ID,  self.mdot_r / self.Ncircuits, self.psat_r);
        
        #Compute Fins Efficiency based on FinsType 
        if self.FinsType == 'WavyLouveredFins':
            WavyLouveredFins(self.Fins)
        elif self.FinsType == 'HerringboneFins':
            HerringboneFins(self.Fins)
        elif self.FinsType == 'PlainFins':
            PlainFins(self.Fins)
            
        self.mdot_da=self.Fins.mdot_da
        
        UA_overall = 1 / (1 / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a) + 1 / (self.h_r_supercritical * self.A_r_wetted));
        epsilon_supercritical=(Tout_r_cr-self.Tin_r)/(self.Tin_a-self.Tin_r)
        Ntu=UA_overall/(self.mdot_da*self.Fins.cp_da)
        if epsilon_supercritical>1.0:
            epsilon_supercritical=1.0-1e-12
        self.w_supercritical=-log(1-epsilon_supercritical)*self.mdot_r*cp_r/((1-exp(-Ntu))*self.mdot_da*self.Fins.cp_da)
        
        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
        # Heat is removed here from the refrigerant since it is being cooled
        self.Q_supercritical = self.mdot_r * cp_r * (Tout_r_cr-self.Tin_r)

        self.Tout_r_cr=Tout_r_cr
        
        #Pressure drop calculations for supercritical refrigerant
        v_r=1./rho_supercritical
        #Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_supercritical*v_r*self.G_r**2/(2*self.ID) #Pressure gradient
        self.DP_r_supercritical=dpdz_r*self.Lcircuit*self.w_supercritical
        #charge for the supercritical portion
        self.Charge_supercritical = self.w_supercritical * self.V_r * rho_supercritical    
        
        if self.Verbosity>7:
            print('supercritical cond resid', 1-self.w_supercritical)
            print('h_r_supercritical',self.h_r_suopercritical)
        
        
        return 1-self.w_supercritical
    
    
    def _Subcool_Forward(self):
        self.w_subcool=1-self.w_supercritical
        
        if self.w_subcool<0:
            raise ValueError('w_subcool in Condenser cannot be less than zero')
        
        #AbstractState
        AS = self.AS
        ## Bubble and dew temperatures (same for fluids without glide)
        Tcr=self.Tcr
        
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
          self.mdot_r / self.Ncircuits, self.ID, Tcr-1.0, self.psat_r, self.AS,
          "Single")
        
        AS.update(CP.PT_INPUTS, self.psat_r, Tcr-1)
        cp_r = AS.cpmass() #[J/kg-K]
    
        # Cross-flow in the subcooled region.
        R_a=1. / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a)
        R_r=1. / (self.h_r_subcool * self.A_r_wetted)
        UA_subcool = self.w_subcool / (R_a + R_r)
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
            
        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
        # Heat is removed here from the refrigerant since it is condensing
        self.Q_subcool=-epsilon_subcool*Cmin*(Tcr-self.Tin_a)
        self.DT_sc=-self.Q_subcool/(self.mdot_r*cp_r)
        self.Tout_r=Tcr-self.DT_sc
        
        AS.update(CP.PT_INPUTS, self.psat_r, (Tcr + self.Tout_r) / 2)
        rho_subcool=AS.rhomass() #[kg/m^3]
        self.Charge_subcool = self.w_subcool * self.V_r * rho_subcool
    
        #Pressure drop calculations for subcooled refrigerant
        v_r=1/rho_subcool
        #Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_subcool*v_r*self.G_r**2/(2*self.ID)  #Pressure gradient
        self.DP_r_subcool=dpdz_r*self.Lcircuit*self.w_subcool
        
def SampleCondenser():
    Fins=FinInputs()
    Fins.Tubes.NTubes_per_bank=18       #number of tubes per bank or row
    Fins.Tubes.Nbank=3                  #number of banks or rows
    Fins.Tubes.Ncircuits=1              #number of circuits
    Fins.Tubes.Ltube=0.61           #one tube length
    Fins.Tubes.OD=7.9/1000
    Fins.Tubes.ID=7.5/1000
    Fins.Tubes.Pl=19/1000                #distance between center of tubes in flow direction                                                
    Fins.Tubes.Pt=25/1000                #distance between center of tubes orthogonal to flow direction
    
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
        'Ref': 'R744',
        'mdot_r': 0.076,
        'Tin_r': 110+273.15,
        'psat_r': 11000000,
        'Fins': Fins,
        'FinsType': 'WavyLouveredFins',  #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'Verbosity':0,
        'Backend':'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    }
    Cond=CondenserClass(**params)
    Cond.Calculate()
    return Cond
    
if __name__=='__main__':
    #This runs if you run this file directly
    Cond=SampleCondenser()
    print(Cond.OutputList())
    
    print('Heat transfer rate in gas cooler is', Cond.Q,'W')
    print('Heat transfer rate in gas cooler (supercritical section) is',Cond.Q_supercritical,'W')
    print('Heat transfer rate in gas cooler (supercritical_liquid section) is',Cond.Q_subcool,'W')
    print('Fraction of circuit length in supercritical section is',Cond.w_supercritical)
    print('Fraction of circuit length in supercritical_liquid section is',Cond.w_subcool)
    print('Refrigerant outlet temperature is',Cond.Tout_r-273.15, 'C')
    print('Air outlet temperature is',Cond.Tout_a-273.15, 'C')