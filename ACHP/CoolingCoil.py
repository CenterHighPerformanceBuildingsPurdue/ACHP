from __future__ import division #Make integer 3/2 give 1.5 in python 2.x
from math import pi
from CoolProp.CoolProp import PropsSI
from Correlations import f_h_1phase_Tube
from FinCorrelations import WavyLouveredFins,HerringboneFins, PlainFins, FinInputs 
from matplotlib import docstring
from BaseDoc import BaseDocClass
from DryWetSegment import DWSVals, DryWetSegment

class CoolingCoilClass():
    """
    The module that implements a cooling coil.  See documentation for further information
    """
    
    def __init__(self,**kwargs):
        """Load the parameters passed in using the dictionary"""
        self.__dict__.update(kwargs)
    
    @docstring.copy_dedent(BaseDocClass.Update) #Use docs from Base Class
    def Update(self,**kwargs):
        """Update the parameters passed in using the dictionary"""
        self.__dict__.update(kwargs)
    
    @docstring.copy_dedent(BaseDocClass.OutputList) #Use docs from Base Class
    def OutputList(self):
        return [
            ('Volumetric flow rate','m^3/s',self.Fins.Air.Vdot_ha),
            ('Inlet Dry bulb temp','K',self.Tin_a),
            ('Inlet Air pressure','kPa',self.Fins.Air.p),
            ('Inlet Air Relative Humidity','-',self.Fins.Air.RH),
            ('Tubes per bank','-',self.Fins.Tubes.NTubes_per_bank),
            ('Number of banks','-',self.Fins.Tubes.Nbank),
            ('Number circuits','-',self.Fins.Tubes.Ncircuits),
            ('Length of tube','m',self.Fins.Tubes.Ltube),
            ('Tube OD','m',self.Fins.Tubes.OD),
            ('Tube ID','m',self.Fins.Tubes.ID),
            ('Tube Long. Pitch','m',self.Fins.Tubes.Pl),
            ('Tube Transverse Pitch','m',self.Fins.Tubes.Pt),
            ('Fins per inch','1/in',self.Fins.Fins.FPI),
            ('Fin waviness pd','m',self.Fins.Fins.Pd),
            ('Fin waviness xf','m',self.Fins.Fins.xf),
            ('Fin thickness','m',self.Fins.Fins.t),
            ('Fin Conductivity','W/m-K',self.Fins.Fins.k_fin),
            ('Q Total','W',self.Q),
            ('Mean HTC glycol','W/m^2-K',self.h_g),
            ('Reynolds # glycol','-',self.Re_g),
            ('Pressure drop glycol','Pa',self.DP_g),
            ('Inlet glycol temp','K',self.Tin_g),
            ('Outlet glyol temp','K',self.Tout_g),
            ('Outlet air temp','K',self.Tout_a),
            ('Mean Air HTC','W/m^2-K',self.Fins.h_a),
            ('Surface Effectiveness','-',self.Fins.eta_a),
            ('Air-side area (fin+tubes)','m^2',self.Fins.A_a),
            ('Mass Flow rate Air','kg/s',self.Fins.mdot_da),
            ('Pressure Drop Air-side','Pa',self.Fins.dP_a),
            ('Wetted area fraction','-',self.f_dry),
            ('Sensible Heat Ratio','-',self.SHR)]
     
    def Initialize(self):
        self.Update()
        
        # Retrieve some parameters from nested structures 
        # for code compactness
        self.ID=self.Fins.Tubes.ID
        self.OD=self.Fins.Tubes.OD
        self.Ltube=self.Fins.Tubes.Ltube
        self.NTubes_per_bank=self.Fins.Tubes.NTubes_per_bank
        self.Nbank=self.Fins.Tubes.Nbank
        self.Ncircuits=self.Fins.Tubes.Ncircuits
        self.Tin_a=self.Fins.Air.Tdb
        self.pin_a=self.Fins.Air.p
        self.RHin_a=self.Fins.Air.RH
        
        # Calculate an effective length of circuit if circuits are 
        # not all the same length
        TotalLength=self.Ltube*self.NTubes_per_bank*self.Nbank
        self.Lcircuit=TotalLength/self.Ncircuits
        # Wetted area on the glycol side
        self.A_g_wetted=self.Ncircuits*pi*self.ID*self.Lcircuit
        
        # Evaluate the air-side heat transfer and pressure drop
        if self.FinsType == 'WavyLouveredFins':
            WavyLouveredFins(self.Fins)
        elif self.FinsType == 'HerringboneFins':
            HerringboneFins(self.Fins)
        elif self.FinsType == 'PlainFins':
            PlainFins(self.Fins)
        
    def Calculate(self):
        """
        This function is now simply a wrapper around the DryWetSegment() 
        function in order to decrease the amount of code replication
        """
        
        self.Initialize()
    
        DWS=DWSVals() #DryWetSegment structure
    
        # Store temporary values to be passed to DryWetSegment
        DWS.A_a=self.Fins.A_a
        DWS.cp_da=self.Fins.cp_da
        DWS.eta_a=self.Fins.eta_a
        DWS.h_a=self.Fins.h_a  #Heat transfer coefficient
        DWS.mdot_da=self.Fins.mdot_da
        DWS.pin_a=self.Fins.Air.p
        DWS.Tin_a=self.Tin_a
        DWS.RHin_a=self.Fins.Air.RH
        DWS.Fins=self.Fins
    
        DWS.Tin_r=self.Tin_g
        DWS.A_r=self.A_g_wetted
        DWS.cp_r=PropsSI('C','T',(self.Tin_g+DWS.Tin_a)/2.0, 'P', self.pin_g, self.Ref_g)#*1000 #Use a guess value of 6K superheat to calculate cp 
        DWS.pin_r=self.pin_g
        DWS.mdot_r=self.mdot_g
        DWS.IsTwoPhase=False
        
        #Use a guess value of 6K superheat to calculate the properties
        self.f_g, self.h_g, self.Re_g=f_h_1phase_Tube(self.mdot_g / self.Ncircuits, self.ID, 
            (self.Tin_g+DWS.Tin_a)/2.0, self.pin_g, self.Ref_g, "Single");
        
        # Average Refrigerant heat transfer coefficient
        DWS.h_r=self.h_g
        
        #Run DryWetSegment
        DryWetSegment(DWS)
        
        #Average mass flux of glycol in circuit
        self.G_g = self.mdot_g/(self.Ncircuits*pi*self.ID**2/4.0) #[kg/m^2-s]
    
        #Pressure drop calculations for glycol (water)
        Dh_g=self.ID
        v_g=1/PropsSI('D','T',self.Tin_g, 'P',self.pin_g,self.Ref_g)
        #Pressure gradient using Darcy friction factor
        dp_dz_g=-self.f_g*v_g*self.G_g**2/(2*Dh_g)
        DP_g=dp_dz_g*self.Lcircuit
    
        self.f_dry=DWS.f_dry
        self.DP_g=DP_g
        self.Q=DWS.Q
        self.Tout_g=DWS.Tout_r
        self.Tout_a=DWS.Tout_a
        self.hout_a=DWS.hout_a
        self.hin_a=DWS.hin_a
        self.SHR=self.Fins.cp_da*(DWS.Tout_a-DWS.Tin_a)/(DWS.hout_a-DWS.hin_a)
        self.Capacity=DWS.Q-self.Fins.Air.FanPower

def TestCase():
    CC=CoolingCoilClass()
    FinsTubes=FinInputs()
    FinsTubes.Tubes.NTubes_per_bank=32
    FinsTubes.Tubes.Nbank=3
    FinsTubes.Tubes.Ncircuits=5
    FinsTubes.Tubes.Ltube=0.452
    FinsTubes.Tubes.OD=0.009525
    FinsTubes.Tubes.ID=0.0089154
    FinsTubes.Tubes.Pl=0.0254
    FinsTubes.Tubes.Pt=0.0219964
    
    FinsTubes.Fins.FPI=14.5
    FinsTubes.Fins.Pd=0.001
    FinsTubes.Fins.xf=0.001
    FinsTubes.Fins.t=0.00011
    FinsTubes.Fins.k_fin=237
    
    FinsTubes.Air.Vdot_ha=0.5663
    FinsTubes.Air.Tmean=299.8
    FinsTubes.Air.Tdb= 299.8
    FinsTubes.Air.p=101325              #Air pressure in Pa
    FinsTubes.Air.RH=0.51
    FinsTubes.Air.RHmean=0.51
    FinsTubes.Air.FanPower=438
        
    CC.Fins = FinsTubes
    CC.FinsType = 'WavyLouveredFins'    #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
    CC.Ref_g = 'Water'
    CC.mdot_g = 0.15
    CC.Tin_g = 278
    CC.pin_g = 300000                    #Refrigerant vapor pressure in Pa
    CC.Verbosity = 3
    
    CC.Calculate()
    print CC.OutputList()
        
if __name__=='__main__':
    TestCase()