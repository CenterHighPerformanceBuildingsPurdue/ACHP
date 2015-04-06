from __future__ import division #Make integer 3/2 give 1.5 in python 2.x
from math import pi,log,exp
from CoolProp.CoolProp import PropsSI #,UseSaturationLUT
from Correlations import f_h_1phase_Tube,ShahEvaporation_Average, LockhartMartinelli,LMPressureGradientAvg,AccelPressureDrop,TwoPhaseDensity
from scipy.optimize import brentq #solver to find roots (zero points) of functions
from scipy.interpolate import interp1d
from FinCorrelations import WavyLouveredFins,FinInputs,IsFinsClass, HerringboneFins, PlainFins
from DryWetSegment import DWSVals, DryWetSegment
from ACHPTools import ValidateFields
import numpy as np
# Turn on saturation curve lookup for CoolProp
#UseSaturationLUT(1)

class EvaporatorClass():
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)
    def Update(self,**kwargs):
        self.__dict__.update(kwargs)
        
    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        Output_List=[]
        #append optional parameters, if applicable
        if hasattr(self,'TestName'):
            Output_List.append(('Name','N/A',self.TestName)) 
        if hasattr(self,'TestDescription'):
            Output_List.append(('Description','N/A',self.TestDescription))
        if hasattr(self,'TestDetails'):
            Output_List.append(('Details','N/A',self.TestDetails))
        Output_List_default=[                                                   #default output list
            ('Volumetric flow rate','m^3/s',self.Fins.Air.Vdot_ha),
            ('Inlet Dry bulb temp','K',self.Tin_a),
            ('Inlet Air pressure','Pa',self.Fins.Air.p),
            ('Inlet Air Relative Humidity','-',self.Fins.Air.RH),
            ('Tubes per bank','-',self.Fins.Tubes.NTubes_per_bank),
            ('Number of banks','-',self.Fins.Tubes.Nbank),
            ('Number circuits','-',self.Fins.Tubes.Ncircuits),
            ('Length of tube','m',self.Fins.Tubes.Ltube),
            ('Tube OD','m',self.Fins.Tubes.OD),
            ('Tube ID','m',self.Fins.Tubes.ID),
            ('Tube Long. Pitch','m',self.Fins.Tubes.Pl),
            ('Tube Transverse Pitch','m',self.Fins.Tubes.Pt),
            ('Outlet superheat','K',self.Tout_r-self.Tdew_r),
            ('Fins per inch','1/in',self.Fins.Fins.FPI),
            ('Fin waviness pd','m',self.Fins.Fins.Pd),
            ('Fin waviness xf','m',self.Fins.Fins.xf),
            ('Fin thickness','m',self.Fins.Fins.t),
            ('Fin Conductivity','W/m-K',self.Fins.Fins.k_fin),
            ('Fins Type','-',self.FinsType),
            ('Q Total','W',self.Q),
            ('Q Superheat','W',self.Q_superheat),
            ('Q Two-Phase','W',self.Q_2phase),
            ('Inlet ref. temp','K',self.Tin_r),
            ('Outlet ref. temp','K',self.Tout_r),
            ('Outlet air temp','K',self.Tout_a),
            ('Evaporator P_sat in','Pa',self.psat_r),
            ('Evaporator inlet quality','-',self.xin_r),
            ('Evaporator ref. flowrate','kg/s',self.mdot_r),
            ('Pressure Drop Total','Pa',self.DP_r),
            ('Pressure Drop Superheat','Pa',self.DP_r_superheat),
            ('Pressure Drop Two-Phase','Pa',self.DP_r_2phase),
            ('Charge Total','kg',self.Charge),
            ('Charge Superheat','kg',self.Charge_superheat),
            ('Charge Two-Phase','kg',self.Charge_2phase),
            ('Mean HTC Superheat','W/m^2-K',self.h_r_superheat),
            ('Mean HTC Two-phase','W/m^2-K',self.h_r_2phase),
            ('Wetted Area Fraction Superheat','-',self.w_superheat),
            ('Wetted Area Fraction Two-phase','-',self.w_2phase),
            ('Mean Air HTC','W/m^2-K',self.Fins.h_a),
            ('Surface Effectiveness','-',self.Fins.eta_a),
            ('Air-side area (fin+tubes)','m^2',self.Fins.A_a),
            ('Mass Flow rate dry Air','kg/s',self.Fins.mdot_da),
            ('Mass Flow rate humid Air','kg/s',self.Fins.mdot_ha),
            ('Pressure Drop Air-side','Pa',self.Fins.dP_a),
            ('Sensible Heat Ratio','-',self.SHR),
            ('Bend Temperature profile',self.Tbends)
        ]
        for i in range(0,len(Output_List_default)):                             #append default parameters to output list
            Output_List.append(Output_List_default[i])
        return Output_List
        
    def AirSideCalcs(self):
        #Update with user FinType
        if self.FinsType == 'WavyLouveredFins':
            WavyLouveredFins(self.Fins)
        elif self.FinsType == 'HerringboneFins':
            HerringboneFins(self.Fins)
        elif self.FinsType == 'PlainFins':
            PlainFins(self.Fins)
    
    def Initialize(self):
        #Input validation the first call of Initialize
        if False:#not hasattr(self,'IsValidated'):
            self.Fins.Validate()
            reqFields=[
                       ('Ref',str,None,None),
                       ('psat_r',float,0.001,100000000),                        #0.000001,100000 is updated to 0.001,100000000
                       ('Fins',IsFinsClass,None,None),
                       ('FinsType',str,None,None),
                       ('hin_r',float,-100000,10000000),
                       ('mdot_r',float,0.000001,10),
                       ]
            optFields=['Verbosity']
            d=self.__dict__ #Current fields in model
            ValidateFields(d,reqFields,optFields)
            self.IsValidated=True
            
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
        # Wetted area on the refrigerant side
        self.A_r_wetted=self.Ncircuits*pi*self.ID*self.Lcircuit
        self.V_r=self.Ncircuits*self.Lcircuit*pi*self.ID**2/4.0
        #Average mass flux of refrigerant in circuit
        self.G_r = self.mdot_r/(self.Ncircuits*pi*self.ID**2/4.0) #[kg/m^2-s]
        
        """
        Tsat() is a relatively slow function since it does a Dekker solve
        over the full two-phase region.  So store the value in order to cut
        down on computational work. 
        """
        ## Bubble and dew temperatures (same for fluids without glide) 
        self.Tbubble_r=PropsSI('T','P',self.psat_r,'Q',0,self.Ref)
        self.Tdew_r=PropsSI('T','P',self.psat_r,'Q',1,self.Ref)
        ## Mean temperature for use in HT relationships
        self.Tsat_r=(self.Tbubble_r+self.Tdew_r)/2
        # Latent heat
        self.h_fg=(PropsSI('H','T',self.Tdew_r,'Q',1.0,self.Ref)-PropsSI('H','T',self.Tbubble_r,'Q',0.0,self.Ref))*1. #*1000. #[J/kg]
        
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
        
        self.Initialize()
        
        # Input and output thermodynamic properties
        ssatL=PropsSI('S','T',self.Tbubble_r,'Q',0,self.Ref) #*1000
        ssatV=PropsSI('S','T',self.Tdew_r,'Q',1,self.Ref) #*1000
        hsatL=PropsSI('H','T',self.Tbubble_r,'Q',0,self.Ref) #*1000
        hsatV=PropsSI('H','T',self.Tdew_r,'Q',1,self.Ref) #*1000
        
        #Must give enthalpy and pressure as inputs
        self.xin_r=(self.hin_r-hsatL)/(hsatV-hsatL)
        self.sin_r=self.xin_r*ssatV+(1-self.xin_r)*ssatL
        self.hin_r=self.xin_r*hsatV+(1-self.xin_r)*hsatL
        self.Tin_r=self.xin_r*self.Tdew_r+(1-self.xin_r)*self.Tbubble_r
        
        if self.xin_r>1.0:
            raise ValueError
        #Begin by assuming that you go all the way to saturated vapor at least
        self.xout_2phase=1.0
        if self._TwoPhase_Forward(1.0)<0:
            # Evaporator outlet is in two-phase region, use all the area and find the outlet quality
            existsSuperheat=False
            self.w_2phase=1.0
            def OBJECTIVE(xout):
                self.xout_2phase=xout
                Q_target=self.mdot_r*(xout-self.xin_r)*(hsatV-hsatL)
                self._TwoPhase_Forward(self.w_2phase)
                return self.Q_2phase-Q_target
            #Use a solver to find outlet quality
            brentq(OBJECTIVE,self.xin_r,1.0)
            self.w_superheat=0.0
            self.Q_superheat=0.0
            self.h_r_superheat=0.0
            self.Re_r_superheat=0.0
            self.Charge_superheat=0.0
            self.Q_sensible_superheat=0.0
            self.Tout_a_superheat=0.0
            self.DP_r_superheat=0.0
            self.fdry_superheat=0.0
        else:
            existsSuperheat=True
            # Evaporator outlet is in superheated region, everything is ok
            self.w_2phase=brentq(self._TwoPhase_Forward,0.00000000001,0.9999999999)
            self._Superheat_Forward(1-self.w_2phase)
        
        self.Q=self.Q_superheat+self.Q_2phase
        self.Charge=self.Charge_superheat+self.Charge_2phase
        if self.Verbosity>4: 
            print self.Q,"Evaporator.Q"
        self.Capacity=self.Q-self.Fins.Air.FanPower
        
        #Sensible heat ratio [-]
        self.SHR=(self.Q_sensible_2phase+self.Q_sensible_superheat)/self.Q
        #Average air outlet temperature (area fraction weighted average) [K]
        self.Tout_a=self.w_superheat*self.Tout_a_superheat+self.w_2phase*self.Tout_a_2phase
        self.DP_r=self.DP_r_superheat+self.DP_r_2phase
        
        #Outlet enthalpy obtained from energy balance
        self.hout_r=self.hin_r+self.Q/self.mdot_r
        
        #Outlet entropy
        if existsSuperheat==True:
            self.sout_r=PropsSI('S','T',self.Tout_r,'P',self.psat_r,self.Ref)#*1000
        else:
            xout_r=(self.hout_r-hsatL)/(hsatV-hsatL)
            self.sout_r=ssatV*xout_r+(1-xout_r)*ssatL
        
        #Outlet superheat an temperature (in case of two phase)
        if existsSuperheat:
            self.DT_sh_calc=self.Tout_r-self.Tdew_r
        else:
            self.DT_sh_calc=(self.hout_r-hsatV)/(PropsSI('C','T',self.Tdew_r,'Q',1,self.Ref)) #*1000
            self.Tout_r=PropsSI('T','P',self.psat_r+self.DP_r/1.0,'Q',xout_r,self.Ref) #self.DP_r/1000.0 is updated by removing /1000.0 #saturated temperature at outlet quality
        self.hmean_r=self.w_2phase*self.h_r_2phase+self.w_superheat*self.h_r_superheat
        self.UA_r=self.hmean_r*self.A_r_wetted
        self.UA_a=self.Fins.h_a*self.Fins.A_a*self.Fins.eta_a
        
        #Build a vector of temperatures at each point where there is a phase transition along the averaged circuit
        if existsSuperheat:
            #Insert the shoulder point
            Tv=[self.Tin_r,self.Tdew_r,self.Tout_r]
            x=[0,self.w_2phase,1]
        else:
            Tv=[self.Tin_r,xout_r*self.Tdew_r+(1-xout_r)*self.Tbubble_r]
            x=[0,1]
        
        #Determine each bend temperature by interpolation
        #------------------------------------------------
        #Number of bends (including inlet and outlet of coil)
        Nbends=1+self.Lcircuit/self.Ltube
        #x-position of each point
        xv=np.linspace(0,1,Nbends)
        
        self.Tbends=interp1d(x,Tv)(xv)
        
        
        
    def _TwoPhase_Forward(self,w_2phase):
    
        DWS=DWSVals() #DryWetSegment structure
    
        # Store temporary values to be passed to DryWetSegment
        DWS.Fins=self.Fins
        DWS.FinsType = self.FinsType                                            #Added to pass FinsType to DryWetSegment
        DWS.A_a=self.Fins.A_a*w_2phase
        DWS.cp_da=self.Fins.cp_da
        DWS.eta_a=self.Fins.eta_a
        DWS.h_a=self.Fins.h_a  #Heat transfer coefficient, not enthalpy
        DWS.mdot_da=self.mdot_da*w_2phase
        DWS.pin_a=self.Fins.Air.p
        DWS.Tdew_r=self.Tdew_r
        DWS.Tbubble_r=self.Tbubble_r

        DWS.Tin_a=self.Tin_a
        DWS.RHin_a=self.Fins.Air.RH
    
        DWS.Tin_r=self.Tsat_r
        DWS.A_r=self.A_r_wetted*w_2phase
        DWS.cp_r=1.0e15 #In the two-phase region the cp is infinite, use 1e15 as a big number;
        DWS.pin_r=self.psat_r
        DWS.mdot_r=self.mdot_r
        DWS.IsTwoPhase=True
        
        #Target heat transfer to go from inlet quality to saturated vapor
        Q_target=self.mdot_r*(self.xout_2phase-self.xin_r)*self.h_fg
        
        if Q_target<0:
            raise ValueError('Q_target in Evaporator must be positive')
        
        # Average Refrigerant heat transfer coefficient
        DWS.h_r=ShahEvaporation_Average(self.xin_r,self.xout_2phase,self.Ref,self.G_r,self.ID,self.psat_r,Q_target/DWS.A_r,self.Tbubble_r,self.Tdew_r)
        
        #Run the DryWetSegment to carry out the heat and mass transfer analysis
        DryWetSegment(DWS)
        
        self.Q_2phase=DWS.Q
        self.Q_sensible_2phase=DWS.Q_sensible
        self.h_r_2phase=DWS.h_r
        self.fdry_2phase=DWS.f_dry
        self.Tout_a_2phase=DWS.Tout_a
        
        rho_average=TwoPhaseDensity(self.Ref,self.xin_r,self.xout_2phase,self.Tdew_r,self.Tbubble_r,slipModel='Zivi')
        self.Charge_2phase = rho_average * w_2phase * self.V_r        
        
        #Frictional pressure drop component
        DP_frict=LMPressureGradientAvg(self.xin_r,self.xout_2phase,self.Ref,self.G_r,self.ID,self.Tbubble_r,self.Tdew_r)*self.Lcircuit*w_2phase
        #Accelerational pressure drop component    
        DP_accel=AccelPressureDrop(self.xin_r,self.xout_2phase,self.Ref,self.G_r,self.Tbubble_r,self.Tdew_r)*self.Lcircuit*w_2phase
        self.DP_r_2phase=DP_frict+DP_accel;
        
        if self.Verbosity>7:
            print w_2phase,DWS.Q,Q_target,self.xin_r,"w_2phase,DWS.Q,Q_target,self.xin_r"
        return DWS.Q-Q_target
    
    def _Superheat_Forward(self,w_superheat):
        self.w_superheat=w_superheat
        DWS=DWSVals() #DryWetSegment structure
    
        # Store temporary values to be passed to DryWetSegment
        DWS.A_a=self.Fins.A_a*w_superheat
        DWS.cp_da=self.Fins.cp_da
        DWS.eta_a=self.Fins.eta_a
        DWS.h_a=self.Fins.h_a  #Heat transfer coefficient
        DWS.mdot_da=self.mdot_da*w_superheat
        DWS.pin_a=self.Fins.Air.p
        DWS.Fins=self.Fins
        DWS.FinsType = self.FinsType           #Added to pass FinsType to DryWetSegment
    
        # Inputs on the air side to two phase region are inlet air again
        DWS.Tin_a=self.Tin_a
        DWS.RHin_a=self.Fins.Air.RH
    
        DWS.Tin_r=self.Tdew_r
        DWS.A_r=self.A_r_wetted*w_superheat
        DWS.cp_r=PropsSI('C','T',self.Tdew_r+2.5, 'P', self.psat_r, self.Ref)#*1000 #Use a guess value of 6K superheat to calculate cp 
        DWS.pin_r=self.psat_r
        DWS.mdot_r=self.mdot_r
        DWS.IsTwoPhase=False
        
        #Use a guess value of 6K superheat to calculate the properties
        self.f_r_superheat, self.h_r_superheat, self.Re_r_superheat=f_h_1phase_Tube(self.mdot_r / self.Ncircuits, self.ID, 
            self.Tdew_r+3, self.psat_r, self.Ref, "Single");
        
        # Average Refrigerant heat transfer coefficient
        DWS.h_r=self.h_r_superheat
        
        #Run DryWetSegment
        DryWetSegment(DWS)
        
        rho_superheat=PropsSI('D','T',(DWS.Tout_r+self.Tdew_r)/2.0, 'P', self.psat_r, self.Ref)
        self.Charge_superheat = w_superheat * self.V_r * rho_superheat
        
        #Pressure drop calculations for superheated refrigerant
        v_r=1/rho_superheat
        #Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_superheat*v_r*self.G_r**2/(2*self.ID)  #Pressure gradient
        self.DP_r_superheat=dpdz_r*self.Lcircuit*self.w_superheat
        
        #Set values
        self.Q_superheat=DWS.Q
        self.Q_sensible_superheat=DWS.Q_sensible
        self.fdry_superheat=DWS.f_dry
        self.Tout_a_superheat=DWS.Tout_a
        self.Tout_r=DWS.Tout_r
    
if __name__=='__main__':
    #This code runs if this file is run by itself, but otherwise doesn't run
    TT=[]
    Q=[]
    ff=[]
    QQ=[]
    hh=[]
    import numpy as np,pylab
    for Tdew in np.linspace(270,299.7,101):  
        FinsTubes=FinInputs()
    
        FinsTubes.Tubes.NTubes_per_bank=32
        FinsTubes.Tubes.Ncircuits=5
        FinsTubes.Tubes.Nbank=3
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
        FinsTubes.Air.Tdb=299.8
        FinsTubes.Air.p=101325                                                  #updated from 101.325kPa to 101325Pa
        FinsTubes.Air.RH=0.51
        FinsTubes.Air.RHmean=0.51
        FinsTubes.Air.FanPower=438
            
        kwargs={'Ref': 'R410A',
                'mdot_r':  0.0708,
                'psat_r':  PropsSI('P','T',Tdew,'Q',1.0,'R410A'),
                'Fins': FinsTubes,
                'FinsType': 'WavyLouveredFins',                                  #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
                'hin_r': PropsSI('H','T',Tdew,'Q',0.15,'R410A'), #*1000
                'Verbosity': 0
        }
        
        Evap=EvaporatorClass(**kwargs)
        Evap.Update(**kwargs) #update not necessary here, but left for illu
        Evap.Calculate()
        
        #print Evap.OutputList()
        Q.append(Evap.Q)
        QQ.append(Evap.Q_2phase)
        TT.append(Tdew)
        ff.append(Evap.w_2phase)
        hh.append(Evap.h_r_2phase)
        
    pylab.plot(TT,QQ,TT,Q)
    pylab.show()
    pylab.plot(TT,hh)
    pylab.show()
    
        
    