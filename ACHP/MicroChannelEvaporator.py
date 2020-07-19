from __future__ import division, print_function, absolute_import
from math import pi,log,exp
from scipy.optimize import brentq #solver to find roots (zero points) of functions
from scipy.interpolate import interp1d
import numpy as np

import CoolProp as CP

from ACHP.Correlations import f_h_1phase_MicroTube,ShahEvaporation_Average,LMPressureGradientAvg,AccelPressureDrop,TwoPhaseDensity,Bertsch_MC_Average,KM_Evap_Average,KandlikarEvaporation_average
from ACHP.MicroFinCorrelations import MultiLouveredMicroFins, MicroFinInputs, IsFinsClass
from ACHP.DryWetSegment import DWSVals, DryWetSegment
from ACHP.ACHPTools import ValidateFields


class MicroChannelEvaporatorClass():
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
            ('Outlet superheat','K',self.Tout_r-self.Tdew_r),
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
            ('Mean Air HTC','W/m^2-K',self.Fins.h_a*self.h_a_tuning),
            ('Surface Effectiveness','-',self.Fins.eta_a),
            ('Air-side area (fin+tubes)','m^2',self.Fins.A_a),
            ('Mass Flow rate dry Air','kg/s',self.Fins.mdot_da),
            ('Mass Flow rate humid Air','kg/s',self.Fins.mdot_ha),
            ('Pressure Drop Air-side','Pa',self.Fins.dP_a),
            ('Sensible Heat Ratio','-',self.SHR),
            ('Bend Temperature profile','K',self.Tbends)
        ]
        for i in range(0,len(Output_List_default)):         #append default parameters to output list
            Output_List.append(Output_List_default[i])
        return Output_List
        
    def AirSideCalcs(self):
        #Update with user FinType
        if self.FinsType == 'MultiLouveredMicroFins':
            MultiLouveredMicroFins(self.Fins)
    
    def Initialize(self):
        #Input validation the first call of Initialize
        if False:#not hasattr(self,'IsValidated'):
            self.Fins.Validate()
            reqFields=[
                       ('psat_r',float,0.001,100000000),                        
                       ('Fins',IsFinsClass,None,None),
                       ('FinsType',str,None,None),
                       ('hin_r',float,-100000,10000000),
                       ('mdot_r',float,0.000001,10),
                       ]
            optFields=['Verbosity','AS','h_a_tuning','h_tp_tuning','DP_tuning']
            d=self.__dict__ #Current fields in model
            ValidateFields(d,reqFields,optFields)
            self.IsValidated=True
        
        #set tuning factors to 1 in case not given by user
        if not hasattr(self,'h_a_tuning'):
            self.h_a_tuning = 1
        if not hasattr(self,'h_tp_tuning'):
            self.h_tp_tuning = 1
        if not hasattr(self,'DP_tuning'):
            self.DP_tuning = 1
            
        # Make sure AbstractState is passed
        assert hasattr(self,'AS'), 'Please specify the Abstract State'
            
        # Retrieve some parameters from nested structures 
        # for code compactness
        self.Ltube=self.Fins.Tubes.Ltube    #tube length
        self.NTubes=self.Fins.Tubes.NTubes  #number of tube (per bank)
        self.Nbank=self.Fins.Tubes.Nbank    #number of banks
        self.Tin_a=self.Fins.Air.Tdb        #inlet air temperature
        self.Pin_a=self.Fins.Air.p         #inlet air pressure
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
        
        # Ratio of PH over PF where PH: heated perimeter of channel, PF: wetted perimeter of channel
        # In this application, heated perimeter = wetted perimeter, times the number of ports
        self.PH_PF = 1*self.Nports
        
        ## Bubble and dew temperatures (same for fluids without glide)
        self.AS.update(CP.PQ_INPUTS, self.psat_r, 0.0)
        self.Tbubble_r=self.AS.T() #[K]
        h_l = self.AS.hmass() #[J/kg]
        self.AS.update(CP.PQ_INPUTS, self.psat_r, 1.0)
        self.Tdew_r=self.AS.T() #[K]
        h_v = self.AS.hmass() #[J/kg]
        ## Mean temperature for use in HT relationships
        self.Tsat_r=(self.Tbubble_r+self.Tdew_r)/2
        # Latent heat
        self.h_fg=h_v - h_l #[J/kg]
        
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
        
        # Input and output thermodynamic properties
        AS.update(CP.QT_INPUTS, 0.0, self.Tbubble_r)
        ssatL=AS.smass() #[J/kg-K]
        hsatL=AS.hmass() #[J/kg]
        AS.update(CP.QT_INPUTS, 1.0, self.Tdew_r)
        ssatV=AS.smass() #[J/kg-K]
        hsatV=AS.hmass() #[J/kg]
        
        #if give enthalpy and pressure as inputs
        if hasattr(self,'hin_r'):
            self.xin_r=(self.hin_r-hsatL)/(hsatV-hsatL)
            self.sin_r=self.xin_r*ssatV+(1-self.xin_r)*ssatL
            self.Tin_r=self.xin_r*self.Tdew_r+(1-self.xin_r)*self.Tbubble_r
        elif hasattr(self,'xin_r'): #if given quality and pressure as inputs
            self.hin_r=self.xin_r*hsatV+(1-self.xin_r)*hsatL
            self.sin_r=self.xin_r*ssatV+(1-self.xin_r)*ssatL
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
            print (self.Q,"MicroChannelEvaporator.Q")
        self.Capacity=self.Q-self.Fins.Air.FanPower
        
        #Sensible heat ratio [-]
        self.SHR=(self.Q_sensible_2phase+self.Q_sensible_superheat)/self.Q
        #Average air outlet temperature (area fraction weighted average) [K]
        self.Tout_a=self.w_superheat*self.Tout_a_superheat+self.w_2phase*self.Tout_a_2phase
        self.DP_r=(self.DP_r_superheat+self.DP_r_2phase)*self.DP_tuning #correcting the total pressure drop
        
        #Outlet enthalpy obtained from energy balance
        self.hout_r=self.hin_r+self.Q/self.mdot_r
        
        #Outlet entropy
        if existsSuperheat==True:
            AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r)
            self.sout_r=AS.smass() #[J/kg-K]
        else:
            xout_r=(self.hout_r-hsatL)/(hsatV-hsatL)
            self.sout_r=ssatV*xout_r+(1-xout_r)*ssatL
        
        #Outlet superheat an temperature (in case of two phase)
        if existsSuperheat:
            self.DT_sh_calc=self.Tout_r-self.Tdew_r
        else:
            AS.update(CP.QT_INPUTS, 1.0, self.Tdew_r)
            cp_sh = AS.cpmass() #[J/kg-K]
            self.DT_sh_calc=(self.hout_r-hsatV)/cp_sh #Effective superheat
            AS.update(CP.PQ_INPUTS, self.psat_r+self.DP_r, xout_r)
            self.Tout_r=AS.T() #saturated temperature at outlet quality [K]
        self.hmean_r=self.w_2phase*self.h_r_2phase+self.w_superheat*self.h_r_superheat
        self.UA_r=self.hmean_r*self.A_r_wetted
        self.UA_a=(self.Fins.h_a*self.h_a_tuning)*self.Fins.A_a*self.Fins.eta_a
        self.UA_w=1/self.Rw
        
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
        DWS.FinsType = self.FinsType                                            
        DWS.A_a=self.Fins.A_a*w_2phase
        DWS.cp_da=self.Fins.cp_da
        DWS.eta_a=self.Fins.eta_a
        DWS.h_a=self.Fins.h_a*self.h_a_tuning  #Heat transfer coefficient, not enthalpy
        DWS.mdot_da=self.mdot_da*w_2phase
        DWS.pin_a=self.Fins.Air.p
        DWS.Tdew_r=self.Tdew_r
        DWS.Tbubble_r=self.Tbubble_r

        DWS.Tin_a=self.Tin_a
        DWS.RHin_a=self.Fins.Air.RH
    
        DWS.Tin_r=self.Tsat_r
        DWS.A_r=self.A_r_wetted*w_2phase
        DWS.Rw=self.Rw/w_2phase
        DWS.cp_r=1.0e15 #In the two-phase region the cp is infinite, use 1e15 as a big number;
        DWS.pin_r=self.psat_r
        DWS.mdot_r=self.mdot_r
        DWS.IsTwoPhase=True
        
        #Target heat transfer to go from inlet quality to saturated vapor
        Q_target=self.mdot_r*(self.xout_2phase-self.xin_r)*self.h_fg
        
        if Q_target<0:
            raise ValueError('Q_target in Evaporator must be positive')
        
        # Average Refrigerant heat transfer coefficient
        #Shah correlation
        #h_r=ShahEvaporation_Average(self.xin_r,self.xout_2phase,self.AS,self.G_r,self.Dh,self.psat_r,Q_target/DWS.A_r,self.Tbubble_r,self.Tdew_r)
        #Bertsch correlation
        #h_r=Bertsch_MC_Average(self.xin_r,self.xout_2phase,self.AS,self.G_r,self.Dh,Q_target/DWS.A_r,self.Lcircuit,self.Tbubble_r,self.Tdew_r)
        #Kandlikar correlation
        #h_r=KandlikarEvaporation_average(self.xin_r,self.xout_2phase,self.AS,self.G_r,self.Dh,self.psat_r,Q_target/DWS.A_r,self.Tbubble_r,self.Tdew_r)
        #Kim and Mudawar (2013)
        DPDZ_frict_2phase, h_r = KM_Evap_Average(self.xin_r,self.xout_2phase,self.AS,self.G_r,self.Dh,self.Tbubble_r,self.Tdew_r,self.psat_r,self.beta,Q_target/DWS.A_r,self.PH_PF)
        DWS.h_r=h_r*self.h_tp_tuning #correct refrigerant side convection heat transfer
        
        #Run the DryWetSegment to carry out the heat and mass transfer analysis
        DryWetSegment(DWS)
        
        self.Q_2phase=DWS.Q
        self.Q_sensible_2phase=DWS.Q_sensible
        self.h_r_2phase=DWS.h_r
        self.fdry_2phase=DWS.f_dry
        self.Tout_a_2phase=DWS.Tout_a
        
        rho_average=TwoPhaseDensity(self.AS,self.xin_r,self.xout_2phase,self.Tdew_r,self.Tbubble_r,slipModel='Zivi')
        self.Charge_2phase = rho_average * w_2phase * self.V_r        
        
        #Frictional pressure drop component
        #Lockhart-Martinelli correlation
        #DP_frict=LMPressureGradientAvg(self.xin_r,self.xout_2phase,self.AS,self.G_r,self.Dh,self.Tbubble_r,self.Tdew_r)*self.Lcircuit*w_2phase
        #using the pressure gradient of Kim & Mudawar (2013) 
        DP_frict=DPDZ_frict_2phase*self.Lcircuit*w_2phase
        
        #Accelerational pressure drop component    
        if self.Ref in 'CarbonDioxide' or 'CO2' or 'R744': #TODO: self.AS.name() is not available yet for all backends, using self.Ref for now
            DP_accel=AccelPressureDrop(self.xin_r,self.xout_2phase,self.AS,self.G_r,self.Tbubble_r,self.Tdew_r,D=self.Dh,slipModel='Premoli')*self.Lcircuit*w_2phase
        else:
            DP_accel=AccelPressureDrop(self.xin_r,self.xout_2phase,self.AS,self.G_r,self.Tbubble_r,self.Tdew_r,slipModel='Zivi')*self.Lcircuit*w_2phase
        self.DP_r_2phase=DP_frict+DP_accel;
        
        if self.Verbosity>7:
            print (w_2phase,DWS.Q,Q_target,self.xin_r,"w_2phase,DWS.Q,Q_target,self.xin_r")
        return DWS.Q-Q_target
    
    def _Superheat_Forward(self,w_superheat):
        self.w_superheat=w_superheat
        DWS=DWSVals() #DryWetSegment structure
        AS=self.AS #AbstractState
    
        # Store temporary values to be passed to DryWetSegment
        DWS.A_a=self.Fins.A_a*w_superheat
        DWS.cp_da=self.Fins.cp_da
        DWS.eta_a=self.Fins.eta_a
        DWS.h_a=self.Fins.h_a*self.h_a_tuning  #Heat transfer coefficient
        DWS.mdot_da=self.mdot_da*w_superheat
        DWS.pin_a=self.Fins.Air.p
        DWS.Fins=self.Fins
        DWS.FinsType = self.FinsType           
        
        # Inputs on the air side to two phase region are inlet air again
        DWS.Tin_a=self.Tin_a
        DWS.RHin_a=self.Fins.Air.RH
    
        DWS.Tin_r=self.Tdew_r
        DWS.A_r=self.A_r_wetted*w_superheat
        DWS.Rw=self.Rw/w_superheat
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tdew_r+2.5)
        DWS.cp_r=AS.cpmass() #Use a guess value of 6K superheat to calculate cp [J/kg-K] 
        DWS.pin_r=self.psat_r
        DWS.mdot_r=self.mdot_r
        DWS.IsTwoPhase=False
        
        #Use a guess value of 6K superheat to calculate the properties
        self.f_r_superheat, self.h_r_superheat, self.Re_r_superheat=f_h_1phase_MicroTube(self.G_r, self.Dh, 
            self.Tdew_r+3, self.psat_r, self.AS, "Single");

        # Average Refrigerant heat transfer coefficient
        DWS.h_r=self.h_r_superheat
        
        #Run DryWetSegment
        DryWetSegment(DWS)
        
        AS.update(CP.PT_INPUTS, self.psat_r, (DWS.Tout_r+self.Tdew_r)/2.0)
        rho_superheat=AS.rhomass() #[kg/m^3]
        self.Charge_superheat = w_superheat * self.V_r * rho_superheat
        
        #Pressure drop calculations for superheated refrigerant
        v_r=1/rho_superheat
        #Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_superheat*v_r*self.G_r**2/(2*self.Dh)  #Pressure gradient
        self.DP_r_superheat=dpdz_r*self.Lcircuit*self.w_superheat
        
        #Set values
        self.Q_superheat=DWS.Q
        self.Q_sensible_superheat=DWS.Q_sensible
        self.fdry_superheat=DWS.f_dry
        self.Tout_a_superheat=DWS.Tout_a
        self.Tout_r=DWS.Tout_r
    
if __name__=='__main__':
    #Example usage for a parametric study
    from CoolProp.CoolProp import PropsSI
    import pylab
    
    num_points= 101
    T_dews= np.linspace(265,292,num_points)
    TT= np.empty(num_points)
    Q_2p= np.empty(num_points)
    w_2p= np.empty(num_points)
    w_sh= np.empty(num_points)
    Q_tot= np.empty(num_points)
    h_2p= np.empty(num_points)
    h_sh= np.empty(num_points)

    FinsTubes=MicroFinInputs()
    
    FinsTubes.Tubes.NTubes=61.354           #Number of tubes (per bank for now!)
    FinsTubes.Tubes.Nbank=1                 #Number of banks (set to 1 for now!)
    FinsTubes.Tubes.Npass=3                 #Number of passes (per bank-averaged)
    FinsTubes.Tubes.Nports=1                #Number of rectangular ports
    FinsTubes.Tubes.Ltube=0.30213           #length of a single tube
    FinsTubes.Tubes.Td=0.0333               #Tube outside width (depth)
    FinsTubes.Tubes.Ht= 0.002               #Tube outside height (major diameter)
    FinsTubes.Tubes.b=0.00635               #Tube spacing     
    FinsTubes.Tubes.tw=0.0003               #Tube wall thickness     
    FinsTubes.Tubes.twp=0.0003              #Port (channel) wall thickness     
    FinsTubes.Tubes.beta=1                  #Port (channel) aspect ratio (=width/height)
    FinsTubes.Tubes.kw=117                  #wall thermal conductivity
        
    FinsTubes.Fins.FPI=11.0998              #Fin per inch
    FinsTubes.Fins.Lf=0.0333                #Fin length
    FinsTubes.Fins.t=0.000152               #Fin thickness
    FinsTubes.Fins.k_fin=117                #Fin thermal conductivity
    
    FinsTubes.Air.Vdot_ha=0.5663            #Air volume flow rate in m^3/s
    FinsTubes.Air.Tmean=299.9 
    FinsTubes.Air.Tdb=299.9                   #Air inlet temperature, K
    FinsTubes.Air.p=101325                  #Air pressure in Pa
    FinsTubes.Air.RHmean=0.51
    FinsTubes.Air.RH=0.51                   #Air inlet relative humidity
    FinsTubes.Air.FanPower=438              #Fan power, Watts
    
    FinsTubes.Louvers.Lalpha=20             #Louver angle, in degree
    FinsTubes.Louvers.lp=0.001              #Louver pitch
    FinsTubes.Louvers.Llouv=0.005737        #Louver cut length
    
    #Abstract State
    Ref = 'R744'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref)
    
    kwargs={'AS': AS,
            'Ref':'R744', 
            'mdot_r':  0.0708,
            'psat_r':  PropsSI('P','T',T_dews[0],'Q',1.0,Ref),
            'Fins': FinsTubes,
            'FinsType': 'MultiLouveredMicroFins',
            'hin_r': PropsSI('H','P',PropsSI('P','T',282,'Q',1.0,Ref),'Q',0.15,Ref), 
            'Verbosity': 0,
            'h_a_tuning':1,
            'h_tp_tuning':1,
            'DP_tuning':1,
        }
        
    MicroEvap=MicroChannelEvaporatorClass(**kwargs) #generate new micro-channel evaporator instance and update kwargs
    
    for i in range(0, len(T_dews)):
        kwargs={'psat_r':  PropsSI('P','T',T_dews[i],'Q',1.0,Ref)}
        MicroEvap.Update(**kwargs)
        MicroEvap.Calculate()
        Q_tot[i] = MicroEvap.Q
        Q_2p[i]= MicroEvap.Q_2phase
        w_2p[i]= MicroEvap.w_2phase
        w_sh[i]= MicroEvap.w_superheat
        h_2p[i]= MicroEvap.h_r_2phase
        h_sh[i]= MicroEvap.h_r_superheat

    print ("Demonstrate output list")
    #print (Evap.OutputList())
    for id, unit, value in MicroEvap.OutputList():                
        print (str(id) + ' = ' + str(value) + ' ' + str(unit))

    pylab.plot(T_dews,Q_2p,T_dews,Q_tot)
    pylab.title('Parametric Study With Fixed flowrates - Capacity')
    pylab.legend(['two-phase','total'],loc='best')
    pylab.title('Parametric Study With Fixed flowrates - Capacity')
    pylab.xlabel('Evaporation Dew Temperature in Kelvin')
    pylab.ylabel('Capacity in Watt')
    #pylab.savefig('Evaporator_py_capacity.pdf')
    pylab.show()
    pylab.plot(T_dews,h_2p,T_dews, h_sh)
    pylab.title('Parametric Study with fixed Flowrates - Heat Transfer Coefficients')
    pylab.legend(['two-phase','superheat'],loc='best')
    pylab.xlabel('Evaporation Dew Temperature in Kelvin')
    pylab.ylabel('Heat Transfer Coefficient in W/m2-K')
    #pylab.savefig('Evaporator_py_HTC.pdf')
    pylab.show()
    pylab.plot(T_dews,w_2p, T_dews, w_sh)
    pylab.title('Parametric Study with fixed Flowrates - Area Fraction')
    pylab.legend(['two-phase', 'superheat'],loc='best')
    pylab.xlabel('Evaporation Dew Temperature in Kelvin')
    pylab.ylabel('Two-phase Wetted Area Fraction')
    pylab.ylim(-0.01,1.01)
    #pylab.savefig('Evaporator_py_wetted_area.pdf')
    pylab.show()