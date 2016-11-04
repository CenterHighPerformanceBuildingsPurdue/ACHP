from __future__ import division #Make integer 3/2 give 1.5 in python 2.x
from math import pi,log,exp
from CoolProp.CoolProp import Props#,UseSaturationLUT
from CoolProp.HumidAirProp import HAProps
from Correlations import f_h_1phase_Tube,ShahEvaporation_Average, LockhartMartinelli,LMPressureGradientAvg,AccelPressureDrop,TwoPhaseDensity
from scipy.optimize import brentq #solver to find roots (zero points) of functions
from scipy.interpolate import interp1d
#import numpy as np
from FinCorrelations import WavyLouveredFins,FinInputs,IsFinsClass
from DryWetSegment import DWSVals, DryWetSegment
from ACHPTools import ValidateFields
import numpy as np
from scipy.optimize import newton, fsolve
# Turn on saturation curve lookup for CoolProp
#UseSaturationLUT(1)

class EvaporatorClass():
    def __init__(self,**kwargs):
        self.h_tp_tuning=1.0
        self.cp_r_iter=False   #use iteration to find correct value of CP? -> default if not updated using kwargs
        self.__dict__.update(kwargs)
        print "h_evap_tuning is manually set for 2-phase in Evaporator. Correct as necessary!"
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
        Output_List_default=[                                                                           #default output list
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
            ('Outlet superheat','K',self.Tout_r-self.Tdew_r),
            ('Fins per inch','1/in',self.Fins.Fins.FPI),
            ('Fin waviness pd','m',self.Fins.Fins.Pd),
            ('Fin waviness xf','m',self.Fins.Fins.xf),
            ('Fin thickness','m',self.Fins.Fins.t),
            ('Fin Conductivity','W/m-K',self.Fins.Fins.k_fin),
            ('Q Total','W',self.Q),
            ('Q Superheat','W',self.Q_superheat),
            ('Q Two-Phase','W',self.Q_2phase),
            ('Inlet ref. temp','K',self.Tin_r),
            ('Outlet ref. temp','K',self.Tout_r),
            ('Outlet air temp','K',self.Tout_a),
            ('Evaporator P_sat in','kPa',self.psat_r),
            ('Evaporator inlet quality','-',self.xin_r),
            ('Evaporator ref. flowrate','g/s',self.mdot_r*1000.0),
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
            ('Sensible Heat Ratio','-',self.SHR)
        ]
        for i in range(0,len(Output_List_default)):                             #append default parameters to output list
            Output_List.append(Output_List_default[i])
        return Output_List
        
    def AirSideCalcs(self):
        self.Fins.h_a_tuning= self.h_a_tuning #pass tuning factor
        WavyLouveredFins(self.Fins)
    
    def Initialize(self):
        #Input validation the first call of Initialize
        if False:#not hasattr(self,'IsValidated'):
            self.Fins.Validate()
            reqFields=[
                       ('Ref',str,None,None),
                       ('psat_r',float,1e-6,100000),
                       ('Fins',IsFinsClass,None,None),
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
        self.Tbubble_r=Props('T','P',self.psat_r,'Q',0,self.Ref)
        self.Tdew_r=Props('T','P',self.psat_r,'Q',1,self.Ref)
        ## Mean temperature for use in HT relationships
        self.Tsat_r=(self.Tbubble_r+self.Tdew_r)/2
        # Latent heat
        self.h_fg=(Props('H','T',self.Tdew_r,'Q',1.0,self.Ref)-Props('H','T',self.Tbubble_r,'Q',0.0,self.Ref))*1000. #[J/kg]
        
        self.Fins.Air.RHmean=self.Fins.Air.RH
        self.Fins.h_tp_tuning= self.h_tp_tuning #pass tuning factor
        WavyLouveredFins(self.Fins)
        self.mdot_ha=self.Fins.mdot_ha #[kg_ha/s]
        self.mdot_da=self.Fins.mdot_da #[kg_da/s]

    def Calculate_PD(self):
        printDP_iteration=False
        printDP=True
        #calculate evaporator with consideration of pressure drop
        self.pin_r=self.psat_r*1.0  #store saturation pressure (needed for iteration)
        
        #first treat 2-phase part
        def objective_evaporator2p(DP_guess):
            #wrapper for Evap.Calculate() to determine real 2-phase pressure drop
            #input is estimate of two phase pressure drop
            #DP_guess=max(-0.1,DP_guess)  #sometimes fsolve uses positive pressure drop; prevent from crashing Calculate()
            #DP_guess=DP_guess[0]
            self.psat_r=self.pin_r+DP_guess/2.0  #assume pressure drop in evaporator is linear
            self.Calculate()
            if printDP_iteration: print "x2 guess",DP_guess,
            return self.DP_r_2phase/1000.0+DP_guess
        self.Calculate()  #calculate normal evaporator to obtain guess for pressure drop
        DP_guess=self.DP_r_2phase/1000.0
        from scipy.optimize import fminbound
        lower_fminbound=max(DP_guess*1.5,50-self.pin_r)  #ensure that no negative pressures occur during solving process
        if printDP: print "DP_2Pguess is",DP_guess," with 2-phase fraction ",self.w_2phase,"bounds",lower_fminbound,(DP_guess*1.5,50-self.pin_r),self.mdot_r
        fminbound(objective_evaporator2p,lower_fminbound,0.0)
        #fsolve(objective_evaporator2p, DP_guess)
        #brentq(objective_evaporator2p,DP_guess,0.0)
        if printDP: print "DP_2Pacrtual is",self.DP_r_2phase/1000.0," with 2-phase fraction ",self.w_2phase
        
        #second treat superheated section, if applicable
        if self.existsSuperheat==True:
            def objective_evaporator1p(DP_guess):
                #function for solving heat transfer of single phase section correctly
                if printDP_iteration: print "x1",DP_guess,
                
                self.psat_r=self.pin_r+self.DP_r_2phase/1000.0+DP_guess/2.0   #assuming linear profile of pressure drop
    
                #need to repeat some part of the initialize functions
                ## Bubble and dew temperatures (same for fluids without glide) 
                self.Tbubble_r=Props('T','P',self.psat_r,'Q',0,self.Ref)
                self.Tdew_r=Props('T','P',self.psat_r,'Q',1,self.Ref)
                ## Mean temperature for use in HT relationships
                self.Tsat_r=(self.Tbubble_r+self.Tdew_r)/2
                # Latent heat
                self.h_fg=(Props('H','T',self.Tdew_r,'Q',1.0,self.Ref)-Props('H','T',self.Tbubble_r,'Q',0.0,self.Ref))*1000. #[J/kg]
                
                #need to repeat part of the calculate function; already know length of two phase section...
                if self.w_2phase>=0.0:  #already know outlet quality from finding the 2-phase pressure drop
                    T_inlet_r=self.Tdew_r #saturated inlet to sensible section
                else: #have superheated conditions, two phase section doesn't exist
                    self.Q_2phase=0
                    self.Charge_2phase=0
                    self.Q_sensible_2phase=0
                    self.h_r_2phase=0.0
                    self.DP_r_2phase=0.0
                    self.xout_2phase=99  #EES-style
                    self.Tout_a_2phase=self.Tin_a  #since we have no two-phase capacity
                    T_inlet_r=Props('T','H',self.hin_r/1000.0, 'P', self.psat_r, self.Ref) #calculated inlet temperature based on enthalpy
                    if self.Tin_a<T_inlet_r: 
                        print "Evap: inlet air temperature smaller than refrigerant temperature",self.Tin_a,"<",T_inlet_r,"self.hin_r/1000.0",self.hin_r/1000.0,"self.psat_r",self.psat_r
                        raise ValueError
                if self.cp_r_iter: #solve for proper cp
                    #print "using iteration to find better estimate for cp_r"
                    try:
                        def objective_cp(x0):
                            self.cp_r=1000.0*(Props('C','T',self.Tdew_r+x0*0.3333, 'P', self.psat_r, self.Ref)+Props('C','T',self.Tdew_r+x0, 'P', self.psat_r, self.Ref)+Props('C','T',self.Tdew_r+1.66667*x0, 'P', self.psat_r, self.Ref))/3.0  #average value
                            self._Superheat_Forward(1-self.w_2phase)
                            return self.Tdew_r+2*x0-self.Tout_r
                        brentq(objective_cp,0,50)
                    except:
                            #if above approach does not work, do not iterate.
                            self.cp_r=1000.0*Props('C','T',self.Tdew_r+3, 'P', self.psat_r, self.Ref)
                            print "warning-cp_r iter in Evaporator did not work. using non-iteratively solved value for superheat instead", self.cp_r,self.w_2phase
                            self._Superheat_Forward(1-self.w_2phase,T_inlet_r)
                else: #use fixed point cp assuming 6K superheat
                    self.cp_r=1000.0*Props('C','T',self.Tdew_r+3, 'P', self.psat_r, self.Ref)
                    self._Superheat_Forward(1-self.w_2phase,T_inlet_r)
                return DP_guess-self.DP_r_superheat/1000.0
            #####actually solve objective function
            DP_guess=self.DP_r_superheat/1000.0
            lower_fminbound=max(DP_guess*1.1,50-self.pin_r-self.DP_r_2phase/1000.0)  #ensure that no negative pressures occur during solving process
            if printDP: print "DP_1Pguess is",DP_guess, lower_fminbound,(DP_guess*1.5,10-self.pin_r-self.DP_r_2phase/1000.0)
            try:
                fminbound(objective_evaporator1p,lower_fminbound,0.0)
            except:
                print "fminbound failed, trying fsolve instead"
                fsolve(objective_evaporator1p, DP_guess/3.0)
            if printDP: print "DP_1Pacrtual is",self.DP_r_superheat/1000.0,"resulting in an outlet pressure of",self.pin_r+(self.DP_r_superheat+self.DP_r_2phase)/1000.0
        
        #update values
        self.psat_r=self.pin_r #back to original values
        #need to repeat some part of the initialize functions
        ## Bubble and dew temperatures (same for fluids without glide) 
        self.Tbubble_r=Props('T','P',self.psat_r,'Q',0,self.Ref)
        self.Tdew_r=Props('T','P',self.psat_r+self.DP_r_2phase/1000.0,'Q',1.0,self.Ref)  #consider pressure drop
        ## Mean temperature for use in HT relationships
        self.Tsat_r=(self.Tbubble_r+self.Tdew_r)/2
        # Input and output thermodynamic properties
        ssatL=Props('S','T',self.Tbubble_r,'Q',0,self.Ref)*1000
        ssatV=Props('S','T',self.Tdew_r,'Q',1,self.Ref)*1000
        hsatL=Props('H','T',self.Tbubble_r,'Q',0,self.Ref)*1000
        hsatV=Props('H','T',self.Tdew_r,'Q',1,self.Ref)*1000
        self.h_fg=hsatV-hsatL  #not considered in heat transfer calculations to prevent instability
        
        ##below as in normal Calculate function, but considering pressure drop
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
        self.pout_r=max(self.pin_r+self.DP_r/1000.0,10.0)  #pressure drop negative and in Pa; limit to a positive pressure to avoid problems during solving process
        
        #Outlet enthalpy obtained from energy balance
        self.hout_r=self.hin_r+self.Q/self.mdot_r
        
        #Outlet superheat an temperature (in case of two phase)
        if self.existsSuperheat:  #neglecting change in pressure for superheat calculation...
            try:
                self.Tout_r=newton(lambda T: Props('H','T',T,'P',self.pout_r,self.Ref)-self.hout_r/1000.0,Props('T','P',self.pout_r,'Q',1.0,self.Ref))
                self.sout_r=Props('S','T',self.Tout_r,'P',self.pout_r,self.Ref)*1000.0
                self.DT_sh_calc=self.Tout_r-self.Tdew_r
            except:
                print "Evap.Calculate_PD() self.pout_r",self.pout_r,"self.pin_r",self.pin_r,  "Props('H','Q',1.0,'P',self.pout_r,self.Ref),self.hout_r/1000.0",Props('H','Q',1.0,'P',self.pout_r,self.Ref),self.hout_r/1000.0
                raise()
        else:
            xout_r=(self.hout_r-hsatL)/(hsatV-hsatL)
            self.sout_r=ssatV*xout_r+(1-xout_r)*ssatL
            self.DT_sh_calc=(self.hout_r-hsatV)/(Props('C','T',self.Tdew_r,'Q',1,self.Ref)*1000.0) #continous superheat
            try:  #this statement is necessary, due to inaccuracies in the calculation of the quality
                self.Tout_r=Props('T','P',self.pout_r,'Q',xout_r,self.Ref) #saturated temperature at outlet quality
            except:
                self.Tout_r=Props('T','P',self.pout_r,'Q',1.0,self.Ref) #saturated temperature at outlet quality
                if xout_r>1.0001:
                    print "in Evaporator - outlet quality ",xout_r,"larger than one - is this a mistake?"

        self.hmean_r=self.w_2phase*self.h_r_2phase+self.w_superheat*self.h_r_superheat
        self.UA_r=self.hmean_r*self.A_r_wetted
        self.UA_a=self.Fins.h_a*self.Fins.A_a*self.Fins.eta_a
        
        #Build a vector of temperatures at each point where there is a phase transition along the averaged circuit
            #this is not updated or supported at this point
        
    def Calculate(self):
        #calculate evaporator without consideration of pressure drop for HT
        self.Initialize()
        
        # Input and output thermodynamic properties
        ssatL=Props('S','T',self.Tbubble_r,'Q',0,self.Ref)*1000
        ssatV=Props('S','T',self.Tdew_r,'Q',1,self.Ref)*1000
        hsatL=Props('H','T',self.Tbubble_r,'Q',0,self.Ref)*1000
        hsatV=Props('H','T',self.Tdew_r,'Q',1,self.Ref)*1000
        
        #Must give enthalpy and pressure as inputs
        #print "in evaporator-enthalpies","hin_r,hsatL,hsatV",self.hin_r,hsatL,hsatV,"self.Tbubble_r,self.psat_r",self.Tbubble_r,self.psat_r
        #print "self.__dict__ in Evaporator.py",self.__dict__
        #print "hin_r in Evaporator.py",self.hin_r
        self.xin_r=(self.hin_r-hsatL)/(hsatV-hsatL)
        self.sin_r=self.xin_r*ssatV+(1-self.xin_r)*ssatL
        self.hin_r=self.xin_r*hsatV+(1-self.xin_r)*hsatL
        self.Tin_r=self.xin_r*self.Tdew_r+(1-self.xin_r)*self.Tbubble_r
        
        #Begin by assuming that you go all the way to saturated vapor at least
        self.xout_2phase=1.0
        if self.Tin_a<(self.Tin_r+0.3):
            #too close to saturation temperature - set capacity to 0
            print "!!!!!!!!!!!!!!!!!!!!!warning - setting capacity to 0 and neglecting pressure drop in evaporator!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            if self.xin_r<1.0:
                #since no capacity, everything just goes out as it came in
                self.w_2phase=1.0
                self.w_superheat=0.0
                existsSuperheat=False
            else:
                self.w_2phase=0.0
                self.w_superheat=1.0
                existsSuperheat=True
            self.Q_2phase=0
            self.Charge_2phase=0
            self.Q_sensible_2phase=0
            self.h_r_2phase=0.0
            self.DP_r_2phase=0.0
            self.xout_2phase=99  #EES-style
            self.Tout_a_2phase=self.Tin_a  #since we have no two-phase capacity
            self.omega_out_2phase=0.0
            self.Q_superheat=0.0
            self.h_r_superheat=0.0
            self.Re_r_superheat=0.0
            self.Charge_superheat=0.0
            self.Q_sensible_superheat=0.0
            self.Tout_a_superheat=0.0
            self.DP_r_superheat=0.0
            self.fdry_superheat=0.0
            self.omega_out_superheat=0.0
        elif self.xin_r<1.0 and self._TwoPhase_Forward(1.0)<0:
            #print "assuming two-phase",self.xin_r
            # Evaporator outlet is in two-phase region, use all the area and find the outlet quality
            existsSuperheat=False
            self.w_2phase=1.0
            def OBJECTIVE(xout):
                self.xout_2phase=xout
                Q_target=self.mdot_r*(xout-self.xin_r)*(hsatV-hsatL)
                self._TwoPhase_Forward(self.w_2phase)
                #print "in Evaporator xout_objective",self.xin_r,xout,self.Q_2phase,Q_target
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
            self.omega_out_superheat=0.0
        else:  #have superheated inlet conditions
            existsSuperheat=True
            # Evaporator outlet is in superheated region, everything is ok
            if self.xin_r<1.0:  #need to find outlet quality
                self.w_2phase=brentq(self._TwoPhase_Forward,0.00000000001,0.9999999999)
                if hasattr(self, 'plotit'):
                    if self.plotit==True:
                        #plot how self._TwoPhase_Forward behaves
                        import pylab as plt
                        w_guesses=np.linspace(0.00000000001,0.9999999999,20)
                        w_resids=np.zeros(len(w_guesses))
                        for i in range(len(w_guesses)):
                            w_resids[i]=self._TwoPhase_Forward(w_guesses[i])
                        print "w_resids",w_resids
                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        ax.plot(w_guesses, w_resids)
                        ax.set_title("residual for 2-phase fraction")
                        plt.show()
                T_inlet_r=self.Tdew_r #saturated inlet to sensible section
            else: #have superheated conditions, two phase section doesn't exist
                self.w_2phase=0.0
                self.Q_2phase=0
                self.Charge_2phase=0
                self.Q_sensible_2phase=0
                self.h_r_2phase=0.0
                self.DP_r_2phase=0.0
                self.xout_2phase=99  #EES-style
                self.Tout_a_2phase=self.Tin_a  #since we have no two-phase capacity
                self.omega_out_2phase=0.0
                T_inlet_r=Props('T','H',self.hin_r/1000.0, 'P', self.psat_r, self.Ref) #calculated inlet temperature based on enthalpy
                if self.Tin_a<T_inlet_r: 
                    print "Evap: inlet air temperature smaller than refrigerant temperature, two phase doesn't exist",self.Tin_a,"<",T_inlet_r,"self.hin_r/1000.0",self.hin_r/1000.0,"self.psat_r",self.psat_r
                    raise ValueError
            if self.cp_r_iter: #solve for proper cp
                print "using iteration to find better estimate for cp_r"
                try:
                    def objective_cp(x0):
                        self.cp_r=1000.0*(Props('C','T',self.Tdew_r+x0*0.3333, 'P', self.psat_r, self.Ref)+Props('C','T',self.Tdew_r+x0, 'P', self.psat_r, self.Ref)+Props('C','T',self.Tdew_r+1.66667*x0, 'P', self.psat_r, self.Ref))/3.0  #average value
                        self._Superheat_Forward(1-self.w_2phase)
                        return self.Tdew_r-2*x0-self.Tout_r
                    brentq_result=brentq(objective_cp,0,50)
                    print "accuracy of cp-iteration loop of evaporator -objective_cp(brentq_result)",objective_cp(brentq_result),brentq_result, -(self.Tdew_r-self.Tout_r)/2.0,self.cp_r
                except:
                        #if above approach does not work, do not iterate.
                        self.cp_r=1000.0*Props('C','T',self.Tdew_r+3, 'P', self.psat_r, self.Ref)
                        print "warning-cp_r iter in Evaporator did not work. using non-iteratively solved value for superheat instead", self.cp_r,self.w_2phase
                        self._Superheat_Forward(1-self.w_2phase,T_inlet_r)
            else: #use fixed point cp assuming 6K superheat
                self.cp_r=1000.0*Props('C','T',self.Tdew_r+3, 'P', self.psat_r, self.Ref)
                self._Superheat_Forward(1-self.w_2phase,T_inlet_r)
            
        self.Q=self.Q_superheat+self.Q_2phase
        self.Charge=self.Charge_superheat+self.Charge_2phase
        if self.Verbosity>4: 
            print self.Q,"Evaporator.Q"
        self.Capacity=self.Q-self.Fins.Air.FanPower
        
        #Sensible heat ratio [-]
        if self.Q>0.0:
            self.SHR=(self.Q_sensible_2phase+self.Q_sensible_superheat)/self.Q
        else:
            self.SHR=1.0   #
        #Average air outlet values (area fraction weighted average, neglecting influence of different air mass flowrate due to viscosity/density changes with change of omega or temperature)
        self.Tout_a=self.w_superheat*self.Tout_a_superheat+self.w_2phase*self.Tout_a_2phase
        self.omega_out=self.w_superheat*self.omega_out_superheat+self.w_2phase*self.omega_out_2phase
        self.Fins.Air.RH_out=HAProps('R','T',self.Tout_a,'P',101.325,'W',self.omega_out)  #neglect pressure drop across evaporator
        self.DP_r=self.DP_r_superheat+self.DP_r_2phase
        self.pout_r=self.psat_r+self.DP_r/1000.0  #pressure drop negative and in Pa
        
        #limit minimum outlset pressure
        if self.pout_r<10:
            print "warning - outlet pressure in evaporator too low - limiting to allow for solver to find correct value"
            self.pout_r=30#kPa
        
        #Outlet enthalpy obtained from energy balance
        self.hout_r=self.hin_r+self.Q/self.mdot_r
        
        #print "debug Evaporator - self.hin_r,",self.hin_r/1000.0,self.hout_r/1000.0," self.Tin_a", self.Fins.Air.Tdb,"self.Fins.mdot_ha",self.Fins.mdot_ha
        
        #Outlet superheat an temperature (in case of two phase)
        if existsSuperheat:
            try:
                self.Tout_r=newton(lambda T: Props('H','T',T,'P',self.pout_r,self.Ref)-self.hout_r/1000.0,Props('T','P',self.pout_r,'Q',1.0,self.Ref))
            except:
                print "problem iwith calculating Tout_r in evaporator.py"
                print "self.hout_r/1000.0",self.hout_r/1000.0,"self.hin_r/1000.0",self.hin_r/1000.0,"Props('H','Q',1.0,'P',self.pout_r,self.Ref)",Props('H','Q',1.0,'P',self.pout_r,self.Ref),"self.pout_r",self.pout_r,"self.psat_r",self.psat_r,"self.mdot_r",self.mdot_r,"self.xin_r",self.xin_r
                print "Props('H','Q',0.0,'P',self.pout_r,self.Ref)",Props('H','Q',0.0,'P',self.pout_r,self.Ref),"self.Tin_a,self.Tsat_r",self.Tin_a,self.Tsat_r,"self.DP_r/1000.0",self.DP_r/1000.0,self.Q,self.Q_superheat,self.Q_2phase,existsSuperheat,self.w_2phase
                print self.Fins.cp_da
                print "neglecting pressure drop in sh-section for normal evaproator (no issue for Calcuilate_PD, since result will be overwritten)"
                if self.w_2phase<1e-11:
                    print "neglecting capacity of two phase section, since likely not converged properly due to small fraction of two-phase area (self.w_2phase=",self.w_2phase
                    self.Q_2phase=0.0
                    self.Q=self.Q_superheat+self.Q_2phase
                    #recalculate Outlet enthalpy obtained from energy balance
                    self.hout_r=self.hin_r+self.Q/self.mdot_r
                    print "self.hout_r",self.hout_r
                self.Tout_r=newton(lambda T: Props('H','T',T,'P',self.psat_r,self.Ref)-self.hout_r/1000.0,Props('T','P',self.psat_r,'Q',1.0,self.Ref))
                self.sout_r=Props('S','T',self.Tout_r,'P',self.pin_r,self.Ref)*1000.0
            self.DT_sh_calc=self.Tout_r-self.Tdew_r
        else:
            xout_r=(self.hout_r-hsatL)/(hsatV-hsatL)
            self.sout_r=ssatV*xout_r+(1-xout_r)*ssatL
            self.DT_sh_calc=(self.hout_r-hsatV)/(Props('C','T',self.Tdew_r,'Q',1,self.Ref)*1000.0) #continuous superheat
            self.Tout_r=Props('T','P',self.pout_r,'Q',xout_r,self.Ref) #saturated temperature at outlet quality

        self.hmean_r=self.w_2phase*self.h_r_2phase+self.w_superheat*self.h_r_superheat
        self.UA_r=self.hmean_r*self.A_r_wetted
        self.UA_a=self.Fins.h_a*self.Fins.A_a*self.Fins.eta_a
        
        #Build a vector of temperatures at each point where there is a phase transition along the averaged circuit
        if existsSuperheat:
            #Insert the shoulder point
            Tv=[self.Tin_r,self.Tdew_r,self.Tout_r]
            x=[0,max(self.w_2phase,0.9999999),1]
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
        
        self.existsSuperheat=existsSuperheat  #needed for Calculate_PD()
        
    def _TwoPhase_Forward(self,w_2phase):
    
        DWS=DWSVals() #DryWetSegment structure
    
        # Store temporary values to be passed to DryWetSegment
        DWS.Fins=self.Fins
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
        if self.mdot_r<0:
            print "Warning: Refrigerant mass flowrate in evaporator is negative!"
        Q_target=self.mdot_r*(self.xout_2phase-self.xin_r)*self.h_fg
        
        if Q_target<0:
            raise ValueError('Q_target in Evaporator must be positive')
        
        # Average Refrigerant heat transfer coefficient
        DWS.h_r=ShahEvaporation_Average(self.xin_r,self.xout_2phase,self.Ref,self.G_r,self.ID,self.Tsat_r,Q_target/DWS.A_r,self.Tbubble_r,self.Tdew_r,h_tp_tuning=0.7)
        if DWS.h_r<=0.0001:
            print 'something is wrong with the refrigerant side HT-coefficient as requested by Evaporator.py, inputs are\n',self.xin_r,self.xout_2phase,self.Ref,self.G_r,self.ID,self.Tsat_r,Q_target/DWS.A_r,self.Tbubble_r,self.Tdew_r
        
        #Run the DryWetSegment to carry out the heat and mass transfer analysis
        DryWetSegment(DWS)
        
        self.Q_2phase=DWS.Q
        self.Q_sensible_2phase=DWS.Q_sensible
        self.h_r_2phase=DWS.h_r
        self.fdry_2phase=DWS.f_dry
        self.Tout_a_2phase=DWS.Tout_a
        self.omega_out_2phase=DWS.omega_out
        
        #print "in evaporator",self.xout_2phase, self.xin_r
        rho_average=TwoPhaseDensity(self.Ref,self.xin_r,self.xout_2phase,self.Tdew_r,self.Tbubble_r,slipModel='Zivi')
        self.Charge_2phase = rho_average * w_2phase * self.V_r        
        
        #Frictional pressure drop component
        DP_frict=LMPressureGradientAvg(self.xin_r,self.xout_2phase,self.Ref,self.G_r,self.ID,self.Tbubble_r,self.Tdew_r)*self.Lcircuit*w_2phase
        #Accelerational pressure drop component    
        DP_accel=AccelPressureDrop(self.xin_r,self.xout_2phase,self.Ref,self.G_r,self.Tbubble_r,self.Tdew_r)
        self.DP_r_2phase=DP_frict+DP_accel;
        
        if self.Verbosity>7:
            print w_2phase,DWS.Q,Q_target,self.xin_r,"w_2phase,DWS.Q,Q_target,self.xin_r"
        return DWS.Q-Q_target
    
    def _Superheat_Forward(self,w_superheat,T_inlet_r=-99):
        if T_inlet_r==-99:  #catch if default case is used
            T_inlet_r=self.Tdew_r
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
    
        # Inputs on the air side to two phase region are inlet air again
        DWS.Tin_a=self.Tin_a
        DWS.RHin_a=self.Fins.Air.RH
    
        DWS.Tin_r=T_inlet_r  #default is dew temperature; can change to consider for initial superheat
        DWS.A_r=self.A_r_wetted*w_superheat
        DWS.cp_r=self.cp_r
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
        try:
            rho_superheat=Props('D','T',(DWS.Tout_r+self.Tdew_r)/2.0, 'P', self.psat_r, self.Ref)
            #members = [attr for attr in dir(DWS()) if not callable(attr) and not attr.startswith("__")]
            #print members
        except:
            print "error in Evaporator, unreasonable inputs?","Inputs to calculate density are DWS.Tout_r,self.Tdew_r,self.psat_r",DWS.Tout_r,self.Tdew_r,self.psat_r,"<<"
            print "DWS values are", 'DWS.A_a',DWS.A_a,'DWS.Tin_r',DWS.Tin_r,'DWS.mdot_da',DWS.mdot_da,"DWS.mdot_r",DWS.mdot_r,'DWS.mdot_da',DWS.mdot_da,'DWS.mdot_r',DWS.mdot_r
            print "plot DWS vals to figure out what is going on"
            raise()
        self.Charge_superheat = w_superheat * self.V_r * rho_superheat
        
        #Pressure drop calculations for subcooled refrigerant
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
        self.omega_out_superheat=DWS.omega_out
        
if __name__=='__main__':
    #This code runs if this file is run by itself, but otherwise doesn't run3
    if True:
        Ref='R410a'
        TT=[]
        Q=[]
        ff=[]
        QQ=[]
        hh=[]
        T_outa=[]
        import numpy as np
        import pylab as pylab
        Tdew=260
        lower_value=Props('H','T',Tdew,'Q',0.01,Ref)*1000
        #lower_value=Props('H','P',Props('P','T',Tdew,'Q',1.0,'R410A'),'T',299,Ref)*1000
        upper_value=Props('H','P',Props('P','T',Tdew,'Q',1.0,Ref),'T',299.79999,Ref)*1000
        for hin_r in np.linspace(lower_value,upper_value,20):  
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
            FinsTubes.Air.p=101.325
            FinsTubes.Air.RH=0.51
            FinsTubes.Air.RHmean=0.51
            FinsTubes.Air.FanPower=438
                
            kwargs={'Ref': Ref,
                    'mdot_r':  0.0708,
                    'psat_r':  Props('P','T',Tdew,'Q',1.0,Ref),
                    'Fins': FinsTubes,
                    'hin_r':hin_r,
                    'Verbosity':0
            }
            Evap=EvaporatorClass(**kwargs)
            Evap.Update(**kwargs) #update not necessary here, but left for illustration
            Evap.Calculate()
            #Evap.Calculate_PD()
            
            #print Evap.OutputList()
            Q.append(Evap.Q)
            QQ.append(Evap.Q_2phase)
            TT.append(hin_r)
            ff.append(Evap.w_2phase)
            hh.append(Evap.h_r_2phase)
            T_outa.append(Evap.Tout_a)
        import time
        start_time = time.clock()
        for i in range(100):
            Evap.Calculate()
        print time.clock() - start_time, "seconds"
        pylab.plot(TT,Q,'x-',label='total')
        pylab.plot(TT,QQ,'x-',label='2-phase')
        pylab.legend()
        pylab.figure()
        pylab.plot(TT,hh,'o-',label='HTC 2-phase')
        pylab.legend()
        pylab.figure()
        pylab.plot(TT,T_outa,'o-',label='Tout_a')
        pylab.legend()
        pylab.show()

#########################################################################
    if False: #run a change-in-dewpoint temperature study
        TT=[]
        Q=[]
        ff=[]
        QQ=[]
        hh=[]
        T_outa=[]
        import numpy as np,pylab
        for Tdew in np.linspace(260,290,10):  
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
            FinsTubes.Air.p=101.325
            FinsTubes.Air.RH=0.51
            FinsTubes.Air.RHmean=0.51
            FinsTubes.Air.FanPower=438
                
            kwargs={'Ref': 'R410A',
                    'mdot_r':  0.0708,
                    'psat_r':  Props('P','T',Tdew,'Q',1.0,'R410A'),
                    'Fins': FinsTubes,
                    'hin_r':Props('H','T',Tdew,'Q',0.15,'R410A')*1000,
                    'Verbosity':0
            }
            
            Evap=EvaporatorClass(**kwargs)
            Evap.Update(**kwargs) #update not necessary here, but left for illustration
            Evap.Calculate()
            Evap.Calculate_PD()
            
            #print Evap.OutputList()
            Q.append(Evap.Q)
            QQ.append(Evap.Q_2phase)
            TT.append(Tdew)
            T_outa.append(Evap.Tout_a)
            
            ff.append(Evap.w_2phase)
            hh.append(Evap.h_r_2phase)
        #determine time used for calculation
        import time
        num_runs=100
        Evap.psat_r=Props('P','T',260,'Q',1.0,'R410A')
        start_time = time.clock()
        for i in range(10):
            Evap.Calculate_PD()
        print time.clock() - start_time, "seconds to run calculate_PD(0)"
        print Evap.OutputList()
        start_time = time.clock()
        for i in range(10):
            Evap.Calculate()
        print time.clock() - start_time, "seconds to run normale calculate fct"
        print Evap.OutputList()
        pylab.plot(TT,QQ,TT,Q,'x-',label=['2-phase','total'])
        pylab.legend()
        pylab.show()
        pylab.plot(TT,hh)
        pylab.show()
