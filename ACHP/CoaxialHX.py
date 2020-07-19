from __future__ import division, print_function, absolute_import
from CoolProp.CoolProp import PropsSI
from ACHP.Correlations import f_h_1phase_Annulus,f_h_1phase_Tube,ShahEvaporation_Average
from ACHP.Correlations import TwoPhaseDensity,LMPressureGradientAvg,AccelPressureDrop
from math import pi,exp,log
from scipy.optimize import brentq
import numpy as np
import CoolProp as CP

class CoaxialHXClass():
    def __init__(self,**kwargs):
        #Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def Update(self,**kwargs):
        #Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
        # Wetted area on the refrigerant side
        self.A_r_wetted=pi*self.ID_i*self.L
        # Wetted area of the glycol side (not including outer tube)
        self.A_g_wetted=pi*self.OD_i*self.L
        
        self.V_r=self.L*pi*self.ID_i**2/4.0
        self.V_g=self.L*pi*(self.ID_o**2-self.OD_i**2)/4.0
    
    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        
        return [
            ('Length of tube','m',self.L),
            ('Annulus wetted OD','m',self.ID_o),
            ('Tube wetted OD/Annulus wetted ID','m',self.OD_i),
            ('Tube wetted ID','m',self.ID_i),
            ('Outlet Superheat','K',self.Tin_r-self.Tdew_r),
            ('Q Total','W',self.Q),
            ('Q Superheat','W',self.Q_superheat),
            ('Q Two-Phase','W',self.Q_2phase),
            ('Q Subcooled','W',self.Q_subcool),
            ('Inlet glycol temp','K',self.Tin_g),
            ('Outlet glycol temp','K',self.Tout_g),
            ('Inlet ref. temp','K',self.Tin_r),
            ('Outlet ref. temp','K',self.Tout_r),
            ('Charge Total','kg',self.Charge_r),
            ('Charge Superheat','kg',self.Charge_r_superheat),
            ('Charge Two-Phase','kg',self.Charge_r_2phase),
            ('Charge Subcool','kg',self.Charge_r_subcool),
            ('Mean HTC Ref. Superheat','W/m^2-K',self.h_r_superheat),
            ('Mean HTC Ref. Two-Phase','W/m^2-K',self.h_r_2phase),
            ('Mean HTC Ref. Subcool','W/m^2-K',self.h_r_subcool),
            ('Mean HTC Gly. Superheat','W/m^2-K',self.h_g),
            ('Mean Reynolds # Gly. Superheat','-',self.Re_g),
            ('Pressure Drop Gly.','Pa',self.DP_g),
            ('Pressure Drop Ref.','Pa',self.DP_r),
            ('Pressure Drop Ref. Superheat','Pa',self.DP_r_superheat),
            ('Pressure Drop Ref. Two-Phase','Pa',self.DP_r_2phase),
            ('Pressure Drop Ref. Subcool','Pa',self.DP_r_subcool),                                                                           
            ('Area fraction Superheat','-',self.w_superheat),
            ('Area fraction Two-Phase','-',self.w_2phase),
            ('Area fraction Subcooled','-',self.w_subcool)
         ]
        
    def Calculate(self):
        #AbstractState (Ref)
        AS_r = self.AS_r
        #AbstractState (Glycol)
        AS_g = self.AS_g
        if hasattr(self,'MassFrac_g'):
            AS_g.set_mass_fractions([self.MassFrac_g])
        elif hasattr(self, 'VoluFrac_g'):
            AS_g.set_volu_fractions([self.VoluFrac_g])
                
        #set tuning factors to 1 in case not given by user
        if not hasattr(self,'h_g_tuning'):
            self.h_g_tuning = 1
        if not hasattr(self,'h_tp_tuning'):
            self.h_tp_tuning = 1
        if not hasattr(self,'DP_g_tuning'):
            self.DP_g_tuning = 1
        if not hasattr(self,'DP_r_tuning'):
            self.DP_r_tuning = 1
            
        #Update the parameters
        self.Update()
        
        #Average mass flux of refrigerant [kg/m^2-s]
        self.G_r = self.mdot_r/(pi*self.ID_i**2/4.0) 
        #Average mass flux of glycol [kg/m^2-s]
        self.G_g = self.mdot_g/(pi*(self.ID_o**2-self.OD_i**2)/4.0) 
        #Hydraulic diameter
        self.Dh_g=self.ID_o-self.OD_i
        #Evaporation hydraulic diameter [m]
        self.Dh_r=self.ID_i
        #Thermal conductivity of the intermediate wall pipe
        self.k =self.Conductivity
        
        #Thermal Conduction Resistance of the intermediate wall
        self.Rw=log(self.OD_i/self.ID_i) / (2*pi*self.k*self.L)
        
        AS_r.update(CP.PQ_INPUTS,self.pin_r,0.0)
        self.Tbubble_r=AS_r.T() #[K]
        hsatL=AS_r.hmass() #[J/kg]
        ssatL=AS_r.smass() #[J/kg-K]
        AS_r.update(CP.PQ_INPUTS,self.pin_r,1.0)
        self.Tdew_r=AS_r.T() #[K]
        hsatV=AS_r.hmass() #[J/kg]
        ssatV=AS_r.smass() #[J/kg-K]
        
        #Saturation temperature
        self.Tsat_r=(self.Tbubble_r+self.Tdew_r)/2.0
        
        #Inlet quality        
        self.xin_r=(self.hin_r-hsatL)/(hsatV-hsatL)
        
        #Change in enthalpy through two-phase region [J/kg]
        self.h_fg=hsatV - hsatL
        self.Tin_r=self.xin_r*self.Tdew_r+(1-self.xin_r)*self.Tbubble_r
        #Inlet entropy
        self.sin_r=self.xin_r*ssatV+(1-self.xin_r)*ssatL
        
        #Mean values for the glycol side based on average of inlet temperatures
        Tavg_g=(self.Tsat_r+self.Tin_g)/2.0
        self.f_g,self.h_g,self.Re_g=f_h_1phase_Annulus(self.mdot_g, self.ID_o, self.OD_i, Tavg_g, self.pin_g, self.AS_g)
        self.h_g = self.h_g * self.h_g_tuning #correct h_g with tuning factor
        
        AS_g.update(CP.PT_INPUTS,self.pin_g,Tavg_g)
        self.cp_g=AS_g.cpmass() #[J/kg-K]
        v_g=1/AS_g.rhomass() #[m^3/kg]
        
        #Glycol pressure drop
        dpdz_g=-self.f_g*v_g*self.G_g**2/(2.*self.Dh_g) #Pressure gradient
        self.DP_g=dpdz_g*self.L*self.DP_g_tuning
        
        def OBJECTIVE(w_superheat):
            """Nested function for driving the Brent's method solver"""
            #Run the superheated portion
            self._Superheat_Forward(w_superheat)
            #Run the two-phase portion and return residual
            return self._TwoPhase_Forward(1-w_superheat)
        
        def OBJECTIVE_2phase(xout_2phase):
            """Nested function for finding outlet quality"""
            #Need to pass in outlet quality but still maintain the full w_2phase 
            return self._TwoPhase_Forward(1.0,xout_2phase)
        
        # First see if you have a superheated portion.  Try to use the entire HX
        # for the two-phase portion
        # ---------------------
        # Intermediate glycol temp between superheated and 2phase sections [K] 
        self.T_g_x=self.Tin_g
        #Call two-phase forward method
        error=self._TwoPhase_Forward(1.0)
        # If HT greater than required
        if error>0:
            #Too much HT if all is 2phase, there is a superheated section
            existsSuperheat=True
            #Solve for the break between 2phase and superheated parts
            w_superheat=brentq(OBJECTIVE,0.00001,0.99999)
            self.w_2phase=1-w_superheat
            self.w_superheat=w_superheat
        else:
            existsSuperheat=False
            # Solve for outlet quality in 2phase section, lowest possible outlet 
            # quality is the inlet quality
            self.xout_2phase=brentq(OBJECTIVE_2phase,self.xin_r,0.99999)
            #Dummy variables for the superheated section which doesn't exist
            self.Q_superheat=0.0
            self.Charge_r_superheat=0.0
            self.h_r_superheat=0.0
            self.DP_r_superheat=0.0
            self.w_superheat=0.0
            self.w_2phase=1.0
            
        self.Charge_r=self.Charge_r_2phase+self.Charge_r_superheat
        self.Q=self.Q_2phase+self.Q_superheat
        self.Tout_g=self.Tin_g-self.Q/(self.cp_g*self.mdot_g)
        
        self.DP_r=(self.DP_r_2phase+self.DP_r_superheat)*self.DP_r_tuning
        
        if existsSuperheat==True:
            AS_r.update(CP.PT_INPUTS,self.pin_r,self.Tout_r)
            self.hout_r=AS_r.hmass() #[J/kg]
            self.sout_r=AS_r.smass() #[J/kg-K]
        else:
            self.Tout_r=self.xout_2phase*self.Tdew_r+(1-self.xout_2phase)*self.Tbubble_r
            AS_r.update(CP.QT_INPUTS,self.xout_2phase,self.Tout_r)
            self.hout_r=AS_r.hmass() #[J/kg]
            self.sout_r=AS_r.smass() #[J/kg-K]
        
#         #Dummy variables for the subcooled section which doesn't exist
#         self.Q_subcool=0.0
#         self.DP_r_subcool=0.0
#         self.h_r_subcool=0.0
#         self.Re_r_subcool=0.0
#         self.Charge_r_subcool=0.0
#         self.w_subcool=0.0
        
    def _Superheat_Forward(self,w_superheat):
        # Superheated portion
        # Mean temperature for superheated part can be taken to be average
        # of dew and glycol inlet temps
        Tavg_sh_r=(self.Tdew_r+self.Tin_g)/2.0
        self.f_r_superheat,self.h_r_superheat,self.Re_r_superheat=f_h_1phase_Tube(self.mdot_r, self.ID_i, Tavg_sh_r, self.pin_r, self.AS_r)
        # Refrigerant specific heat
        self.AS_r.update(CP.PT_INPUTS,self.pin_r,Tavg_sh_r)
        cp_r_superheat=self.AS_r.cpmass() #[J/kg-K]
        # Overall conductance of heat transfer surface in superheated
        # portion
        UA_superheat=w_superheat/(1/(self.h_g*self.A_g_wetted)+1/(self.h_r_superheat*self.A_r_wetted)+self.Rw)
        #List of capacitance rates [W/K]
        C=[cp_r_superheat*self.mdot_r,self.cp_g*self.mdot_g]
        Cmin=min(C)
        Cr=Cmin/max(C)
        Ntu_superheat=UA_superheat/Cmin
        
        #for Cr<1 and pure counter flow (Incropera. Table 11.3)
        epsilon_superheat = ((1 - exp(-Ntu_superheat * (1 - Cr))) / 
            (1 - Cr * exp(-Ntu_superheat * (1 - Cr))))
        
        self.Q_superheat=epsilon_superheat*Cmin*(self.Tin_g-self.Tdew_r)
        
        self.Tout_r=self.Tdew_r+self.Q_superheat/(self.mdot_r*cp_r_superheat)
        # Refrigerant density (supeheated)
        self.AS_r.update(CP.PT_INPUTS,self.pin_r,(self.Tin_g+self.Tdew_r)/2.0) 
        rho_superheat=self.AS_r.rhomass() #[kg/m^3]
        # Regrigerant charge (supeheated)
        self.Charge_r_superheat = w_superheat * self.V_r * rho_superheat
        
        #Pressure drop calculations for superheated refrigerant
        v_r=1./rho_superheat
        #Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_superheat*v_r*self.G_r**2/(2.*self.Dh_r) #Pressure gradient
        self.DP_r_superheat=dpdz_r*self.L*w_superheat
        
        # Temperature of "glycol" at the point where the refrigerant is at
        # a quality of 1.0 [K]
        self.T_g_x=self.Tin_g-self.Q_superheat/(self.cp_g*self.mdot_g)
        
    def _TwoPhase_Forward(self,w_2phase=1.0,xout_2phase=1.0):
        # Update outlet quality field [-]
        self.xout_2phase=xout_2phase
        # Heat transfer rate based on inlet quality [W]
        self.Q_2phase=self.mdot_r*(self.xout_2phase-self.xin_r)*self.h_fg
        # Heat flux in 2phase section (for Shah correlation) [W/m^2]
        q_flux=self.Q_2phase/(w_2phase*self.A_r_wetted)
        #
        h_r_2phase=ShahEvaporation_Average(self.xin_r,1.0,self.AS_r,
                    self.G_r,self.Dh_r,self.pin_r,q_flux,self.Tbubble_r,self.Tdew_r)
        self.h_r_2phase = h_r_2phase * self.h_tp_tuning
        UA_2phase=w_2phase/(1/(self.h_g*self.A_g_wetted)+1/(self.h_r_2phase*self.A_r_wetted)+self.Rw)
        C_g=self.cp_g*self.mdot_g
        Ntu_2phase=UA_2phase/(C_g)
        
        #for Cr=0 and counter flow in two-phase:
        epsilon_2phase=1-exp(-Ntu_2phase)
        Q_2phase_eNTU=epsilon_2phase*C_g*(self.T_g_x-self.Tsat_r)
        
        rho_average=TwoPhaseDensity(self.AS_r,self.xin_r,xout_2phase,self.Tdew_r,self.Tbubble_r,slipModel='Zivi')
        self.Charge_r_2phase = rho_average * w_2phase * self.V_r     
        
        # Frictional prssure drop component
        DP_frict=LMPressureGradientAvg(self.xin_r,xout_2phase,self.AS_r,self.G_r,self.Dh_r,self.Tbubble_r,self.Tdew_r)*w_2phase*self.L
        # Accelerational prssure drop component
        DP_accel=AccelPressureDrop(self.xin_r,xout_2phase,self.AS_r,self.G_r,self.Tbubble_r,self.Tdew_r)*w_2phase*self.L
        self.DP_r_2phase=DP_frict+DP_accel
        
        if self.Verbosity>4:
            print (Q_2phase_eNTU-self.Q_2phase)
        return Q_2phase_eNTU-self.Q_2phase
        
if __name__=='__main__':
    
    TT=[]
    QQ=[]
    Q1=[]
    #refigearnt Abstract State
    Ref_r = 'R290'
    Backend_r = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS_r = CP.AbstractState(Backend_r, Ref_r)
    #glycol Abstract State
    Ref_g = 'Water'
    Backend_g = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS_g = CP.AbstractState(Backend_g, Ref_g)
    for Tdew_evap in np.linspace(270,290.4):
        Tdew_cond=317.73
#        Tdew_evap=285.42
        pdew_cond=PropsSI('P','T',Tdew_cond,'Q',1.0,Ref_r)
        h=PropsSI('H','T',Tdew_cond-7,'P',pdew_cond,Ref_r)
        params={
                'ID_i':0.0278,      #inner tube, Internal Diameter (ID)
                'OD_i':0.03415,     #inner tube, Outer Diameter (OD)
                'ID_o':0.045,       #outer tube (annulus), Internal Diameter (ID)
                'L':50,
                'mdot_r':0.040,
                'mdot_g':0.38,
                'hin_r':h,
                'pin_r':PropsSI('P','T',Tdew_evap,'Q',1.0,Ref_r),
                'pin_g':300000,     #pin_g in Pa
                'Tin_g':290.52,
                'AS_r':AS_r, #Abstract state of refigerant
                'AS_g':AS_g, ##Abstract state of glycol
                'Verbosity':0,
                'Conductivity' : 237, #[W/m-K]
                'h_g_tuning':1,
                'h_tp_tuning':1,
                'DP_g_tuning':1,
                'DP_r_tuning':1
                }
        IHX=CoaxialHXClass(**params)
        IHX.Calculate()
        
        TT.append(Tdew_evap)
        QQ.append(IHX.h_r_2phase)#IHX.Q)
        Q1.append(IHX.h_r_superheat)
                  
        print (IHX.Q)