from __future__ import division, absolute_import, print_function
from CoolProp.CoolProp import PropsSI
import CoolProp as CP

from math import pi,exp,log,sqrt,tan,cos,sin
from scipy.optimize import brentq
import numpy as np
import pylab

from ACHP.Correlations import ShahEvaporation_Average,PHE_1phase_hdP,Cooper_PoolBoiling,TwoPhaseDensity,TrhoPhase_ph,Phase_ph,LMPressureGradientAvg,KandlikarPHE,Bertsch_MC,AccelPressureDrop,ShahCondensation_Average,LongoCondensation

class PHEHXClass():
    """
    There are a number of possibilities:
        
        Each fluid can:
        a) Not change phase
        c) Evaporate 
        c) Condense
        
        Possibility matrix
        
                                      Hot stream
         Cold Stream  || Subcooled ||  Two-Phase || Superheated ||
                      --------------------------------------------
         Subcooled    ||           ||            ||             ||
                      --------------------------------------------
         Two-Phase    ||           ||            ||             ||
                      --------------------------------------------
         Superheated  ||           ||            ||             ||
                      --------------------------------------------
    
    Hot stream goes to the left in the matrix, cold stream goes down.  If 
    hot stream comes in subcooled, there are only three combinations that 
    might exist.  
    
    Based on inlet states can figure out what states are possible.  
    """
    
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
        #TODO: fix this list of outputs
        return [
            ('Effective Length','m',self.Lp),
            ('Wetted area','m^2',self.A_h_wetted),
            ('Outlet Superheat','K',self.Tin_c-self.Tdew_c),
            ('Q Total','W',self.Q),
            ('Q Superheat Hot','W',self.Q_superheated_h),
            ('Q Two-Phase Hot','W',self.Q_2phase_h),
            ('Q Subcooled Hot','W',self.Q_subcooled_h),
            ('Q Superheat Cold','W',self.Q_superheated_c),
            ('Q Two-Phase Cold','W',self.Q_2phase_c),
            ('Q Subcooled Cold','W',self.Q_subcooled_c),
            ('Inlet hot stream temp','K',self.Tin_h),
            ('Outlet hot stream temp','K',self.Tout_h),
            ('Inlet cold stream temp','K',self.Tin_c),
            ('Outlet cold stream temp','K',self.Tout_c),
            ('Charge Total','kg',self.Charge_h),
            ('Charge Superheat','kg',self.Charge_superheated_h),
            ('Charge Two-Phase','kg',self.Charge_2phase_h),
            ('Charge Subcool','kg',self.Charge_subcooled_h),
            ('Charge Total','kg',self.Charge_c),
            ('Charge Superheat','kg',self.Charge_superheated_c),
            ('Charge Two-Phase','kg',self.Charge_2phase_c),
            ('Charge Subcool','kg',self.Charge_subcooled_c),
            ('Hot HTC Superheat','W/m^2-K',self.h_superheated_h),
            ('Hot HTC Two-Phase','W/m^2-K',self.h_2phase_h),
            ('Hot HTC Subcool','W/m^2-K',self.h_subcooled_h),
            ('Cold Mean HTC Superheat','W/m^2-K',self.h_superheated_c),
            ('Cold Mean HTC Ref. Two-Phase','W/m^2-K',self.h_2phase_c),
            ('Cold Mean HTC Ref. Subcool','W/m^2-K',self.h_subcooled_c),
            ('Pressure Drop Hot','Pa',self.DP_h),
            ('Pressure Drop Hot superheated','Pa',self.DP_superheated_h),
            ('Pressure Drop Hot 2 phase','Pa',self.DP_2phase_h),
            ('Pressure Drop Hot subcooled','Pa',self.DP_subcooled_h),
            ('Pressure Drop Cold','Pa',self.DP_c),
            ('Pressure Drop Cold superheated','Pa',self.DP_superheated_c),
            ('Pressure Drop Cold 2 phase','Pa',self.DP_2phase_c),
            ('Pressure Drop Cold subcooled','Pa',self.DP_subcooled_c),
            ('Area fraction Superheat Hot','-',self.w_superheated_h),
            ('Area fraction Two-Phase Hot','-',self.w_2phase_h),
            ('Area fraction Subcooled Hot','-',self.w_subcooled_h),
            ('Area fraction Superheat Cold','-',self.w_superheated_c),
            ('Area fraction Two-Phase Cold','-',self.w_2phase_c),
            ('Area fraction Subcooled Cold','-',self.w_subcooled_c)
         ]
        
    def DetermineHTBounds(self):
        # See if each phase could change phase if it were to reach the
        # inlet temperature of the opposite phase 
        
        #AbstractState
        AS_h = self.AS_h
        AS_c = self.AS_c
        
        #Inlet phases
        self.Tin_h,rhoin_h,Phasein_h=TrhoPhase_ph(self.AS_h,self.pin_h,self.hin_h,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)
        self.Tin_c,rhoin_c,Phasein_c=TrhoPhase_ph(self.AS_c,self.pin_c,self.hin_c,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)
        
        # Find the maximum possible rate of heat transfer as the minimum of 
        # taking each stream to the inlet temperature of the other stream
        AS_h.update(CP.PT_INPUTS, self.pin_h, self.Tin_c)
        hout_h=AS_h.hmass() #[J/kg]
        AS_c.update(CP.PT_INPUTS, self.pin_c, self.Tin_h)
        hout_c=AS_c.hmass() #[J/kg]
        Qmax=min([self.mdot_c*(hout_c-self.hin_c),self.mdot_h*(self.hin_h-hout_h)])
        if Qmax<0:
            raise ValueError('Qmax in PHE must be > 0')
        
        # Now we need to check for internal pinch points where the temperature
        # profiles would tend to overlap given the "normal" definitions of 
        # maximum heat transfer of taking each stream to the inlet temperature 
        # of the other stream
        #
        # First we build the same vectors of enthalpies like below
        EnthalpyList_c,EnthalpyList_h=self.BuildEnthalpyLists(Qmax)
        # Then we find the temperature of each stream at each junction
        TList_c=np.zeros_like(EnthalpyList_c)
        TList_h=np.zeros_like(EnthalpyList_h)

        if not len(EnthalpyList_h)==len(EnthalpyList_c):
            raise ValueError('Length of enthalpy lists for both fluids must be the same')
        
        #Make the lists of temperatures of each fluid at each cell boundary
        for i in range(len(EnthalpyList_h)):
            TList_c[i]=TrhoPhase_ph(self.AS_c,self.pin_c,EnthalpyList_c[i],self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0]
            TList_h[i]=TrhoPhase_ph(self.AS_h,self.pin_h,EnthalpyList_h[i],self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0]

#        #Double-check that the edges are not pinched
#        if TList_c[0]-1e-9>TList_h[0] or TList_c[-1]-1e-9>TList_h[-1]:
#            raise ValueError('Outlet or inlet of PHE is pinching.  Why?')
        
        #TODO: could do with more generality if both streams can change phase
        #Check if any internal points are pinched
        if np.sum(TList_c>TList_h)>0:
            #Loop over the internal cell boundaries
            for i in range(1,len(TList_c)-1):
                #If cold stream is hotter than the hot stream
                if TList_c[i]-1e-9>TList_h[i]:
                    #Find new enthalpy of cold stream at the hot stream cell boundary
                    AS_c.update(CP.PT_INPUTS, self.pin_c, TList_h[i])
                    hpinch=AS_c.hmass() #[J/kg]
                    #Find heat transfer of hot stream in right-most cell
                    Qextra=self.mdot_h*(EnthalpyList_h[i+1]-EnthalpyList_h[i])
                    Qmax=self.mdot_c*(hpinch-self.hin_c)+Qextra
        
        return Qmax
        
    def PlateHTDP(self,AS,T,p,mdot_gap):
        """
        This function calls mainly the heat transfer and pressure drop 
        for single phase fluids in PHE. 
        Inputs: T [K] and p [Pa]
        Outputs: h [W/m^2-K] and cp [J/kg-K]
        Note: There are several other output. Check "PHE_1phase_hdP" function
        for more details.
        """
        Inputs={
            'AS':AS, #AS: AbstractState which includes refrigerant or glycol and backend
            'T':T,
            'p':p,
            'mdot_gap' : mdot_gap,                  #mass flow rate per channel
            'PlateAmplitude': self.PlateAmplitude,
            'PlateWavelength' : self.PlateWavelength,
            'InclinationAngle': self.InclinationAngle,
            'Bp': self.Bp,
            'Lp': self.Lp
        }
        Outputs=PHE_1phase_hdP(Inputs)
        return Outputs['h'],Outputs['cp'],Outputs
    
    def BuildEnthalpyLists(self,Q):
        #Start the enthalpy lists with inlet and outlet enthalpies
        #Ordered from lowest to highest enthalpies for both streams
        EnthalpyList_h=[self.hin_h-Q/self.mdot_h, self.hin_h]
        EnthalpyList_c=[self.hin_c,self.hin_c+Q/self.mdot_c]
        
        #Save the value of Q and outlet enthalpies
        self.Q=Q
        self.hout_h=EnthalpyList_h[0]
        self.hout_c=EnthalpyList_c[1]
        
        #Call AbstractState
        AS_h = self.AS_h
        AS_c = self.AS_c
        
        #Find the phase boundaries that exist, and add them to lists
        if 'IncompressibleBackend' in AS_h.backend_name():
            hsatL_h=1e9
            hsatV_h=1e9
        else:
            AS_h.update(CP.DmassT_INPUTS, self.rhosatL_h, self.Tbubble_h)
            hsatL_h=AS_h.hmass() #[J/kg]
            AS_h.update(CP.DmassT_INPUTS, self.rhosatV_h, self.Tdew_h)
            hsatV_h=AS_h.hmass() #[J/kg]
        
        if 'IncompressibleBackend' in AS_c.backend_name():
            hsatL_c=1e9
            hsatV_c=1e9
        else:
            AS_c.update(CP.DmassT_INPUTS, self.rhosatL_c, self.Tbubble_c)
            hsatL_c=AS_c.hmass() #[J/kg]
            AS_c.update(CP.DmassT_INPUTS, self.rhosatV_c, self.Tdew_c)
            hsatV_c=AS_c.hmass() #[J/kg]
        
        # Check whether the enthalpy boundaries are within the bounds set by 
        # the imposed amount of heat transfer
        if hsatV_c<EnthalpyList_c[-1] and hsatV_c>EnthalpyList_c[0]:
            EnthalpyList_c.insert(len(EnthalpyList_c)-1,hsatV_c)
        if hsatL_c<EnthalpyList_c[-1] and hsatL_c>EnthalpyList_c[0]:
            EnthalpyList_c.insert(1,hsatL_c)
            
        if hsatV_h<EnthalpyList_h[-1] and hsatV_h>EnthalpyList_h[0]:
            EnthalpyList_h.insert(len(EnthalpyList_h)-1,hsatV_h)
        if hsatL_h<EnthalpyList_h[-1] and hsatL_h>EnthalpyList_h[0]:
            EnthalpyList_h.insert(1,hsatL_h)
            
        I_h=0
        I_c=0
        while I_h<len(EnthalpyList_h)-1:
            #Try to figure out whether the next phase transition is on the hot or cold side     
            Qbound_h=self.mdot_h*(EnthalpyList_h[I_h+1]-EnthalpyList_h[I_h])
            Qbound_c=self.mdot_c*(EnthalpyList_c[I_c+1]-EnthalpyList_c[I_c])
            if Qbound_h<Qbound_c-1e-9:
                # Minimum amount of heat transfer is on the hot side,
                # add another entry to EnthalpyList_c 
                EnthalpyList_c.insert(I_c+1, EnthalpyList_c[I_c]+Qbound_h/self.mdot_c)
            elif Qbound_h>Qbound_c+1e-9:
                # Minimum amount of heat transfer is on the cold side,
                # add another entry to EnthalpyList_h at the interface
                EnthalpyList_h.insert(I_h+1, EnthalpyList_h[I_h]+Qbound_c/self.mdot_h)
            I_h+=1
            I_c+=1
                    
        self.hsatL_c=hsatL_c
        self.hsatL_h=hsatL_h
        self.hsatV_c=hsatV_c
        self.hsatV_h=hsatV_h

        return EnthalpyList_c,EnthalpyList_h
    
    def PostProcess(self,cellList):
        """
        Combine all the cells to calculate overall parameters like pressure drop
        and fraction of heat exchanger in two-phase on both sides
        """
        #AbstractState
        AS_c = self.AS_c
        AS_h = self.AS_h
        
        def collect(cellList,tag,tagvalue,out):
            collectList=[]
            for cell in cellList:
                if cell[tag]==tagvalue:
                    collectList.append(cell[out])
            return collectList
        
        self.DP_c=0
        self.DP_c_superheat=0
        self.DP_c_2phase=0
        self.DP_c_subcooled=0
        self.DP_h=0
        self.DP_h_superheat=0
        self.DP_h_2phase=0
        self.DP_h_subcooled=0
        self.Charge_c=0
        self.Charge_h=0
        for cell in cellList:
            self.DP_c+=cell['DP_c']
            self.DP_h+=cell['DP_h']
            self.Charge_c+=cell['Charge_c']
            self.Charge_h+=cell['Charge_h']
        self.w_superheated_h=sum(collect(cellList,'Phase_h','Superheated','w'))
        self.w_2phase_h=sum(collect(cellList,'Phase_h','TwoPhase','w'))
        self.w_subcooled_h=sum(collect(cellList,'Phase_h','Subcooled','w'))
        self.w_superheated_c=sum(collect(cellList,'Phase_c','Superheated','w'))
        self.w_2phase_c=sum(collect(cellList,'Phase_c','TwoPhase','w'))
        self.w_subcooled_c=sum(collect(cellList,'Phase_c','Subcooled','w'))
        
        self.DP_superheated_c=sum(collect(cellList,'Phase_c','Superheated','DP_c'))
        self.DP_2phase_c=sum(collect(cellList,'Phase_c','TwoPhase','DP_c'))
        self.DP_subcooled_c=sum(collect(cellList,'Phase_c','Subcooled','DP_c'))
        self.DP_c=self.DP_superheated_c+self.DP_2phase_c+self.DP_subcooled_c
        self.DP_c=self.DP_c*self.DP_cold_tuning #correct the pressure drop on the cold side
        
        self.DP_superheated_h=sum(collect(cellList,'Phase_h','Superheated','DP_h'))
        self.DP_2phase_h=sum(collect(cellList,'Phase_h','TwoPhase','DP_h'))
        self.DP_subcooled_h=sum(collect(cellList,'Phase_h','Subcooled','DP_h'))
        self.DP_h=self.DP_superheated_h+self.DP_2phase_h+self.DP_subcooled_h
        self.DP_h=self.DP_h*self.DP_hot_tuning #correct the pressure drop on the hot side
        
        self.Charge_superheated_c=sum(collect(cellList,'Phase_c','Superheated','Charge_c'))
        self.Charge_2phase_c=sum(collect(cellList,'Phase_c','TwoPhase','Charge_c'))
        self.Charge_subcooled_c=sum(collect(cellList,'Phase_c','Subcooled','Charge_c'))
        self.Charge_c=self.Charge_superheated_c+self.Charge_2phase_c+self.Charge_subcooled_c
        self.Charge_superheated_h=sum(collect(cellList,'Phase_h','Superheated','Charge_h'))
        self.Charge_2phase_h=sum(collect(cellList,'Phase_h','TwoPhase','Charge_h'))
        self.Charge_subcooled_h=sum(collect(cellList,'Phase_h','Subcooled','Charge_h'))
        self.Charge_h=self.Charge_superheated_h+self.Charge_2phase_h+self.Charge_subcooled_h
        
        self.Q_superheated_h=sum(collect(cellList,'Phase_h','Superheated','Q'))
        self.Q_2phase_h=sum(collect(cellList,'Phase_h','TwoPhase','Q'))
        self.Q_subcooled_h=sum(collect(cellList,'Phase_h','Subcooled','Q'))
        self.Q_superheated_c=sum(collect(cellList,'Phase_c','Superheated','Q'))
        self.Q_2phase_c=sum(collect(cellList,'Phase_c','TwoPhase','Q'))
        self.Q_subcooled_c=sum(collect(cellList,'Phase_c','Subcooled','Q'))
        
        w_superheat=collect(cellList,'Phase_c','Superheated','w')
        w_2phase=collect(cellList,'Phase_c','TwoPhase','w')
        h_c_sh=collect(cellList,'Phase_c','Superheated','h_c')
        h_c_2phase=collect(cellList,'Phase_c','TwoPhase','h_c')
        
        if len(w_superheat)>0:
            self.h_superheated_c=float(sum(np.array(h_c_sh)*np.array(w_superheat))/sum(w_superheat))
        else:
            self.h_superheated_c=0
            
        if len(w_2phase)>0:
            self.h_2phase_c=float(sum(np.array(h_c_2phase)*np.array(w_2phase))/sum(w_2phase))
        else:
            self.h_2phase_c=0
            
        ### Collect all the cells on the hot side
        w_subcooled_h=collect(cellList,'Phase_h','Subcooled','w')
        w_superheat_h=collect(cellList,'Phase_h','Superheated','w')
        w_2phase_h=collect(cellList,'Phase_h','TwoPhase','w')
        h_h_sh=collect(cellList,'Phase_h','Superheated','h_h')
        h_h_2phase=collect(cellList,'Phase_h','TwoPhase','h_h')
        h_h_subcool=collect(cellList,'Phase_h','Subcooled','h_h')
        
        w_subcooled_c=collect(cellList,'Phase_c','Subcooled','w')
        w_superheat_c=collect(cellList,'Phase_c','Superheated','w')
        w_2phase_c=collect(cellList,'Phase_c','TwoPhase','w')
        h_c_sh=collect(cellList,'Phase_c','Superheated','h_c')
        h_c_2phase=collect(cellList,'Phase_c','TwoPhase','h_c')
        h_c_subcool=collect(cellList,'Phase_c','Subcooled','h_c')
        
        if len(w_subcooled_h)>0:
            self.h_subcooled_h=float(sum(np.array(h_h_subcool)*np.array(w_subcooled_h))/sum(w_subcooled_h))
        else:
            self.h_subcooled_h=0
        
        if len(w_2phase_h)>0:
            self.h_2phase_h=float(sum(np.array(h_h_2phase)*np.array(w_2phase_h))/sum(w_2phase_h))
        else:
            self.h_2phase_h=0
            
        if len(w_superheat_h)>0:
            self.h_superheated_h=float(sum(np.array(h_h_sh)*np.array(w_superheat_h))/sum(w_superheat_h))
        else:
            self.h_superheated_h=0
            
        if len(w_subcooled_c)>0:
            self.h_subcooled_c=float(sum(np.array(h_c_subcool)*np.array(w_subcooled_c))/sum(w_subcooled_c))
        else:
            self.h_subcooled_c=0
        
        if len(w_2phase_c)>0:
            self.h_2phase_c=float(sum(np.array(h_c_2phase)*np.array(w_2phase_c))/sum(w_2phase_c))
        else:
            self.h_2phase_c=0
            
        if len(w_superheat_c)>0:
            self.h_superheated_c=float(sum(np.array(h_c_sh)*np.array(w_superheat_c))/sum(w_superheat_c))
        else:
            self.h_superheated_c=0
            
        
        
        self.q_flux=collect(cellList,'Phase_c','TwoPhase','q_flux')
            
        self.Tout_h,self.rhoout_h=TrhoPhase_ph(self.AS_h,self.pin_h,self.hout_h,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0:2]
        self.Tout_c,self.rhoout_c=TrhoPhase_ph(self.AS_c,self.pin_c,self.hout_c,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0:2]
        
        if 'IncompressibleBackend' in AS_c.backend_name():
            AS_c.update(CP.PT_INPUTS, self.pin_c, self.Tout_c)
            self.sout_c=AS_c.smass() #[J/kg-K]
            self.DT_sc_c=1e9
        else:
            AS_c.update(CP.DmassT_INPUTS, self.rhoout_c, self.Tout_c)
            self.sout_c=AS_c.smass() #[J/kg-K]
            #Effective subcooling for both streams
            AS_c.update(CP.QT_INPUTS, 0.0, self.Tbubble_c)
            hsatL=AS_c.hmass() #[J/kg]
            cpsatL=AS_c.cpmass() #[J/kg-K]
            if self.hout_c>hsatL:
                #Outlet is at some quality on cold side
                self.DT_sc_c=-(self.hout_c-hsatL)/cpsatL
            else:
                self.DT_sc_c=self.Tbubble_c-self.Tout_c
        
        if 'IncompressibleBackend' in AS_h.backend_name():
            AS_h.update(CP.PT_INPUTS, self.pin_h, self.Tout_h)
            self.sout_h=AS_h.smass() #[J/kg-K]
            self.DT_sc_h=1e9
        else:
            AS_h.update(CP.DmassT_INPUTS, self.rhoout_h, self.Tout_h)
            self.sout_h=AS_c.smass() #[J/kg-K]
            AS_h.update(CP.QT_INPUTS, 0.0, self.Tbubble_h)
            hsatL=AS_c.hmass() #[J/kg]
            cpsatL=AS_c.cpmass() #[J/kg-K]
            if self.hout_h>hsatL:
                #Outlet is at some quality on hot side
                self.DT_sc_h=-(self.hout_h-hsatL)/cpsatL
            else:
                self.DT_sc_h=self.Tbubble_h-self.Tout_h
        
    def eNTU_CounterFlow(self,Cr,Ntu):
        """
        This function returns the effectiveness for counter flow fluids
        """
        return ((1 - exp(-Ntu * (1 - Cr))) / 
            (1 - Cr * exp(-Ntu * (1 - Cr))))
    
    def _OnePhaseH_OnePhaseC_Qimposed(self,Inputs):
        """
        Single phase on both sides (hot and cold)
        Inputs: dictionary of parameters
        Outputs: dictionary of parameters,
        but mainly w, pressure drop and heat transfer coefficient 
        This function calculate the fraction of heat exchanger 
        that would be required for given thermal duty "w" and DP and h 
        """
        #AbstarctState
        AS_h = self.AS_h
        AS_c = self.AS_c
        
        #Calculate the mean temperature
        Tmean_h=Inputs['Tmean_h']
        Tmean_c=Inputs['Tmean_c']
        #Evaluate heat transfer coefficient for both fluids
        h_h,cp_h,PlateOutput_h=self.PlateHTDP(self.AS_h, Tmean_h, Inputs['pin_h'],self.mdot_h/self.NgapsHot)
        h_c,cp_c,PlateOutput_c=self.PlateHTDP(self.AS_c, Tmean_c, Inputs['pin_c'],self.mdot_c/self.NgapsCold)
        
        #Use cp calculated from delta h/delta T
        cp_h=Inputs['cp_h']
        cp_c=Inputs['cp_c']
        #Evaluate UA [W/K] if entire HX was in this section 
        UA_total=1/(1/(h_h*self.A_h_wetted)+1/(h_c*self.A_c_wetted)+self.PlateThickness/(self.PlateConductivity*(self.A_c_wetted+self.A_h_wetted)/2.))
        #Get Ntu [-]
        C=[cp_c*self.mdot_c,cp_h*self.mdot_h]
        Cmin=min(C)
        Cr=Cmin/max(C)
        
        #Effectiveness [-]
        Q=Inputs['Q']
        Qmax=Cmin*(Inputs['Tin_h']-Inputs['Tin_c'])
        epsilon = Q/Qmax
        
        #Pure counterflow with Cr<1 (Incropera Table 11.4)
        NTU=1/(Cr-1)*log((epsilon-1)/(epsilon*Cr-1))
        
        #Required UA value
        UA_req=Cmin*NTU
        
        #w is required part of heat exchanger for this duty
        w=UA_req/UA_total
        
        #Determine both charge components
        AS_h.update(CP.PT_INPUTS, self.pin_h,Tmean_h)
        rho_h=AS_h.rhomass()#[kg/m^3]
        Charge_h = w * self.V_h * rho_h
        AS_c.update(CP.PT_INPUTS, self.pin_c,Tmean_c)
        rho_c=AS_c.rhomass()#[kg/m^3]
        Charge_c = w * self.V_c * rho_c
        
        
        #Pack outputs
        Outputs={
            'w': w,
            'Tout_h': Inputs['Tin_h']-Q/(self.mdot_h*cp_h),
            'Tout_c': Inputs['Tin_c']+Q/(self.mdot_c*cp_c),
            'Charge_c': Charge_c,
            'Charge_h': Charge_h,
            'DP_h': -PlateOutput_h['DELTAP'],
            'DP_c': -PlateOutput_c['DELTAP'],
            'h_h':h_h,
            'h_c':h_c,
            
        }
        o = Inputs
        o.update(**Outputs)
        return o
        

    def _OnePhaseH_TwoPhaseC_Qimposed(self,Inputs):
        """
        The hot stream is single phase, and the cold stream is evaporating (two phase)
        Inputs: dictionary of parameters
        Outputs: dictionary of parameters, 
        but mainly w, pressure drop and heat transfer coefficient 
        This function calculate the fraction of heat exchanger 
        that would be required for given thermal duty "w" and DP and h
        """
        #AbstarctState
        AS_h = self.AS_h
        AS_c = self.AS_c
        
        #Calculate the mean temperature for the hot single-phase fluid
        h_h,cp_h,PlateOutput_h=self.PlateHTDP(self.AS_h, Inputs['Tmean_h'], Inputs['pin_h'],self.mdot_h/self.NgapsHot)
        #Use cp calculated from delta h/delta T
        cp_h=Inputs['cp_h']
        #Mole mass of refrigerant for Cooper correlation
        M=AS_c.molar_mass()*1000 #[kg/kmol]
        #Reduced pressure for Cooper Correlation
        pcrit_c = AS_c.p_critical() #critical pressure of Ref_c [Pa]
        pstar=Inputs['pin_c']/pcrit_c
        change=999
        w=1
        Q=Inputs['Q']
        """
        The Cooper Pool boiling relationship is a function of the heat flux, 
        therefore the heat flux must be iteratively determined
        
        According to Cleasson J. PhD Thesis "Thermal and Hydraulic Performance
        of Compact Brazed Plate Heat Exchangers Operating a Evaporators in Domestic
        Heat Pumps", KTH, 2004, pp. 98: the saturated nucleate pool boiling 
        correlation by Cooper (1984) works rather well at varying conditions, 
        if multiplied by a factor C=1.5.
        
        To this end, a tuning coefficient, i.e. h_tp_cold_tuninig, is added to the
        Cooper pool boiling correlation.
        
        """
        while abs(change)>1e-6:
            q_flux=Q/(w*self.A_c_wetted)

            if hasattr(self,'Rp'): #check if surface roughness is given
                Rp = self.Rp
            else: #otherwise, use the default surface roughness
                Rp = 1.0

            #Heat transfer coefficient from Cooper Pool Boiling with
            #correction for the two-phase zone of the cold side
            h_c_2phase=self.h_tp_cold_tuning*Cooper_PoolBoiling(pstar,Rp,q_flux,M) 
             
            G=self.mdot_c/self.A_c_flow
            Dh=self.Dh_c
            x=(Inputs['xin_c']+Inputs['xout_c'])/2

            UA_total=1/(1/(h_h*self.A_h_wetted)+1/(h_c_2phase*self.A_c_wetted)+self.PlateThickness/(self.PlateConductivity*self.A_c_wetted))
            C_h=cp_h*self.mdot_h
            
            Qmax=C_h*(Inputs['Tin_h']-Inputs['Tsat_c'])
            epsilon=Q/Qmax
            
            if epsilon>=1.0:
                epsilon=1.0-1e-12 
            NTU=-log(1-epsilon)
            UA_req=NTU*C_h
            
            change=UA_req/UA_total-w
            w=UA_req/UA_total
        
        #Refrigerant charge
        AS_h.update(CP.PT_INPUTS,self.pin_h,Inputs['Tmean_h'])
        rho_h=AS_h.rhomass() #[kg/m^3]
        Charge_h = w * self.V_h * rho_h
        rho_c=TwoPhaseDensity(self.AS_c,Inputs['xin_c'],Inputs['xout_c'],self.Tdew_c,self.Tbubble_c,slipModel='Zivi')
        Charge_c = rho_c * w * self.V_c
        
        #Use Lockhart Martinelli to calculate the pressure drop.  Claesson found good agreement using C parameter of 4.67
        DP_frict_c=LMPressureGradientAvg(Inputs['xin_c'],Inputs['xout_c'],self.AS_c,self.mdot_c/self.A_c_flow,self.Dh_c,self.Tbubble_c,self.Tdew_c,C=4.67)*w*self.Lp
        #Accelerational pressure drop component    
        DP_accel_c=AccelPressureDrop(Inputs['xin_c'],Inputs['xout_c'],self.AS_c,self.mdot_c/self.A_c_flow,self.Tbubble_c,self.Tdew_c)*w*self.Lp
        
        #Pack outputs
        Outputs={
            'w':w,
            'Tout_h': Inputs['Tin_h']-Q/(self.mdot_h*cp_h),
            'Charge_c': Charge_c,
            'Charge_h': Charge_h,
            'DP_h': -PlateOutput_h['DELTAP'],
            'DP_c': DP_frict_c+DP_accel_c,
            'h_h':h_h,
            'h_c':h_c_2phase,
            'q_flux':q_flux
        }
        o = Inputs
        o.update(**Outputs)
        return o
    
    
    def _TwoPhaseH_OnePhaseC_Qimposed(self,Inputs):
        """
        Hot stream is condensing (two phase), cold stream is single phase
        Inputs: dictionary of parameters
        Outputs: dictionary of parameters, 
        but mainly w, pressure drop and heat transfer coefficient 
        This function calculate the fraction of heat exchanger 
        that would be required for given thermal duty "w" and DP and h
        """
        #AbstarctState
        AS_h = self.AS_h
        AS_c = self.AS_c
        
        h_h_2phase=LongoCondensation((Inputs['xout_h']+Inputs['xin_h'])/2,self.mdot_h/self.A_h_flow,self.Dh_h,self.AS_h,self.Tbubble_h,self.Tdew_h);
        h_h_2phase = h_h_2phase*self.h_tp_hot_tuning #correct the convection heat transfer of the two-phase of the hot side
        
        h_c,cp_c,PlateOutput_c=self.PlateHTDP(self.AS_c, Inputs['Tmean_c'], Inputs['pin_c'],self.mdot_c/self.NgapsCold)
        #Use cp calculated from delta h/delta T
        cp_c=Inputs['cp_c']
        UA_total=1/(1/(h_c*self.A_c_wetted)+1/(h_h_2phase*self.A_h_wetted)+self.PlateThickness/(self.PlateConductivity*(self.A_c_wetted+self.A_h_wetted)/2.))
        C_c=cp_c*self.mdot_c
        
        Q=Inputs['Q']
        Qmax=C_c*(Inputs['Tsat_h']-Inputs['Tin_c'])
        epsilon=Q/Qmax
        
        #Cr = 0, so NTU is simply
        try:
            NTU=-log(1-epsilon)
        except:
            pass
        
        UA_req=NTU*C_c
        w=UA_req/UA_total
        
        AS_c.update(CP.PT_INPUTS,self.pin_c,Inputs['Tmean_c'])
        rho_c=AS_c.rhomass() #[kg/m^3]
        Charge_c = w * self.V_c * rho_c
        rho_h=TwoPhaseDensity(self.AS_h,Inputs['xout_h'],Inputs['xin_h'],self.Tdew_h,self.Tbubble_h,slipModel='Zivi')
        Charge_h = w * self.V_h * rho_h
        
        #Use Lockhart Martinelli to calculate the pressure drop.  Claesson found good agreement using C parameter of 4.67
        DP_frict_h=LMPressureGradientAvg(Inputs['xout_h'],Inputs['xin_h'],self.AS_h,self.mdot_h/self.A_h_flow,self.Dh_h,self.Tbubble_h,self.Tdew_h,C=4.67)*w*self.Lp
        #Accelerational pressure drop component    
        DP_accel_h=-AccelPressureDrop(Inputs['xout_h'],Inputs['xin_h'],self.AS_h,self.mdot_h/self.A_h_flow,self.Tbubble_h,self.Tdew_h)*w*self.Lp
        
        #Pack outputs
        Outputs={
            'w': w,
            'Tout_c': Inputs['Tin_c']-Q/(self.mdot_c*cp_c),
            'DP_c': -PlateOutput_c['DELTAP'],
            'DP_h': DP_accel_h+DP_frict_h,
            'Charge_c':Charge_c,
            'Charge_h':Charge_h,
            'h_h':h_h_2phase,
            'h_c':h_c,
        }
        
        return dict(Inputs.items()+Outputs.items())
        
    def Calculate(self):
        """
        Calculate the PHE
        
        For now, only evaporation against glycol is possible
        Cold: Ref
        Hot: Glycol
        """
        #AbstractState
        AS_c = self.AS_c
        if hasattr(self,'MassFrac_c'):
            AS_c.set_mass_fractions([self.MassFrac_c])
        elif hasattr(self, 'VoluFrac_c'):
            AS_c.set_volu_fractions([self.VoluFrac_c])
        
        AS_h = self.AS_h
        if hasattr(self,'MassFrac_h'):
            AS_h.set_mass_fractions([self.MassFrac_h])
        elif hasattr(self, 'VoluFrac_h'):
            AS_h.set_volu_fractions([self.VoluFrac_h])
        
        #set tuning factors to 1 in case not given by user
        if not hasattr(self,'h_tp_cold_tuning'):
            self.h_tp_cold_tuning = 1
        if not hasattr(self,'h_tp_hot_tuning'):
            self.h_tp_hot_tuning = 1
        if not hasattr(self,'DP_hot_tuning'):
            self.DP_hot_tuning = 1
        if not hasattr(self,'DP_cold_tuning'):
            self.DP_cold_tuning = 1    
        
        # Allocate channels between hot and cold streams
        if not hasattr(self,'MoreChannels') or self.MoreChannels not in ['Hot','Cold']:
            raise KeyError("MoreChannels not found, options are 'Hot' or 'Cold'")
        #There are (Nplates - 1) gaps between the plates
        if self.MoreChannels=='Hot':
            #Hot stream gets the extra channel
            self.NgapsHot=(self.Nplates-1)//2+1
            self.NgapsCold=self.Nplates-1-self.NgapsHot
        else:
            #Cold stream gets the extra channel
            self.NgapsCold=(self.Nplates-1)//2+1
            self.NgapsHot=self.Nplates-1-self.NgapsCold
        
        #Saturation temperatures for cold fluid
        if 'IncompressibleBackend' in AS_c.backend_name():
            self.rhosatL_c= None
            self.rhosatV_c= None
            self.Tbubble_c = None
            self.Tdew_c = None
            self.Tsat_c = None
        else:
            AS_c.update(CP.PQ_INPUTS, self.pin_c, 0.0)
            self.Tbubble_c=AS_c.T() #[K]
            self.rhosatL_c=AS_c.rhomass() #[kg/m^3]
            AS_c.update(CP.PQ_INPUTS, self.pin_c, 1.0)
            self.Tdew_c=AS_c.T() #[K]
            self.rhosatV_c=AS_c.rhomass() #[kg/m^3]
            self.Tsat_c=(self.Tbubble_c+self.Tdew_c)/2.0
        
        if 'IncompressibleBackend' in AS_h.backend_name():
            self.Tbubble_h=None
            self.Tdew_h=None
            self.Tsat_h=None
            self.rhosatL_h=None
            self.rhosatV_h=None
        else:
            #Saturation temperatures for hot fluid
            AS_h.update(CP.PQ_INPUTS, self.pin_h, 0.0)
            self.Tbubble_h=AS_h.T() #[K]
            self.rhosatL_h=AS_h.rhomass() #[kg/m^3]
            AS_h.update(CP.PQ_INPUTS, self.pin_h, 1.0)
            self.Tdew_h=AS_h.T() #[K]
            self.rhosatV_h=AS_h.rhomass() #[kg/m^3]
            self.Tsat_h=(self.Tbubble_h+self.Tdew_h)/2.0
        
        #The rest of the inlet states
        self.Tin_h,self.rhoin_h=TrhoPhase_ph(self.AS_h,self.pin_h,self.hin_h,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0:2]
        self.Tin_c,self.rhoin_c=TrhoPhase_ph(self.AS_c,self.pin_c,self.hin_c,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0:2]
        
        if 'IncompressibleBackend' in AS_c.backend_name():
            AS_c.update(CP.PT_INPUTS, self.pin_c, self.Tin_c)
            self.sin_c=AS_c.smass() #[J/kg-K]
        else:
            AS_c.update(CP.DmassT_INPUTS, self.rhoin_c,self.Tin_c)
            self.sin_c=AS_c.smass() #[J/kg-K]
            
        if 'IncompressibleBackend' in AS_h.backend_name():
            AS_h.update(CP.PT_INPUTS, self.pin_h, self.Tin_h)
            self.sin_h=AS_h.smass() #[J/kg-K]
        else:
            AS_h.update(CP.DmassT_INPUTS, self.rhoin_h, self.Tin_h)
            self.sin_h=AS_h.smass() #[J/kg-K]
            
        # Find HT and Delta P on the hot side
        #---------------
        #Mean values for the hot side based on average of inlet temperatures
        HotPlateInputs={
            'PlateAmplitude': self.PlateAmplitude,
            'PlateWavelength' : self.PlateWavelength,
            'InclinationAngle': self.InclinationAngle,
            'Bp': self.Bp,
            'Lp': self.Lp
        }
        HotPlateOutputs=PHE_1phase_hdP(HotPlateInputs,JustGeo=True)
        #There are (Nplates-2) active plates (outer ones don't do anything)
        self.A_h_wetted=HotPlateOutputs['Ap']*(self.Nplates-2)
        self.V_h=HotPlateOutputs['Vchannel']*self.NgapsHot
        self.A_h_flow=HotPlateOutputs['Aflow']*self.NgapsHot
        self.Dh_h=HotPlateOutputs['Dh']
        
        # Find geometric parameters for cold side of plates
        ColdPlateInputs={
            'PlateAmplitude': self.PlateAmplitude,
            'PlateWavelength' : self.PlateWavelength,
            'InclinationAngle': self.InclinationAngle,
            'Bp': self.Bp,
            'Lp': self.Lp
        }
        ColdPlateOutputs=PHE_1phase_hdP(ColdPlateInputs,JustGeo=True)
        #There are (Nplates-2) active plates (outer ones don't do anything)
        self.A_c_wetted=ColdPlateOutputs['Ap']*(self.Nplates-2)
        self.V_c=ColdPlateOutputs['Vchannel']*self.NgapsCold
        self.A_c_flow=ColdPlateOutputs['Aflow']*self.NgapsCold
        self.Dh_c=ColdPlateOutputs['Dh']
        
        #Figure out the limiting rate of heat transfer
        self.Qmax=self.DetermineHTBounds()
        
        def GivenQ(Q):
            """
            In this function, the heat transfer rate is imposed.  Therefore the
            outlet states for both fluids are known, and each element can be solved
            analytically in one shot without any iteration.
            """
            
            EnthalpyList_c,EnthalpyList_h=self.BuildEnthalpyLists(Q)
                
#            #Plot temperature versus enthalpy profiles
#            for i in range(len(EnthalpyList_c)-1):
#                hc=np.linspace(EnthalpyList_c[i],EnthalpyList_c[i+1])
#                Tc=np.zeros_like(hc)
#                for j in range(len(hc)):
#                    Tc[j],r,Ph=TrhoPhase_ph(self.Ref_c,self.pin_c,hc[j],self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)
#                pylab.plot(self.mdot_c*(hc-EnthalpyList_c[0])/1000,Tc,'b')
#                
#            for i in range(len(EnthalpyList_h)-1):
#                hh=np.linspace(EnthalpyList_h[i],EnthalpyList_h[i+1])
#                Th=np.zeros_like(hh)
#                for j in range(len(hh)):
#                    Th[j],r,Ph=TrhoPhase_ph(self.Ref_h,self.pin_h,hh[j],self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)
#                pylab.plot(self.mdot_h*(hh-EnthalpyList_h[0])/1000,Th,'r')
#            pylab.show()
                
#            Ph(self.Ref_h)
#            pylab.plot(np.array(EnthalpyList_h)/1000,self.pin_h/1000*np.ones_like(EnthalpyList_h))
#            pylab.show()
            
            I_h=0
            I_c=0
            wList=[]
            cellList=[]
            while I_h<len(EnthalpyList_h)-1:
                #Try to figure out whether the next phase transition is on the hot or cold side     
                Qbound_h=self.mdot_h*(EnthalpyList_h[I_h+1]-EnthalpyList_h[I_h])
                Qbound_c=self.mdot_c*(EnthalpyList_c[I_c+1]-EnthalpyList_c[I_c])
                if Qbound_h<Qbound_c-1e-9:
                    # Minimum amount of heat transfer is on the hot side,
                    # add another entry to EnthalpyList_c 
                    Qbound=Qbound_h
                    EnthalpyList_c.insert(I_c+1, EnthalpyList_c[I_c]+Qbound/self.mdot_c)
                elif Qbound_h>Qbound_c+1e-9:
                    # Minimum amount of heat transfer is on the cold side,
                    # add another entry to EnthalpyList_h at the interface
                    Qbound=Qbound_c
                    EnthalpyList_h.insert(I_h+1, EnthalpyList_h[I_h]+Qbound/self.mdot_h)
                else:
                    Qbound=Qbound_h
                
                #Figure out the inlet and outlet enthalpy for this cell
                hout_h=EnthalpyList_h[I_h]
                hin_h=hout_h+Qbound/self.mdot_h
                hin_c=EnthalpyList_c[I_c]
                hout_c=hin_c+Qbound/self.mdot_c
                
                # Figure out what combination of phases you have:
                # -------------------------------------------------
                # Hot stream is either single phase or condensing (two phase)
                # Cold stream is either single phase or evaporating (two phase)
                
                #Use midpoint enthalpies to figure out the phase in the cell
                Phase_h=Phase_ph(self.AS_h,self.pin_h,(hin_h+hout_h)/2,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)
                Phase_c=Phase_ph(self.AS_c,self.pin_c,(hin_c+hout_c)/2,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)
                #Determine inlet and outlet temperatures to the cell ([0] gives the first element of the tuple which is temeperature)
                Tin_h=TrhoPhase_ph(self.AS_h,self.pin_h,hin_h,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0]
                Tin_c=TrhoPhase_ph(self.AS_c,self.pin_c,hin_c,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0]
                Tout_h=TrhoPhase_ph(self.AS_h,self.pin_h,hout_h,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0]
                Tout_c=TrhoPhase_ph(self.AS_c,self.pin_c,hout_c,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0]
                
                if Phase_h in ['Subcooled','Superheated'] and Phase_c in ['Subcooled','Superheated']:
                    # Both are single-phase
                    Inputs={
                        'Q':Qbound,
                        'cp_h':(hin_h-hout_h)/(Tin_h-Tout_h),
                        'cp_c':(hin_c-hout_c)/(Tin_c-Tout_c),
                        'Tmean_h':(Tin_h+Tout_h)/2,
                        'Tmean_c':(Tin_c+Tout_c)/2,
                        'Tin_h':Tin_h,
                        'Tin_c':Tin_c,
                        'pin_h':self.pin_h,
                        'pin_c':self.pin_c,
                        'Phase_c':Phase_c,
                        'Phase_h':Phase_h
                    }
                    Outputs=self._OnePhaseH_OnePhaseC_Qimposed(Inputs)
                    wList.append(Outputs['w'])
                    cellList.append(Outputs)
                    if self.Verbosity>6:
                        print('w[1-1]: ', Outputs['w'])
                elif Phase_h=='TwoPhase' and Phase_c in ['Subcooled','Superheated']:
                    # Hot stream is condensing, and cold stream is single-phase (SH or SC)
                    # TODO: bounding state can be saturated state if hot stream is condensing
                    #Must be two-phase so quality is defined
                    xin_h=(hin_h-self.hsatL_h)/(self.hsatV_h-self.hsatL_h)
                    xout_h=(hout_h-self.hsatL_h)/(self.hsatV_h-self.hsatL_h)
                    Inputs={
                        'Q':Qbound,
                        'xin_h':xin_h,
                        'xout_h':xout_h,
                        'Tsat_h':self.Tsat_h,
                        'Tmean_c':(Tin_c+Tout_c)/2,
                        'cp_c':(hin_c-hout_c)/(Tin_c-Tout_c),
                        'Tin_c':Tin_c,
                        'pin_h':self.pin_h,
                        'pin_c':self.pin_c,
                        'Phase_c':Phase_c,
                        'Phase_h':Phase_h
                    }
                    Outputs=self._TwoPhaseH_OnePhaseC_Qimposed(Inputs)
                    if self.Verbosity>6:
                        print('w[2-1]: ', Outputs['w'])
                    wList.append(Outputs['w'])
                    cellList.append(Outputs)
                elif Phase_c=='TwoPhase' and Phase_h in ['Subcooled','Superheated']:
                    # Cold stream is evaporating, and hot stream is single-phase (SH or SC)
                    
                    #Must be two-phase so quality is defined
                    xin_c=(hin_c-self.hsatL_c)/(self.hsatV_c-self.hsatL_c)
                    xout_c=(hout_c-self.hsatL_c)/(self.hsatV_c-self.hsatL_c)
                    
                    Inputs={
                        'Q':Qbound,
                        'xin_c':xin_c,
                        'xout_c':xout_c,
                        'Tsat_c':self.Tsat_c,
                        'cp_h':(hin_h-hout_h)/(Tin_h-Tout_h),
                        'Tmean_h':(Tin_h+Tout_h)/2,
                        'Tin_h':Tin_h,
                        'pin_h':self.pin_h,
                        'pin_c':self.pin_c,
                        'Phase_c':Phase_c,
                        'Phase_h':Phase_h
                    }
                    Outputs=self._OnePhaseH_TwoPhaseC_Qimposed(Inputs)
                    if self.Verbosity>6:
                        print('w[1-2]: ', Outputs['w'])
                    wList.append(Outputs['w'])
                    cellList.append(Outputs)
                    
                I_h+=1
                I_c+=1
                
            self.cellList=cellList
            if self.Verbosity>6:
                print('wsum:', np.sum(wList))
            return np.sum(wList)-1.0
        try:
            brentq(GivenQ,0.0000001*self.Qmax,0.999999999*self.Qmax,xtol=0.000001*self.Qmax)#,xtol=0.000001*self.Qmax) is commented
        except ValueError:
            print(GivenQ(0.0000001*self.Qmax),GivenQ(0.999999999*self.Qmax))
            raise
        # Collect parameters from all the pieces
        self.PostProcess(self.cellList)
        
def WyattPHEHX():
    #Abstract State        
    Ref_c = 'R134a'
    Backend_c = 'HEOS'
    AS_c = CP.AbstractState(Backend_c, Ref_c)
    Ref_h = 'Water'
    Backend_h = 'HEOS'
    AS_h = CP.AbstractState(Backend_h, Ref_h)
    
    Tdew=PropsSI('T','P',962833,'Q',1.0,Ref_c)
    params={
        'AS_c':AS_c,
        'mdot_c':0.073,
        'pin_c':962833,
        'hin_c':PropsSI('H','T',Tdew,'Q',0.0,Ref_c), #[J/kg-K]
        'xin_c':0.0,
        
        'AS_h':AS_h,
        'mdot_h':100.017,
        'pin_h':PropsSI('P','T',115.5+273.15,'Q',1,Ref_h),
        'hin_h':PropsSI('H','T',115.5+273.15,'Q',1,Ref_h), #[J/kg-K]
        
        #Geometric parameters
        'Bp' : 0.119,
        'Lp' : 0.526, #Center-to-center distance between ports
        'Nplates' : 110,
        'PlateAmplitude' : 0.00102, #[m]
        'PlateThickness' : 0.0003, #[m]
        'PlateWavelength' : 0.0066, #[m]
        'InclinationAngle' : pi/3,#[rad]
        'PlateConductivity' : 15.0, #[W/m-K]
        'Rp': 1.0, #[microns] Surface roughness
        'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
    
        'Verbosity':10,
        
        'h_tp_cold_tuning':1,
        'h_tp_hot_tuning':1,
        'DP_hot_tuning':1,
        'DP_cold_tuning':1
    }
    PHE=PHEHXClass(**params)
    PHE.Calculate()

#    pylab.plot(TT,QQ)
#    pylab.show()

def SWEPVariedmdot():
    #Abstract State        
    Ref_c = 'R290'
    Backend_c = 'HEOS'
    AS_c = CP.AbstractState(Backend_c, Ref_c)
    Ref_h = 'Water'
    Backend_h = 'HEOS'
    AS_h = CP.AbstractState(Backend_h, Ref_h)
    
    Tin=8+273.15
    for mdot_h in [0.4176,0.5013,0.6267,0.8357,1.254,2.508]:
        params={
            'AS_c':AS_c,
            'mdot_c':0.03312,
            'pin_c':PropsSI('P','T',Tin,'Q',1.0,Ref_c),
            'hin_c':PropsSI('H','T',Tin,'Q',0.15,Ref_c), #[J/kg-K]
            
            'AS_h':AS_h,
            'mdot_h':mdot_h,
            'pin_h':200000,
            'hin_h':PropsSI('H','T',15+273.15,'P',200000,Ref_h), #[J/kg-K]
            
            #Geometric parameters
            'Bp' : 0.101,
            'Lp' : 0.455, #Center-to-center distance between ports
            'Nplates' : 46,
            'PlateAmplitude' : 0.00102, #[m]
            'PlateThickness' : 0.0003, #[m]
            'PlateWavelength' : 0.00626, #[m]
            'InclinationAngle' : 65/180*pi,#[rad]
            'PlateConductivity' : 15.0, #[W/m-K]
            'Rp': 1.0, #[microns] Surface roughness
            'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
        
            'Verbosity':0,
            
            'h_tp_cold_tuning':1,
            'h_tp_hot_tuning':1,
            'DP_hot_tuning':1,
            'DP_cold_tuning':1
        }
        PHE=PHEHXClass(**params)
        PHE.Calculate()
        print(PHE.Q,',',PHE.h_subcooled_h,',',-PHE.DP_h/1000)

def SamplePHEHX():
    
#    TT=[]
#    QQ=[]
#    for Tin in np.linspace(270,295,50):
#        psat_h=PropsSI('P','T',315.42,'Q',1.0,'R410A')
#        params={
#                'Ref_h':'R410A',
#                'mdot_h':0.075,
#                'pin_h':psat_h,
#                'hin_h':PropsSI('H','T',335.42,'P',psat_h,'R410A'), #*1000
#                
#                'Ref_c':'Water',
#                'mdot_c':0.38,
#                'pin_c':300000,        #pin_c in Pa
#                'hin_c':PropsSI('H','T',282.52,'P',300000,'Water'), #*1000
#                
#                #Geometric parameters
#                'Bp' : 0.11,
#                'Lp' : 0.300, #Center-to-center distance between ports
#                'Nplates' : 16,
#                'PlateAmplitude' : 0.002, #[m]
#                'PlateThickness' : 0.0003, #[m]
#                'PlateConductivity' : 15.0, #[W/m-K]
#                'PlateWavelength' : 0.0126, #[m]
#                'InclinationAngle' : pi/3,#[rad]
#                'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
#                
#                'Verbosity':1
#        }
#        PHE=PHEHXClass(**params)
#        PHE.Calculate()
#        QQ.append(PHE.Q)
#        TT.append(Tin)
#        print(PHE.Q, PHE.Q/PHE.Qmax)
#    pylab.plot(TT,QQ)
#    pylab.show()
    
    #Abstract State        
    Ref_c = 'R290'
    Backend_c = 'HEOS'
    AS_c = CP.AbstractState(Backend_c, Ref_c)
    Ref_h = 'Water'
    Backend_h = 'HEOS'
    AS_h = CP.AbstractState(Backend_h, Ref_h)
    
    TT=[]
    QQ=[]
    Q1=[]
    Tin=8+273.15
#    for Tin in np.linspace(275,287,101):##
    for mdot_h in [0.4176,0.5013,0.6267,0.8357,1.254,2.508]:
        params={
            'AS_c':AS_c,
            'mdot_c':0.03312,
            'pin_c':PropsSI('P','T',Tin,'Q',1.0,Ref_c),
            'hin_c':PropsSI('H','T',Tin,'Q',0.15,Ref_c), #[J/kg-K]
            
            'AS_h':AS_h,
            'mdot_h':mdot_h,
            'pin_h':200000,
            'hin_h':PropsSI('H','T',15+273.15,'P',200000,Ref_h), #[J/kg-K]
            
            #Geometric parameters
            'Bp' : 0.101,
            'Lp' : 0.455, #Center-to-center distance between ports
            'Nplates' : 46,
            'PlateAmplitude' : 0.00102, #[m]
            'PlateThickness' : 0.0003, #[m]
            'PlateWavelength' : 0.00626, #[m]
            'InclinationAngle' : 65/180*pi,#[rad]
            'PlateConductivity' : 15.0, #[W/m-K]
            'Rp': 1.0, #[microns] Surface roughness
            'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
            'Verbosity':0,
            
            'h_tp_cold_tuning':1,
            'h_tp_hot_tuning':1,
            'DP_hot_tuning':1,
            'DP_cold_tuning':1
        }
        PHE=PHEHXClass(**params)
        PHE.Calculate()
        TT.append(Tin)
        QQ.append(PHE.h_2phase_c)#PHE.Q/PHE.Qmax)
        Q1.append(PHE.q_flux)#w_2phase_c)#PHE.Q/PHE.Qmax)
#        print(PHE.Q/PHE.Qmax,PHE.Q)
        print(PHE.Q,',',PHE.h_subcooled_h,',',-PHE.DP_h/1000)
        print(PHE.OutputList())
    print(TT)
    print(QQ)
    print(Q1)
#    pylab.plot(TT,QQ)
#    pylab.show()
        
if __name__=='__main__':
    #SamplePHEHX()
    #WyattPHEHX()
    
    SWEPVariedmdot()