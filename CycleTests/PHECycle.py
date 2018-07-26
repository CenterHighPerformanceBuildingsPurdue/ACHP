from __future__ import division, print_function, absolute_import
from ACHP.Compressor import CompressorClass  #Compressor
from ACHP.Condenser import CondenserClass    #Condenser
from ACHP.MicroChannelCondenser import MicroCondenserClass
from ACHP.Evaporator import EvaporatorClass  #Evaporator
from ACHP.CoolingCoil import CoolingCoilClass #Cooling Coil
from ACHP.CoaxialHX import CoaxialHXClass #Coaxial internal heat exchanger
from ACHP.Pump import PumpClass # Secondary loop pump class
from ACHP.PHEHX import PHEHXClass
from scipy.optimize import brent, fsolve 
#^^ fsolve - roots (multiple variables); brent - root of one variable fct
#from CoolProp.CoolProp import PropsSI #,Tsat        #refrigerant properties
from ACHP.FinCorrelations import WavyLouveredFins,FinInputs     #fin correlations
from ACHP.MicroFinCorrelations import MicroFinInputs
import numpy as np                  #NumPy is fundamental scientific package                                   
import CoolProp as CP

class SecondaryCycleClass():
    def __init__(self):
        """
        Load up the necessary sub-structures to be filled with
        the code that follows
        """
        self.Compressor=CompressorClass()
        self.CoolingCoil=CoolingCoilClass()
        self.CoolingCoil.Fins=FinInputs()
        self.PHEHX=PHEHXClass()
        self.Pump=PumpClass()
    
    def Update(self):
        '''
        Update cycle class with selected HX type
        Update cyle class with Abstract State
        '''
        if self.EvapSolver == 'Moving-Boundary':
            if self.EvapType == 'Fin-tube':
                self.Evaporator=EvaporatorClass()
                self.Evaporator.Fins=FinInputs()
            elif self.EvapType == 'Micro-channel':
                raise
            else:
                raise
        elif self.EvapSolver == 'Finite-Element':
            raise
        else:
            raise    
        
        if self.CondSolver == 'Moving-Boundary':
            if self.CondType == 'Fin-tube':
                self.Condenser=CondenserClass()
                self.Condenser.Fins=FinInputs()
            elif self.CondType == 'Micro-channel':
                self.Condenser=MicroCondenserClass()
                self.Condenser.Fins=MicroFinInputs()
            else:
                raise
        elif self.CondSolver == 'Finite-Element':
            raise
        else:
            raise
        
        #Abstract State   
        self.AS = CP.AbstractState(self.Backend, self.Ref)
        if hasattr(self,'MassFrac'):
            self.AS.set_mass_fractions([self.MassFrac])
        elif hasattr(self,'VoluFrac'):
            self.AS.set_volu_fractions([self.VoluFrac])
        #Abstract State for SecLoopFluid  
        self.AS_SLF = CP.AbstractState(self.Backend_SLF, self.SecLoopFluid)
        if hasattr(self,'MassFrac_SLF'):
            self.AS_SLF.set_mass_fractions([self.MassFrac_SLF])
        elif hasattr(self,'VoluFrac_SLF'):
            self.AS_SLF.set_volu_fractions([self.VoluFrac_SLF])
                        
    def Calculate(self,DT_evap,DT_cond,Tin_IHX):
        """
        Inputs are differences in temperature [K] between HX air inlet temperature 
        and the dew temperature for the heat exchanger.
        
        Required Inputs:
            DT_evap: 
                Difference in temperature [K] between cooling coil air inlet temperature and refrigerant dew temperature
            DT_cond:
                Difference in temperature [K] between condenser air inlet temperature and refrigeant dew temperature
            Tin_IHX:
                Inlet "glycol" temperature to IHX   
        """
        if self.Verbosity>1:
            print ('Inputs: DTevap %7.4f DTcond %7.4f fT_IHX %7.4f'%(DT_evap,DT_cond,Tin_IHX))
        
        #AbstractState
        AS = self.AS
        #AbstractState for SecLoopFluid
        AS_SLF = self.AS_SLF
        
        """
        The coldest the glycol entering the cooling coil could be would be the 
        """
        self.Tdew_cond=self.Condenser.Fins.Air.Tdb+DT_cond
        self.Tdew_evap=self.CoolingCoil.Fins.Air.Tdb-DT_evap
        AS.update(CP.QT_INPUTS,1.0,self.Tdew_cond)
        psat_cond=AS.p() #[Pa]
        AS.update(CP.QT_INPUTS,1.0,self.Tdew_evap)
        psat_evap=AS.p() #[Pa]
        AS.update(CP.PQ_INPUTS,psat_evap,0.0)
        self.Tbubble_evap=AS.T() #[K]    
        
        params={               #dictionary -> key:value, e.g. 'key':2345,
            'pin_r': psat_evap,   
            'pout_r': psat_cond,
            'Tin_r': self.Tdew_evap+self.Compressor.DT_sh,
            'AS': AS,
        }
        self.Compressor.Update(**params)
        self.Compressor.Calculate()
        
        params={
            'mdot_r': self.Compressor.mdot_r,
            'Tin_r': self.Compressor.Tout_r,
            'psat_r': psat_cond,
            'AS': AS,
        }
        self.Condenser.Update(**params)
        self.Condenser.Calculate()
        
        AS.update(CP.QT_INPUTS,0.0,self.Tbubble_evap)
        hL=AS.hmass() #[J/kg]
        AS.update(CP.QT_INPUTS,1.0,self.Tdew_evap)
        hV=AS.hmass() #[J/kg]
        xin_r=(self.Condenser.hout_r-hL)/(hV-hL)
        AS_SLF.update(CP.PT_INPUTS,300000,Tin_IHX)
        h_in = AS_SLF.hmass() #[J/kg]
        params={
            'mdot_h': self.Pump.mdot_g,
            'hin_h': h_in,
            'hin_c': self.Condenser.hout_r,
            'mdot_c': self.Compressor.mdot_r,
            'pin_c': psat_evap,
            'xin_c': xin_r,
            'AS_c': AS,
            'AS_h': AS_SLF,
        }
        self.PHEHX.Update(**params)
        self.PHEHX.Calculate()
        
        #Now run CoolingCoil to predict inlet glycol temperature to IHX
        params={
            'mdot_g': self.Pump.mdot_g,
            'Tin_g': self.PHEHX.Tout_h,
            'AS_g': AS_SLF,
        }
        self.CoolingCoil.Update(**params)
        self.CoolingCoil.Calculate()
        
        params={
            'DP_g': self.PHEHX.DP_h+self.CoolingCoil.DP_g,
            'Tin_g': self.CoolingCoil.Tout_g,
            'AS_g': AS_SLF,
        }
        self.Pump.Update(**params)
        self.Pump.Calculate()
        
        self.Charge=self.Condenser.Charge+self.PHEHX.Charge_c
        self.EnergyBalance=self.Compressor.CycleEnergyIn+self.Condenser.Q+self.PHEHX.Q
        
        resid=np.zeros((3))
        resid[0]=self.EnergyBalance
        
        if self.ImposedVariable=='Subcooling':
            resid[1]=self.Condenser.DT_sc-self.DT_sc_target    
        elif self.ImposedVariable=='Charge':
            resid[1]=self.Charge-self.Charge_target
        resid[2]=self.PHEHX.Q-self.CoolingCoil.Q
        
        if self.Verbosity>1:
            print ('Qres % 12.6e Resid2: % 12.6e ResSL %10.4f Charge %10.4f SC: %8.4f' %(resid[0],resid[1],resid[2],self.Charge,self.Condenser.DT_sc))
            
        self.Capacity=self.CoolingCoil.Capacity
        self.COP=self.CoolingCoil.Q/self.Compressor.W
        self.COSP=self.CoolingCoil.Capacity/(self.Compressor.W+self.Pump.W+self.CoolingCoil.Fins.Air.FanPower+self.Condenser.Fins.Air.FanPower)
        self.SHR=self.CoolingCoil.SHR
        return resid

    def PreconditionedSolve(self):
        # Takes the place of a lambda function since lambda functions do not
        # bubble the error up properly
        def OBJECTIVE(x):
            return self.Calculate(x[0],x[1],x[2])
        
        # Increase DT_evap in increments of 0.5 K until the system solves, and
        # two steps cause a sign change in cycle energy balance
        # Should give a starting point quite close to the "right" soluton
        oldDT_evap=None
        oldQresid=None
        Tin_IHX=None
        for DT_evap in np.linspace(12,25,21):
            try:
                #Run the cooler to get a good starting inlet temperature for IHX
                self.CoolingCoil.mdot_g=self.Pump.mdot_g
                self.CoolingCoil.Tin_g=self.CoolingCoil.Fins.Air.Tdb-DT_evap
                self.CoolingCoil.AS_g = self.AS_SLF
                self.CoolingCoil.Calculate()
                Tin_IHX=self.CoolingCoil.Tout_g
                resid=self.Calculate(DT_evap,8,Tin_IHX)
            except Exception,e:
                if self.Verbosity>1:
                    print ('Failed: ',e.__str__())
                raise
                pass
            else:
                #Not set yet, first one that works
                if oldQresid==None:
                    oldQresid=resid[0]
                    oldDT_evap=DT_evap
                #Has been set, and sign changes, use average of old and new DT
                elif oldQresid*resid[0]<0:
                    DT_evap=(DT_evap+oldDT_evap)/2
                    break
        
        #Run the Newton-Raphson solver to solve the system
        x=fsolve(OBJECTIVE,[DT_evap,20,Tin_IHX])
        
        if self.Verbosity>1:
            print ('Capacity: ', self.Capacity)
            print ('COP: ',self.COP)
            print ('COP (w/ both fans): ',self.COSP)
            print ('SHR: ',self.SHR)
        
if __name__=='__main__': 
    for i in range(1):
        #########################################################################
        #################     SECONDARY CYCLE INITIALIZATION    #################
        #########################################################################
        
        ## Here we load parameters that are not a function of operating conditions
        ## They are primarily geometric parameters
        
        Cycle=SecondaryCycleClass()
        
        #--------------------------------------
        #--------------------------------------
        #       Cycle parameters
        #--------------------------------------
        #--------------------------------------
        Cycle.Verbosity = 0 #the idea here is to have different levels of debug output
        Cycle.ImposedVariable = 'Subcooling'                                                           
        Cycle.DT_sc_target = 7.0
        Cycle.Charge_target = 3.3
        Cycle.Ref='R410A'
        Cycle.Backend='TTSE&HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
        Cycle.Oil = 'POE32'
        Cycle.shell_pressure = 'low-pressure'
        Cycle.SecLoopFluid = 'Water'
        Cycle.Backend_SLF = 'INCOMP' #backend of SecLoopFluid
        Cycle.EvapSolver = 'Moving-Boundary' #choose the type of Evaporator solver scheme (for now only 'Moving-Boundary')
        Cycle.EvapType = 'Fin-tube' #if EvapSolver = 'Moving-Boundary', choose the type of evaporator (for now only 'Fin-tube')
        Cycle.CondSolver = 'Moving-Boundary' #choose the type of Condenser solver scheme (for now only 'Moving-Boundary')
        Cycle.CondType = 'Fin-tube' #if CondSolver = 'Moving-Boundary', choose the type of condenser ('Fin-tube' or 'Micro-channel')
        Cycle.Update()
        
        #--------------------------------------
        #--------------------------------------
        #       Compressor parameters
        #--------------------------------------
        #--------------------------------------
    
        #A few 3 ton cooling capacity compressor maps 
        if Cycle.Ref=='R410A':
            M=[217.3163128,5.094492028,-0.593170311,4.38E-02,-2.14E-02,1.04E-02,7.90E-05,-5.73E-05,1.79E-04,-8.08E-05]#compressor map coefficients
            P=[-561.3615705,-15.62601841,46.92506685,-0.217949552,0.435062616,-0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03]
        
        params={
            'M':M,
            'P':P,
            'Ref':Cycle.Ref,                                                              #refrigerant
            'Oil':Cycle.Oil, #Compressor lubricant oil
            'shell_pressure':Cycle.shell_pressure, #Compressor shell pressure
            'V_oil_sump':0, #Volume of oil in the sump
            'fp':0.15, #Fraction of electrical power lost as heat to ambient            #shell heat loss
            'Vdot_ratio': 1.0, #Displacement Scale factor                               #up- or downsize compressor (1=original)
            'Verbosity': 0, # How verbose should the debugging statements be [0 to 10]
            'DT_sh': 5. #Superheat at inlet to compressor [K]                           
          }
        Cycle.Compressor.Update(**params)
        
        #--------------------------------------
        #--------------------------------------
        #      Condenser parameters 
        #      -> see GUI for illustration/units
        #--------------------------------------
        #--------------------------------------                                                                                                     
        Cycle.Condenser.Fins.Tubes.NTubes_per_bank=41  #number of tubes per bank=row
        Cycle.Condenser.Fins.Tubes.Nbank=1             #number of baks/rows
        Cycle.Condenser.Fins.Tubes.Ncircuits=5             #number of baks/rows
        Cycle.Condenser.Fins.Tubes.Ltube=2.286
        Cycle.Condenser.Fins.Tubes.OD=0.007
        Cycle.Condenser.Fins.Tubes.ID=0.0063904
        Cycle.Condenser.Fins.Tubes.Pl=0.0191  #distance between center of tubes in flow direction                                                
        Cycle.Condenser.Fins.Tubes.Pt=0.0222  #distance between center of tubes orthogonal to flow direction
        Cycle.Condenser.Fins.Tubes.kw=237     #wall thermal conductivity (i.e pipe material)
        
        Cycle.Condenser.Fins.Fins.FPI=25      #Number of fins per inch
        Cycle.Condenser.Fins.Fins.Pd=0.001    #2* amplitude of wavy fin
        Cycle.Condenser.Fins.Fins.xf=0.001    #1/2 period of fin
        Cycle.Condenser.Fins.Fins.t=0.00011   #Thickness of fin material
        Cycle.Condenser.Fins.Fins.k_fin=237   #Thermal conductivity of fin material
        
        Cycle.Condenser.Fins.Air.Vdot_ha=1.7934 #rated volumetric flowrate
        Cycle.Condenser.Fins.Air.Tmean=308.15   
        Cycle.Condenser.Fins.Air.Tdb=308.15     #Dry Bulb Temperature
        Cycle.Condenser.Fins.Air.p=101325       #Air pressure in Pa
        Cycle.Condenser.Fins.Air.RH=0.51        #Relative Humidity
        Cycle.Condenser.Fins.Air.RHmean=0.51    
        Cycle.Condenser.Fins.Air.FanPower=160
        
        params={
            'FinsType': 'WavyLouveredFins',                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
            'Verbosity': 0
        }
        Cycle.Condenser.Update(**params)
            
        #--------------------------------------
        #--------------------------------------
        #           Cooling Coil
        #           -> see Condenser and GUI for explanations
        #--------------------------------------
        #--------------------------------------
        Cycle.CoolingCoil.Fins.Tubes.NTubes_per_bank=32
        Cycle.CoolingCoil.Fins.Tubes.Nbank=4
        Cycle.CoolingCoil.Fins.Tubes.Ncircuits=6
        Cycle.CoolingCoil.Fins.Tubes.Ltube=0.452
        Cycle.CoolingCoil.Fins.Tubes.OD=0.009525
        Cycle.CoolingCoil.Fins.Tubes.ID=0.0089154
        Cycle.CoolingCoil.Fins.Tubes.Pl=0.0254
        Cycle.CoolingCoil.Fins.Tubes.Pt=0.0219964
        Cycle.CoolingCoil.Fins.Tubes.kw=237     #wall thermal conductivity (i.e pipe material)
        
        Cycle.CoolingCoil.Fins.Fins.FPI=14.5
        Cycle.CoolingCoil.Fins.Fins.Pd=0.001
        Cycle.CoolingCoil.Fins.Fins.xf=0.001
        Cycle.CoolingCoil.Fins.Fins.t=0.00011
        Cycle.CoolingCoil.Fins.Fins.k_fin=237
        
        Cycle.CoolingCoil.Fins.Air.Vdot_ha=0.5663
        Cycle.CoolingCoil.Fins.Air.Tmean=299.8
        Cycle.CoolingCoil.Fins.Air.Tdb=299.8
        Cycle.CoolingCoil.Fins.Air.p=101325
        Cycle.CoolingCoil.Fins.Air.RH=0.51
        Cycle.CoolingCoil.Fins.Air.RHmean=0.51
        Cycle.CoolingCoil.Fins.Air.FanPower=438
               
        params={
            'pin_g': 200000,
            'FinsType': 'WavyLouveredFins',                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'                                                    #pin_g in Pa
            'Verbosity':0
        }
        Cycle.CoolingCoil.Update(**params)
        
        params={
            'pin_h':300000,
            'Verbosity':0,
            
            #Geometric parameters
            'Bp' : 0.117,
            'Lp' : 0.300, #Center-to-center distance between ports
            'Nplates' : 46,
            'PlateAmplitude' : 0.001, #[m]
            'PlateThickness' : 0.0003, #[m]
            'PlateConductivity' : 15.0, #[W/m-K]
            'Rp': 1.0, #[microns] Surface roughness
            'PlateWavelength' : 0.0126, #[m]
            'InclinationAngle' : 3.14159/2.3,#[rad]
            'MoreChannels' : 'Hot' #Which stream gets the extra channel, 'Hot' or 'Cold'
            }
        Cycle.PHEHX.Update(**params)
        
        params={
            'eta':0.5,  #Pump+motor efficiency
            'mdot_g':0.38, #Flow Rate kg/s
            'pin_g':300000,
            'Verbosity':0,
            }
        Cycle.Pump.Update(**params)
        
        #Now solve
        Cycle.PreconditionedSolve()
        
        print (Cycle.Pump.DP_g,Cycle.Pump.W)
