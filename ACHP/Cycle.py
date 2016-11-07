from __future__ import division
import sys
from Compressor import CompressorClass  #Compressor
from Condenser import CondenserClass    #Condenser
from Evaporator import EvaporatorClass  #Evaporator
from CoolingCoil import CoolingCoilClass #Cooling Coil
from MultiCircuitEvaporator import MultiCircuitEvaporatorClass
from CoaxialHX import CoaxialHXClass #Coaxial internal heat exchanger
from PHEHX import PHEHXClass #Plate-Heat-Exchanger 
from LineSet import LineSetClass #Line set class
from Pump import PumpClass # Secondary loop pump class
from scipy.optimize import brentq, fsolve,newton 
#^^ fsolve - roots (multiple variables); brent - root of one variable fct

import CoolProp as CP
from CoolProp.Plots import PropertyPlot
from FinCorrelations import FinInputs     #fin correlations
import numpy as np                  #NumPy is fundamental scientific package
from Correlations import TrhoPhase_ph            
from Solvers import MultiDimNewtRaph, Broyden
from Preconditioners import DXPreconditioner,SecondaryLoopPreconditioner
        
class SecondaryCycleClass():
    def __init__(self):
        """
        Load up the necessary sub-structures to be filled with
        the code that follows
        """
        self.Compressor=CompressorClass()
        
        #Outdoor coil is a Condenser in cooling mode and evaporator in heating mode
        self.Condenser=CondenserClass()
        self.Condenser.Fins=FinInputs()
        self.Evaporator=EvaporatorClass()
        self.Evaporator.Fins=FinInputs()
        
        self.CoolingCoil=CoolingCoilClass()
        self.CoolingCoil.Fins=FinInputs()
        self.Pump=PumpClass()
        #Add both types of internal heat exchangers
        self.CoaxialIHX=CoaxialHXClass()
        self.PHEIHX=PHEHXClass()
        self.LineSetSupply=LineSetClass()
        self.LineSetReturn=LineSetClass()
        #Make IHX an empty class for holding parameters common to PHE and Coaxial IHX
        class struct: pass
        self.IHX=struct()
        
    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        return [
            ('Charge','kg',self.Charge),
            ('Condenser Subcooling','K',self.DT_sc),
            ('Primary Ref.','-',self.Ref),
            ('Secondary Ref.','-',self.SecLoopFluid),
            ('Imposed Variable','-',self.ImposedVariable),
            ('IHX Type','-',self.IHXType),
            ('COP','-',self.COP),
            ('COSP','-',self.COSP),
            ('Net Capacity','W',self.CoolingCoil.Capacity),
            ('Net Power','W',self.Power),
            ('SHR','-',self.SHR),
            ('Condensation temp (dew)','K',self.Tdew_cond),
            ('Evaporation temp (dew)','K',self.Tdew_evap),
         ]

    def Calculate(self,DT_evap,DT_cond,Tin_CC):
        """
        Inputs are differences in temperature [K] between HX air inlet temperature 
        and the dew temperature for the heat exchanger.
        
        Required Inputs:
            DT_evap: 
                Difference in temperature [K] between cooling coil air inlet temperature and refrigerant dew temperature
            DT_cond:
                Difference in temperature [K] between condenser air inlet temperature and refrigerant dew temperature
            Tin_CC:
                Inlet "glycol" temperature to line set feeding cooling coil   
        """
        if self.Verbosity>1:
            print 'Inputs: DTevap %7.4f DTcond %7.4f fT_IHX %7.4f'%(DT_evap,DT_cond,Tin_CC) 
        
        #AbstractState
        if hasattr(self,'Backend'): #check if backend is given
            AS = CP.AbstractState(self.Backend, self.Ref)
            if hasattr(self,'MassFrac'):
                AS.set_mass_fractions([self.MassFrac])
        else: #otherwise, use the defualt backend
            AS = CP.AbstractState('HEOS', self.Ref)
            self.Backend = 'HEOS'
        self.AS = AS
        #AbstractState for SecLoopFluid
        if hasattr(self,'Backend_SLF'): #check if backend_SLF is given
            AS_SLF = CP.AbstractState(self.Backend_SLF, self.SecLoopFluid)
            if hasattr(self,'MassFrac_SLF'):
                AS_SLF.set_mass_fractions([self.MassFrac_SLF])
        else: #otherwise, use the defualt backend
            AS_SLF = CP.AbstractState('HEOS', self.SecLoopFluid)
            self.Backend_SLF = 'HEOS'
        self.AS_SLF = AS_SLF
        
        #Store the values to save on computation for later
        self.DT_evap=DT_evap
        self.DT_cond=DT_cond
        self.Tin_CC=Tin_CC
        
        #If the user doesn't include the Mode, set it to Air Conditioning 
        if not hasattr(self,'Mode'):
            self.Mode='AC'
            
        if self.Mode=='AC':
            self.Tdew_cond=self.Condenser.Fins.Air.Tdb+DT_cond
            self.Tdew_evap=self.CoolingCoil.Fins.Air.Tdb-DT_evap
        elif self.Mode=='HP':
            self.Tdew_cond=Tin_CC+DT_cond
            self.Tdew_evap=self.Evaporator.Fins.Air.Tdb-DT_evap
        else:
            raise ValueError('Mode must be AC or HP')
        
        #Evaporator and condeser saturation pressures
        AS.update(CP.QT_INPUTS,1.0,self.Tdew_cond)
        psat_cond=AS.p() #[Pa]
        AS.update(CP.QT_INPUTS,1.0,self.Tdew_evap)
        psat_evap=AS.p() #[Pa]
        
        #Evaporator and condeser bubble temepratures
        AS.update(CP.PQ_INPUTS,psat_evap,0.0)
        self.Tbubble_evap=AS.T() #[K]
        AS.update(CP.PQ_INPUTS,psat_cond,0.0)
        self.Tbubble_cond=AS.T() #[K]
        
        #Cycle solver for 'AC' mode
        if self.Mode=='AC':
            
            params={               #dictionary -> key:value, e.g. 'key':2345,
                'pin_r': psat_evap+self.DP_low,   
                'pout_r': psat_cond-self.DP_high,
                'Tin_r': self.Tdew_evap+self.PHEIHX.DT_sh, # TrhoPhase_ph(self.Ref,psat_evap,self.LineSetReturn.hout,self.Tbubble_evap,self.Tdew_evap)[0],
                'Ref':  self.Ref,
                'Backend': self.Backend
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()
            
            params={
                'mdot_r': self.Compressor.mdot_r,
                'Tin_r': self.Compressor.Tout_r,
                'psat_r': psat_cond,
                'Ref': self.Ref,
                'Backend': self.Backend
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()
            
            #Inlet enthalpy to LineSetSupply
            AS_SLF.update(CP.PT_INPUTS,self.Pump.pin_g,Tin_CC)
            h_in_LineSetSupply = AS_SLF.hmass() #[J/kg]
            params={
                'mdot': self.Pump.mdot_g,
                'hin': h_in_LineSetSupply,
            }
            self.LineSetSupply.Update(**params)
            self.LineSetSupply.Calculate()
            
            #Now run CoolingCoil to predict inlet glycol temperature to IHX
            params={
                'mdot_g': self.Pump.mdot_g,
                'Tin_g': self.LineSetSupply.Tout,
            }
            self.CoolingCoil.Update(**params)
            self.CoolingCoil.Calculate()
            
            #Inlet enthalpy to LineSetReturn
            AS_SLF.update(CP.PT_INPUTS,self.Pump.pin_g,self.CoolingCoil.Tout_g)
            h_in_LineSetReturn = AS_SLF.hmass() #[J/kg]
            params={
                'mdot': self.Pump.mdot_g,
                'hin': h_in_LineSetReturn
            }
            self.LineSetReturn.Update(**params)
            self.LineSetReturn.Calculate()
            
            if self.IHXType=='Coaxial':
                params={
                    'mdot_g': self.Pump.mdot_g,
                    'Tin_g': self.CoolingCoil.Tout_g,
                    'pin_r': psat_evap,
                    'hin_r': self.Condenser.hout_r,
                    'Ref_r': self.Ref,
                    'Backend_r': self.Backend,
                    'mdot_r': self.Compressor.mdot_r,
                }
                self.CoaxialIHX.Update(**params)
                self.CoaxialIHX.Calculate()
                self.IHX.Charge_r=self.CoaxialIHX.Charge_r
                self.IHX.Q=self.CoaxialIHX.Q
                self.IHX.Tout_g=self.CoaxialIHX.Tout_g
                self.IHX.DP_g=self.CoaxialIHX.DP_g
                self.IHX.hout_r=self.CoaxialIHX.hout_r
                self.IHX.DP_r=self.CoaxialIHX.DP_r
                if hasattr(self,'PHEIHX'):
                    del self.PHEIHX
            elif self.IHXType=='PHE':
                #Inlet enthalpy to PHEIHX
                AS_SLF.update(CP.PT_INPUTS,self.PHEIHX.pin_h,self.CoolingCoil.Tout_g)
                h_in_PHEIHX = AS_SLF.hmass() #[J/kg]
                params={
                    'mdot_h': self.Pump.mdot_g,
                    'hin_h': h_in_PHEIHX,
                    'mdot_c': self.Compressor.mdot_r,
                    'pin_c': psat_evap,
                    'hin_c': self.Condenser.hout_r,
                }
                self.PHEIHX.Update(**params)
                self.PHEIHX.Calculate()
                self.IHX.Charge_r=self.PHEIHX.Charge_c
                self.IHX.Q=self.PHEIHX.Q
                self.IHX.Tout_g=self.PHEIHX.Tout_h
                self.IHX.DP_g=self.PHEIHX.DP_h
                self.IHX.DP_r=self.PHEIHX.DP_c
                self.IHX.hout_r=self.PHEIHX.hout_c
                if hasattr(self,'CoaxialIHX'):
                    del self.CoaxialIHX
            
            params={
                'DP_g': self.IHX.DP_g+self.CoolingCoil.DP_g+self.LineSetSupply.DP+self.LineSetReturn.DP,
                'Tin_g': self.CoolingCoil.Tout_g
            }
            self.Pump.Update(**params)
            self.Pump.Calculate()
            
            self.Charge=self.Condenser.Charge+self.IHX.Charge_r
            self.EnergyBalance=self.Compressor.CycleEnergyIn+self.Condenser.Q+self.IHX.Q
            #Calculate properties:
            AS.update(CP.QT_INPUTS,0.0,self.Tbubble_cond)
            h_L = AS.hmass() #[J/kg]
            cp_L = AS.cpmass() #[J/kg-K]
            AS.update(CP.PT_INPUTS,psat_cond,self.Tbubble_cond-self.DT_sc_target)
            h_target = AS.hmass() #[J/kg]
            self.DT_sc=(h_L - self.Condenser.hout_r)/cp_L
            deltaH_sc=self.Compressor.mdot_r*(h_L-h_target)
            
#            ## Plot a p-h plot
#            plot = PropertyPlot('HEOS::' + self.Ref, 'PH', unit_system='KSI')
#            plot.plot([self.Compressor.hin_r/1000,self.Compressor.hout_r/1000,self.Condenser.hout_r/1000,self.PHEIHX.hin_c/1000,self.Compressor.hin_r/1000],[psat_evap/1000,psat_cond/1000,psat_cond/1000,psat_evap/1000,psat_evap/1000])
#            plot.show()
            
            resid=np.zeros((3))
            resid[0]=self.Compressor.mdot_r*(self.Compressor.hin_r-self.IHX.hout_r)
            
            if self.ImposedVariable=='Subcooling':
                resid[1]=self.Condenser.DT_sc-self.DT_sc_target    
            elif self.ImposedVariable=='Charge':
                resid[1]=self.Charge-self.Charge_target
#            resid[2]=self.IHX.Q-self.CoolingCoil.Q+self.Pump.W

            self.residSL=self.IHX.Q-self.CoolingCoil.Q+self.Pump.W+self.LineSetSupply.Q+self.LineSetReturn.Q
            resid[2]=self.residSL
            
            
            
            if self.Verbosity>7:
                print 'Wcomp % 12.6e Qcond: % 12.6e QPHE %10.4f ' %(self.Compressor.W,self.Condenser.Q,self.IHX.Q)
            if self.Verbosity>1:
                print 'Qres % 12.6e Resid2: % 12.6e ResSL %10.4f  Charge %10.4f SC: %8.4f' %(resid[0],resid[1],self.residSL,self.Charge,self.Condenser.DT_sc)
                
            self.Capacity=self.CoolingCoil.Capacity
            self.COP=self.CoolingCoil.Q/self.Compressor.W
            self.COSP=self.CoolingCoil.Capacity/(self.Compressor.W+self.Pump.W+self.CoolingCoil.Fins.Air.FanPower+self.Condenser.Fins.Air.FanPower)
            self.SHR=self.CoolingCoil.SHR
            self.Power=self.Compressor.W+self.Pump.W+self.CoolingCoil.Fins.Air.FanPower+self.Condenser.Fins.Air.FanPower
            self.DP_high_Model=self.Condenser.DP_r #[Pa]
            self.DP_low_Model=self.IHX.DP_r #[Pa]
        
        #Cycle solver for 'HP' mode
        elif self.Mode=='HP':
            if psat_evap+self.DP_low<0:
                raise ValueError('Compressor inlet pressure less than zero ['+str(psat_evap+self.DP_low)+' Pa] - is low side pressure drop too high?')
            params={               #dictionary -> key:value, e.g. 'key':2345,
                'pin_r': psat_evap+self.DP_low,   
                'pout_r': psat_cond-self.DP_high,
                'Tin_r': self.Tdew_evap+self.Evaporator.DT_sh, # TrhoPhase_ph(self.Ref,psat_evap,self.LineSetReturn.hout,self.Tbubble_evap,self.Tdew_evap)[0],
                'Ref':  self.Ref,
                'Backend': self.Backend
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()
            
            #Inlet enthalpy to LineSetSupply
            AS_SLF.update(CP.PT_INPUTS,self.Pump.pin_g,Tin_CC)
            h_in_LineSetSupply = AS_SLF.hmass() #[J/kg]
            params={
                'mdot': self.Pump.mdot_g,
                'hin': h_in_LineSetSupply
            }
            self.LineSetSupply.Update(**params)
            self.LineSetSupply.Calculate()
            
            #Now run CoolingCoil to predict inlet glycol temperature to IHX
            params={
                'mdot_g': self.Pump.mdot_g,
                'Tin_g': self.LineSetSupply.Tout,
            }
            self.CoolingCoil.Update(**params)
            self.CoolingCoil.Calculate()
            
            #Inlet enthalpy to LineSetReturn
            AS_SLF.update(CP.PT_INPUTS,self.Pump.pin_g,self.CoolingCoil.Tout_g)
            h_in_LineSetReturn = AS_SLF.hmass() #[J/kg]
            params={
                'mdot': self.Pump.mdot_g,
                'hin': h_in_LineSetReturn
            }
            self.LineSetReturn.Update(**params)
            self.LineSetReturn.Calculate()
            
            #Inlet enthalpy to PHEIHX
            AS_SLF.update(CP.PT_INPUTS,self.PHEIHX.pin_c,self.LineSetReturn.Tout)
            h_in_PHEIHX = AS_SLF.hmass() #[J/kg]
            params={
                'mdot_h': self.Compressor.mdot_r,
                'hin_h': self.Compressor.hout_r,
                'pin_h': psat_cond,
                'mdot_c': self.Pump.mdot_g,
                'hin_c': h_in_PHEIHX
            }
            self.PHEIHX.Update(**params)
            self.PHEIHX.Calculate()
            
            params={
                'mdot_r': self.Compressor.mdot_r,
                'psat_r': psat_evap,
                'hin_r': self.PHEIHX.hout_h,
                'Ref': self.Ref,
                'Backend': self.Backend
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()
            
            params={
                'DP_g': self.PHEIHX.DP_c+self.CoolingCoil.DP_g,
                'Tin_g': self.CoolingCoil.Tout_g
            }
            self.Pump.Update(**params)
            self.Pump.Calculate()
            
            self.Charge=self.Evaporator.Charge+self.PHEIHX.Charge_h
            
            #Calculate properties:
            AS.update(CP.QT_INPUTS,0.0,self.Tbubble_cond)
            h_L = AS.hmass() #[J/kg]
            cp_L = AS.cpmass() #[J/kg-K]
            AS.update(CP.PT_INPUTS,psat_cond,self.Tbubble_cond-self.DT_sc_target)
            h_target = AS.hmass() #[J/kg]
            #Calculate an effective subcooling amount by deltah/cp_satL
            #Can be positive or negative (negative is quality at outlet
            self.DT_sc=self.PHEIHX.DT_sc_h#(PropsSI('H','T',self.Tbubble_cond,'Q',0,self.Ref)-self.PHEIHX.hout_h)/(PropsSI('C','T',self.Tbubble_cond,'Q',0,self.Ref)) #*1000 #*1000
            deltaH_sc=self.Compressor.mdot_r*(h_L-h_target)
            
            resid=np.zeros((3))
            resid[0]=self.Compressor.mdot_r*(self.Compressor.hin_r-self.Evaporator.hout_r)
            
            if self.ImposedVariable=='Subcooling':
                resid[1]=self.DT_sc-self.DT_sc_target
            elif self.ImposedVariable=='Charge':
                resid[1]=self.Charge-self.Charge_target
            self.residSL=self.PHEIHX.Q+self.CoolingCoil.Q+self.Pump.W+self.LineSetSupply.Q+self.LineSetReturn.Q
            resid[2]=self.residSL
            
            if self.Verbosity>1:
                print 'Qres % 12.6e Resid2: % 12.6e ResSL %10.4f Charge %10.4f SC: %8.4f' %(resid[0],resid[1],self.residSL,self.Charge,self.DT_sc)
                
            self.Capacity=-self.CoolingCoil.Q+self.CoolingCoil.Fins.Air.FanPower
            self.COP=-self.CoolingCoil.Q/self.Compressor.W
            self.COSP=self.Capacity/(self.Compressor.W+self.Pump.W+self.CoolingCoil.Fins.Air.FanPower+self.Evaporator.Fins.Air.FanPower)
            self.Power=self.Compressor.W+self.Pump.W+self.CoolingCoil.Fins.Air.FanPower+self.Evaporator.Fins.Air.FanPower
            self.SHR=-self.CoolingCoil.SHR
            self.DP_high_Model=self.PHEIHX.DP_h #[Pa]
            self.DP_low_Model=self.Evaporator.DP_r #[Pa]
        self.DT_evap=DT_evap
        self.DT_cond=DT_cond   
        return resid

    def PreconditionedSolve(self,PrecondValues=None):
        '''
        PrecondValues = dictionary of values DT_evap, DT_cond and Tin_CC
        '''
        
        def OBJECTIVE(x):
            """
            Takes the place of a lambda function since lambda functions do not bubble error properly
            """
            return self.Calculate(x[0],x[1],x[2])
        def OBJECTIVE2(x,Tin):
            """
            Takes the place of a lambda function since lambda functions do not bubble error properly
            """
            return self.Calculate(x[0],x[1],Tin)
        def OBJECTIVE_SL(Tin_CC):
            """
            Objective function for the inner loop of the vapor compression system
            
            Using the MultiDimNewtRaph function will re-evaluate the Jacobian at 
            every step.  Slower, but more robust since the solution surfaces aren't
            smooth enough
            
            Note: This function is not currently used!
            """
            x=MultiDimNewtRaph(OBJECTIVE2,[self.DT_evap,self.DT_cond],args=(Tin_CC,))
   
            # Update the guess values for Delta Ts starting 
            # at the third step (after at least one update away 
            # from the boundaries)
            if self.OBJ_SL_counter>=0:
                self.DT_evap=x[0]
                self.DT_cond=x[1]
                pass
            self.OBJ_SL_counter+=1
            return self.residSL
            
        def PrintDPs():
            print 'DP_LP :: Input:',self.DP_low,'Pa / Model calc:',self.DP_low_Model,'Pa'
            print 'DP_HP :: Input:',self.DP_high,'Pa / Model calc:',self.DP_high_Model,'Pa'   
        
        #Some variables need to be initialized
        self.DP_low=0 #The actual low-side pressure drop to be used in Pa
        self.DP_high=0 #The actual low-side pressure drop to be used in Pa
        self.OBJ_SL_counter=0
        
        #Run the preconditioner to get guess values for the temperatures
        if PrecondValues is None:
            self.DT_evap,self.DT_cond,Tin_CC=SecondaryLoopPreconditioner(self)
        else:
            self.DT_evap=PrecondValues['DT_evap']
            self.DT_cond=PrecondValues['DT_cond']
            Tin_CC=PrecondValues['Tin_CC']
            
        
        #Remove the other, non-used IHX class if found
        if self.IHXType=='PHE':
            if hasattr(self,'CoaxialIHX'):
                del self.CoaxialIHX
        else:
            if hasattr(self,'PHEIHX'):
                del self.PHEIHX
            
        #Remove the condenser if in heating mode and condenser found
        if self.Mode=='HP':
            if hasattr(self,'Condenser'):
                del self.Condenser
        iter=1
        max_error_DP=999
        #Outer loop with a more relaxed convergence criterion
        while max_error_DP>0.5:
            iter_inner=1
            #Inner loop to determine pressure drop for high and low sides
            while max_error_DP>0.05 and iter_inner<10:
                
                #Run to calculate the pressure drop as starting point
                OBJECTIVE([self.DT_evap,self.DT_cond,Tin_CC])
                
                #Calculate the max error
                max_error_DP=max([abs(self.DP_low_Model-self.DP_low),abs(self.DP_high_Model-self.DP_high)])
                
                if self.Verbosity>0:
                    PrintDPs()
                    print 'Max pressure drop error [inner loop] is',max_error_DP,'Pa'
                    
                #Update the pressure drop terms
                self.DP_low=self.DP_low_Model #/1000
                self.DP_high=self.DP_high_Model #/1000
                
                iter_inner+=1
                
            if self.Verbosity > 0:
                print "Done with the inner loop on pressure drop"
                
            # Use Newton-Raphson solver
            (self.DT_evap,self.DT_cond,Tin_CC)=MultiDimNewtRaph(OBJECTIVE,[self.DT_evap,self.DT_cond,Tin_CC],dx=0.1)
            
            #Calculate the error
            max_error_DP=max([abs(self.DP_low_Model-self.DP_low),abs(self.DP_high_Model-self.DP_high)])
            
            if self.Verbosity>0:
                PrintDPs()    
                print 'Max pressure drop error [outer loop] is',max_error_DP,'Pa'
        
        if self.Verbosity>1:
            print 'Capacity: ', self.Capacity
            print 'COP: ',self.COP
            print 'COP (w/ both fans): ',self.COSP
            print 'SHR: ',self.SHR
        return
        
class DXCycleClass():
    def __init__(self):
        """
        Load up the necessary sub-structures to be filled with
        the code that follows
        """
        self.Compressor=CompressorClass()
        self.Condenser=CondenserClass()
        self.Condenser.Fins=FinInputs()
        self.Evaporator=EvaporatorClass()
        self.Evaporator.Fins=FinInputs()
        self.LineSetSupply=LineSetClass()
        self.LineSetReturn=LineSetClass()
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
            ('Charge','kg',self.Charge),
            ('Condensation temp (dew)','K',self.Tdew_cond),
            ('Evaporation temp (dew)','K',self.Tdew_evap),
            ('Condenser Subcooling','K',self.DT_sc),
            ('Primary Ref.','-',self.Ref),
            ('COP','-',self.COP),
            ('COSP','-',self.COSP),
            ('Net Capacity','W',self.Capacity),
            ('Net Power','W',self.Power),
            ('SHR','-',self.SHR),
            ('Imposed Variable','-',self.ImposedVariable),
         ]
        for i in range(0,len(Output_List_default)):                             #append default parameters to output list
            Output_List.append(Output_List_default[i])
        return Output_List

    def Calculate(self,DT_evap,DT_cond):
        """
        Inputs are differences in temperature [K] between HX air inlet temperature 
        and the dew temperature for the heat exchanger.
        
        Required Inputs:
            DT_evap: 
                Difference in temperature [K] between evaporator air inlet temperature and refrigerant dew temperature
            DT_cond:
                Difference in temperature [K] between condenser air inlet temperature and refrigeant dew temperature 
        """
        if self.Verbosity>1:
            print 'DTevap %7.4f DTcond %7.4f,' %(DT_evap,DT_cond)
        
        #AbstractState
        if hasattr(self,'Backend'): #check if backend is given
            AS = CP.AbstractState(self.Backend, self.Ref)
        else: #otherwise, use the defualt backend
            AS = CP.AbstractState('HEOS', self.Ref)
            self.Backend = 'HEOS'
        self.AS = AS
        
        #Condenser and evaporator dew temperature (guess)
        Tdew_cond=self.Condenser.Fins.Air.Tdb+DT_cond#the values (Tin_a,..) come from line 128ff
        Tdew_evap=self.Evaporator.Fins.Air.Tdb-DT_evap
        #Condenser and evaporator saturation pressures
        AS.update(CP.QT_INPUTS,1.0,Tdew_cond)
        psat_cond=AS.p() #[Pa]
        AS.update(CP.QT_INPUTS,1.0,Tdew_evap)
        psat_evap=AS.p() #[Pa]
        #evaporator bubble temparture
        AS.update(CP.PQ_INPUTS,psat_evap,0.0)
        Tbubble_evap=AS.T() #[T]
        
        self.Tdew_cond=Tdew_cond
        self.Tdew_evap=Tdew_evap
        
        #If the user doesn't include the Mode, fail
        assert hasattr(self,'Mode')
        
        #Cycle Solver in 'AC' model
        if self.Mode=='AC':
            if not hasattr(self.Compressor,'mdot_r') or self.Compressor.mdot_r<0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value 
                # for the mass flow rate for the line set 
                params={               #dictionary -> key:value, e.g. 'key':2345,
                    'pin_r': psat_evap,   
                    'pout_r': psat_cond,
                    'Tin_r': Tdew_evap+self.Evaporator.DT_sh,
                    'Ref':  self.Ref,
                    'Backend': self.Backend,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()
            
            #Calculate inlet enthalpy 
            AS.update(CP.PT_INPUTS,psat_evap,Tdew_evap+self.Evaporator.DT_sh)
            h_in = AS.hmass() #[J/kg]
            params={
                'pin': psat_evap,
                'hin': h_in,
                'mdot': self.Compressor.mdot_r,
                'Ref':  self.Ref,
                'Backend': self.Backend
            }
            self.LineSetReturn.Update(**params)
            self.LineSetReturn.Calculate()
            
            params={               #dictionary -> key:value, e.g. 'key':2345,
                'pin_r': psat_evap-self.DP_low,   
                'pout_r': psat_cond+self.DP_high,
                'Tin_r': TrhoPhase_ph(self.AS,psat_evap,self.LineSetReturn.hout,Tbubble_evap,Tdew_evap)[0],
                'Ref':  self.Ref,
                'Backend': self.Backend
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()
            if self.Verbosity>1:
                print 'Comp DP L H',self.DP_low,self.DP_high
            
            params={
                'mdot_r': self.Compressor.mdot_r,
                'Tin_r': self.Compressor.Tout_r,
                'psat_r': psat_cond,
                'Ref': self.Ref,
                'Backend': self.Backend
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()
            
            params={
                'pin':psat_cond,
                'hin':self.Condenser.hout_r,
                'mdot':self.Compressor.mdot_r,
                'Ref':self.Ref,
                'Backend': self.Backend
            }
            self.LineSetSupply.Update(**params)
            self.LineSetSupply.Calculate()
            
            params={
                'mdot_r': self.Compressor.mdot_r,
                'psat_r': psat_evap,
                'hin_r': self.LineSetSupply.hout,
                'Ref': self.Ref,
                'Backend': self.Backend
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()
            
            self.Charge=self.Condenser.Charge+self.Evaporator.Charge+self.LineSetSupply.Charge+self.LineSetReturn.Charge
            self.EnergyBalance=self.Compressor.CycleEnergyIn+self.Condenser.Q+self.Evaporator.Q
            
            resid=np.zeros((2))
            self.DP_HighPressure=self.Condenser.DP_r+self.LineSetSupply.DP
            self.DP_LowPressure=self.Evaporator.DP_r+self.LineSetReturn.DP
            resid[0]=self.Compressor.mdot_r*(self.LineSetReturn.hin-self.Evaporator.hout_r)
            
            if self.ImposedVariable=='Subcooling':
                resid[1]=self.Condenser.DT_sc-self.DT_sc_target    
            elif self.ImposedVariable=='Charge':
                resid[1]=self.Charge-self.Charge_target
            
            if self.Verbosity>1:
                print resid
            
            self.Capacity=self.Evaporator.Capacity
            self.Power=self.Compressor.W+self.Evaporator.Fins.Air.FanPower+self.Condenser.Fins.Air.FanPower
            self.COP=self.Evaporator.Q/self.Compressor.W
            self.COSP=self.Evaporator.Capacity/self.Power
            self.SHR=self.Evaporator.SHR
            self.DT_sc=self.Condenser.DT_sc
        
        #Cycle Solver in 'HP' model
        elif self.Mode=='HP':            
            params={               #dictionary -> key:value, e.g. 'key':2345,
                'pin_r': psat_evap-self.DP_low,   
                'pout_r': psat_cond+self.DP_high,
                'Tin_r': Tdew_evap+self.Evaporator.DT_sh,
                'Ref':  self.Ref,
                'Backend': self.Backend
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()
            
            params={
                'pin': psat_cond,
                'hin': self.Compressor.hout_r,
                'mdot': self.Compressor.mdot_r,
                'Ref':  self.Ref,
                'Backend': self.Backend
            }
            self.LineSetSupply.Update(**params)
            self.LineSetSupply.Calculate()
            
            params={
                'mdot_r': self.Compressor.mdot_r,
                'Tin_r': self.Compressor.Tout_r,
                'psat_r': psat_cond,
                'Ref': self.Ref,
                'Backend': self.Backend
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()
            
            params={
                'pin': psat_cond,
                'hin': self.Condenser.hout_r,
                'mdot': self.Compressor.mdot_r,
                'Ref': self.Ref,
                'Backend': self.Backend
            }
            self.LineSetReturn.Update(**params)
            self.LineSetReturn.Calculate()
            
            params={
                'mdot_r': self.Compressor.mdot_r,
                'psat_r': psat_evap,
                'hin_r': self.LineSetReturn.hout,
                'Ref': self.Ref,
                'Backend': self.Backend
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()
            
            self.Charge=self.Condenser.Charge+self.Evaporator.Charge+self.LineSetSupply.Charge+self.LineSetReturn.Charge
            self.EnergyBalance=self.Compressor.CycleEnergyIn+self.Condenser.Q+self.Evaporator.Q
            
            resid=np.zeros((2))
            resid[0]=self.Compressor.mdot_r*(self.Compressor.hin_r-self.Evaporator.hout_r)
            
            if self.ImposedVariable=='Subcooling':
                resid[1]=self.Condenser.DT_sc-self.DT_sc_target    
            elif self.ImposedVariable=='Charge':
                resid[1]=self.Charge-self.Charge_target
            
            self.Capacity=-self.Condenser.Q+self.Condenser.Fins.Air.FanPower
            self.DT_sc=self.Condenser.DT_sc
            self.Power=self.Compressor.W+self.Evaporator.Fins.Air.FanPower+self.Condenser.Fins.Air.FanPower
            self.COP=-self.Condenser.Q/self.Compressor.W
            self.COSP=self.Capacity/self.Power
            self.SHR=self.Evaporator.SHR
            self.DP_HighPressure=self.Condenser.DP_r+self.LineSetSupply.DP
            self.DP_LowPressure=self.Evaporator.DP_r+self.LineSetReturn.DP
        else:
            ValueError("DX Cycle mode must be 'AC', or 'HP'")
        if self.Verbosity>1:
            print 'DTevap %7.4f DTcond %7.4f Qres % 12.6e DTsc: % 12.6e Charge %10.4f SC: %8.4f' %(DT_evap,DT_cond,resid[0],resid[1],self.Charge,self.Condenser.DT_sc)
        self.DT_evap=DT_evap
        self.DT_cond=DT_cond
        return resid
    
    def PreconditionedSolve(self):
        """
        Solver that will precondition by trying a range of DeltaT until the model
        can solve, then will kick into 2-D Newton Raphson solve
        
        The two input variables for the system solver are the differences in 
        temperature between the inlet air temperature of the heat exchanger and the
        dew temperature of the refrigerant.  This is important for refrigerant blends
        with temperature glide during constant-pressure evaporation or condensation.
        Good examples of common working fluid with glide would be R404A or R410A.
        """
        def OBJECTIVE_DXCycle(x):
            """
            A wrapper function to convert input vector for fsolve to the proper form for the solver
            """
            try:
                resids=self.Calculate(DT_evap=float(x[0]),DT_cond=float(x[1]))#,DP_low=float(x[2]),DP_high=float(x[3]))
            except ValueError:
                raise
            return resids
        
        # Use the preconditioner to determine a reasonably good starting guess
        DT_evap_init,DT_cond_init=DXPreconditioner(self)

        GoodRun=False
        while GoodRun==False:
            try:
                self.DP_low=0
                self.DP_high=0
                DP_converged=False        
                while DP_converged==False:
                    #Actually run the Newton-Raphson solver to get the solution
                    x=Broyden(OBJECTIVE_DXCycle,[DT_evap_init,DT_cond_init])
                    delta_low=abs(self.DP_low-abs(self.DP_LowPressure))
                    delta_high=abs(self.DP_high-abs(self.DP_HighPressure))
                    self.DP_low=abs(self.DP_LowPressure)
                    self.DP_high=abs(self.DP_HighPressure)
                    #Update the guess values based on last converged values
                    DT_evap_init=self.DT_evap
                    DT_cond_init=self.DT_cond
                    if delta_low<1 and delta_high<1:
                        DP_converged=True
                    if self.Verbosity>4:
                        print self.DP_HighPressure,self.DP_LowPressure,'DPHP'
                    GoodRun=True
            except AttributeError:
                # This will be a fatal error !! Should never have attribute error
                raise 
            except:
                print "--------------  Exception Caught ---------------- " 
                print "Error of type",sys.exc_info()[0]," is: " + sys.exc_info()[1].message
                raise
        
        if self.Verbosity>0:
            print 'Capacity: ', self.Capacity
            print 'COP: ',self.COP
            print 'COP (w/ both fans): ',self.COSP
            print 'SHR: ',self.SHR
            print 'UA_r_evap',self.Evaporator.UA_r
            print 'UA_a_evap',self.Evaporator.UA_a
            print 'UA_r_cond',self.Condenser.UA_r
            print 'UA_a_cond',self.Condenser.UA_a