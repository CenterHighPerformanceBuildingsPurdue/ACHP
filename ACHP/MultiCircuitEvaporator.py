from __future__ import division, print_function, absolute_import
from math import floor,ceil
from CoolProp.CoolProp import PropsSI
from ACHP.FinCorrelations import FinInputs
from ACHP.Evaporator import EvaporatorClass
import numpy as np
import matplotlib.pyplot as plt
import copy
import CoolProp
from scipy.optimize import newton
from ACHP.ACHPTools import Write2CSV
import CoolProp as CP

#MultiCircuitEvaporator inherits things from the Evaporator base class
class MultiCircuitEvaporatorClass(EvaporatorClass):
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
        Output_list=[]
        num_evaps= int(self.Fins.Tubes.Ncircuits)
        Output_list_i=[]
        for i in range(num_evaps):  #Generate output list with all evaporators
            Output_list_i.append([('Volumetric flow rate'+' '+str(i),'m^3/s',self.Evaps[i].Fins.Air.Vdot_ha),
                                    ('Inlet Dry bulb temp'+' '+str(i),'K',self.Evaps[i].Fins.Air.Tdb),
                                        ('Inlet Air pressure'+' '+str(i),'Pa',self.Evaps[i].Fins.Air.p),
                                        ('Inlet Air Relative Humidity'+' '+str(i),'-',self.Evaps[i].Fins.Air.RH),
                                        ('Tubes per bank'+' '+str(i),'-',self.Evaps[i].Fins.Tubes.NTubes_per_bank),
                                        ('Number of banks'+' '+str(i),'-',self.Evaps[i].Fins.Tubes.Nbank),
                                        ('Number circuits'+' '+str(i),'-',self.Evaps[i].Fins.Tubes.Ncircuits),
                                        ('Length of tube'+' '+str(i),'m',self.Evaps[i].Fins.Tubes.Ltube),
                                        ('Tube OD'+' '+str(i),'m',self.Evaps[i].Fins.Tubes.OD),
                                        ('Tube ID'+' '+str(i),'m',self.Evaps[i].Fins.Tubes.ID),
                                        ('Tube Long. Pitch'+' '+str(i),'m',self.Evaps[i].Fins.Tubes.Pl),
                                        ('Tube Transverse Pitch'+' '+str(i),'m',self.Evaps[i].Fins.Tubes.Pt),
                                        ('Outlet superheat'+' '+str(i),'K',self.Evaps[i].Tout_r-self.Evaps[i].Tdew_r),
                                        ('Fins per inch'+' '+str(i),'1/in',self.Evaps[i].Fins.Fins.FPI),
                                        ('Fin waviness pd'+' '+str(i),'m',self.Evaps[i].Fins.Fins.Pd),
                                        ('Fin waviness xf'+' '+str(i),'m',self.Evaps[i].Fins.Fins.xf),
                                        ('Fin thickness'+' '+str(i),'m',self.Evaps[i].Fins.Fins.t),
                                        ('Fin Conductivity'+' '+str(i),'W/m-K',self.Evaps[i].Fins.Fins.k_fin),
                                        ('Fins Type'+' '+str(i),'-',self.Evaps[i].FinsType),
                                        ('Refrigerant flowrate'+' '+str(i),'kg/s',self.Evaps[i].mdot_r),
                                        ('Q Total'+' '+str(i),'W',self.Evaps[i].Q),
                                        ('Q Superheat'+' '+str(i),'W',self.Evaps[i].Q_superheat),
                                        ('Q Two-Phase'+' '+str(i),'W',self.Evaps[i].Q_2phase),
                                        ('Inlet ref. temp'+' '+str(i),'K',self.Evaps[i].Tin_r),
                                        ('Outlet ref. temp'+' '+str(i),'K',self.Evaps[i].Tout_r),
                                        ('Outlet air temp'+' '+str(i),'K',self.Evaps[i].Tout_a),
                                        ('Pressure Drop Total'+' '+str(i),'Pa',self.Evaps[i].DP_r),
                                        ('Pressure Drop Superheat'+' '+str(i),'Pa',self.Evaps[i].DP_r_superheat),
                                        ('Pressure Drop Two-Phase'+' '+str(i),'Pa',self.Evaps[i].DP_r_2phase),
                                        ('Charge Total'+' '+str(i),'kg',self.Evaps[i].Charge),
                                        ('Charge Superheat'+' '+str(i),'kg',self.Evaps[i].Charge_superheat),
                                        ('Charge Two-Phase'+' '+str(i),'kg',self.Evaps[i].Charge_2phase),
                                        ('Mean HTC Superheat'+' '+str(i),'W/m^2-K',self.Evaps[i].h_r_superheat),
                                        ('Mean HTC Two-phase'+' '+str(i),'W/m^2-K',self.Evaps[i].h_r_2phase),
                                        ('Wetted Area Fraction Superheat'+' '+str(i),'-',self.Evaps[i].w_superheat),
                                        ('Wetted Area Fraction Two-phase'+' '+str(i),'-',self.Evaps[i].w_2phase),
                                        ('Mean Air HTC'+' '+str(i),'W/m^2-K',self.Evaps[i].Fins.h_a),
                                        ('Surface Effectiveness'+' '+str(i),'-',self.Evaps[i].Fins.eta_a),
                                        ('Air-side area (fin+tubes)'+' '+str(i),'m^2',self.Evaps[i].Fins.A_a),
                                        ('Mass Flow rate dry Air'+' '+str(i),'kg/s',self.Evaps[i].Fins.mdot_da),
                                        ('Mass Flow rate humid Air'+' '+str(i),'kg/s',self.Evaps[i].Fins.mdot_ha),
                                        ('Pressure Drop Air-side'+' '+str(i),'Pa',self.Evaps[i].Fins.dP_a),
                                        ('Sensible Heat Ratio'+' '+str(i),'-',self.Evaps[i].SHR)])
        for i in range(len(Output_list_i[0])):  #sort output list, such that corresponding values are next to each other
            sumsi=0    #need sums and avgs
            for n in range(num_evaps):
                Output_list.append(Output_list_i[n][i])
                if type(Output_list_i[n][i][2]) is not type("string"):
                    sumsi+=Output_list_i[n][i][2] #sum up values
                    if n == num_evaps-1:
                        Output_list.append((Output_list_i[n][i][0][:-2]+" sum",Output_list_i[n][i][1],sumsi))
                        Output_list.append((Output_list_i[n][i][0][:-2]+" avg",Output_list_i[n][i][1],sumsi/num_evaps))
        Output_List_tot=[]
        #append optional parameters, if applicable
        if hasattr(self,'TestName'):
            Output_List_tot.append(('Name','N/A',self.TestName)) 
        if hasattr(self,'TestDescription'):
            Output_List_tot.append(('Description','N/A',self.TestDescription))
        if hasattr(self,'TestDetails'):
            Output_List_tot.append(('Details','N/A',self.TestDetails))
        for i in range(0,len(Output_list)):      #append default parameters to output list
            Output_List_tot.append(Output_list[i])
        return Output_List_tot

    def Calculate(self):
        #AbstractState
        AS = self.AS
        
        #Check that the length of lists of mdot_r and FinsTubes.Air.Vdot_ha 
        #match the number of circuits or are all equal to 1 (standard evap)
        Ncircuits=int(self.Fins.Tubes.Ncircuits)

        # Make Ncircuits copies of evaporator classes defined 
        #  by the inputs to the MCE superclass
        EvapDict=self.__dict__
        #Delete AbstractState class from dictionary to avoid problem with deepcopy. However, AS can still be used later (see line 218)
        del EvapDict['AS'] 
        
        self.Evaps=[]
        for i in range(Ncircuits):
            # Make a deep copy to break all the links between the Fins structs
            # of each of the evaporator instances
            ED=copy.deepcopy(EvapDict)
            #Create new evaporator instanciated with new deep copied dictionary
            E=EvaporatorClass(**ED)
            #Add to list of evaporators
            self.Evaps.append(E)
        
        #Upcast single values to lists, and convert numpy arrays to lists
        self.Fins.Air.Vdot_ha=np.atleast_1d(self.Fins.Air.Vdot_ha).tolist()
        self.mdot_r=np.atleast_1d(self.mdot_r).tolist()
        
        if Ncircuits != len(self.mdot_r) and len(self.mdot_r)>1:
            print ("Problem with length of vector for mdot_r for MCE")
        else:
            if len(self.mdot_r)==1: #Single value passed in for mdot_r
                if hasattr(self,'mdot_r_coeffs'):
                    if len(self.mdot_r_coeffs)!=Ncircuits:
                        raise AttributeError("Size of array mdot_r_coeffs: "+str(len(self.mdot_r_coeffs))+" does not equal Ncircuits: "+str(Ncircuits))
                    elif abs(np.sum(self.mdot_r_coeffs)-1)>=100*np.finfo(float).eps:
                        raise AttributeError("mdot_r_coeffs must sum to 1.0.  Sum *100000 is: "+str(100000*np.sum(self.mdot_r_coeffs)))
                    else:
                        # A vector of weighting factors multiplying the total mass flow rate is provided with the right length
                        for i in range(Ncircuits):
                            self.Evaps[i].mdot_r=self.mdot_r[-1]*self.mdot_r_coeffs[i]
                else:
                    # Refrigerant flow is evenly distributed between circuits,
                    # give each evaporator an equal portion of the refrigerant
                    for i in range(Ncircuits):
                        self.Evaps[i].mdot_r=self.mdot_r[-1]/Ncircuits
            else:
                for i in range(Ncircuits):
                    self.Evaps[i].mdot_r=self.mdot_r[i]
                   
        # Deal with the possibility that the quality might be varied among circuits
        if hasattr(self,'mdot_v_coeffs'):
            if len(self.mdot_v_coeffs)!=Ncircuits:
                raise AttributeError("Size of array mdot_v_coeffs: "+str(len(self.mdot_v_coeffs))+" does not equal Ncircuits: "+str(Ncircuits))
            elif abs(np.sum(self.mdot_v_coeffs)-1)>=10*np.finfo(float).eps:
                raise AttributeError("mdot_v_coeffs must sum to 1.0.  Sum is: "+str(np.sum(self.mdot_v_coeffs)))
            else:
                AS.update(CP.PQ_INPUTS, self.psat_r, 0.0) 
                hsatL=AS.hmass() #[J/kg]
                AS.update(CP.PQ_INPUTS, self.psat_r, 1.0)
                hsatV=AS.hmass() #[J/kg]
                x_inlet=(self.hin_r-hsatL)/(hsatV-hsatL)
                mdot_v=x_inlet*sum(self.mdot_r)
                for i in range(Ncircuits):
                    mdot_v_i=self.mdot_v_coeffs[i]*mdot_v
                    x_i=mdot_v_i/self.Evaps[i].mdot_r
                    AS.update(CP.PQ_INPUTS, self.psat_r, x_i)
                    self.Evaps[i].hin_r=AS.hmass() #[J/kg]
                 
        #For backwards compatibility, if the coefficients are provided in the FinInputs class, copy them to the base class
        if hasattr(self.Fins.Air,'Vdot_ha_coeffs'):
            self.Vdot_ha_coeffs=self.Fins.Air.Vdot_ha_coeffs
            print ("Warning: please put the vector Vdot_ha_coeffs in the base MCE class, accesssed as MCE.Vdot_ha_coeffs")
            
        if Ncircuits !=len(self.Fins.Air.Vdot_ha) and len(self.Fins.Air.Vdot_ha)>1:
            print ("Problem with length of vector for Vdot_ha for MCE")
        else:
            if len(self.Fins.Air.Vdot_ha)==1:
                if hasattr(self,'Vdot_ha_coeffs'):
                    if len(self.Vdot_ha_coeffs)!=Ncircuits:
                        raise AttributeError("Size of array Vdot_ha_coeffs: "+str(len(self.Vdot_ha_coeffs))+" does not equal Ncircuits: "+str(Ncircuits))
                    elif abs(np.sum(self.Vdot_ha_coeffs)-1)>=10*np.finfo(float).eps:
                        raise AttributeError("Vdot_ha_coeffs does not sum to 1.0! Sum is: "+str(np.sum(self.Vdot_ha_coeffs)))                        
                    else:                        
                        # A vector of factors multiplying the total volume flow rate is provided
                        for i in range(Ncircuits):
                            self.Evaps[i].Fins.Air.Vdot_ha=self.Fins.Air.Vdot_ha[-1]*self.Vdot_ha_coeffs[i]
                else:
                    # Air flow is evenly distributed between circuits,
                    # give each circuit an equal portion of the air flow
                    for i in range(Ncircuits):
                        self.Evaps[i].Fins.Air.Vdot_ha=self.Fins.Air.Vdot_ha[-1]/Ncircuits
            else:
                for i in range(Ncircuits):
                    self.Evaps[i].Fins.Air.Vdot_ha=self.Fins.Air.Vdot_ha[i]
                    
        # Distribute the tubes of the bank among the different circuits
        # If Tubes per bank is divisible by the number of circuits, all the 
        # circuits have the same number of tubes per bank
        
        # The circuits are ordered from fewer to more if they are not evenly distributed
        NTubes_min=int(floor(self.Fins.Tubes.NTubes_per_bank/Ncircuits))
        NTubes_max=int(ceil(self.Fins.Tubes.NTubes_per_bank/Ncircuits))
        
        if NTubes_min==NTubes_max:
            #If evenly divisible, use the tubes per circuit from the division
            A=Ncircuits
        else:
            # Total number of tubes per bank is given by
            A=(self.Fins.Tubes.NTubes_per_bank-Ncircuits*NTubes_max)/(NTubes_min-NTubes_max)
        
        for i in range(Ncircuits):
            if i+1<=A:
                self.Evaps[i].Fins.Tubes.NTubes_per_bank=NTubes_min
            else:
                self.Evaps[i].Fins.Tubes.NTubes_per_bank=NTubes_max
            
        for i in range(Ncircuits):
            self.Evaps[i].Fins.Tubes.Ncircuits=1
            #Actually run each Evaporator
            self.Evaps[i].AS = AS
            self.Evaps[i].Calculate()
            
        #Collect the outputs from each of the evaporators individually
        #Try to mirror the outputs of each of the evaporators
        self.Q=np.sum([self.Evaps[i].Q for i in range(Ncircuits)])
        self.Charge=np.sum([self.Evaps[i].Charge for i in range(Ncircuits)])
        self.Charge_superheat=np.sum([self.Evaps[i].Charge_superheat for i in range(Ncircuits)])
        self.Charge_2phase=np.sum([self.Evaps[i].Charge_2phase for i in range(Ncircuits)])
        self.DP_r=np.mean([self.Evaps[i].DP_r for i in range(Ncircuits)])  #simplified
        self.DP_r_superheat=np.mean([self.Evaps[i].DP_r_superheat for i in range(Ncircuits)])  #simplified
        self.DP_r_2phase=np.mean([self.Evaps[i].DP_r_2phase for i in range(Ncircuits)])  #simplified
        self.Tin_r=np.mean([self.Evaps[i].Tin_r for i in range(Ncircuits)])  #simplified
        self.h_r_superheat=np.mean([self.Evaps[i].h_r_superheat for i in range(Ncircuits)])  #simplified, really should consider flowrate
        self.h_r_2phase=np.mean([self.Evaps[i].h_r_2phase for i in range(Ncircuits)])  #simplified, really should consider flowrate
        self.w_superheat=np.sum([self.Evaps[i].w_superheat for i in range(Ncircuits)])/float(Ncircuits)
        self.w_2phase=np.sum([self.Evaps[i].w_2phase for i in range(Ncircuits)])/float(Ncircuits)
        self.hout_r=0.0
        for i in range(Ncircuits): self.hout_r+=self.Evaps[i].hout_r*self.Evaps[i].mdot_r 
        self.hout_r=(self.hout_r/sum(self.mdot_r))
        self.Tin_a=self.Evaps[0].Fins.Air.Tdb #assuming equal temperature for all circuits
        self.Tout_a=0.0
        for i in range(Ncircuits): self.Tout_a+=self.Evaps[i].Tout_a*self.Evaps[i].Fins.Air.Vdot_ha
        self.Tout_a=(self.Tout_a/sum(self.Fins.Air.Vdot_ha))
        Pout_r=self.psat_r+self.DP_r/1.0
        AS.update(CP.PQ_INPUTS, Pout_r, 1.0)
        hsatV_out=AS.hmass() #[J/kg]
        AS.update(CP.PQ_INPUTS, Pout_r, 0.0)
        hsatL_out=AS.hmass() #[J/kg]
        if self.hout_r>hsatV_out:
            AS.update(CP.HmassP_INPUTS, self.hout_r, Pout_r)
            self.Tout_r= AS.T() #superheated temperature at outlet [K]
        else:
            xout_r=((self.hout_r-hsatL_out)/(hsatV_out-hsatL_out))
            AS.update(CP.PQ_INPUTS, Pout_r, xout_r)
            self.Tout_r=AS.T() #saturated temperature at outlet quality [K]
        self.Capacity=np.sum([self.Evaps[i].Q for i in range(Ncircuits)])-self.Fins.Air.FanPower
        self.SHR=np.mean([self.Evaps[i].SHR for i in range(Ncircuits)])
        self.UA_a=np.sum([self.Evaps[i].UA_a for i in range(Ncircuits)])
        self.UA_r=np.sum([self.Evaps[i].UA_r for i in range(Ncircuits)])
        self.Q_superheat=np.sum([self.Evaps[i].Q_superheat for i in range(Ncircuits)])
        self.Q_2phase=np.sum([self.Evaps[i].Q_2phase for i in range(Ncircuits)])
        #Convert back to a single value for the overall evaporator
        self.Fins.Air.Vdot_ha=float(self.Fins.Air.Vdot_ha[-1])
        if self.Verbosity>0:
            print (chr(127)) #progress bar


if __name__=='__main__':
    # This code runs if this file is run by itself, but otherwise doesn't run
    # Exampe usage of multi circuit evaporator
    
    FinsTubes=FinInputs()

    FinsTubes.Tubes.NTubes_per_bank=32
    FinsTubes.Tubes.Ncircuits=5
    FinsTubes.Tubes.Nbank=3
    FinsTubes.Tubes.Ltube=0.452
    FinsTubes.Tubes.OD=0.009525
    FinsTubes.Tubes.ID=0.0089154
    FinsTubes.Tubes.Pl=0.0254
    FinsTubes.Tubes.Pt=0.0219964
    FinsTubes.Tubes.kw=237
    
    FinsTubes.Fins.FPI=14.5
    FinsTubes.Fins.Pd=0.001
    FinsTubes.Fins.xf=0.001
    FinsTubes.Fins.t=0.00011
    FinsTubes.Fins.k_fin=237
    
    FinsTubes.Air.Vdot_ha=0.5663
    FinsTubes.Air.Tmean=299.8
    FinsTubes.Air.Tdb=299.8
    FinsTubes.Air.p=101325
    FinsTubes.Air.RH=0.51
    FinsTubes.Air.RHmean=0.51
    FinsTubes.Air.FanPower=438
    
    #Abstract State
    Ref = 'R410A'
    Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref)
    
    Tdew=282.0
    kwargs={'AS': AS,
            'mdot_r': 0.0708,
            'psat_r': PropsSI('P','T',Tdew,'Q',1.0,Ref),
            'Fins': FinsTubes,
            'FinsType': 'WavyLouveredFins', #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
            'hin_r': PropsSI('H','P',PropsSI('P','T',Tdew,'Q',1.0,Ref),'Q',0.15,Ref),
            'Verbosity': 0,
            'h_a_tuning':1,
            'h_tp_tuning':1,
            'DP_tuning':1,
    }
     
    Evap=EvaporatorClass(**kwargs)
    Evap.Update(**kwargs)
    Evap.Calculate()
    
    Tdew=282.0
    kwargs={'AS': AS,
            'mdot_r':  0.0708,
            'mdot_r_coeffs':[0.1,0.1,0.3,0.31,0.19],   #Mass flow distribution at distributor
            'mdot_v_coeffs':[0.2,0.2,0.2,0.2,0.2],   #Quality distribution at distributor
            'Vdot_ha_coeffs':[0.2,0.2,0.2,0.2,0.2],  #airside flow distribution
            'psat_r': PropsSI('P','T',Tdew,'Q',1.0,Ref),
            'Fins': FinsTubes,
            'FinsType': 'WavyLouveredFins',                                  
            'hin_r': PropsSI('H','P',PropsSI('P','T',Tdew,'Q',1.0,Ref),'Q',0.15,Ref),
            'Verbosity':0,
            'TestName':'MCE-0014',
            'TestDescription':'shows application of MCE',
            'TestDetails':'This is the sample multi circuited evaporator',
            'h_a_tuning':1,
            'h_tp_tuning':1,
            'DP_tuning':1,
    }
    
    MCE=MultiCircuitEvaporatorClass(**kwargs)
    MCE.Update(**kwargs)
    MCE.Calculate()
    print ('Q_standard_evap = ' + str(Evap.Q) + ' W')
    print ('Q_MCE = '+str(MCE.Q)+' W')
    print ('Heat transfer reduction due to maldisribution'+str((MCE.Q-Evap.Q)*100./Evap.Q)+' %')

    print ("demonstrating output list\n",'='*20,'\n')
    print (MCE.OutputList(), '\n', '='*20)

    csv_file= "Evaporator_MCE.csv"
    print ("\n demonstrating write to csv",csv_file,'\n')
    #Write2CSV(MCE,open(csv_file,'w'),append=False)
      
    print ("check MCE capacity summation", MCE.hin_r*MCE.mdot_r[-1])
    print (np.sum([MCE.Evaps[i].hin_r*MCE.Evaps[i].mdot_r for i in range(MCE.Fins.Tubes.Ncircuits)]))

    # plot maldistribution
    # h = np.array([MCE.Evaps[i].hout_r for i in range(len(MCE.Evaps))])
    # p = MCE.psat_r*(1+0*h)
    # plot = PropertyPlot('HEOS::R410A', 'PH', unit_system='KSI')
    # plot.calc_isolines(CoolProp.iQ, num=2)
    # plot.axis.plot(MCE.hin_r/1000,MCE.psat_r/1000,'>', label= 'Inlet')
    # plot.axis.plot(h/1000,p/1000,'x', label= 'Individual circuit exits')
    # plot.axis.plot(MCE.hout_r/1000,MCE.psat_r/1000,'o', label= 'Overall exit')
    # plt.legend(loc= 'best', numpoints=1)
    # plot.savefig('MultiCircuitEvaporator_py_example.pdf')
    # plot.figure.show()
    # plt.close(plot.figure)
    
    print  ("outlet enthalpies", [MCE.Evaps[i].hout_r for i in range(len(MCE.Evaps))], MCE.hout_r, PropsSI('H','T',MCE.Evaps[-1].Tdew_r,'Q',1,Ref))
    print  ("outlet superheats", [MCE.Evaps[i].DT_sh_calc for i in range(len(MCE.Evaps))])

    print ("\nrerunning with smaller mass flowrate")

#     kwargs={'AS': AS,
#             'mdot_r':  0.05,
#             'TestDescription':'shows application of MCE',
#             'TestDetails':'Reduced mass flowrate'
#     }
#     
#     MCE.Update(**kwargs)
#     MCE.Calculate()
#     Write2CSV(MCE,open(csv_file,'a'),append=True) #does

    # plot maldistribution
    # h = np.array([MCE.Evaps[i].hout_r for i in range(len(MCE.Evaps))])
    # p = MCE.psat_r*np.ones_like(h)
    # plot = PropertyPlot('HEOS::R410A', 'PH', unit_system='KSI')
    # plot.calc_isolines(CoolProp.iQ, num=2)
    # plot.axis.plot(MCE.hin_r/1000,MCE.psat_r/1000,'>', label= 'Inlet')
    # plot.axis.plot(h/1000,p/1000,'x', label= 'Individual circuit exits')
    # plot.axis.plot(MCE.hout_r/1000,MCE.psat_r/1000,'o', label= 'Overall exit')
    # plt.legend(loc= 'best', numpoints=1)
    # plot.savefig('MultiCircuitEvaporator_py_example_reduced_flowrate.pdf')
    # plt.show()
    # plt.close(plot.figure)