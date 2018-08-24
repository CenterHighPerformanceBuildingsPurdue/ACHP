from __future__ import division, print_function, absolute_import
from ACHP.Cycle import DXCycleClass,SecondaryCycleClass
from ACHP.convert_units import F2K
from ACHP.ACHPTools import Write2CSV

def SampleSecondaryLoopSystem():
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
    Cycle.Verbosity = 10 #the idea here is to have different levels of debug output
    Cycle.ImposedVariable = 'Subcooling'                                                  
    Cycle.DT_sc_target = 7.0
    Cycle.Charge_target = 2.4
    Cycle.DT_sh=5
    Cycle.Ref='R410A'
    Cycle.Backend='HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    Cycle.Oil = 'POE32'
    Cycle.shell_pressure = 'low-pressure'
    Cycle.SecLoopFluid = 'MEG'
    Cycle.Backend_SLF = 'INCOMP'
    Cycle.IHXType = 'PHE'# or could be 'Coaxial'
    Cycle.Mode='AC'
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
        'V_oil_sump':0, #Volume of oil in the sump
        'shell_pressure':Cycle.shell_pressure, #Compressor shell pressure
        'fp':0.15, #Fraction of electrical power lost as heat to ambient            #shell heat loss
        'Vdot_ratio': 1.0, #Displacement Scale factor                               #up- or downsize compressor (1=original)
        'Verbosity': 0, # How verbose should the debugging statements be [0 to 10]
        'DT_sh': 5. #Superheat at inlet to compressor [K]                           
      }
    Cycle.Compressor.Update(**params)
    
    #--------------------------------------
    #--------------------------------------
    #      Evaporator parameters 
    #      -> see GUI for illustration/units
    #--------------------------------------
    #--------------------------------------                                                                                                     
    Cycle.Condenser.Fins.Tubes.NTubes_per_bank=41  #number of tubes per bank=row
    Cycle.Condenser.Fins.Tubes.Nbank=1             #number of baks/rows
    Cycle.Condenser.Fins.Tubes.Ncircuits=5         #number of baks/rows
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
    Cycle.Condenser.Fins.Air.FanPower=260
    
    params={
        'Verbosity':0,
        'FinsType':'WavyLouveredFins'                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'

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
    Cycle.CoolingCoil.Fins.Tubes.Ncircuits=4
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
        'Ref_g': Cycle.SecLoopFluid,
        'Backend_g': Cycle.Backend_SLF,
        'pin_g': 200000,
        'Verbosity':0,
        'mdot_g':0.38,
        'FinsType':'WavyLouveredFins'                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        
    }
    Cycle.CoolingCoil.Update(**params)
    
    params={
        'ID_i':0.0278,
        'OD_i':0.03415,
        'ID_o':0.045,
        'L':50,
        'pin_g':300000,
        'Ref_r':Cycle.Ref,
        'Ref_g':Cycle.SecLoopFluid,
        'Backend_g': Cycle.Backend_SLF,
        'Verbosity':0
        }
    Cycle.CoaxialIHX.Update(**params)
    
    params={
        'pin_h':300000,
        'Ref_h':Cycle.SecLoopFluid,
        'Backend_h': Cycle.Backend_SLF,
        'Ref_c':Cycle.Ref,
        
        #Geometric parameters
        'Bp' : 0.117,
        'Lp' : 0.300, #Center-to-center distance between ports
        'Nplates' : 46,
        'PlateAmplitude' : 0.001, #[m]
        'PlateThickness' : 0.0003, #[m]
        'PlateConductivity' : 15.0, #[W/m-K]
        'Rp': 1.0, #[microns] Surface roughness
        'PlateWavelength' : 0.00628, #[m]
        'InclinationAngle' : 3.14159/3,#[rad]
        'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
        'Verbosity':0,
    }
    Cycle.PHEIHX.Update(**params)
    
    params={
        'eta':0.5,  #Pump+motor efficiency
        'mdot_g':0.38, #Flow Rate kg/s
        'pin_g':300000,
        'Ref_g':Cycle.SecLoopFluid,
        'Backend_g': Cycle.Backend_SLF,
        'Verbosity':0,
        }
    Cycle.Pump.Update(**params)
    
    params={
        'L':5,
        'k_tube':0.19,
        't_insul':0.02,
        'k_insul':0.036,
        'T_air':297,
        'Ref': Cycle.SecLoopFluid,
        'Backend': Cycle.Backend_SLF,
        'pin': 300000,
        'h_air':0.0000000001,
    }
    
    Cycle.LineSetSupply.Update(**params)
    Cycle.LineSetReturn.Update(**params)
    Cycle.LineSetSupply.OD=0.009525
    Cycle.LineSetSupply.ID=0.007986
    Cycle.LineSetReturn.OD=0.01905
    Cycle.LineSetReturn.ID=0.017526
    
    #Now solve
    Cycle.PreconditionedSolve()
    
    print (Cycle.OutputList())

def SampleSecondaryLoopHPSystem():
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
    Cycle.Verbosity = 5 #the idea here is to have different levels of debug output
    Cycle.ImposedVariable = 'Subcooling'                                                           
    Cycle.DT_sc_target = 7.0
    Cycle.Charge_target = 3.3
    Cycle.DT_sh=5
    Cycle.Ref='R410A'
    Cycle.Backend='HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    Cycle.Oil = 'POE32'
    Cycle.shell_pressure = 'low-pressure'
    Cycle.SecLoopFluid = 'Water'
    Cycle.Backend_SLF = 'INCOMP'
    Cycle.IHXType = 'PHE'
    Cycle.Mode = 'HP'
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
        'V_oil_sump':0, #Volume of oil in the sump
        'shell_pressure':Cycle.shell_pressure, #Compressor shell pressure
        'fp':0.15, #Fraction of electrical power lost as heat to ambient            #shell heat loss
        'Vdot_ratio': 1.0, #Displacement Scale factor                               #up- or downsize compressor (1=original)
        'Verbosity': 0, # How verbose should the debugging statements be [0 to 10]
                           
      }
    Cycle.Compressor.Update(**params)
    
    #--------------------------------------
    #--------------------------------------
    #      Condenser parameters 
    #      -> see GUI for illustration/units
    #--------------------------------------
    #--------------------------------------                                                                                                     
    Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=41  #number of tubes per bank=row
    Cycle.Evaporator.Fins.Tubes.Nbank=1             #number of baks/rows
    Cycle.Evaporator.Fins.Tubes.Ncircuits=5             #number of baks/rows
    Cycle.Evaporator.Fins.Tubes.Ltube=2.286
    Cycle.Evaporator.Fins.Tubes.OD=0.007
    Cycle.Evaporator.Fins.Tubes.ID=0.0063904
    Cycle.Evaporator.Fins.Tubes.Pl=0.0191  #distance between center of tubes in flow direction                                                
    Cycle.Evaporator.Fins.Tubes.Pt=0.0222  #distance between center of tubes orthogonal to flow direction
    Cycle.Evaporator.Fins.Tubes.kw=237     #wall thermal conductivity (i.e pipe material)
    
    Cycle.Evaporator.Fins.Fins.FPI=25      #Number of fins per inch
    Cycle.Evaporator.Fins.Fins.Pd=0.001    #2* amplitude of wavy fin
    Cycle.Evaporator.Fins.Fins.xf=0.001    #1/2 period of fin
    Cycle.Evaporator.Fins.Fins.t=0.00011   #Thickness of fin material
    Cycle.Evaporator.Fins.Fins.k_fin=237   #Thermal conductivity of fin material
    
    Cycle.Evaporator.Fins.Air.Vdot_ha=1.7934 #rated volumetric flowrate
    Cycle.Evaporator.Fins.Air.Tmean=F2K(40)   
    Cycle.Evaporator.Fins.Air.Tdb=F2K(40)     #Dry Bulb Temperature
    Cycle.Evaporator.Fins.Air.p=101325       #Air pressure in Pa
    Cycle.Evaporator.Fins.Air.RH=0.51        #Relative Humidity
    Cycle.Evaporator.Fins.Air.RHmean=0.51    
    Cycle.Evaporator.Fins.Air.FanPower=160
    
    Cycle.Evaporator.DT_sh= 5. #Superheat at inlet to compressor [K]
    params={
        'Verbosity':0,
        'FinsType':'WavyLouveredFins'                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'

    }
    Cycle.Evaporator.Update(**params)
        
    #--------------------------------------
    #--------------------------------------
    #           Cooling Coil
    #           -> see Condenser and GUI for explanations
    #--------------------------------------
    #--------------------------------------
    
    Cycle.CoolingCoil.Fins.Tubes.NTubes_per_bank=32
    Cycle.CoolingCoil.Fins.Tubes.Nbank=4
    Cycle.CoolingCoil.Fins.Tubes.Ncircuits=4
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
    Cycle.CoolingCoil.Fins.Air.Tmean=F2K(70)
    Cycle.CoolingCoil.Fins.Air.Tdb=F2K(70)
    Cycle.CoolingCoil.Fins.Air.p=101325
    Cycle.CoolingCoil.Fins.Air.RH=0.51
    Cycle.CoolingCoil.Fins.Air.RHmean=0.51
    Cycle.CoolingCoil.Fins.Air.FanPower=438
           
    params={
        'Ref_g': Cycle.SecLoopFluid,
        'pin_g': 200000,
        'Verbosity':0,
        'FinsType':'WavyLouveredFins'                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
    }
    Cycle.CoolingCoil.Update(**params)
    
    params={
        'pin_c':300000,
        'Ref_c':Cycle.SecLoopFluid,
        'Backend_c': Cycle.Backend_SLF,
        'Ref_h':Cycle.Ref,
        
        #Geometric parameters
        'Bp' : 0.117,
        'Lp' : 0.300, #Center-to-center distance between ports
        'Nplates' : 46,
        'PlateAmplitude' : 0.001, #[m]
        'PlateThickness' : 0.0003, #[m]
        'PlateConductivity' : 15.0, #[W/m-K]
        'Rp': 1.0, #[microns] Surface roughness
        'PlateWavelength' : 0.00628, #[m]
        'InclinationAngle' : 3.14159/3,#[rad]
        'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
        'Verbosity':0
    }
    Cycle.PHEIHX.Update(**params)
    
    params={
        'eta':0.5,  #Pump+motor efficiency
        'mdot_g':0.38, #Flow Rate kg/s
        'pin_g':300000,
        'Ref_g':Cycle.SecLoopFluid,
        'Backend_g': Cycle.Backend_SLF,
        'Verbosity':0,
        }
    Cycle.Pump.Update(**params)
    
    params={
        'L':5,
        'k_tube':0.19,
        't_insul':0.02,
        'k_insul':0.036,
        'T_air':297,
        'Ref': Cycle.SecLoopFluid,
        'Backend': Cycle.Backend_SLF,
        'h_air':0.0000000001,
    }
    
    Cycle.LineSetSupply.pin=300000
    Cycle.LineSetReturn.pin=300000
    Cycle.LineSetSupply.Update(**params)
    Cycle.LineSetReturn.Update(**params)
    Cycle.LineSetSupply.OD=0.009525
    Cycle.LineSetSupply.ID=0.007986
    Cycle.LineSetReturn.OD=0.01905
    Cycle.LineSetReturn.ID=0.017526
    
    #Now solve
    Cycle.PreconditionedSolve()
    
    print (Cycle.Pump.DP_g,Cycle.Pump.W)
        
def SampleDXACSystem(Calculate=True):    
    """
    A sample DX Air Conditioning system.  This is based on the work of Bo Shen.
    """
    #########################################################################
    ######################     CYCLE INITIALIZATION    ######################
    #########################################################################
    
    ## Here we load parameters that are not a function of operating conditions
    ## They are primarily geometric parameters
    
    Cycle=DXCycleClass()
    
    #--------------------------------------
    #--------------------------------------
    #       Cycle parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Verbosity = 10 #the idea here is to have different levels of debug output
    Cycle.ImposedVariable = 'Subcooling'                                                           
    Cycle.DT_sc_target = 7.0
    Cycle.Charge_target = 2.8
    Cycle.Mode='AC' 
    Cycle.Ref='R410A'
    Cycle.Backend='HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    Cycle.Oil = 'POE32'
    Cycle.shell_pressure = 'low-pressure'
    Cycle.TestName='DXAC-0014'  #this and the two next lines can be used to specify exact test conditions
    Cycle.TestDescription='shows application of DXAC system'
    Cycle.TestDetails='This is the sample cycle for the DXAC system which can be modified for other systems'
    Cycle.EvapSolver = 'Moving-Boundary' #choose the type of Evaporator solver scheme (for now only 'Moving-Boundary'')
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
        'V_oil_sump':0, #Volume of oil in the sump
        'shell_pressure':Cycle.shell_pressure, #Compressor shell pressure
        'fp':0.0, #Fraction of electrical power lost as heat to ambient            #shell heat loss
        'Vdot_ratio': 1.0, #Displacement Scale factor                               #up- or downsize compressor (1=original)
        'Verbosity': 0, # How verbose should the debugging statements be [0 to 10]                           
      }
    Cycle.Compressor.Update(**params)                                          
    
    #--------------------------------------
    #--------------------------------------
    #      Condenser parameters 
    #      -> see GUI for illustration/units
    #--------------------------------------
    #--------------------------------------
    Cycle.Condenser.Fins.Tubes.NTubes_per_bank=24  #number of tubes per bank=row
    Cycle.Condenser.Fins.Tubes.Nbank=1             #number of banks/rows
    Cycle.Condenser.Fins.Tubes.Ncircuits=3
    Cycle.Condenser.Fins.Tubes.Ltube=2.252
    Cycle.Condenser.Fins.Tubes.OD=0.00913
    Cycle.Condenser.Fins.Tubes.ID=0.00849
    Cycle.Condenser.Fins.Tubes.Pl=0.0191  #distance between center of tubes in flow direction                                                
    Cycle.Condenser.Fins.Tubes.Pt=0.0254  #distance between center of tubes orthogonal to flow direction
    Cycle.Condenser.Fins.Tubes.kw=237     #wall thermal conductivity (i.e pipe material)
    
    Cycle.Condenser.Fins.Fins.FPI=25      #Number of fins per inch
    Cycle.Condenser.Fins.Fins.Pd=0.001    #2* amplitude of wavy fin
    Cycle.Condenser.Fins.Fins.xf=0.001    #1/2 period of fin
    Cycle.Condenser.Fins.Fins.t=0.00011   #Thickness of fin material
    Cycle.Condenser.Fins.Fins.k_fin=237   #Thermal conductivity of fin material
    
    Cycle.Condenser.Fins.Air.Vdot_ha=1.7934 #rated volumetric flowrate
    Cycle.Condenser.Fins.Air.Tmean=308.15   
    Cycle.Condenser.Fins.Air.Tdb=Cycle.Condenser.Fins.Air.Tmean     #Dry Bulb Temperature
    Cycle.Condenser.Fins.Air.p=101325       #Air pressure in Pa
    Cycle.Condenser.Fins.Air.RH=0.51        #Relative Humidity
    Cycle.Condenser.Fins.Air.RHmean=0.51
    Cycle.Condenser.Fins.Air.FanPower=260
    
    params={
        'Verbosity':0,
        'FinsType':'WavyLouveredFins'                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'

    }
    Cycle.Condenser.Update(**params)
    
    # ----------------------------------
    #       Expanison device Parameters
    # ----------------------------------
    params={
            'ExpType':'Ideal',     #expansion device type
        }
    Cycle.ExpDev.Update(**params)

    #--------------------------------------
    #--------------------------------------
    #           Evaporator
    #           -> see Condenser and GUI for explanations
    #--------------------------------------
    #--------------------------------------
    if hasattr(Cycle,'TestName'):                                           #update parameters for output list, if applicable
        Cycle.Evaporator.TestName=Cycle.TestName 
    if hasattr(Cycle,'TestDescription'):
        Cycle.Evaporator.TestDescription=Cycle.TestDescription
    if hasattr(Cycle,'TestDetails'):
        Cycle.Evaporator.TestDetails=Cycle.TestDetails
    Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=32  #dimensional parameters
    Cycle.Evaporator.Fins.Tubes.Nbank=3
    Cycle.Evaporator.Fins.Tubes.Ltube=0.452
    Cycle.Evaporator.Fins.Tubes.OD=0.00913
    Cycle.Evaporator.Fins.Tubes.ID=0.00849
    Cycle.Evaporator.Fins.Tubes.Pl=0.0191
    Cycle.Evaporator.Fins.Tubes.Pt=0.0254
    Cycle.Evaporator.Fins.Tubes.Ncircuits=5
    Cycle.Evaporator.Fins.Tubes.kw=237     #wall thermal conductivity (i.e pipe material)
    
    Cycle.Evaporator.Fins.Fins.FPI=14.5
    Cycle.Evaporator.Fins.Fins.Pd=0.001
    Cycle.Evaporator.Fins.Fins.xf=0.001
    Cycle.Evaporator.Fins.Fins.t=0.00011
    Cycle.Evaporator.Fins.Fins.k_fin=237
    
    Cycle.Evaporator.Fins.Air.Vdot_ha=0.56319
    Cycle.Evaporator.Fins.Air.Tmean=297.039
    Cycle.Evaporator.Fins.Air.Tdb=297.039
    Cycle.Evaporator.Fins.Air.p=101325
    Cycle.Evaporator.Fins.Air.RH=0.5
    Cycle.Evaporator.Fins.Air.RHmean=0.5
    Cycle.Evaporator.Fins.Air.FanPower=438
           
    params={
        'Ref': Cycle.Ref,
        'Verbosity':0,
        'DT_sh':5,
        'FinsType':'WavyLouveredFins'                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
    }
    Cycle.Evaporator.Update(**params)
    
    params={
        'L':7.6,
        'k_tube':0.19,
        't_insul':0.02,
        'k_insul':0.036,
        'T_air':297,
        'Ref': Cycle.Ref,
        'h_air':0.0000000001
    }
    
    Cycle.LineSetLiquid.Update(**params)
    Cycle.LineSetSuction.Update(**params)
    Cycle.LineSetLiquid.OD=0.009525
    Cycle.LineSetLiquid.ID=0.007986
    Cycle.LineSetSuction.OD=0.01905
    Cycle.LineSetSuction.ID=0.017526
    
    # ----------------------------------
    # ----------------------------------
    #       Line Set Discharge Parameters
    # ----------------------------------
    # ----------------------------------
    params={
            'L':0.3,                #tube length in m
            'k_tube':0.19,
            't_insul':0, #no insulation
            'k_insul':0.036,
            'T_air':297,
            'Ref': Cycle.Ref,
            'h_air':0.0000000001,
            'LineSetOption': 'Off'
            }
      
    Cycle.LineSetDischarge.Update(**params)
    Cycle.LineSetDischarge.OD=0.009525
    Cycle.LineSetDischarge.ID=0.007986
    
    #Now solve if Calculate has not been set to False
    if Calculate==True:
        Cycle.PreconditionedSolve()
    
    return Cycle
    

        
def SampleDXHPSystem():    
    #########################################################################
    ######################     CYCLE INITIALIZATION    ######################
    #########################################################################
    
    ## Here we load parameters that are not a function of operating conditions
    ## They are primarily geometric parameters
    
    Cycle=DXCycleClass()
    
    #--------------------------------------
    #--------------------------------------
    #       Cycle parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Verbosity = 5 #the idea here is to have different levels of debug output
    Cycle.ImposedVariable = 'Subcooling'                                                           
    Cycle.DT_sc_target = 7.0
    Cycle.Charge_target = 3.3
    Cycle.Mode='HP' 
    Cycle.Ref='R410A'
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
        'fp':0.0, #Fraction of electrical power lost as heat to ambient            #shell heat loss
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

    Cycle.Condenser.Fins.Tubes.NTubes_per_bank=32
    Cycle.Condenser.Fins.Tubes.Nbank=3
    Cycle.Condenser.Fins.Tubes.Ncircuits=6
    Cycle.Condenser.Fins.Tubes.Ltube=0.452
    Cycle.Condenser.Fins.Tubes.OD=0.009525
    Cycle.Condenser.Fins.Tubes.ID=0.0089154
    Cycle.Condenser.Fins.Tubes.Pl=0.0254
    Cycle.Condenser.Fins.Tubes.Pt=0.0219964
    Cycle.Condenser.Fins.Tubes.kw=237                   #wall thermal conductivity (i.e pipe material)
    
    Cycle.Condenser.Fins.Fins.FPI=14.5
    Cycle.Condenser.Fins.Fins.Pd=0.001
    Cycle.Condenser.Fins.Fins.xf=0.001
    Cycle.Condenser.Fins.Fins.t=0.00011
    Cycle.Condenser.Fins.Fins.k_fin=237
    
    Cycle.Condenser.Fins.Air.Vdot_ha=0.5663
    Cycle.Condenser.Fins.Air.Tmean=F2K(70)
    Cycle.Condenser.Fins.Air.Tdb=F2K(70)
    Cycle.Condenser.Fins.Air.p=101325
    Cycle.Condenser.Fins.Air.RH=0.51
    Cycle.Condenser.Fins.Air.RHmean=0.51
    Cycle.Condenser.Fins.Air.FanPower=438
    
    params={
        'Verbosity':0,
        'FinsType':'WavyLouveredFins'                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
    }
    Cycle.Condenser.Update(**params)
        
    #--------------------------------------
    #--------------------------------------
    #           Evaporator
    #           -> see Condenser and GUI for explanations
    #--------------------------------------
    #--------------------------------------
    Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=41  #number of tubes per bank=row
    Cycle.Evaporator.Fins.Tubes.Nbank=1             #number of banks/rows
    Cycle.Evaporator.Fins.Tubes.Ncircuits=5
    Cycle.Evaporator.Fins.Tubes.Ltube=2.286
    Cycle.Evaporator.Fins.Tubes.OD=0.007
    Cycle.Evaporator.Fins.Tubes.ID=0.0063904
    Cycle.Evaporator.Fins.Tubes.Pl=0.0191  #distance between center of tubes in flow direction                                                
    Cycle.Evaporator.Fins.Tubes.Pt=0.0222  #distance between center of tubes orthogonal to flow direction
    
    Cycle.Evaporator.Fins.Fins.FPI=25      #Number of fins per inch
    Cycle.Evaporator.Fins.Fins.Pd=0.001    #2* amplitude of wavy fin
    Cycle.Evaporator.Fins.Fins.xf=0.001    #1/2 period of fin
    Cycle.Evaporator.Fins.Fins.t=0.00011   #Thickness of fin material
    Cycle.Evaporator.Fins.Fins.k_fin=237   #Thermal conductivity of fin material
    
    Cycle.Evaporator.Fins.Air.Vdot_ha=1.7934 #rated volumetric flowrate
    Cycle.Evaporator.Fins.Air.Tmean=F2K(47)   
    Cycle.Evaporator.Fins.Air.Tdb=F2K(47)     #Dry Bulb Temperature
    Cycle.Evaporator.Fins.Air.p=101325       #Air pressure in Pa
    Cycle.Evaporator.Fins.Air.RH=0.51        #Relative Humidity
    Cycle.Evaporator.Fins.Air.FanPower=160
    Cycle.Evaporator.Fins.Air.RHmean=0.51
           
    params={
        'Ref': Cycle.Ref,
        'Verbosity':0,
        'DT_sh':5,
        'FinsType':'WavyLouveredFins'                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
    }
    Cycle.Evaporator.Update(**params)
    
    params={
        'L':7.6,
        'k_tube':0.19,
        't_insul':0.02,
        'k_insul':0.036,
        'T_air':297,
        'Ref': Cycle.Ref,
        'h_air':6,
        'LineSetOption': 'Off'
    }
    
    Cycle.LineSetSupply.Update(**params)
    Cycle.LineSetReturn.Update(**params)
    Cycle.LineSetSupply.OD=0.01905 
    Cycle.LineSetSupply.ID=0.017526
    Cycle.LineSetReturn.OD=0.009525
    Cycle.LineSetReturn.ID=0.007986
    
    #Now solve
    Cycle.PreconditionedSolve()
    
if __name__=='__main__':
    
    print ('Running SampleDXACSystem')
    cycle=SampleDXACSystem()
    #Write the outputs to file
    Write2CSV(cycle.Evaporator,open('Evaporator.csv','w'),append=False)
    Write2CSV(cycle,open('Cycle.csv','w'),append=False)
    #append a second run with different temperauture
    new_outdoor_temp=290
    params={
        'Tin_a':new_outdoor_temp,
    }
    cycle.Condenser.Update(**params)
    cycle.Condenser.Fins.Air.Tdb=new_outdoor_temp
    cycle.TestName='DXAC-0018'  #this and the two next lines can be used to specify exact test conditions
    cycle.TestDescription='another point'
    cycle.TestDetails='Here we changed the air temperature'
    cycle.Calculate(5,7.0,5)  #Calculate(DT_evap,DT_cond,DT_sh), as defined in the DXCycleClass
    Write2CSV(cycle.Evaporator,open('Evaporator.csv','a'),append=True)
    Write2CSV(cycle,open('Cycle.csv','a'),append=True)
    #SampleDXHPSystem()
    print ()
    print ('Running SampleSecondaryLoopSystem')
    SampleSecondaryLoopSystem()
    print ()
    #print ('Running SampleSecondaryHPLoopSystem')
    #SampleSecondaryLoopHPSystem()
    
