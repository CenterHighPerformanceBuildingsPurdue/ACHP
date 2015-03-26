from Cycle import SecondaryCycleClass

#Instantiate the class
Cycle=SecondaryCycleClass()

#--------------------------------------
#       Cycle parameters
#--------------------------------------
Cycle.Verbosity = 0 #the idea here is to have different levels of debug output
Cycle.ImposedVariable = 'Subcooling' #or this could be 'Charge' for imposed charge
Cycle.DT_sc_target = 7.0
#Cycle.Charge_target = 2.4  #Needed if charge is imposed, not otherwise
Cycle.Ref='R410A'
Cycle.SecLoopFluid = 'EG-20%'
Cycle.IHXType = 'PHE'# or could be 'Coaxial'
Cycle.Mode='AC'

#--------------------------------------
#--------------------------------------
#       Compressor parameters
#--------------------------------------
#--------------------------------------

#A 3 ton cooling capacity compressor map
if Cycle.Ref=='R410A':
    M=[217.3163128,5.094492028,-0.593170311,4.38E-02,-2.14E-02,1.04E-02,
    7.90E-05,-5.73E-05,1.79E-04,-8.08E-05]
    P=[-561.3615705,-15.62601841,46.92506685,-0.217949552,0.435062616,
    -0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03]

params={
    'M':M,
    'P':P,
    'Ref':Cycle.Ref,    #refrigerant
    'fp':0.15, #Fraction of electrical power lost as heat to ambient            
    'Vdot_ratio': 1.0, #Displacement Scale factor to up- or downsize compressor (1=original)
    'Verbosity': 0, # How verbose should the debugging statements be [0 to 10]
  }
Cycle.Compressor.Update(**params)

#--------------------------------------
#      Condenser parameters 
#--------------------------------------
Cycle.Condenser.Fins.Tubes.NTubes_per_bank=24  #number of tubes per bank=row
Cycle.Condenser.Fins.Tubes.Nbank=1             #number of banks/rows
Cycle.Condenser.Fins.Tubes.Ncircuits=3
Cycle.Condenser.Fins.Tubes.Ltube=2.252
Cycle.Condenser.Fins.Tubes.OD=0.00913
Cycle.Condenser.Fins.Tubes.ID=0.00849
Cycle.Condenser.Fins.Tubes.Pl=0.0191  #distance between center of tubes in flow direction                                                
Cycle.Condenser.Fins.Tubes.Pt=0.0254  #distance between center of tubes orthogonal to flow direction

Cycle.Condenser.Fins.Fins.FPI=25      #Number of fins per inch
Cycle.Condenser.Fins.Fins.Pd=0.001    #2* amplitude of wavy fin
Cycle.Condenser.Fins.Fins.xf=0.001    #1/2 period of fin
Cycle.Condenser.Fins.Fins.t=0.00011   #Thickness of fin material
Cycle.Condenser.Fins.Fins.k_fin=237   #Thermal conductivity of fin material

Cycle.Condenser.Fins.Air.Vdot_ha=1.7934 #rated volumetric flowrate
Cycle.Condenser.Fins.Air.Tmean=308.15   
Cycle.Condenser.Fins.Air.Tdb=308.15     #Dry Bulb Temperature
Cycle.Condenser.Fins.Air.p=101.325      #Air pressure
Cycle.Condenser.Fins.Air.RH=0.51        #Relative Humidity
Cycle.Condenser.Fins.Air.RHmean=0.51
Cycle.Condenser.Fins.Air.FanPower=260

Cycle.Condenser.Ref=Cycle.Ref
Cycle.Condenser.Verbosity=0
    
#--------------------------------------
#      Cooling Coil parameters
#--------------------------------------
Cycle.CoolingCoil.Fins.Tubes.NTubes_per_bank=32
Cycle.CoolingCoil.Fins.Tubes.Nbank=3
Cycle.CoolingCoil.Fins.Tubes.Ncircuits=5
Cycle.CoolingCoil.Fins.Tubes.Ltube=0.452
Cycle.CoolingCoil.Fins.Tubes.OD=0.00913
Cycle.CoolingCoil.Fins.Tubes.ID=0.00849
Cycle.CoolingCoil.Fins.Tubes.Pl=0.0191
Cycle.CoolingCoil.Fins.Tubes.Pt=0.0254

Cycle.CoolingCoil.Fins.Fins.FPI=14.5
Cycle.CoolingCoil.Fins.Fins.Pd=0.001
Cycle.CoolingCoil.Fins.Fins.xf=0.001
Cycle.CoolingCoil.Fins.Fins.t=0.00011
Cycle.CoolingCoil.Fins.Fins.k_fin=237

Cycle.CoolingCoil.Fins.Air.Vdot_ha=0.56319
Cycle.CoolingCoil.Fins.Air.Tmean=297.039
Cycle.CoolingCoil.Fins.Air.Tdb=297.039
Cycle.CoolingCoil.Fins.Air.p=101.325
Cycle.CoolingCoil.Fins.Air.RH=0.5
Cycle.CoolingCoil.Fins.Air.RHmean=0.5
Cycle.CoolingCoil.Fins.Air.FanPower=438
       
params={
    'Ref_g': Cycle.SecLoopFluid,
    'pin_g': 200,
    'Verbosity':0,
    'mdot_g':0.38
}
Cycle.CoolingCoil.Update(**params)

params={
    'ID_i':0.0278,
    'OD_i':0.03415,
    'ID_o':0.045,
    'L':50,
    'pin_g':300,
    'Ref_r':Cycle.Ref,
    'Ref_g':Cycle.SecLoopFluid,
    'Verbosity':0
    }
Cycle.CoaxialIHX.Update(**params)

params={
    'pin_h':300,
    'Ref_h':Cycle.SecLoopFluid,
    'Ref_c':Cycle.Ref,
    
    #Geometric parameters
    'Bp' : 0.117,
    'Lp' : 0.300, #Center-to-center distance between ports
    'Nplates' : 46,
    'PlateAmplitude' : 0.001, #[m]
    'PlateThickness' : 0.0003, #[m]
    'PlateConductivity' : 15.0, #[W/m-K]
    'PlateWavelength' : 0.00628, #[m]
    'InclinationAngle' : 3.14159/3,#[rad]
    'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
    'Verbosity':0,
    'DT_sh':5
}
Cycle.PHEIHX.Update(**params)

params={
    'eta':0.5,  #Pump+motor efficiency
    'mdot_g':0.38, #Flow Rate kg/s
    'pin_g':300,
    'Ref_g':Cycle.SecLoopFluid,
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
    'pin': 300,
    'h_air':0.0000000001,
}

Cycle.LineSetSupply.Update(**params)
Cycle.LineSetReturn.Update(**params)
Cycle.LineSetSupply.OD=0.009525
Cycle.LineSetSupply.ID=0.007986
Cycle.LineSetReturn.OD=0.01905
Cycle.LineSetReturn.ID=0.017526

#Now solve
from time import time
t1=time()
Cycle.PreconditionedSolve()

#Outputs
print 'Took '+str(time()-t1)+' seconds to run Cycle model'
print 'Cycle coefficient of system performance is '+str(Cycle.COSP)
print 'Cycle refrigerant charge is '+str(Cycle.Charge)+' kg'