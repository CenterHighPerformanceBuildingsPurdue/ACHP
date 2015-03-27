from Cycle import DXCycleClass

#Instantiate the cycle class
Cycle=DXCycleClass()

#--------------------------------------
#       Cycle parameters
#--------------------------------------
Cycle.Verbosity = 0 #the idea here is to have different levels of debug output
Cycle.ImposedVariable = 'Subcooling'                                                           
Cycle.DT_sc_target = 7.0
Cycle.Mode='AC' 
Cycle.Ref='R410A'

#--------------------------------------
#       Compressor parameters
#--------------------------------------

#A 3 ton cooling capacity compressor map 
M=[217.3163128,5.094492028,-0.593170311,4.38E-02,-2.14E-02,
    1.04E-02,7.90E-05,-5.73E-05,1.79E-04,-8.08E-05]
P=[-561.3615705,-15.62601841,46.92506685,-0.217949552,
    0.435062616,-0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03]

params={
    'M':M,
    'P':P,
    'Ref':Cycle.Ref,  #Refrigerant
    'fp':0.0, #Fraction of electrical power lost as heat to ambient
    'Vdot_ratio': 1.0, #Displacement Scale factor                  
    'Verbosity': 0, # How verbose should the debugging be [0-10]
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
#     Evaporator Parameters
#--------------------------------------
Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=32
Cycle.Evaporator.Fins.Tubes.Nbank=3
Cycle.Evaporator.Fins.Tubes.Ltube=0.452
Cycle.Evaporator.Fins.Tubes.OD=0.00913
Cycle.Evaporator.Fins.Tubes.ID=0.00849
Cycle.Evaporator.Fins.Tubes.Pl=0.0191
Cycle.Evaporator.Fins.Tubes.Pt=0.0254
Cycle.Evaporator.Fins.Tubes.Ncircuits=5

Cycle.Evaporator.Fins.Fins.FPI=14.5
Cycle.Evaporator.Fins.Fins.Pd=0.001
Cycle.Evaporator.Fins.Fins.xf=0.001
Cycle.Evaporator.Fins.Fins.t=0.00011
Cycle.Evaporator.Fins.Fins.k_fin=237

Cycle.Evaporator.Fins.Air.Vdot_ha=0.56319
Cycle.Evaporator.Fins.Air.Tmean=297.039
Cycle.Evaporator.Fins.Air.Tdb=297.039
Cycle.Evaporator.Fins.Air.p=101.325
Cycle.Evaporator.Fins.Air.RH=0.5
Cycle.Evaporator.Fins.Air.RHmean=0.5
Cycle.Evaporator.Fins.Air.FanPower=438
       
Cycle.Evaporator.Ref=Cycle.Ref
Cycle.Evaporator.Verbosity=0
Cycle.Evaporator.DT_sh=5

# ----------------------------------
#       Line Set Parameters
# ----------------------------------
params={
    'L':7.6,
    'k_tube':0.19,
    't_insul':0.02,
    'k_insul':0.036,
    'T_air':297,
    'Ref': Cycle.Ref,
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
print 'Took '+str(time()-t1)+' seconds to run Cycle model'
print 'Cycle coefficient of system performance is '+str(Cycle.COSP)
print 'Cycle refrigerant charge is '+str(Cycle.Charge)+' kg'
