from Cycle import DXCycleClass,F2K

#Instantiate the class
Cycle=DXCycleClass()

#--------------------------------------
#       Cycle parameters
#--------------------------------------
Cycle.Verbosity = 0 #the idea here is to have different levels of debug output
Cycle.ImposedVariable = 'Subcooling' #or it could be 'Charge'
Cycle.DT_sc_target = 7.0
#Cycle.Charge_target = 3.3 #uncomment for use with imposed charge
Cycle.Mode='HP' 
Cycle.Ref='R410A'

#--------------------------------------
#       Compressor parameters
#--------------------------------------

#A few 3 ton cooling capacity compressor maps 
M=[217.3163128,5.094492028,-0.593170311,4.38E-02,
    -2.14E-02,1.04E-02,7.90E-05,-5.73E-05,1.79E-04,-8.08E-05]
P=[-561.3615705,-15.62601841,46.92506685,-0.217949552,
    0.435062616,-0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03]

params={
'M':M,
'P':P,
'Ref':Cycle.Ref,                                                              #refrigerant
'fp':0.0, #Fraction of electrical power lost as heat to ambient            #shell heat loss
'Vdot_ratio': 1.0, #Displacement Scale factor                               #up- or downsize compressor (1=original)
'Verbosity': 0, # How verbose should the debugging statements be [0 to 10]                        
}
Cycle.Compressor.Update(**params)                                          

#--------------------------------------
#      Condenser parameters 
#--------------------------------------
Cycle.Condenser.Fins.Tubes.NTubes_per_bank=32
Cycle.Condenser.Fins.Tubes.Nbank=3
Cycle.Condenser.Fins.Tubes.Ncircuits=6
Cycle.Condenser.Fins.Tubes.Ltube=0.452
Cycle.Condenser.Fins.Tubes.OD=0.009525
Cycle.Condenser.Fins.Tubes.ID=0.0089154
Cycle.Condenser.Fins.Tubes.Pl=0.0254
Cycle.Condenser.Fins.Tubes.Pt=0.0219964

Cycle.Condenser.Fins.Fins.FPI=14.5
Cycle.Condenser.Fins.Fins.Pd=0.001
Cycle.Condenser.Fins.Fins.xf=0.001
Cycle.Condenser.Fins.Fins.t=0.00011
Cycle.Condenser.Fins.Fins.k_fin=237

Cycle.Condenser.Fins.Air.Vdot_ha=0.5663
Cycle.Condenser.Fins.Air.Tmean=F2K(70)
Cycle.Condenser.Fins.Air.Tdb=F2K(70)
Cycle.Condenser.Fins.Air.p=101.325
Cycle.Condenser.Fins.Air.RH=0.51
Cycle.Condenser.Fins.Air.RHmean=0.51
Cycle.Condenser.Fins.Air.FanPower=438

Cycle.Condenser.Ref=Cycle.Ref
Cycle.Condenser.Verbosity=0

#--------------------------------------
#           Evaporator parameters
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
Cycle.Evaporator.Fins.Air.p=101.325      #Air pressure
Cycle.Evaporator.Fins.Air.RH=0.51        #Relative Humidity
Cycle.Evaporator.Fins.Air.RHmean=0.51
Cycle.Evaporator.Fins.Air.FanPower=160

Cycle.Evaporator.Ref=Cycle.Ref
Cycle.Evaporator.Verbosity=0
Cycle.Evaporator.DT_sh=5

# ----------------------------------
#       Line Set parameters
# ----------------------------------
params={
'L':7.6,
'k_tube':0.19,
't_insul':0.02,
'k_insul':0.036,
'T_air':297,
'Ref': Cycle.Ref,
'h_air':6,
}

Cycle.LineSetSupply.Update(**params)
Cycle.LineSetReturn.Update(**params)
Cycle.LineSetSupply.OD=0.01905 
Cycle.LineSetSupply.ID=0.017526
Cycle.LineSetReturn.OD=0.009525
Cycle.LineSetReturn.ID=0.007986

#Now solve
from time import time
t1=time()
Cycle.PreconditionedSolve()

#Outputs
print 'Took '+str(time()-t1)+' seconds to run Cycle model'
print 'Cycle coefficient of system performance is '+str(Cycle.COSP)
print 'Cycle refrigerant charge is '+str(Cycle.Charge)+' kg'