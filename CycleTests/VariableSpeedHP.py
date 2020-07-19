'''
This example is for Variable Speed 5 tons Heat Pump unit 
with a rotary compressor in Cooling Mode
'''
from __future__ import division, absolute_import, print_function
from ACHP.Cycle import VariableSpeedHPClass 
from ACHP.Plots import PlotsClass
from ACHP.convert_units import *

#Instantiate the cycle class
Cycle=VariableSpeedHPClass()

#--------------------------------------
#         Cycle parameters
#--------------------------------------
Cycle.Verbosity = 1 #the idea here is to have different levels of debug output 
Cycle.ImposedVariable = 'Subcooling'
Cycle.CycleType = 'DX'
Cycle.DT_sc_target = 5.2
Cycle.Mode='AC'
Cycle.Ref='R410A'
Cycle.Backend='HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
Cycle.Oil = 'POE32'
Cycle.shell_pressure = 'low-pressure'
Cycle.CompModel = 'Miranda-Mendoza'
Cycle.CompType = 'Hitachi-RollingPiston'
Cycle.EvapSolver = 'Moving-Boundary' #choose the type of Evaporator solver scheme ('Moving-Boundary' or 'Finite-Element')
Cycle.EvapType = 'Fin-tube' #if EvapSolver = 'Moving-Boundary', choose the type of evaporator ('Fin-tube' or 'Micro-channel')
Cycle.CondSolver = 'Moving-Boundary' #choose the type of Condenser solver scheme ('Moving-Boundary' or 'Finite-Element')
Cycle.CondType = 'Fin-tube' #if CondSolver = 'Moving-Boundary', choose the type of condenser ('Fin-tube' or 'Micro-channel')
Cycle.Update()

#--------------------------------------
#     Charge correction parameters (activate by setting Cycle.ImposedVariable to 'Charge' and Cycle.ChargeMethod to either 'One-point' or 'Two-point')
#--------------------------------------
Cycle.C = 0 #[kg]
Cycle.K = 0 #[kg]
Cycle.w_ref = 0 #[-]

#--------------------------------------
#       Compressor parameters
#--------------------------------------
#A 4 ton cooling capacity compressor map of rolling piston compressor
params={
        'Ref': Cycle.Ref,
        'Tamb': 35 + 273.15,
        'fp':0, #Fraction of electrical power lost as heat to ambient
        'Vdisp': 47e-6, #Displacement volume, m^3
        'Vdot_ratio': 1.0, #Displacement Scale factor
        'N': 3600, #Rotational speed, rpm
        'a_etav':[1.,-0.32094114,0.00170576,-0.03695206], #Volumetric eff. coeff.
        'a_etais':[1.,0.24148598,0.37491376,0.07996134,-0.03366503], #Isentropic eff. coeff.
        'a_etaoi':[1.,-0.01360404,-0.07435644,0.18584579,0.63589724], #Overall eff. coeff.
        'shell_pressure': Cycle.shell_pressure,
        'Oil': Cycle.Oil,
        'V_oil_sump': 0.0,
        'Verbosity': 0.0
        }

Cycle.Compressor.Update(**params)

#--------------------------------------
#      Condenser parameters
#--------------------------------------
Cycle.Condenser.Fins.Tubes.NTubes_per_bank=24       #number of tubes per bank=row 
Cycle.Condenser.Fins.Tubes.Nbank=2                  #number of banks/rows 
Cycle.Condenser.Fins.Tubes.Ncircuits=8 
Cycle.Condenser.Fins.Tubes.Ltube=2.252
Cycle.Condenser.Fins.Tubes.OD=0.00913
Cycle.Condenser.Fins.Tubes.ID=0.00849
Cycle.Condenser.Fins.Tubes.Pl=0.0191                #distance between center of tubes in flow direction 
Cycle.Condenser.Fins.Tubes.Pt=0.0254                #distance between center of tubes orthogonal to flow direction
Cycle.Condenser.Fins.Tubes.kw=237                   #wall thermal conductivity (i.e pipe material)

Cycle.Condenser.Fins.Fins.FPI=25                    #Number of fins per inch
Cycle.Condenser.Fins.Fins.Pd=0.001                  #2* amplitude of wavy fin
Cycle.Condenser.Fins.Fins.xf=0.001                  #1/2 period of fin
Cycle.Condenser.Fins.Fins.t=0.00011                 #Thickness of fin material
Cycle.Condenser.Fins.Fins.k_fin=237                 #Thermal conductivity of fin material

Cycle.Condenser.Fins.Air.Vdot_ha=1.7934             #rated volumetric flowrate
Cycle.Condenser.Fins.Air.Tmean=35.0+273.1
Cycle.Condenser.Fins.Air.Tdb=35.0+273.15            # AHRI Standard 210/240 - Cooling test variable speed A2
Cycle.Condenser.Fins.Air.p=101325                   #Condenser Air pressure in Pa
Cycle.Condenser.Fins.Air.RH=0.51
Cycle.Condenser.Fins.Air.RHmean=0.51
Cycle.Condenser.Fins.Air.FanPower=260

Cycle.Condenser.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
Cycle.Condenser.Verbosity=0

#--------------------------------------
# Evaporator Parameters 
#--------------------------------------
Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=32
Cycle.Evaporator.Fins.Tubes.Nbank=3
Cycle.Evaporator.Fins.Tubes.Ltube=0.452
Cycle.Evaporator.Fins.Tubes.OD=0.00913
Cycle.Evaporator.Fins.Tubes.ID=0.00849
Cycle.Evaporator.Fins.Tubes.Pl=0.0191
Cycle.Evaporator.Fins.Tubes.Pt=0.0254
Cycle.Evaporator.Fins.Tubes.Ncircuits=5
Cycle.Evaporator.Fins.Tubes.kw=237                   #wall thermal conductivity (i.e pipe material)

Cycle.Evaporator.Fins.Fins.FPI=14.5
Cycle.Evaporator.Fins.Fins.Pd=0.001
Cycle.Evaporator.Fins.Fins.xf=0.001
Cycle.Evaporator.Fins.Fins.t=0.00011
Cycle.Evaporator.Fins.Fins.k_fin=237

Cycle.Evaporator.Fins.Air.Vdot_ha=0.56319
Cycle.Evaporator.Fins.Air.Tmean=26.7+273.15
Cycle.Evaporator.Fins.Air.Tdb=26.7+273.15           # AHRI Standard 210/240 - Cooling test variable speed A2
Cycle.Evaporator.Fins.Air.p=101325                                              #Evaporator Air pressure in Pa
Cycle.Evaporator.Fins.Air.RH=0.5
Cycle.Evaporator.Fins.Air.RHmean=0.5
Cycle.Evaporator.Fins.Air.FanPower=438

Cycle.Evaporator.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
Cycle.Evaporator.Verbosity=0
Cycle.Evaporator.DT_sh=6.7                    #target superheat

# ----------------------------------
#       Expansion device Parameters
# ----------------------------------
params={
        'ExpType':'Linear-TXV',     #expansion device type
        'Tsh_static':4,             #static superheat
        'Tsh_max':8,                #maximum superheat
        'D':0.0006604,              #inside diameter [m]
        'C':1.2656e-6,              #constant from manufacturer [m^2/K]
        'Adj':0.7630,               #Adjust the diameter (tuning factor)
    }
Cycle.ExpDev.Update(**params)

# ----------------------------------
#       Liquid Line Set Parameters
# ----------------------------------
params={
        'L':7.6,
        'k_tube':0.19,
        't_insul':0.02,
        'k_insul':0.036,
        'T_air':297,
        'h_air':0.0000000001,
        'LineSetOption': 'On'
        }

Cycle.LineSetLiquid.Update(**params)
Cycle.LineSetLiquid.OD=in2m(3/8)    #3/8"
Cycle.LineSetLiquid.ID=in2m(0.25)   #1/4"

# ----------------------------------
#       Suction Line Set Parameters
# ----------------------------------
params={
        'L':7.6,
        'k_tube':0.19,
        't_insul':0.02,
        'k_insul':0.036,
        'T_air':297,
        'h_air':0.0000000001,
        'LineSetOption': 'On'
        }

Cycle.LineSetSuction.Update(**params)
Cycle.LineSetSuction.OD=in2m(1.125)    #1-1/8"
Cycle.LineSetSuction.ID=in2m(1.025)

# ----------------------------------
#       Line Set Discharge Parameters
# ----------------------------------
params={
        'L':0.3,                #tube length in m
        'k_tube':0.19,
        't_insul':0, #no insulation
        'k_insul':0.036,
        'T_air':297,
        'h_air':0.0000000001,
        'LineSetOption': 'On'
        }
  
Cycle.LineSetDischarge.Update(**params)
Cycle.LineSetDischarge.OD=in2m(1.125)    #1-1/8"
Cycle.LineSetDischarge.ID=in2m(1.025)

# ----------------------------------
#       Liquid Receiver Parameters
# ----------------------------------
params={
        'Tamb':35+273.15,
        'ID':0.115,
        'OD': 0.1236,
        'h_tank': 0.274,
        }
  
Cycle.SuctionAccumulator.Update(**params)


#Now solve
from time import time
t1=time()
Cycle.PreconditionedSolve()
print ('Took '+str(time()-t1)+' seconds to run Cycle model')
print ('Cycle COP is '+str(Cycle.COSP))
print ('Cycle refrigerant charge is '+str(Cycle.Charge)+' kg')
print (Cycle.OutputList())
#Now do cycle plotting
plot = PlotsClass()
plot.TSOverlay(Cycle)
plot.PHOverlay(Cycle)