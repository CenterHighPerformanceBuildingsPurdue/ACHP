'''This code is for secondary loop cycle in HP Mode'''
from __future__ import division, absolute_import, print_function
from ACHP.Cycle import SecondaryCycleClass 
from ACHP.convert_units import F2K 

#Instantiate the class
Cycle=SecondaryCycleClass()

#--------------------------------------
# Cycle parameters
#--------------------------------------

Cycle.Verbosity = 0 #the idea here is to have different levels of debug output
Cycle.ImposedVariable = 'Subcooling' #or this could be 'Charge' for imposed charge
Cycle.DT_sc_target = 7.0
#Cycle.Charge_target = 2.4  #Needed if charge is imposed, not otherwise
Cycle.Ref='R410A'
Cycle.Backend='HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
Cycle.Oil='POE32'
Cycle.SecLoopFluid = 'MEG'
Cycle.MassFrac_SLF = 0.21 #Mass fraction of incompressible SecLoopFluid [i.e MEG-20%]
Cycle.Backend_SLF = 'INCOMP' #backend of SecLoopFluid
Cycle.IHXType = 'PHE'# or could be 'Coaxial'
Cycle.shell_pressure='low-pressure'
Cycle.Mode='HP'
Cycle.DT_sh = 5 #superheat
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
#--------------------------------------
#       Compressor parameters
#--------------------------------------
#--------------------------------------

#A 3 ton cooling capacity compressor map
if Cycle.Ref=='R410A' or Cycle.Ref == 'R410A.mix':
    M=[217.3163128,5.094492028,-0.593170311,4.38E-02,-2.14E-02,1.04E-02, 7.90E-05,-5.73E-05,1.79E-04,-8.08E-05]
    P=[-561.3615705,-15.62601841,46.92506685,-0.217949552,0.435062616, -0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03]

params={
        'M':M,
        'P':P,
        'Ref':Cycle.Ref, #refrigerant
        'Oil':Cycle.Oil, #Compressor lubricant oil
        'shell_pressure':Cycle.shell_pressure, #Compressor shell pressure
        'fp':0.15, #Fraction of electrical power lost as heat to ambient
        'Vdot_ratio': 1.0, #Displacement Scale factor to up- or downsize compressor (1=original)
        'V_oil_sump':0, #Volume of oil in the sump
        'Verbosity': 0, # How verbose should the debugging statements be [0 to 10]
        }

Cycle.Compressor.Update(**params)

#--------------------------------------
# Evaportor parameters
#--------------------------------------
Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=24   #number of tubes per bank=row
Cycle.Evaporator.Fins.Tubes.Nbank=1              #number of banks/rows
Cycle.Evaporator.Fins.Tubes.Ncircuits=3
Cycle.Evaporator.Fins.Tubes.Ltube=2.252
Cycle.Evaporator.Fins.Tubes.OD=0.00913
Cycle.Evaporator.Fins.Tubes.ID=0.00849
Cycle.Evaporator.Fins.Tubes.Pl=0.0191            #distance between center of tubes in flow direction
Cycle.Evaporator.Fins.Tubes.Pt=0.0254            #distance between center of tubes orthogonal to flow direction
Cycle.Evaporator.Fins.Tubes.kw=237               #wall thermal conductivity (i.e pipe material)

Cycle.Evaporator.Fins.Fins.FPI=25                #Number of fins per inch
Cycle.Evaporator.Fins.Fins.Pd=0.001              #2* amplitude of wavy fin
Cycle.Evaporator.Fins.Fins.xf=0.001              #1/2 period of fin
Cycle.Evaporator.Fins.Fins.t=0.00011             #Thickness of fin material
Cycle.Evaporator.Fins.Fins.k_fin=237             #Thermal conductivity of fin material

Cycle.Evaporator.Fins.Air.Vdot_ha=0.5663
Cycle.Evaporator.Fins.Air.Tmean=F2K(40)
Cycle.Evaporator.Fins.Air.Tdb=F2K(40)
Cycle.Evaporator.Fins.Air.p=101325               #Evaporator Air Pressure in Pa
Cycle.Evaporator.Fins.Air.RH=0.51
Cycle.Evaporator.Fins.Air.RHmean=0.51
Cycle.Evaporator.Fins.Air.FanPower=438

Cycle.Evaporator.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
Cycle.Evaporator.Verbosity=0
Cycle.Evaporator.DT_sh=5            #target superheat

#--------------------------------------
# Cooling Coil parameters
#--------------------------------------
Cycle.CoolingCoil.Fins.Tubes.NTubes_per_bank=32
Cycle.CoolingCoil.Fins.Tubes.Nbank=3
Cycle.CoolingCoil.Fins.Tubes.Ncircuits=5
Cycle.CoolingCoil.Fins.Tubes.Ltube=0.452
Cycle.CoolingCoil.Fins.Tubes.OD=0.00913
Cycle.CoolingCoil.Fins.Tubes.ID=0.00849
Cycle.CoolingCoil.Fins.Tubes.Pl=0.0191
Cycle.CoolingCoil.Fins.Tubes.Pt=0.0254
Cycle.CoolingCoil.Fins.Tubes.kw=237               #wall thermal conductivity (i.e pipe material)

Cycle.CoolingCoil.Fins.Fins.FPI=14.5
Cycle.CoolingCoil.Fins.Fins.Pd=0.001
Cycle.CoolingCoil.Fins.Fins.xf=0.001
Cycle.CoolingCoil.Fins.Fins.t=0.00011
Cycle.CoolingCoil.Fins.Fins.k_fin=237

Cycle.CoolingCoil.Fins.Air.Vdot_ha=0.56319
Cycle.CoolingCoil.Fins.Air.Tmean=F2K(70)
Cycle.CoolingCoil.Fins.Air.Tdb=F2K(70)
Cycle.CoolingCoil.Fins.Air.p=101325
Cycle.CoolingCoil.Fins.Air.RH=0.5
Cycle.CoolingCoil.Fins.Air.RHmean=0.5

Cycle.CoolingCoil.Fins.Air.FanPower=438

params={
        'pin_g': 300000,
        'Verbosity':0,
        'mdot_g':0.38,
        'FinsType': 'WavyLouveredFins',                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        }
Cycle.CoolingCoil.Update(**params)

# ----------------------------------
#       IHX Plate Heat Exchanger
# ----------------------------------
params={
        'pin_c':300000,
    
        #Geometric parameters
        'HXType':'Plate-HX',
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

# ----------------------------------
# ----------------------------------
#       Seconday Loop Pump
# ----------------------------------
# ----------------------------------
params={
        'eta':0.5, #Pump+motor efficiency
        'mdot_g':0.38, #Flow Rate kg/s
        'pin_g':300000,
        'Verbosity':0,
        }
Cycle.Pump.Update(**params)

# ----------------------------------
#       Line Set Supply and Return Parameters
# ----------------------------------
params={
        'L':5,
        'k_tube':0.19,
        't_insul':0.02,
        'k_insul':0.036,
        'T_air':297,
        'pin': 300000,
        'h_air':0.0000000001,
        'LineSetOption': 'Off'
        }
Cycle.LineSetSupply.Update(**params)
Cycle.LineSetReturn.Update(**params)
Cycle.LineSetSupply.OD=0.01905
Cycle.LineSetSupply.ID=0.017526
Cycle.LineSetReturn.OD=0.01905
Cycle.LineSetReturn.ID=0.017526

# --------------------------------------------
#       Line Set Parameters Suction and Liquid
# --------------------------------------------
params={
        'L':7.6,
        'k_tube':0.19,
        't_insul':0.02,
        'k_insul':0.036,
        'T_air':297,
        'h_air':0.0000000001,
        'LineSetOption': 'Off'
        }

Cycle.LineSetLiquid.Update(**params)
Cycle.LineSetSuction.Update(**params)
Cycle.LineSetLiquid.OD=0.009525
Cycle.LineSetLiquid.ID=0.007986
Cycle.LineSetSuction.OD=0.01905
Cycle.LineSetSuction.ID=0.017526

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
        'LineSetOption': 'Off'
        }
  
Cycle.LineSetDischarge.Update(**params)
Cycle.LineSetDischarge.OD=0.009525
Cycle.LineSetDischarge.ID=0.007986


#Now solve
from time import time 
t1=time()
Cycle.PreconditionedSolve()

#Outputs
print ('Took '+str(time()-t1)+' seconds to run Cycle model')
print ('Cycle coefficient of system performance is '+str(Cycle.COSP))
print ('Cycle refrigerant charge is '+str(Cycle.Charge)+' kg')