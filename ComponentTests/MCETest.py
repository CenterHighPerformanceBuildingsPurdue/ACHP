'''
Created on Mar 22, 2015

@author: AmmarBahman
'''

from FinCorrelations import FinInputs
from Evaporator import EvaporatorClass
from MultiCircuitEvaporator import MultiCircuitEvaporatorClass
from CoolProp.CoolProp import PropsSI

FinsTubes=FinInputs()

FinsTubes.Tubes.NTubes_per_bank=32
FinsTubes.Tubes.Ncircuits=5
FinsTubes.Tubes.Nbank=3
FinsTubes.Tubes.Ltube=0.452
FinsTubes.Tubes.OD=0.009525
FinsTubes.Tubes.ID=0.0089154
FinsTubes.Tubes.Pl=0.0254
FinsTubes.Tubes.Pt=0.0219964

FinsTubes.Fins.FPI=14.5
FinsTubes.Fins.Pd=0.001
FinsTubes.Fins.xf=0.001
FinsTubes.Fins.t=0.00011
FinsTubes.Fins.k_fin=237

FinsTubes.Air.Vdot_ha=0.5663
FinsTubes.Air.Tmean=299.8
FinsTubes.Air.Tdb=299.8
FinsTubes.Air.p=101325                  #Air pressure in Pa
FinsTubes.Air.RH=0.51
FinsTubes.Air.RHmean=0.51
FinsTubes.Air.FanPower=438

#This uses the normal baseline evaporator model
kwargs={'Ref': 'R410A',
        'mdot_r': 0.0708,
        'psat_r': PropsSI('P','T',282.0,'Q',1.0,'R410A'),
        'Fins': FinsTubes,
        'FinsType': 'WavyLouveredFins',                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'hin_r': PropsSI('H','T',282.0,'Q',0.15,'R410A'), #*1000
        'Verbosity':0
        }
Evap=EvaporatorClass(**kwargs)
Evap.Calculate()
print 'Evap Q=' + str(Evap.Q) + ' W'

#This uses the multi-circuited evaporator model but with no mal-distribution
kwargs={'Ref': 'R410A',
        'mdot_r': 0.0708,
        'psat_r': PropsSI('P','T',282.0,'Q',1.0,'R410A'),
        'Fins': FinsTubes,
        'FinsType': 'WavyLouveredFins',                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'hin_r': PropsSI('H','T',282.0,'Q',0.15,'R410A'), #*1000
        'Verbosity':0
        }
MCE=MultiCircuitEvaporatorClass(**kwargs)
MCE.Calculate()
print 'MCE Q='+str(MCE.Q)+' W w/o mal-distribution' 

#Not exactly the same since
# This uses the multi-circuited evaporator model with mal-distribution of 
# refrigerant, refrigerant quality, and air volumetric flow rate 
kwargs={'Ref': 'R410A',
        'mdot_r': 0.0708,
        'psat_r': PropsSI('P','T',282.0,'Q',1.0,'R410A'),
        'mdot_r_coeffs': [0.3,0.2,0.1,0.2,0.2],
        'mdot_v_coeffs': [0.4,0.2,0.1,0.2,0.1],
        'Vdot_ha_coeffs': [0.3,0.2,0.2,0.2,0.1],
        'Fins': FinsTubes,
        'FinsType': 'WavyLouveredFins',                   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'hin_r': PropsSI('H','T',282.0,'Q',0.15,'R410A'), #*1000
        'Verbosity':0
        }
MCE=MultiCircuitEvaporatorClass(**kwargs)
MCE.Calculate()
print 'MCE Q='+str(MCE.Q)+' W w/ mal-distribution'
