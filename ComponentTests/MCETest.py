from __future__ import division, print_function, absolute_import
from ACHP.FinCorrelations import FinInputs
from ACHP.Evaporator import EvaporatorClass
from ACHP.MultiCircuitEvaporator import MultiCircuitEvaporatorClass
import CoolProp as CP
from CoolProp.CoolProp import PropsSI
import numpy as np

FinsTubes=FinInputs()

FinsTubes.Tubes.NTubes_per_bank=32
FinsTubes.Tubes.Ncircuits=5
FinsTubes.Tubes.Nbank=3
FinsTubes.Tubes.Ltube=0.452
FinsTubes.Tubes.OD=0.009525
FinsTubes.Tubes.ID=0.0089154
FinsTubes.Tubes.Pl=0.0254
FinsTubes.Tubes.Pt=0.0219964
FinsTubes.Tubes.kw=237                   #Wall thermal conductivity

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

Ref = 'R410A'
Backend = 'TTSE&HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
AS = CP.AbstractState(Backend, Ref)

# This uses the normal baseline evaporator model
kwargs={'AS': AS,
        'mdot_r': 0.0708,
        'psat_r': PropsSI('P','T',282.0,'Q',1.0,Ref),
        'Fins': FinsTubes,
        'FinsType': 'WavyLouveredFins', #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'hin_r': PropsSI('H','P',PropsSI('P','T',282.0,'Q',1.0,Ref),'Q',0.15,Ref),
        'Verbosity':0,
        }
Evap=EvaporatorClass(**kwargs)
Evap.Calculate()
print ('Evap Q=' + str(Evap.Q) + ' W')

# This uses the multi-circuited evaporator model but with no mal-distribution
kwargs={'AS': AS,
        'mdot_r': 0.0708,
        'psat_r': PropsSI('P','T',282.0,'Q',1.0,Ref),
        'Fins': FinsTubes,
        'FinsType': 'WavyLouveredFins', #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'hin_r': PropsSI('H','P',PropsSI('P','T',282.0,'Q',1.0,Ref),'Q',0.15,Ref),
        'Verbosity':0,
        }
MCE=MultiCircuitEvaporatorClass(**kwargs)
MCE.Calculate()
print ('MCE Q='+str(MCE.Q)+' W w/o mal-distribution')

# Not exactly the same since
# This uses the multi-circuited evaporator model with mal-distribution of 
# refrigerant, refrigerant quality, and air volumetric flow rate 
kwargs={'AS': AS,
        'mdot_r': 0.0708,
        'psat_r': PropsSI('P','T',282.0,'Q',1.0,Ref),
        'mdot_r_coeffs': [0.3,0.2,0.1,0.2,0.2],
        'mdot_v_coeffs': [0.4,0.2,0.1,0.2,0.1],
        'Vdot_ha_coeffs': [0.3,0.2,0.2,0.2,0.1],
        'Fins': FinsTubes,
        'FinsType': 'WavyLouveredFins', #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'hin_r': PropsSI('H','P',PropsSI('P','T',282.0,'Q',1.0,Ref),'Q',0.15,Ref),
        'Verbosity':0,
        }
MCE=MultiCircuitEvaporatorClass(**kwargs)
MCE.Calculate()
print ('MCE Q='+str(MCE.Q)+' W w/ mal-distribution')


# Another way to express air maldirtribution for the last example, 
# using a vector of Vdot_ha and m_dot_r instead of vector of coefficients
FinsTubes.Air.Vdot_ha=0.5663*np.array([0.3,0.2,0.2,0.2,0.1])
kwargs={'AS': AS,
        'mdot_r': 0.0708*np.array([0.3,0.2,0.1,0.2,0.2]),
        'mdot_v_coeffs': [0.4,0.2,0.1,0.2,0.1],
        'psat_r': PropsSI('P','T',282.0,'Q',1.0,Ref),
        'Fins': FinsTubes,
        'FinsType': 'WavyLouveredFins', #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'hin_r': PropsSI('H','P',PropsSI('P','T',282.0,'Q',1.0,Ref),'Q',0.15,Ref),
        'Verbosity':0,
        }
MCE=MultiCircuitEvaporatorClass(**kwargs)
MCE.Calculate()
print ('MCE Q='+str(MCE.Q)+' W w/ mal-distribution')
