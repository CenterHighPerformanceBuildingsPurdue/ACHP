from __future__ import division
from Evaporator import EvaporatorClass
from FinCorrelations import FinInputs
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
FinsTubes.Air.p=101325             #updated from 101.325kPa to 101325Pa
FinsTubes.Air.RH=0.51
FinsTubes.Air.RHmean=0.51
FinsTubes.Air.FanPower=438

# Temporary value for inlet quality 
# (Quality is passed in from cycle model)
xin_r=0.15
    
kwargs={'Ref': 'R410A',
        'mdot_r':  0.0708,
        'psat_r':  PropsSI('P','T',282,'Q',1.0,'R410A'),
        'Vdot_a':  0.5663,
        'pin_a':   101325,      #Air pressure in Pa 
        'Tin_a':   299.8,
        'mdot_a': 0.662,
        'RHin_a': 0.51,
        'OD':0.009525,
        'ID':0.0089154,
        'Ncircuits':5,
        'Ltube':0.452,
        'Nbank':3,
        'NTubes_per_bank':32,
        'FanPower': 438,
        'Fins': FinsTubes,
        'xin_r': xin_r,
        'Verbosity': 0
}

Evap=EvaporatorClass(**kwargs)
Evap.Update(**kwargs)
Evap.Calculate()

print 'Evaporator charge total is', Evap.Charge, 'kg'
print 'Evaporator heat transfer rate is',Evap.Q,'W'
print 'Evaporator capacity (less fan power) is',Evap.Capacity,'W'
print 'Evaporator fraction of length in two-phase section',Evap.w_2phase,'W'
print 'Evaporator sensible heat ratio',Evap.SHR