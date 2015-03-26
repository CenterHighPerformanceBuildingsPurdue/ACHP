from CoolProp.CoolProp import Props
from FinCorrelations import FinInputs
from Evaporator import EvaporatorClass

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
FinsTubes.Air.Tdb=299.8
FinsTubes.Air.p=101.325
FinsTubes.Air.RH=0.51
FinsTubes.Air.FanPower=438

kwargs={'Ref': 'R410A',
    'mdot_r':  0.0708,
    'psat_r':  Props('P','T',282,'Q',1.0,'R410A'),
    'Fins': FinsTubes,
    'hin_r':Props('H','T',282,'Q',0.15,'R410A')*1000, #[J/kg]
    'Verbosity':0
}

Evap=EvaporatorClass(**kwargs)
Evap.Update(**kwargs)
Evap.Calculate()

print 'Evaporator heat trasfer rate is',Evap.Q,'W'
print 'Evaporator capacity (less fan power) is',Evap.Capacity,'W'
print 'Evaporator fraction of length in two-phase section',Evap.w_2phase,'W'
print 'Evaporator sensible heat ratio',Evap.SHR