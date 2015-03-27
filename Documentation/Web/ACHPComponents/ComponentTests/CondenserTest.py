from CoolProp.CoolProp import Props
from Condenser import CondenserClass
from FinCorrelations import FinInputs

Fins=FinInputs()
Fins.Tubes.NTubes_per_bank=41  #number of tubes per bank=row
Fins.Tubes.Nbank=1             #number of baks/rows
Fins.Tubes.Ncircuits=5             #number of baks/rows
Fins.Tubes.Ltube=2.286
Fins.Tubes.OD=0.007
Fins.Tubes.ID=0.0063904
Fins.Tubes.Pl=0.0191  #distance between center of tubes in flow direction                                                
Fins.Tubes.Pt=0.0222  #distance between center of tubes orthogonal to flow direction

Fins.Fins.FPI=25      #Number of fins per inch
Fins.Fins.Pd=0.001    #2* amplitude of wavy fin
Fins.Fins.xf=0.001    #1/2 period of fin
Fins.Fins.t=0.00011   #Thickness of fin material
Fins.Fins.k_fin=237   #Thermal conductivity of fin material

Fins.Air.Vdot_ha=1.7934 #rated volumetric flowrate
Fins.Air.Tdb=308.15     #Dry Bulb Temperature
Fins.Air.p=101.325      #Air pressure
Fins.Air.RH=0.51        #Relative Humidity
Fins.Air.FanPower=160    

params={
    'Ref': 'R410A',
    'mdot_r':  0.0708,
    'Tin_r':   333.15,
    'psat_r':  Props('P','T',323.15,'Q',1.0,'R410A'),
    'Fins':Fins,
    'Verbosity':0
}
Cond=CondenserClass(**params)
Cond.Calculate()

print 'Heat transfer rate in condenser is', Cond.Q,'W'
print 'Fraction of circuit length in subcooled section is',Cond.w_subcool
print 'Fraction of circuit length in twophase section is',Cond.w_2phase
print 'Fraction of circuit length in superheated section is',Cond.w_superheat