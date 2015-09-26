from ACHP.CoolingCoil import CoolingCoilClass
from ACHP.FinCorrelations import FinInputs

FinsTubes=FinInputs()
FinsTubes.Tubes.NTubes_per_bank=32
FinsTubes.Tubes.Nbank=3
FinsTubes.Tubes.Ncircuits=5
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
FinsTubes.Air.Tdb= 299.8
FinsTubes.Air.p=101325              #Air pressure in Pa
FinsTubes.Air.RH=0.51
FinsTubes.Air.RHmean=0.51
FinsTubes.Air.FanPower=438          #fan power in Watts

# Here are two equivalent methods for setting parameters
# 1. Create an empty instance of the class, then set parameters CC=CoolingCoilClass()
CC=CoolingCoilClass()
CC.Fins = FinsTubes
CC.FinsType = 'WavyLouveredFins'    #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
CC.Ref_g = 'Water'
CC.mdot_g = 0.15
CC.Tin_g = 278
CC.pin_g = 300000                   #Refrigerant vapor pressure in Pa
CC.Verbosity = 3
CC.Calculate()

print "Method 1:"
print "Cooling Coil Q: " + str(CC.Q) + " W"
print "Cooling Coil SHR: " + str(CC.SHR) + " "

# 2. Build a dictionary of values, then use that to initialize the class
kwds={'Fins':FinsTubes,
      'FinsType': 'WavyLouveredFins',   #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
      'Ref_g': 'Water',
      'mdot_g': 0.15,
      'Tin_g': 278,
      'pin_g': 300000,               #Refrigerant vapor pressure in Pa
      'Verbosity':3}
CC2=CoolingCoilClass(**kwds)
CC2.Calculate()

print "Method 2:"
print "Cooling Coil Q: " + str(CC2.Q) + " W"
print "Cooling Coil SHR: " + str(CC2.SHR) + " "
