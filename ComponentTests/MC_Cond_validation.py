from __future__ import division, print_function, absolute_import
import CoolProp as CP
from CoolProp.CoolProp import PropsSI
from ACHP.MicroChannelCondenser import MicroCondenserClass
from ACHP.MicroFinCorrelations import MicroFinInputs
from ACHP.convert_units import in2m, mm2m, cfm2cms, F2K, kPa2Pa, C2K

Fins=MicroFinInputs()
Fins.Tubes.NTubes=30               #Number of tubes (per bank for now!)
Fins.Tubes.Nbank=2                 #Number of banks (set to 1 for now!)
Fins.Tubes.Npass=2                 #Number of passes (per bank) #averaged if not even
Fins.Tubes.Nports=11               #Number of rectangular ports
Fins.Tubes.Ltube=in2m(18.24)       #length of a single tube
Fins.Tubes.Td=in2m(1)              #Tube outside width (depth)
Fins.Tubes.Ht=in2m(0.072)          #Tube outside height (major diameter)
Fins.Tubes.b=in2m(0.488)           #Tube spacing   
Fins.Tubes.tw=in2m(0.015)          #Tube wall thickness
Fins.Tubes.twp=in2m(0.016)         #Port (channel) wall thickness     
Fins.Tubes.beta=1.7675             #Port (channel) aspect ratio (=width/height)
Fins.Tubes.kw=237                  #Wall thermal conductivity

Fins.Fins.FPI=13                   #Fin per inch
Fins.Fins.Lf=in2m(1)               #Fin length = tube outside width in this HX
Fins.Fins.t=in2m(0.0045)           ##measured## #Fin thickness
Fins.Fins.k_fin=117                #Fin thermal conductivity for pure Aluminum
    
Fins.Air.Vdot_ha=cfm2cms(1500)     #Air volume flow rate in m^3/s
Fins.Air.Tdb=F2K(125)              #Air inlet temperature, K
Fins.Air.p=101325                  #Air inlet pressure in Pa
Fins.Air.RH=0.199                  #Air inlet relative humidity
Fins.Air.FanPower=855              #Fan power, Watts
    
Fins.Louvers.Lalpha=25             ##estimated## #Louver angle, in degree
Fins.Louvers.lp=mm2m(1.12)         ##measured## #Louver pitch

Ref = 'R407C'
Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
AS = CP.AbstractState(Backend, Ref)

params={
        'AS': AS,
        'mdot_r': 0.04472,
        'Tin_r': C2K(110),
        'psat_r': kPa2Pa(3108), 
        'Fins': Fins,
        'FinsType': 'MultiLouveredMicroFins',
        'Verbosity':0,
        }


Cond=MicroCondenserClass(**params)
Cond.Calculate()

print ('Heat transfer rate in condenser is', Cond.Q,'W')
print ('Heat transfer rate in condenser (superheat section) is',Cond.Q_superheat,'W')
print ('Heat transfer rate in condenser (twophase section) is',Cond.Q_2phase,'W')
print ('Heat transfer rate in condenser (subcooled section) is',Cond.Q_subcool,'W')
print ('Fraction of circuit length in superheated section is',Cond.w_superheat)
print ('Fraction of circuit length in twophase section is',Cond.w_2phase)
print ('Fraction of circuit length in subcooled section is',Cond.w_subcool)
# for id, unit, value in Cond.OutputList():
#     print (str(id) + ' = ' + str(value) + ' ' + str(unit))