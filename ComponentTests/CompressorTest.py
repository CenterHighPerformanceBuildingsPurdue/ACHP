from __future__ import division, print_function, absolute_import
from ACHP.Compressor import CompressorClass
import CoolProp as CP
from CoolProp.CoolProp import PropsSI

Ref = 'R134a'
Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
AS = CP.AbstractState(Backend, Ref)
    
kwds={
      'M':[217.3163128,5.094492028,-0.593170311,4.38E-02,-2.14E-02,1.04E-02,7.90E-05,-5.73E-05,1.79E-04,-8.08E-05],
      'P':[-561.3615705,-15.62601841,46.92506685,-0.217949552,0.435062616,-0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03],
      'Ref':Ref,
      'Tin_r':280,
      'pin_r':PropsSI('P','T',279,'Q',1,Ref),
      'pout_r':PropsSI('P','T',315,'Q',1,Ref),
      'fp':0.15, #Fraction of electrical power lost as heat to ambient
      'Vdot_ratio': 1.0, #Displacement Scale factor
      'shell_pressure': 'low-pressure',
      'Oil': 'POE32',
      'V_oil_sump': 0.0,
      'AS':AS
      }
Comp=CompressorClass(**kwds)
Comp.Calculate()

print ('Electrical power is: ' + str(Comp.W) + ' W')
print ('Actual mass flow rate is: ' + str(Comp.mdot_r) + ' kg/s')
print ('Isentropic Efficiency is: ' + str(Comp.eta_oi))
print ('Discharge Refrigerant Temperature is: ' + str(Comp.Tout_r) + ' K')
print (' ')

'''to print all the output, uncomment the next 2 lines'''
# for id, unit, value in Comp.OutputList():                
#     print (str(id) + ' = ' + str(value) + ' ' + str(unit))