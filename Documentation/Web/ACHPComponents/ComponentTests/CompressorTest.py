from Compressor import CompressorClass
from CoolProp.CoolProp import Props
kwds={
      'M':[217.3163128,5.094492028,-0.593170311,4.38E-02,
        -2.14E-02,1.04E-02,7.90E-05,-5.73E-05,1.79E-04,-8.08E-05],
      'P':[-561.3615705,-15.62601841,46.92506685,-0.217949552,
        0.435062616,-0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03],
      'Ref':'R134a',
      'Tin_r':280,
      'pin_r':Props('P','T',279,'Q',1,'R134a'),
      'pout_r':Props('P','T',315,'Q',1,'R134a'),
      'fp':0.15, #Fraction of electrical power lost as heat to ambient
      'Vdot_ratio': 1.0 #Displacement Scale factor
      }
Comp=CompressorClass(**kwds)
Comp.Calculate()

print 'Electrical power is: ' + str(Comp.W) + ' W'
print 'Actual mass flow rate is: ' + str(Comp.mdot_r) + ' kg/s'