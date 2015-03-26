
from CoolProp.CoolProp import Props
import pylab,numpy as np

import sys
sys.path.append('../../..')
from PyACHP.Correlations import ShahCondensation_Average

x=np.linspace(0,1,1000)
h=np.zeros_like(x)
TsatL,TsatV=300.,300.
p=Props('P','T',TsatL,'Q',0,'R134a')
for i in range(len(x)):
    h[i]=ShahCondensation_Average(x[i],x[i],'R134a',200,0.01,p,TsatL,TsatV)
    
havg=np.trapz(h,x=x)
pylab.figure(figsize=(7,5))
pylab.axhline(havg,ls='--')
pylab.text(0.2,havg,r'$\alpha_{avg}$',ha='center',va='bottom')
pylab.text(0.7,1500,'R134a\nG=200 kg/m$^2$\nD=0.01 m\np=%0.1f kPa' %(p) )
pylab.gca().set_xlabel('x [-]')
pylab.gca().set_ylabel(r'$\alpha$ [W/m$^2$/K]')
pylab.plot(x,h)
pylab.title('Shah Condensation HTC as a function of quality')
pylab.show()