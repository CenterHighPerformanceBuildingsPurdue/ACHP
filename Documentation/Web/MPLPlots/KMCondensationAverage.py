
import CoolProp as CP
import pylab,numpy as np

from ACHP.Correlations import KM_Cond_Average

x=np.linspace(0,1,1000)
h=np.zeros_like(x)
TsatL,TsatV=300.,300.
AS= CP.AbstractState('HEOS','R134a')
AS.update(CP.QT_INPUTS,0.0,TsatL)
p = AS.p() #[Pa]
beta = 1 #channel aspect ratio (=width/height)
for i in range(len(x)):
    h[i]=KM_Cond_Average(x[i],x[i],AS,200,0.01,TsatL,TsatV,p,beta)[1]
    
havg=np.trapz(h,x=x)
pylab.figure(figsize=(7,5))
pylab.axhline(havg,ls='--')
pylab.text(0.2,havg,r'$\alpha_{avg}$',ha='center',va='bottom')
pylab.text(0.7,1500,'R134a\nG=200 kg/m$^2$\nD=0.01 m\np=%0.1f kPa' %(p/1000) )
pylab.gca().set_xlabel('x [-]')
pylab.gca().set_ylabel(r'$\alpha$ [W/m$^2$/K]')
pylab.plot(x,h)
pylab.title('KM Condensation HTC as a function of quality')
pylab.show()
