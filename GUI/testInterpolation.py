import numpy as np, pylab
from scipy.interpolate import interp1d
from CoolProp.CoolProp import Props

N=100
T=np.linspace(180,350,N)
p=np.zeros_like(T)

for i in range(len(T)):
    p[i]=Props('S','T',T[i],'Q',1,'R290')

quadddd=interp1d(T,p,kind='quadratic')

Tint=(T[0:N-1]+T[1:N])/2
pint=np.zeros_like(Tint)
pval=np.zeros_like(Tint)
for i in range(len(Tint)):
    pint[i]=quadddd(Tint[i])
    pval[i]=Props('S','T',Tint[i],'Q',1,'R290')
    
pylab.plot(Tint,(pint/pval-1)*100)
pylab.show()

pylab.plot(T,p,Tint,pint)
pylab.show()