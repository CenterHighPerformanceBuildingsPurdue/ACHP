import pylab, numpy as np

r1=1
r2=0.6
t=np.linspace(0,2*np.pi)
#Tube insulation
x=r1*np.cos(t)
y=r1*np.sin(t)
pylab.fill(x,y,'k')

#Tube
x=r2*np.cos(t)
y=r2*np.sin(t)
pylab.fill(x,y,'w',color='r',fc='white',lw=8)

pylab.gca().text(0,(r1+r2)/2,'Insulation',color='w',ha='center',va='center',size=20)
pylab.gca().text(0,0.8*(r2),'Tube',color='k',ha='center',va='center',size=20)

pylab.gca().text(0.9*r1,0.9*(r1),r'$T_{\infty}$',color='k',ha='center',va='center',size=20)
pylab.plot(0,0,'k.',ms=10)
pylab.gca().text(-0.1,0,r'$T$',color='k',ha='center',va='center',size=20)
pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.show()