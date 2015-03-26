import pylab,numpy as np
from numpy import sin
from matplotlib.patches import FancyArrowPatch

fig=pylab.figure()
w=1
h=1
th=3.14159/25.
x=np.r_[0,0,w,w,0]
y=np.r_[0,h,h-w*sin(th),0-w*sin(th),0]
pylab.plot(x,y)

x=np.r_[0,0,w/2.0,w/2.0,0]
y=np.r_[0,h/6.0,h/6.0-w/2.0*sin(th),0-w/2.0*sin(th),0]
pylab.plot(x,y,'--')
pylab.text(w/4.0,h/12.0-w/4.0*sin(th)-h/30.,'$A_{a,subcool}$',ha='center',va='center')

h0=h-w/2.0*sin(th)-h/6.0
x=np.r_[w/2.0,w/2.0,w,w,w/2.0]
y=np.r_[0+h0,h/6.0+h0,h/6.0-w/2.0*sin(th)+h0,0-w/2.0*sin(th)+h0,0+h0]
pylab.plot(x,y,'--')
pylab.text(0.75*w,h-h/12.0-0.75*w*sin(th)-h/30.,'$A_{a,superheat}$',ha='center',va='center')

pylab.text(0.5*w,h/2.0-0.5*w*sin(th),'$A_{a,two-phase}$',ha='center',va='center')

##Add the circuits
for y0 in [h/12.,h/12.+h/6.,h/12.+2*h/6.,h/12.+3*h/6.,h/12.+4*h/6.,h/12.+5*h/6.]:
    pylab.plot(np.r_[0,w],np.r_[y0,y0-w*sin(th)],'k',lw=4)
    
pylab.gca().add_patch(FancyArrowPatch((w+w/10.,h-h/12.0-(w+w/10.)*sin(th)),(w,h-h/12.0-w*sin(th)),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.gca().add_patch(FancyArrowPatch((0,h/12.0),(-w/10.,h/12.0-(-w/10.)*sin(th)),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
    
pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.show()