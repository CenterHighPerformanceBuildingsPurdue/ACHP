import matplotlib,pylab, numpy as np
from matplotlib.patches import FancyArrowPatch

def circle(x,y,r):
    """ This function makes lines for a circle given center and radius"""
    t=np.linspace(0,2*np.pi)
    xv=x+r*np.cos(t)
    yv=y+r*np.sin(t)
    return xv,yv

fig=pylab.figure()

## Actually make the set of tubes
nR=2
nC=3
offset = 0.5
r=0.25
for i in range(nR):
    for j in range(nC):
        if j%2==0:
            x=j
            y=i
        else:
            x=j
            y=i+offset
        xv,yv=circle(x,y,r)
        pylab.plot(xv,yv,'b')


##Dimension lines
pylab.plot(np.r_[0,0],np.r_[1,2],'k')
pylab.plot(np.r_[1,1],np.r_[1+offset,2],'k')
pylab.text(0.5,2,'Longitudinal\nSpacing\n$p_l$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((0,1.8),(1,1.8),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

pylab.plot(np.r_[-0.5,0],np.r_[1,1],'k')
pylab.plot(np.r_[-0.5,0],np.r_[0,0],'k')
pylab.text(-0.6,0.5,'Transverse\nSpacing\n$p_t$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((-0.3,0),(-0.3,1),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

#Airflow arrow
pylab.gca().add_patch(FancyArrowPatch((1.75,0.5),(3,0.5),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.text(2.25,0.5,'Airflow\nDirection',ha='center',va='center')

pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.show()