from __future__ import division
import matplotlib,pylab, numpy as np
from matplotlib.patches import FancyArrowPatch

wC=0.1
hC=1.0
h=1
Nc=6

fig=pylab.figure(figsize=(6.5,6))
ax=fig.add_axes((0,0,1,1))
for i in range(Nc):
    #The cold channel
    x0=wC*(2*i)
    y0=hC/2
    x=[x0+wC/2,x0+wC/2,x0-wC/2,x0-wC/2,x0+wC/2]
    y=[y0-hC/2,y0+hC/2,y0+hC/2,y0-hC/2,y0-hC/2]
    pylab.fill(x,y,'c')
    pylab.gca().add_patch(FancyArrowPatch((x0,y0-0.3*hC),(x0,y0+0.3*hC),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
    
    #The cold channel
    x0=wC*(2*i+1)
    y0=hC/2
    x=[x0+wC/2,x0+wC/2,x0-wC/2,x0-wC/2,x0+wC/2]
    y=[y0-hC/2,y0+hC/2,y0+hC/2,y0-hC/2,y0-hC/2]
    pylab.fill(x,y,'r')
    pylab.gca().add_patch(FancyArrowPatch((x0,y0+0.3*hC),(x0,y0-0.3*hC),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
    
ax.axis('equal')    
ax.axis('off')

pylab.show()