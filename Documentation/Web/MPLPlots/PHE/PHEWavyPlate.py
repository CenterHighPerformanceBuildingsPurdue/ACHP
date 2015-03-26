from __future__ import division
import matplotlib,pylab, numpy as np
from matplotlib.patches import FancyArrowPatch

fig=pylab.figure(figsize=(6,2))
ax=fig.add_axes((0,0,1,1))

t=np.linspace(0,4*np.pi,100)
pylab.plot(t,np.cos(t))
    

pylab.gca().add_patch(FancyArrowPatch((4*np.pi,-1),(4*np.pi,1),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.text(4*np.pi,0,r' $2\hat a$',ha='left',va='center')
pylab.gca().add_patch(FancyArrowPatch((2*np.pi,1),(4*np.pi,1),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.text(3*np.pi,1.01,r' $\Lambda$',ha='center',va='bottom')

ax.axis('equal')
ax.axis('off')

pylab.show()