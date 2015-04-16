import pylab,numpy as np
from matplotlib.patches import FancyArrowPatch

fig=pylab.figure(figsize=(8,4))

pylab.fill(np.r_[0,5,5,0,0],np.r_[-1,-1,0,0,-1],'b',alpha=0.1)
pylab.plot(np.r_[0,5],np.r_[1,1],'k')
pylab.plot(np.r_[0,5],np.r_[0,0],'k')
pylab.plot(np.r_[0,5],np.r_[-1,-1],'k')

pylab.plot(np.r_[3,3],np.r_[-1,1],'k--')

pylab.plot(3,0.5,'ko')
pylab.plot(3,-0.5,'ko')
pylab.text(3.05,0.5,'$T_{a,x}$',ha='left',va='center')
pylab.text(3.05,-0.5,'$T_{g,x}$',ha='left',va='center')

pylab.plot(0,0.5,'ko')
pylab.plot(0,-0.5,'ko')
pylab.text(0.05,0.5,'$T_{a,i}$',ha='left',va='center')
pylab.text(0.05,-0.5,'$T_{g,o}$',ha='left',va='center')

pylab.plot(5,0.5,'ko')
pylab.plot(5,-0.5,'ko')
pylab.text(5.05,0.5,'$T_{a,o}$',ha='left',va='center')
pylab.text(5.05,-0.5,'$T_{g,i}$',ha='left',va='center')

pylab.plot(5,0.5,'ko')
pylab.plot(5,-0.5,'ko')
pylab.text(5.05,0.5,'$T_{a,o}$',ha='left',va='center')
pylab.text(5.05,-0.5,'$T_{g,i}$',ha='left',va='center')
pylab.gca().add_patch(FancyArrowPatch((1,0.5),(0.5,0.5),arrowstyle='<|-',fc='k',ec='k',mutation_scale=20,lw=0.8))

#Wet line
pylab.plot(np.r_[3,5],np.r_[0.01,0.01],'b',lw=4)

pylab.plot(0,0,'ro')
pylab.text(-0.05,0,'$T_{i,s}$',ha='right',va='center')

pylab.plot(5,0,'ro')
pylab.text(5.05,0,'$T_{o,s}$',ha='left',va='center')

pylab.text(4,-1.1,'Wet Part',ha='center',va='top')
pylab.text(1.5,-1.1,'Dry Part',ha='center',va='top')
pylab.gca().add_patch(FancyArrowPatch((4,-0.5),(4.5,-0.5),arrowstyle='<|-',fc='k',ec='k',mutation_scale=20,lw=0.8))

pylab.gca().set_xlim(-0.1,5.1)
pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.show()