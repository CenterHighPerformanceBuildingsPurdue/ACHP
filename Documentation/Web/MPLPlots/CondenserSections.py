import pylab, numpy as np

fig=pylab.figure()
pylab.plot(np.r_[0,4,4,0,0],np.r_[0,0,1,1,0],'k')

pylab.plot(np.r_[1,1],np.r_[0,1],'k-.')
pylab.plot(np.r_[3,3],np.r_[0,1],'k-.')

pylab.gca().text(2,0.5,'Two-Phase',ha='center',va='center')
pylab.gca().text(3.5,0.5,'Superheated',ha='center',va='center')
pylab.gca().text(0.5,0.5,'Subcooled',ha='center',va='center')
pylab.gca().text(2,1.1,'Superheated, Subcooled and Two-Phase Sections',ha='center',va='bottom')

pylab.gca().text(0.5,-0.1,'$w_{subcool}$',ha='center',va='center')
pylab.gca().text(2,-0.1,'$w_{two-phase}$',ha='center',va='center')
pylab.gca().text(3.5,-0.1,'$w_{superheat}$',ha='center',va='center')

y0=1.7 #shift of second axis
pylab.plot(np.r_[0,4,4,0,0],np.r_[-y0,-y0,-y0+1,-y0+1,-y0],'k')
pylab.plot(np.r_[3,3],np.r_[-y0,-y0+1],'k-.')
pylab.gca().text(1.5,-y0+0.5,'Two-Phase',ha='center',va='center')
pylab.gca().text(3.5,-y0+0.5,'Superheated',ha='center',va='center')
pylab.gca().text(2,-y0+1.1,'Superheated and Two-Phase Sections',ha='center',va='bottom')

pylab.gca().text(3.5,-y0-0.1,'$w_{superheat}$',ha='center',va='center')
pylab.gca().text(1.5,-y0-0.1,'$w_{two-phase}$',ha='center',va='center')

pylab.xlim(-0.1,4.3)
pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.show()