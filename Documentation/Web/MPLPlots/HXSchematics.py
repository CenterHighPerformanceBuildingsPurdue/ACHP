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
nR=4
nC=6
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
pylab.gca().axis('equal')
pylab.gca().axis('off')


## Connect up each circuit(assumed to be number of rows)
for i in range(nR):
    for j in range(nC-1):
        if j%2==0:
            x1=j
            y1=i
        else:
            x1=j
            y1=i+offset
            
        if (j+1)%2==0:
            x2=j+1
            y2=i
        else:
            x2=j+1
            y2=i+offset
            
        pylab.plot(np.r_[x1,x2],np.r_[y1,y2],'k')

## Inlet and outlet arrows for air stream
for i in range(2*nR):
    pylab.gca().add_patch(FancyArrowPatch((-1,i*0.5),(-0.5,i*0.5),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
    pylab.gca().add_patch(FancyArrowPatch((nC-0.5,i*0.5),(nC,i*0.5),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
    pylab.gca().text(-1.25,nR/2,'Air Inlet',rotation=90)
    pylab.gca().text(nC+0.25,nR/2,'Air Outlet',rotation=90)
    
## "Header" for inlet and outlet of refrigerant
pylab.plot(np.r_[nC-1.4,nC-0.6,nC-0.6,nC-1.4,nC-1.4],np.r_[-0.25,-0.25,nR,nR,-0.25],'k--')
pylab.plot(np.r_[1-1.4,1-0.6,1-0.6,1-1.4,1-1.4],np.r_[-0.25,-0.25,nR,nR,-0.25],'k--')
pylab.gca().text(nC-1,nR+0.1,'Fluid Inlet\nBank',ha='center',va='bottom')
pylab.gca().text(0,nR+0.1,'Fluid Outlet\nBank',ha='center',va='bottom')
pylab.show()