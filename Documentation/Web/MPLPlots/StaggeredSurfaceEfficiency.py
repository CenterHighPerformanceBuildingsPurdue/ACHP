import matplotlib,pylab, numpy as np
from matplotlib.patches import FancyArrowPatch

from numpy import arcsin,cos,sin,sqrt
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
pylab.gca().axis('equal')
pylab.gca().axis('off')

#Each cell gets its boundaries
pl=1.0
pt=1.0
X_D=sqrt(pl**2+pt**2/4)/2.0
X_T = pt / 2.0
th= arcsin(X_T/(2*X_D)) #Angle to midpoint of cell

#pylab.plot(X_D*cos(th),X_D*sin(th),'o')

for i in range(nR):
    for j in range(nC):
        if j%2==0:
            x=j
            y=i
        else:
            x=j
            y=i+offset
            
        #Plot the tube
        xv,yv=circle(x,y,r)
        pylab.plot(xv,yv,'b')
            
        #Plot the boundaries of the unit cell
        x0=x
        y0=y
        L=X_D/cos(th)
        x1,y1=(L+x0,0+y0)
        x2,y2=(L-2*(L-X_D*cos(th))+x0,0+2*X_D*sin(th)+y0)
        x3,y3=(-(L-2*(L-X_D*cos(th)))+x0,0+2*X_D*sin(th)+y0)
        x4,y4=(-L+x0,0+y0)
        x5,y5=(-(L-2*(L-X_D*cos(th)))+x0,-(0+2*X_D*sin(th))+y0)
        x6,y6=(L-2*(L-X_D*cos(th))+x0,-(0+2*X_D*sin(th))+y0)
        x7,y7=(L+x0,0+y0)
        pylab.plot(np.r_[x1,x2,x3,x4,x5,x6,x7],np.r_[y1,y2,y3,y4,y5,y6,y7],'k:')
        
pylab.gca().add_patch(FancyArrowPatch((0,0),(1,0.5),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.gca().add_patch(FancyArrowPatch((0,0),(0,1),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.gca().text(0,0.6,'$2X_T=P_t$',ha='right',va='bottom')
pylab.gca().text(X_D*cos(th)+0.15,X_D*sin(th)+0.05,'$X_D$',ha='right',va='bottom')

pylab.show()