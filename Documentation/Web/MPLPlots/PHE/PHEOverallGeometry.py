import pylab,numpy as np
from math import cos,sin
from matplotlib.patches import FancyArrowPatch

fig=pylab.figure(figsize=(4,7))
ax=fig.add_axes((0,0,1,1))

# Slanted lines
N=7
w=0.75
phi=np.pi/3.0
for y0 in np.linspace(0,1,N):
    x=[-w/2.0,0,w/2.0]
    y=[-w/2.0*sin(phi)+y0,y0,-w/2.0*sin(phi)+y0]
    pylab.plot(x,y,'k')

# The frame
x=[w/2.0,w/2.0,-w/2.0,-w/2.0,w/2.0]
y=[-w/2.0*sin(phi),1,1,-w/2.0*sin(phi),-w/2.0*sin(phi)]
pylab.plot(x,y,'k--')

# The inclination angle
pylab.text(0.02,-0.08,r'$\varphi$',ha='left',va='top',size=20)
pylab.plot([0,0],[0,-w/2.0*sin(phi)],'k:')

# The ports
t=np.linspace(0,np.pi*2.0,300)
rP=0.125
portsxy=[(-w/2.0+rP,1+rP),(-w/2.0+rP,-w/2.0*sin(phi)-rP),(w/2.0-rP,1+rP),(w/2.0-rP,-w/2.0*sin(phi)-rP)]

for x0,y0 in portsxy:
    pylab.plot(x0+rP*np.cos(t),y0+rP*np.sin(t),'k')
    pylab.fill(x0+rP*np.cos(t),y0+rP*np.sin(t),'grey')

#rounded edges
rE=0.15
y0E=(1.0-w/2.0*sin(phi))/2
wE=1.2*w
hE=2.9*(1.0-w/2.0*sin(phi))

t=np.linspace(0,np.pi/2,100)
x1=wE/2.0-rE+rE*np.cos(t)
y1=hE/2.0-rE+rE*np.sin(t)+y0E
t=np.linspace(np.pi/2,np.pi,100)
x2=-wE/2.0+rE+rE*np.cos(t)
y2=hE/2.0-rE+rE*np.sin(t)+y0E
t=np.linspace(np.pi,3.0*np.pi/2.0,100)
x3=-wE/2.0+rE+rE*np.cos(t)
y3=-hE/2.0+rE+rE*np.sin(t)+y0E
t=np.linspace(3.0*np.pi/2.0,2*np.pi,100)
x4=wE/2.0-rE+rE*np.cos(t)
y4=-hE/2.0+rE+rE*np.sin(t)+y0E
pylab.plot(np.r_[x1,x2,x3,x4,x1[0]],np.r_[y1,y2,y3,y4,y1[0]],'k')

#Length labels
pylab.gca().add_patch(FancyArrowPatch((wE/2.0+0.08,-w/2.0*sin(phi)),(wE/2.0+0.08,1),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.text(wE/2.0+0.1,y0E,'$L$',ha='left',va='center',size=20)
pylab.gca().add_patch(FancyArrowPatch((-w/2.0,y0E+hE/2.0+0.06),(w/2.0,y0E+hE/2.0+0.06),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.text(0,y0E+hE/2.0+0.08,'$B$',ha='center',va='bottom',size=20)

pylab.gca().add_patch(FancyArrowPatch((-wE/2.0-0.08,portsxy[2][1]),(-wE/2.0-0.08,portsxy[3][1]),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.text(-wE/2.0-0.1,y0E,'$L_p$',ha='right',va='center',size=20)
pylab.gca().add_patch(FancyArrowPatch((portsxy[3][0],y0E-hE/2.0-0.06),(portsxy[0][0],y0E-hE/2.0-0.06),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.text(0,y0E-hE/2.0-0.08,'$B_p$',ha='center',va='top',size=20)

pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.show()