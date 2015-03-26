# A file to generate a cycle schematic for the preconditioner

#Imports
from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import numpy as np
from matplotlib.patches import FancyArrowPatch

def rect (x0,y0,w,h):
    return (np.r_[x0-w/2,x0+w/2,x0+w/2,x0-w/2,x0-w/2],np.r_[y0-h/2,y0-h/2,y0+h/2,y0+h/2,y0-h/2])

def TwoTriangles(H,W,x,y):
    X=np.r_[x-W/2.0,x+W/2.0,x+W/2.0,x-W/2.0,x-W/2.0]
    Y=np.r_[y-H/2.0,y+H/2.0,y-H/2.0,y+H/2.0,y-H/2.0]
    return (X,Y)
    
def circle(x0,y0,r,N):
    t=np.linspace(0,2*np.pi,N)
    x=x0+r*np.cos(t)
    y=y0+r*np.sin(t)
    return (x,y)
    
fig=plt.figure(figsize=(8,4))
ax=fig.add_axes((0,0,1,1))

#################################################
############# The components ####################
#################################################

# Condenser
x,y=rect(3,0,1,2); fig.gca().fill(x,y,'red') 
fig.gca().text(3,0,'Condenser',ha='center',va='center')

# Compressor
x,y=circle(1.5,-1.5,0.25,100); fig.gca().plot(x,y,'k') 
fig.gca().text(1.5,-1.9,'Compressor',ha='center',va='center')

#Internal Heat Exchanger
x,y=rect(0,0,1,2); fig.gca().fill(x,y,'cyan')
fig.gca().text(0,0,'Internal\nHeat\nExchanger',ha='center',va='center')

# Expansion Device
x,y=TwoTriangles(0.5,0.5,1.5,1.5); fig.gca().plot(x,y,'k')
fig.gca().text(1.5,1.1,'Expansion Device',ha='center',va='center')

# Cooling Coil
x,y=rect(-3,0,1,2); fig.gca().fill(x,y,'lightblue')
fig.gca().text(-3,0,'Cooling\nCoil',ha='center',va='center')

# Pump
x,y=circle(-0.85,-1.5,0.25,100); fig.gca().plot(x,y,'k') 
fig.gca().text(-0.85,-1.9,'Pump',ha='center',va='center')

# Line Sets
fig.gca().plot( np.r_[-2.8,-1.7],np.r_[1.5,1.5],'k',lw=6)
fig.gca().text(-2.25,1.7,'Supply Line',ha='center',va='center')
fig.gca().plot( np.r_[-2.8,-1.7],np.r_[-1.5,-1.5],'k',lw=6)
fig.gca().text(-2.25,-1.7,'Return Line',ha='center',va='center')

#################################################
##### Connect the components with lines #########
#################################################

# Condenser to XV
fig.gca().plot(np.r_[3,3],np.r_[1.0,1.5],'k')
fig.gca().add_patch(FancyArrowPatch((3,1.5),(1.75,1.5),arrowstyle='-|>',mutation_scale=20,fc='k',ec='k',lw=0.5))

# XV to IHX
fig.gca().plot(np.r_[1.25,0.33],np.r_[1.5,1.5],'k')
fig.gca().add_patch(FancyArrowPatch((0.33,1.5),(0.33,1),arrowstyle='-|>',mutation_scale=20,fc='k',ec='k',lw=0.5))

# IHX to Compressor
fig.gca().plot(np.r_[0.33,0.33],np.r_[-1,-1.5],'k')
fig.gca().add_patch(FancyArrowPatch((0.33,-1.5),(1.25,-1.5),arrowstyle='-|>',mutation_scale=20,fc='k',ec='k',lw=0.5))

# Compressor to Condenser
fig.gca().plot(np.r_[1.75,3],np.r_[-1.5,-1.5],'k')
fig.gca().add_patch(FancyArrowPatch((3,-1.5),(3,-1),arrowstyle='-|>',mutation_scale=20,fc='k',ec='k',lw=0.5))

# Cooling Coil to Pump
fig.gca().plot(np.r_[-3,-3],np.r_[-1,-1.5],'k')
fig.gca().add_patch(FancyArrowPatch((-3,-1.5),(-1.1,-1.5),arrowstyle='-|>',mutation_scale=20,fc='k',ec='k',lw=0.5))

# Pump to IHX
fig.gca().plot(np.r_[-0.6,-0.33],np.r_[-1.5,-1.5],'k')
fig.gca().add_patch(FancyArrowPatch((-0.33,-1.5),(-0.33,-1),arrowstyle='-|>',mutation_scale=20,fc='k',ec='k',lw=0.5))

# Pump to Cooling Coil
fig.gca().plot(np.r_[-0.33,-0.33,-3,-3],np.r_[1,1.5,1.5,1],'k')

#################################################
################## Cycle Inputs #################
#################################################

fig.gca().text(-0.65,0,'$T_{evap}$',ha='right',bbox=dict(boxstyle='rarrow',fc='w',ec='k'))
fig.gca().text(2.35,0,'$T_{cond}$',ha='right',bbox=dict(boxstyle='rarrow',fc='w',ec='k'))
fig.gca().text(-3.2,1.2,'$T_{g,i,cc}$',ha='right',bbox=dict(boxstyle='rarrow',fc='w',ec='k'))

fig.gca().text(0.33,-1.5,'1',ha='right',va='center',bbox=dict(boxstyle='round',fc='w',ec='k'))
fig.gca().text(3,-1.5,'2',ha='right',va='center',bbox=dict(boxstyle='round',fc='w',ec='k'))
fig.gca().text(3,1.5,'3',ha='right',va='center',bbox=dict(boxstyle='round',fc='w',ec='k'))
fig.gca().text(0.33,1.5,'4',ha='right',va='center',bbox=dict(boxstyle='round',fc='w',ec='k'))

fig.gca().text(-3,1.5,'5',ha='right',va='center',bbox=dict(boxstyle='round',fc='w',ec='k'))
fig.gca().text(-3,-1.5,'6',ha='right',va='center',bbox=dict(boxstyle='round',fc='w',ec='k'))
fig.gca().text(-1.4,-1.5,'7',ha='right',va='center',bbox=dict(boxstyle='round',fc='w',ec='k'))
fig.gca().text(-0.33,-1.5,'8',ha='right',va='center',bbox=dict(boxstyle='round',fc='w',ec='k'))
fig.gca().text(-0.33,1.5,'9',ha='right',va='center',bbox=dict(boxstyle='round',fc='w',ec='k'))

ax.axis('equal')
ax.axis('off')
plt.show()