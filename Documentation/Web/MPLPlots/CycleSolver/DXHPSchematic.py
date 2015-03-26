# A file to generate a cycle schematic for the DX-mode AC system

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
    
fig=plt.figure()
ax=fig.add_axes((0,0,1,1))

#################################################
############# The components ####################
#################################################

# Condenser
x,y=rect(3,0,1,2); fig.gca().fill(x,y,'red') 
fig.gca().text(3,0,'Condenser\n(Indoors)',ha='center',va='center')

# Compressor
x,y=circle(1.5,-1.5,0.165,100); fig.gca().plot(x,y,'k') 
fig.gca().text(1.5,-1.8,'Compressor',ha='center',va='center')

# Evaporator
x,y=rect(0,0,1,2); fig.gca().fill(x,y,'lightblue') 
fig.gca().text(0,0,'Evaporator\n(Outdoors)',ha='center',va='center')

# Expansion Device
x,y=TwoTriangles(0.25,0.25,1.5,1.5); fig.gca().plot(x,y,'k') 
fig.gca().text(1.5,1.3,'Expansion Device',ha='center',va='center')

#Liquid Line
fig.gca().plot([2.15,2.85],[1.5,1.5],'k',lw=6)
fig.gca().text(2.5,1.6,'Liquid Line',ha='center',va='center')

#Vapor return Line
fig.gca().plot([2.15,2.85],[-1.5,-1.5],'k',lw=6)
fig.gca().text(2.5,-1.4,'Vapor Line',ha='center',va='center')

#################################################
##### Connect the components with lines #########
#################################################

# Condenser to Expansion Device
fig.gca().plot(np.r_[3,3],np.r_[1.0,1.5],'k')
fig.gca().add_patch(FancyArrowPatch((3,1.5),(1.625,1.5),arrowstyle='-|>',mutation_scale=20,fc='k',ec='k',lw=0.5))

# XV to Evaporator 
fig.gca().plot(np.r_[1.375,0],np.r_[1.5,1.5],'k')
fig.gca().add_patch(FancyArrowPatch((0,1.5),(0,1),arrowstyle='-|>',mutation_scale=20,fc='k',ec='k',lw=0.5))

# Condenser to XV
fig.gca().plot(np.r_[0,0],np.r_[-1,-1.5],'k')
fig.gca().add_patch(FancyArrowPatch((0,-1.5),(1.375,-1.5),arrowstyle='-|>',mutation_scale=20,fc='k',ec='k',lw=0.5))

# Compressor to Condenser
fig.gca().plot(np.r_[1.675,3],np.r_[-1.5,-1.5],'k')
fig.gca().add_patch(FancyArrowPatch((3,-1.5),(3,-1),arrowstyle='-|>',mutation_scale=20,fc='k',ec='k',lw=0.5))

#################################################
################## Cycle Inputs #################
#################################################

fig.gca().text(-0.65,0,'$T_{evap}$',ha='right',bbox=dict(boxstyle='rarrow',fc='w',ec='k'))
fig.gca().text(2.35,0,'$T_{cond}$',ha='right',bbox=dict(boxstyle='rarrow',fc='w',ec='k'))

ms=10
fig.gca().text(0,-1.5,'1',ha='center',va='center',bbox=dict(fc='w',boxstyle='round'))
fig.gca().text(2,-1.5,'2',ha='center',va='center',bbox=dict(fc='w',boxstyle='round'))
fig.gca().text(3,-1.5,'3',ha='center',va='center',bbox=dict(fc='w',boxstyle='round'))
fig.gca().text(3,1.5,'4',ha='center',va='center',bbox=dict(fc='w',boxstyle='round'))
fig.gca().text(2,1.5,'5',ha='center',va='center',bbox=dict(fc='w',boxstyle='round'))
fig.gca().text(0,1.5,'6',ha='center',va='center',bbox=dict(fc='w',boxstyle='round'))

plt.draw()
plt.gca().axis('equal')
plt.gca().axis('off')
plt.show()