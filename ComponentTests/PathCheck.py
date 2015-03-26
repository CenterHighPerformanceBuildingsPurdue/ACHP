'''
Created on Mar 20, 2015

@author: AmmarBahman
'''

''' This code is just for demonstration and has nothing to do with ACHP'''

from __future__ import division

import sys
print sys.path
print ' '

import CoolProp; print CoolProp.__file__
print ' '
import CoolProp; print CoolProp.__version__, CoolProp.__gitrevision__
print ' '

from CoolProp.CoolProp import cair_sat
print cair_sat(300) #Input T in K #kJ/kg-K
print' '

import CoolProp
import CoolProp.CoolProp as CP
print 'Humidity ratio of 50% rel. hum. air at 300 K, 101.325 kPa:', CP.HAProps('W', 'T', 300, 'P', 101.325, 'R', 0.5), 'kg_w/kg_da'
print ' '

print 'Props', CP.Props('C','T',288.209928591,'P',101.325,'R410a')
print 'PropsSI', CP.PropsSI('C','T',288.209928591,'P',101325,'R410a')
print ' '
import CoolProp.HumidAirProp as CPH
print 'HAProps',CPH.HAProps('V','T',300,'P',101.325,'R',0.5)
print 'HAPropsSI',CPH.HAPropsSI('V','T',300,'P',101325,'R',0.5)
print ' '


from CoolProp.CoolProp import FluidsList
print FluidsList()
#print FluidsList().index('Air')
#print FluidsList().count('Air') ==1
print' '

x ='Argon'
if FluidsList().count(x)==1:
    print FluidsList().index(x)
else: 
    print 'NO'



from CoolProp.CoolProp import get_REFPROPname, PropsSI, Props
from CoolProp.CoolProp import State

Fluid = 'INCOMP::TD12'
x = PropsSI('D', 'T', 300, 'P', 300, Fluid)
print 'density = ', x
print' '

if 'INCOMP' in Fluid:
    print 'YES'
else:
    print 'NO'


print(' ')
print('************ USING REFPROP ***************') #This backend provides a clean interface between CoolProp and REFPROP
print(' ')

import CoolProp
print 'Temp of Water at 101.325 kPa is: ',CoolProp.CoolProp.PropsSI('T','P',101325,'Q',0,'REFPROP::Water'), 'K'
print ' '

import sphinx; print sphinx.__version__
print ' '

print PropsSI("T","Q",0,"P",2301391.979,"R410A")
print PropsSI("T","Q",1,"P",2301391.979,"R410A")

#print PropsSI("D","T",310.996631,"P",2301391.979,"REFPROP::R410A")
print PropsSI("T","Q",0,"P",2301391.979,"REFPROP::R410A")
print PropsSI("T","Q",1,"P",2301391.979,"REFPROP::R410A")
