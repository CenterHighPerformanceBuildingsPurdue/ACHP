from __future__ import division, # Make integer 3/2 give 1.5 in python 2.x
                       print_function
from CoolProp.CoolProp import PropsSI
import numpy as np
import pylab as pylab
from math import pi
#from PyACHP.Correlations import AccelPressureDrop, LockhartMartinelli
from ACHP.FinCorrelations import WavyLouveredFins,HerringboneFins,FinInputs,PlainFins

#example evaporator
FinsTubes=FinInputs()

FinsTubes.Tubes.NTubes_per_bank=32
FinsTubes.Tubes.Ncircuits=5
FinsTubes.Tubes.Nbank=6
FinsTubes.Tubes.Ltube=0.452
FinsTubes.Tubes.OD=0.009525
FinsTubes.Tubes.ID=0.0089154
FinsTubes.Tubes.Pl=0.0254
FinsTubes.Tubes.Pt=0.0219964

FinsTubes.Fins.FPI=14.5
FinsTubes.Fins.Pd=0.001
FinsTubes.Fins.xf=0.001
FinsTubes.Fins.t=0.00011
FinsTubes.Fins.k_fin=237

FinsTubes.Air.Vdot_ha=0.5663
FinsTubes.Air.Tmean=299.8
FinsTubes.Air.Tdb=299.8
FinsTubes.Air.p=101325
FinsTubes.Air.RH=0.51
FinsTubes.Air.RHmean=0.51
FinsTubes.Air.FanPower=438

FinsTubes.Validate()
print(FinsTubes)  #just print our inputs

#generate data
V_dots_ha=np.linspace(0.1,1,1000)
#overall fin efficiency
eta_louvered=np.zeros_like(V_dots_ha)  #wavy-louvered fins
eta_wavy=np.zeros_like(V_dots_ha)  #herringbone fins
eta_plain=np.zeros_like(V_dots_ha)  #plain fins
#pressure drop
DP_louvered=np.zeros_like(V_dots_ha)  #wavy-louvered fins
DP_wavy=np.zeros_like(V_dots_ha)  #herringbone fins
DP_plain=np.zeros_like(V_dots_ha)  #plain fins
#Reynolds number
RE_D_louvered=np.zeros_like(V_dots_ha) #Reynolds number
RE_Dc_wavy=np.zeros_like(V_dots_ha) #Reynolds number
RE_Dc_plain=np.zeros_like(V_dots_ha) #Reynolds number

for i in range(len(V_dots_ha)):
    #determine inlet/outlet quality
    FinsTubes.Air.Vdot_ha = V_dots_ha[i]
    #wavy louvered fins
    WavyLouveredFins(FinsTubes)
    DP_louvered[i]=FinsTubes.dP_a #wavy-louvered fins
    eta_louvered[i]=FinsTubes.eta_a  #wavy-louvered fins
    RE_D_louvered[i]=FinsTubes.Re #hReynolds number
    #wavy herringbone fins
    HerringboneFins(FinsTubes) 
    eta_wavy[i]=FinsTubes.eta_a  #herringbone fins
    DP_wavy[i]=FinsTubes.dP_a  #herringbone fins
    RE_Dc_wavy[i]=FinsTubes.Re #hReynolds number
    #plain herringbone fins
    PlainFins(FinsTubes) 
    eta_plain[i]=FinsTubes.eta_a  #herringbone fins
    DP_plain[i]=FinsTubes.dP_a  #herringbone fins
    RE_Dc_plain[i]=FinsTubes.Re #hReynolds number

x_lbl="Air flowrate in m^3/s"

#One function (no inputs) per plot for use with Sphinx+matplotlib
def plot_dp():
    scl='linear';txtx=0.6;txty=2;pos='best';
    pylab.figure(figsize=(7,5))
    pylab.text(txtx,txty,"Data for some example HX")
    pylab.title('Pressure Drop as a function of air flowrate')
    pylab.gca().set_ylabel(r'Pressure drop in [Pa]')
    pylab.gca().set_xlabel(x_lbl)
    pylab.plot(V_dots_ha,DP_wavy,label="Herringbone fins")
    pylab.plot(V_dots_ha,DP_louvered,'--',label="Wavy-louvered fins")
    pylab.plot(V_dots_ha,DP_plain,':',label="Plain fins")
    pylab.gca().set_ylim(0,40)
    leg=pylab.legend(loc=pos, fancybox=True)
    leg.get_frame().set_alpha(0.5)
    pylab.gca().set_yscale(scl)
    pylab.show()

def plot_Re():
    scl='linear';txtx=0.6;txty=100;pos='best';
    pylab.figure(figsize=(7,5))
    pylab.text(txtx,txty,"Data for some example HX")
    pylab.title('Reynolds number as a function of air flowrate')
    pylab.gca().set_ylabel(r'Re in [-]')
    pylab.gca().set_xlabel(x_lbl)
    pylab.plot(V_dots_ha, RE_D_louvered,label="Wavy-louvered fins - base is D_collar")
    pylab.plot(V_dots_ha,RE_Dc_wavy,'--',label="Herringbone fins - base is D")
    pylab.plot(V_dots_ha,RE_Dc_plain,':',label="Plain fins - base is D_collar")
    pylab.gca().set_ylim(0,None)
    leg=pylab.legend(loc=pos, fancybox=True)
    leg.get_frame().set_alpha(0.5)
    pylab.gca().set_yscale(scl)
    pylab.show()

def plot_eta():
    scl='linear';txtx=0.6;txty=0.1;pos='best';
    pylab.figure(figsize=(7,5))
    pylab.text(txtx,txty,"Data for some example HX")
    pylab.title('Fin efficiency as a function of air flowrate')
    pylab.gca().set_ylabel(r'Fin efficiency in [-]')
    pylab.gca().set_xlabel(x_lbl)
    pylab.plot(V_dots_ha,eta_wavy,label="Herringbone fins")
    pylab.plot(V_dots_ha, eta_louvered,'--',label="Wavy-louvered fins")
    pylab.plot(V_dots_ha, eta_plain,':',label="Plain fins")
    pylab.gca().set_ylim(0,None)
    leg=pylab.legend(loc=pos, fancybox=True)
    leg.get_frame().set_alpha(0.5)
    pylab.gca().set_yscale(scl)
    pylab.show()
    
if __name__=='__main__':
    plot_eta()
    plot_dp()
    plot_Re()
