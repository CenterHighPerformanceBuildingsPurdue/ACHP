from __future__ import division, print_function, absolute_import
import sys,os

from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
from math import tan,atan,pi
from CoolProp.CoolProp import PropsSI

class PlotsClass():
    
    def __init__(self, *args, **kwds):
        pass    
        
    def TSOverlay(self,Cycle,**kwargs):
  
        Ref=Cycle.Ref
        
        Tmin=PropsSI('Tmin',Ref)
        Tmax=PropsSI('Tmax',Ref)   
        pmin=PropsSI('pmin',Ref)
        pmax=PropsSI('pmax',Ref)        

        fig = plt.figure(figsize=(6,4))
        self.axes=fig.add_subplot(111)
    
        Tsat = np.linspace(Tmin,PropsSI(Ref,"Tcrit")-0.0000000001,1000)
        (ssatL,ssatV)=(0.0*Tsat,0.0*Tsat)
        hcrit=np.linspace(0,0,2)
        pcrit=np.linspace(0,0,2)
        for i in np.arange(len(Tsat)):
            ssatL[i] = PropsSI('S','T',Tsat[i],'Q',0,Ref)/1000
            ssatV[i] = PropsSI('S','T',Tsat[i],'Q',1,Ref)/1000

        # Create critical iso-therm
        p_critical = Cycle.AS.p_critical()
        T_critical = Cycle.AS.T_critical()
        
        T_lim = np.linspace(Tmin+10,Tmax, 10000)
        s_lim = 0.0*T_lim
        
        for i in np.arange(len(T_lim)):
            s_lim[i] = PropsSI('S','T',T_lim[i],'P',p_critical*0.999,Ref)/1000

        # Check if the cycle is transcritical
        if Cycle.Compressor.Tout_r > T_critical  and Cycle.Compressor.pout_r > p_critical:
            # Average saturated temperatures and pressures
            Tsat_evap = Cycle.Tdew_evap 
            psat_evap = PropsSI('P','T',Tsat_evap,'Q',1,Ref)       

            # Saturated conditions at the evaporator
            sdew_evap = PropsSI('S','P',psat_evap,'Q',1,Ref)
            Tdew_evap = PropsSI('T','P',psat_evap,'Q',1,Ref)
            sbubble_evap = PropsSI('S','P',psat_evap,'Q',0,Ref)
            Tbubble_evap = PropsSI('T','P',psat_evap,'Q',0,Ref)

            Tevap_range1 = np.linspace(Tmin,Tbubble_evap,1000)
            Tevap_range2 = np.linspace(Tdew_evap,Tmax,1000)
            sevap_1 = 0.0*Tevap_range1
            sevap_2 = 0.0*Tevap_range2            

            for i in range(len(Tevap_range1)):
                
                try:
                    sevap_1[i] = PropsSI('S','P',psat_evap,'T',Tevap_range1[i],Ref)/1000
                except:
                    sevap_1[i] = PropsSI('S','P',psat_evap,'Q',0,Ref)/1000
                
                try:
                    sevap_2[i] = PropsSI('S','P',psat_evap,'T',Tevap_range2[i],Ref)/1000
                except:
                    sevap_2[i] = PropsSI('S','P',psat_evap,'Q',1,Ref)/1000
        
        else:
            # Average saturated temperatures and pressures
            Tsat_evap = Cycle.Tdew_evap
            Tsat_cond = Cycle.Tdew_cond
            psat_evap = PropsSI('P','T',Tsat_evap,'Q',1,Ref)
            psat_cond = PropsSI('P','T',Tsat_cond,'Q',1,Ref)
    
            # Saturated conditions at the evaporator
            sdew_evap = PropsSI('S','P',psat_evap,'Q',1,Ref)
            Tdew_evap = PropsSI('T','P',psat_evap,'Q',1,Ref)
            sbubble_evap = PropsSI('S','P',psat_evap,'Q',0,Ref)
            Tbubble_evap = PropsSI('T','P',psat_evap,'Q',0,Ref)
            
            # Saturated conditions at the condenser
            sdew_cond = PropsSI('S','P',psat_cond,'Q',1,Ref)
            Tdew_cond = PropsSI('T','P',psat_cond,'Q',1,Ref)
            sbubble_cond = PropsSI('S','P',psat_cond,'Q',0,Ref)
            Tbubble_cond = PropsSI('T','P',psat_cond,'Q',0,Ref)
            
            Tcond_range1 = np.linspace(Tmin,Tbubble_cond,1000)
            Tcond_range2 = np.linspace(Tdew_cond,Tmax,1000)
            scond_1 = 0.0*Tcond_range1
            scond_2 = 0.0*Tcond_range2
            Tevap_range1 = np.linspace(Tmin,Tbubble_evap,1000)
            Tevap_range2 = np.linspace(Tdew_evap,Tmax,1000)
            sevap_1 = 0.0*Tevap_range1
            sevap_2 = 0.0*Tevap_range2
            
            
            for i in range(len(Tcond_range1)):
                
                try:
                    scond_1[i] = PropsSI('S','P',psat_cond,'T',Tcond_range1[i],Ref)/1000
                except:
                    scond_1[i] = PropsSI('S','P',psat_cond,'Q',0,Ref)/1000
                    
                try:
                    scond_2[i] = PropsSI('S','P',psat_cond,'T',Tcond_range2[i],Ref)/1000
                except:
                    scond_2[i] = PropsSI('S','P',psat_cond,'Q',1,Ref)/1000
                
                try:
                    sevap_1[i] = PropsSI('S','P',psat_evap,'T',Tevap_range1[i],Ref)/1000
                except:
                    sevap_1[i] = PropsSI('S','P',psat_evap,'Q',0,Ref)/1000
                
                try:
                    sevap_2[i] = PropsSI('S','P',psat_evap,'T',Tevap_range2[i],Ref)/1000
                except:
                    sevap_2[i] = PropsSI('S','P',psat_evap,'Q',1,Ref)/1000

            # Isolines condensing temperature
            self.axes.plot(scond_1,Tcond_range1,'k-.',lw=1, alpha = 0.3)
            self.axes.plot(scond_2,Tcond_range2,'k-.',lw=1, alpha = 0.3)
            self.axes.plot([sdew_cond/1000,sbubble_cond/1000],[Tdew_cond,Tbubble_cond],'k-.',lw=1, alpha = 0.3)            

        # Critical isothermal
        self.axes.plot(s_lim,T_lim,'k-.',lw=1, alpha = 0.3)
        
        # Ts saturated lines
        self.axes.plot(ssatL,Tsat,'k',lw=1.5,alpha = 0.5)
        self.axes.plot(ssatV,Tsat,'k',lw=1.5, alpha=0.5)
        
        # Isolines evaporating temperature
        self.axes.plot(sevap_1,Tevap_range1,'k-.',lw=1, alpha = 0.3)
        self.axes.plot(sevap_2,Tevap_range2,'k-.',lw=1, alpha = 0.3)
        self.axes.plot([sdew_evap/1000,sbubble_evap/1000],[Tdew_evap,Tbubble_evap],'k-.',lw=1, alpha = 0.3)
    
        # Axes labels
        self.axes.set_xlabel('Entropy [kJ/(kg-K)]')
        self.axes.set_ylabel('Temperature [K]')

        # Check if the cycle is transcritical
        if Cycle.Compressor.Tout_r > T_critical  and Cycle.Compressor.pout_r > p_critical: 
            # Axes limits
            self.axes.set_xlim(0.9*Cycle.GasCooler.sout_r/1000,1.1*Cycle.Compressor.sout_r/1000)        
            self.axes.set_ylim(0.95*Cycle.Evaporator.Tout_r,1.1*Cycle.Compressor.Tout_r)
        
        else:        
            # Axes limits
            self.axes.set_xlim(0.7*Cycle.Condenser.sout_r/1000,1.2*Cycle.Compressor.sout_r/1000)
            self.axes.set_ylim(0.9*Cycle.Evaporator.Tout_r,1.1*Cycle.Compressor.Tout_r)
        
        # VI cycles compressor
        if Cycle.CycleType=='1INJ_Econ' or Cycle.CycleType=='1INJ_FlashTank':
            s_comp=np.r_[Cycle.Compressor.sin_r/1000.,Cycle.Compressor.s_31/1000,Cycle.Compressor.sinj_r/1000.,Cycle.Compressor.sout_r/1000.]
            T_comp=np.r_[Cycle.Compressor.Tin_r,Cycle.Compressor.T_31,Cycle.Compressor.Tinj_r,Cycle.Compressor.Tout_r]
            self.axes.plot(s_comp,T_comp,'k',lw = 2)
            self.axes.plot(s_comp[0],T_comp[0],'ko', mfc = 'w')
            self.axes.plot(s_comp[1],T_comp[1],'ko', mfc = 'w')
            self.axes.plot(s_comp[2],T_comp[2],'ko', mfc = 'w')
            self.axes.plot(s_comp[3],T_comp[3],'ko', mfc = 'w')
            self.axes.text(s_comp[0],T_comp[0],' 1',ha='left',va='top')
            #self.axes.text(h_comp[1],p_comp[1],' 11',ha='left',va='bottom')
            self.axes.text(s_comp[2],T_comp[2],' 1inj',ha='right',va='top')
            self.axes.text(s_comp[3],T_comp[3],' 2',ha='left',va='bottom')
            
            s_s_comp=np.r_[Cycle.Compressor.sin_r/1000.,Cycle.Compressor.s_31s/1000,Cycle.Compressor.s_41s/1000]
            T_s_comp=np.r_[Cycle.Compressor.Tin_r,Cycle.Compressor.Tinj_r,Cycle.Compressor.Tout_r]
            #self.axes.plot(s_s_comp,T_s_comp,'r--')
            
        else: 
        
            s_comp=np.r_[Cycle.Compressor.sin_r/1000.,Cycle.Compressor.sout_r/1000.]
            T_comp=np.r_[Cycle.Compressor.Tin_r,Cycle.Compressor.Tout_r]
            self.axes.plot(s_comp,T_comp,'k',lw=2)
            self.axes.plot(s_comp[0],T_comp[0],'ko',mfc='w')
            self.axes.plot(s_comp[1],T_comp[1],'ko',mfc='w')
            self.axes.text(s_comp[0],T_comp[0],' 1',ha='left',va='top')
            self.axes.text(s_comp[1],T_comp[1],' 2',ha='left',va='bottom')
           
        if (Cycle.CycleType=='Secondary' and Cycle.Mode=='AC') or Cycle.CycleType=='DX':
            # DX systems and secondary loop in cooling mode have condenser 
            sL=PropsSI('S','T',Cycle.Tdew_cond,'Q',0,Ref)/1000
            sV=PropsSI('S','T',Cycle.Tdew_cond,'Q',1,Ref)/1000
            s_cond=np.r_[Cycle.Condenser.sin_r/1000., sdew_cond/1000,sbubble_cond/1000,Cycle.Condenser.sout_r/1000.]
            T_cond=np.r_[Cycle.Condenser.Tin_r, Tdew_cond,Tbubble_cond,Cycle.Condenser.Tout_r]
            self.axes.plot(s_cond,T_cond,'k',lw=2)
            self.axes.plot(s_cond[3],T_cond[3],'ko',mfc='w')
            self.axes.text(s_cond[3],T_cond[3],'3$\quad\quad$',ha='right',va='bottom')
            
            self.axes.plot([s_cond[3],s_cond[0]],[Cycle.Condenser.Tin_a,Cycle.Condenser.Tout_a],'r-')

            self.axes.text(0.5*s_cond[0]+0.5*s_cond[3],Cycle.Condenser.Tin_a,'Outdoor Air',backgroundcolor='w',ha='center',va='center')
        
        elif Cycle.CycleType == 'CO2-DX' or Cycle.CycleType == 'CO2-LSHX-DX' or Cycle.CycleType == 'Transcrit-DX':
            # Transcritical cycles 
            s_gascool = np.linspace(Cycle.GasCooler.sin_r,Cycle.GasCooler.sout_r,10000,endpoint=True)
            p_gascool = Cycle.GasCooler.psat_r
        
            T_gascool = 0.0*s_gascool
        
            for i in np.arange(len(s_gascool)):
                T_gascool[i] = PropsSI('T','S',s_gascool[i],'P',p_gascool,Ref)
            
            self.axes.plot(s_gascool/1000,T_gascool,'k',lw=2)
    
        elif Cycle.CycleType=='1INJ_Econ' or Cycle.CycleType=='1INJ_FlashTank':
            # VI cycles condenser
            sL=PropsSI('S','T',Cycle.Tdew_cond,'Q',0,Ref)/1000
            sV=PropsSI('S','T',Cycle.Tdew_cond,'Q',1,Ref)/1000
            s_cond=np.r_[Cycle.Condenser.sin_r/1000., sdew_cond/1000,sbubble_cond/1000,Cycle.Condenser.sout_r/1000.]
            T_cond=np.r_[Cycle.Condenser.Tin_r, Tdew_cond,Tbubble_cond,Cycle.Condenser.Tout_r]
            self.axes.plot(s_cond,T_cond,'k',lw=2)
            self.axes.plot(s_cond[3],T_cond[3],'ko',mfc='w')
            #self.axes.text(s_cond[3],T_cond[3],'3$\quad\quad$',ha='right',va='bottom')
            
            self.axes.plot([s_cond[3],s_cond[0]],[Cycle.Condenser.Tin_a,Cycle.Condenser.Tout_a],'r-')
            
            self.axes.text(0.5*s_cond[0]+0.5*s_cond[3],Cycle.Condenser.Tin_a,'Outdoor Air',backgroundcolor='w',ha='center',va='center')       
        
        elif Cycle.CycleType=='Secondary' and Cycle.Mode=='HP':
            sV=PropsSI('S','T',Cycle.Tdew_evap,'Q',1,Ref)/1000
            s_evap=np.r_[Cycle.Evaporator.sin_r/1000., sV,Cycle.Evaporator.sout_r/1000.]
            T_evap=np.r_[Cycle.Evaporator.Tin_r, Cycle.Tdew_evap,Cycle.Evaporator.Tout_r]
            self.axes.plot(s_evap,T_evap,'k',lw=2)
            self.axes.plot(s_evap[0],T_evap[0],'ko',mfc='w')
            self.axes.text(s_evap[0],T_evap[0],'4$\quad\quad$',ha='right',va='bottom')
            
            self.axes.plot([s_evap[2],s_evap[0]],[Cycle.Evaporator.Tin_a,Cycle.Evaporator.Tout_a],'r-')
            
            self.axes.text(0.5*s_evap[0]+0.5*s_evap[2],Cycle.Evaporator.Tin_a,'Outdoor Air',backgroundcolor='w',ha='center',va='center')
                
        if Cycle.CycleType=="Secondary":
            if Cycle.IHXType=='Coaxial':
                IHX=Cycle.CoaxialIHX
                sV=PropsSI('S','T',Cycle.Tdew_evap,'Q',1,Ref)/1000
                s_IHX=np.r_[IHX.sin_r/1000., sV,IHX.sout_r/1000.]
                T_IHX=np.r_[IHX.Tin_r, Cycle.Tdew_evap,IHX.Tout_r]
                s_XV=np.r_[IHX.sin_r/1000.,Cycle.Condenser.sout_r/1000.]
                T_XV=np.r_[IHX.Tin_r,Cycle.Condenser.Tout_r]
                T_IHXg=np.r_[IHX.Tout_g,IHX.Tin_g]
                s_IHXg=np.r_[s_IHX[0],s_IHX[len(s_IHX)-1]]
            elif Cycle.IHXType=='PHE':
                if Cycle.Mode=='AC':
                    IHX=Cycle.PHEIHX
                    sV=PropsSI('S','T',Cycle.Tdew_evap,'Q',1,Ref)/1000
                    s_IHX=np.r_[IHX.sin_c/1000., sV,IHX.sout_c/1000.]
                    T_IHX=np.r_[IHX.Tin_c, Cycle.Tdew_evap,IHX.Tout_c]
                    s_XV=np.r_[IHX.sin_c/1000.,Cycle.Condenser.sout_r/1000.]
                    T_XV=np.r_[IHX.Tin_c,Cycle.Condenser.Tout_r]
                    T_IHXg=np.r_[IHX.Tout_h,IHX.Tin_h]
                    s_IHXg=np.r_[s_IHX[0],s_IHX[len(s_IHX)-1]]
                    self.axes.text(s_IHX[0],T_IHX[0],'4$\quad\quad$',ha='right',va='top')
                    self.axes.plot(s_IHX[0],T_IHX[0],'ko',mfc='w')
                else:
                    IHX=Cycle.PHEIHX
                    sL=PropsSI('S','T',Cycle.Tbubble_cond,'Q',0,Ref)/1000
                    sV=PropsSI('S','T',Cycle.Tdew_cond,'Q',1,Ref)/1000
                    s_IHX=np.r_[IHX.sin_h/1000., sV,sL,IHX.sout_h/1000.]
                    T_IHX=np.r_[IHX.Tin_h, Cycle.Tdew_cond,Cycle.Tbubble_cond,IHX.Tout_h]
                    s_XV=np.r_[Cycle.Evaporator.sin_r/1000.,IHX.sout_h/1000.]
                    T_XV=np.r_[Cycle.Evaporator.Tin_r,IHX.Tout_h]
                    T_IHXg=np.r_[IHX.Tout_c,IHX.Tin_c]
                    s_IHXg=np.r_[s_IHX[0],s_IHX[len(s_IHX)-1]]
                    self.axes.text(s_IHX[3],T_IHX[3],'3$\quad\quad$',ha='right',va='top')
                    self.axes.plot(s_IHX[3],T_IHX[3],'ko',mfc='w')
            else:
                raise ValueError('Secondary loop system must have a coaxial or PHE heat exchanger')
        
            self.axes.plot(s_IHX,T_IHX,'k',lw=2)
            
            
            self.axes.plot(s_XV,T_XV,'k',lw=2)            
            self.axes.plot([s_IHX[2],s_IHX[0]],[Cycle.CoolingCoil.Tin_a,Cycle.CoolingCoil.Tout_a],'b-')
            self.axes.text(0.5*s_IHX[2]+0.5*s_IHX[0],Cycle.CoolingCoil.Tin_a,'Indoor Air',backgroundcolor='w',ha='center',va='center')
            self.axes.text(0.5*s_IHXg[0]+0.5*s_IHXg[1],np.mean(T_IHXg),'Glycol IHX',backgroundcolor='w',ha='center',va='center')            
            self.axes.plot(s_IHXg,T_IHXg,'g-.')
            
        elif Cycle.CycleType=="DX":
            sL=PropsSI('S','T',Cycle.Tdew_evap,'Q',0,Ref)/1000
            sV=PropsSI('S','T',Cycle.Tdew_evap,'Q',1,Ref)/1000
            s_evap=np.r_[Cycle.Evaporator.sin_r/1000., sV,Cycle.Evaporator.sout_r/1000.]
            T_evap=np.r_[Cycle.Evaporator.Tsat_r, Cycle.Tdew_evap,Cycle.Evaporator.Tout_r]
            self.axes.plot(s_evap,T_evap,'k',lw=2)
            self.axes.plot(s_evap[0],T_evap[0],'ko',mfc='w')
            self.axes.text(s_evap[0],T_evap[0],'4$\quad\quad$',ha='right',va='top')
            
            s_XV=np.r_[Cycle.Evaporator.sin_r/1000.,Cycle.Condenser.sout_r/1000.]
            T_XV=np.r_[Cycle.Evaporator.Tsat_r,Cycle.Condenser.Tout_r]
            self.axes.plot(s_XV,T_XV,'k',lw=2)
            
            self.axes.plot([s_evap[2],s_evap[0]],[Cycle.Evaporator.Tin_a,Cycle.Evaporator.Tout_a],'b-')
            self.axes.text(0.5*s_evap[2]+0.5*s_evap[0],Cycle.Evaporator.Tin_a,'Indoor Air',backgroundcolor='w',ha='center',va='center')
       
        elif Cycle.CycleType == 'CO2-DX' or Cycle.CycleType == 'CO2-LSHX-DX' or Cycle.CycleType == 'Transcrit-DX':
            sL=PropsSI('S','T',Cycle.Tdew_evap,'Q',0,Ref)/1000
            sV=PropsSI('S','T',Cycle.Tdew_evap,'Q',1,Ref)/1000
            s_evap=np.r_[Cycle.Evaporator.sin_r/1000., sV,Cycle.Evaporator.sout_r/1000.]
            T_evap=np.r_[Cycle.Evaporator.Tsat_r, Cycle.Tdew_evap,Cycle.Evaporator.Tout_r]
            self.axes.plot(s_evap,T_evap,'k',lw=2)
            self.axes.plot(s_evap[0],T_evap[0],'ko',mfc='w')
            self.axes.text(s_evap[0],T_evap[0],'4$\quad\quad$',ha='right',va='top')
            
            s_XV=np.r_[Cycle.Evaporator.sin_r/1000.,Cycle.GasCooler.sout_r/1000.]
            T_XV=np.r_[Cycle.Evaporator.Tsat_r,Cycle.GasCooler.Tout_r]
            self.axes.plot(s_XV,T_XV,'k',lw=2)
            
            self.axes.plot([s_evap[2],s_evap[0]],[Cycle.Evaporator.Tin_a,Cycle.Evaporator.Tout_a],'b-')
            self.axes.text(0.5*s_evap[2]+0.5*s_evap[0],Cycle.Evaporator.Tin_a,'Indoor Air',backgroundcolor='w',ha='center',va='center')
                        
        elif Cycle.CycleType=='1INJ_Econ':  

            s_evap=np.r_[Cycle.Evaporator.sin_r/1000., sdew_evap/1000,Cycle.Evaporator.sout_r/1000.]
            T_evap=np.r_[Cycle.Evaporator.Tin_r, Cycle.Tdew_evap,Cycle.Evaporator.Tout_r]
            self.axes.plot(s_evap,T_evap,'k',lw=2)
            self.axes.plot(s_evap[0],T_evap[0],'ko',mfc='w')
            self.axes.text(s_evap[0],T_evap[0],'7$\quad\quad$',ha='right',va='top')

            s_PHEHX = np.r_[Cycle.PHEHX.sout_h/1000.,Cycle.Condenser.sout_r/1000.]
            T_PHEHX = np.r_[Cycle.PHEHX.Tout_h,Cycle.Condenser.Tout_r]
            self.axes.plot(s_PHEHX,T_PHEHX,'k', lw = 2)
            self.axes.plot(s_PHEHX[0],T_PHEHX[0],'ko', mfc = 'w')
            self.axes.text(s_PHEHX[0],T_PHEHX[0],'4$\quad\quad$',ha='right',va='bottom')
            s_XV=np.r_[Cycle.Evaporator.sin_r/1000.,Cycle.PHEHX.sout_h/1000.]
            T_XV=np.r_[Cycle.Evaporator.Tin_r,Cycle.PHEHX.Tout_h]
            self.axes.plot(s_XV,T_XV,'k', lw = 2)   

            #h_XV2=np.r_[Cycle.Condenser.hout_r/1000.,Cycle.Condenser.hout_r/1000.]
            #p_XV2=np.r_[Cycle.Condenser.psat_r/1000.,Cycle.PHEHX.pin_c/1000.]
            #self.axes.plot(h_XV[1],p_XV[1],'bo')
            #self.axes.plot(h_XV2,p_XV2,'b')                  
            
            sdew_inj = PropsSI('S','P',Cycle.Compressor.pinj_r,'Q',1,Cycle.Ref)
            Tdew_inj = PropsSI('T','P',Cycle.Compressor.pinj_r,'Q',1,Cycle.Ref)
            
            s_1INJ=np.r_[Cycle.PHEHX.sout_h/1000,sdew_inj/1000,Cycle.PHEHX.sout_c/1000.]
            T_1INJ=np.r_[Cycle.PHEHX.Tin_c,Tdew_inj,Cycle.PHEHX.Tout_c]
            self.axes.plot(s_1INJ,T_1INJ,'k', lw = 2)
            self.axes.plot(s_1INJ[0],T_1INJ[0],'ko', mfc = 'w')
            self.axes.text(s_1INJ[0],T_1INJ[0],'5$\quad\quad$',ha='right',va='bottom')

            self.axes.plot([s_evap[2],s_evap[0]],[Cycle.Evaporator.Tin_a,Cycle.Evaporator.Tout_a],'b-')
            self.axes.text(0.5*s_evap[2]+0.5*s_evap[0],Cycle.Evaporator.Tin_a,'Indoor Air',backgroundcolor='w',ha='center',va='center')


        elif Cycle.CycleType=='1INJ_FlashTank':  

            s_evap=np.r_[Cycle.Evaporator.sin_r/1000., sdew_evap/1000,Cycle.Evaporator.sout_r/1000.]
            T_evap=np.r_[Cycle.Evaporator.Tin_r, Tdew_evap,Cycle.Evaporator.Tout_r]
            self.axes.plot(s_evap,T_evap,'k',lw=2)
            self.axes.plot(s_evap[0],T_evap[0],'ko',mfc='w')
            self.axes.text(s_evap[0],T_evap[0],'9$\quad\quad$',ha='right',va='top')

            s_XV=np.r_[Cycle.FlashTank.sin/1000.,Cycle.Condenser.sout_r/1000.]
            T_XV=np.r_[Cycle.FlashTank.Tin,Cycle.Condenser.Tout_r]
            self.axes.plot(s_XV,T_XV,'k', lw = 2) 

            s_FlashTank_liquid = np.r_[Cycle.FlashTank.sin/1000.,Cycle.FlashTank.sout/1000.]
            T_FlashTank_liquid = np.r_[Cycle.FlashTank.Tin,Cycle.FlashTank.Tout]
            self.axes.plot(s_FlashTank_liquid,T_FlashTank_liquid,'k', lw = 2)
            self.axes.plot(s_FlashTank_liquid[0],T_FlashTank_liquid[0],'ko', mfc = 'w')
            self.axes.text(s_FlashTank_liquid[0],T_FlashTank_liquid[0],'6$\quad\quad$',ha='right',va='top')
            
            
            s_1INJ=np.r_[Cycle.FlashTank.sin/1000,Cycle.Compressor.sinj_r/1000.]
            T_1INJ=np.r_[Cycle.FlashTank.Tin,Cycle.Compressor.Tinj_r]
            self.axes.plot(s_1INJ,T_1INJ,'k', lw = 2)
            self.axes.plot(s_1INJ[0],T_1INJ[0],'ro', mfc = 'w')
            self.axes.text(s_1INJ[0],T_1INJ[0],'8$\quad\quad$',ha='right',va='bottom')

            s_XV2=np.r_[Cycle.FlashTank.sout/1000.,Cycle.Evaporator.sin_r/1000]
            T_XV2=np.r_[Cycle.FlashTank.Tin,Cycle.Evaporator.Tsat_r]
            self.axes.plot(s_XV[1],T_XV[1],'ko', mfc = 'w')
            self.axes.plot(s_XV2,T_XV2,'k', lw = 2) 
            self.axes.plot(s_XV2[0],T_XV2[0],'ko', mfc = 'w')
            self.axes.text(s_XV2[0],T_XV2[0],'8$\quad\quad$',ha='right',va='top')

            self.axes.plot([s_evap[2],s_evap[0]],[Cycle.Evaporator.Tin_a,Cycle.Evaporator.Tout_a],'b-')
            self.axes.text(0.5*s_evap[2]+0.5*s_evap[0],Cycle.Evaporator.Tin_a,'Indoor Air',backgroundcolor='w',ha='center',va='center')
            
        plt.show()
        
    def PHOverlay(self,Cycle,**kwargs):
  
        Ref=Cycle.Ref

        Tmin=PropsSI('Tmin',Ref)
        Tmax=PropsSI('Tmax',Ref)  
        pmin=PropsSI('pmin',Ref)
        pmax=PropsSI('pmax',Ref)
        
        fig = plt.figure(figsize=(6,4))
        self.axes=fig.add_subplot(111)
        
        # Create saturated lines
        Tsat = np.linspace(Tmin,PropsSI(Ref,"Tcrit")-0.0000000001,1000)
        (hsatL,hsatV)=(0.0*Tsat,0.0*Tsat)
        (psatL,psatV)=(0.0*Tsat,0.0*Tsat)
        for i in np.arange(len(Tsat)):
            psatL[i] = PropsSI('P','T',Tsat[i],'Q',0,Ref)/1000
            psatV[i] = PropsSI('P','T',Tsat[i],'Q',1,Ref)/1000
            hsatL[i] = PropsSI('H','T',Tsat[i],'Q',0,Ref)/1000
            hsatV[i] = PropsSI('H','T',Tsat[i],'Q',1,Ref)/1000         
        
        # Create critical iso-therm
        p_critical = Cycle.AS.p_critical()
        T_critical = Cycle.AS.T_critical()
        
        p_lim = np.linspace(pmin,pmax, 10000)
        h_lim = 0.0*p_lim
        
        for i in np.arange(len(p_lim)):
            h_lim[i] = PropsSI('H','P',p_lim[i],'T',T_critical-0.0000000001,Ref)/1000
        
        # Check if the cycle is transcritical
        if Cycle.Compressor.Tout_r > T_critical  and Cycle.Compressor.pout_r > p_critical:
            # Average saturated temperatures and pressures
            Tdew_evap = Cycle.Tdew_evap
            psat_evap = PropsSI('P','T',Tdew_evap,'Q',1,Ref)      

            # Saturated conditions at the evaporator
            hdew_evap = PropsSI('H','P',psat_evap,'Q',1,Ref)
            Tdew_evap = PropsSI('T','P',psat_evap,'Q',1,Ref)
            hbubble_evap = PropsSI('H','P',psat_evap,'Q',0,Ref)
            Tbubble_evap = PropsSI('T','P',psat_evap,'Q',0,Ref)        

            pevap_range1 = np.linspace(pmin,psat_evap,1000)
            pevap_range2 = np.linspace(psat_evap,pmax,1000)
            hevap_1 = 0.0*pevap_range1
            hevap_2 = 0.0*pevap_range2        

            for i in range(len(pevap_range1)):        

                try:    
                    hevap_1[i] = PropsSI('H','T',Tdew_evap,'P',pevap_range1[i],Ref)/1000
                except:
                    hevap_1[i] = PropsSI('H','P',psat_evap,'Q',1,Ref)/1000
                    
                try:
                    hevap_2[i] = PropsSI('H','T',Tbubble_evap,'P',pevap_range2[i],Ref)/1000
                except:
                    hevap_2[i] = PropsSI('H','P',psat_evap,'Q',0,Ref)/1000        
        
        else:
            # Average saturated temperatures and pressures
            Tdew_evap = Cycle.Tdew_evap
            Tdew_cond = Cycle.Tdew_cond
            psat_evap = PropsSI('P','T',Tdew_evap,'Q',1,Ref)
            psat_cond = PropsSI('P','T',Tdew_cond,'Q',1,Ref)
            
            # Saturated conditions at the evaporator
            hdew_evap = PropsSI('H','P',psat_evap,'Q',1,Ref)
            Tdew_evap = PropsSI('T','P',psat_evap,'Q',1,Ref)
            hbubble_evap = PropsSI('H','P',psat_evap,'Q',0,Ref)
            Tbubble_evap = PropsSI('T','P',psat_evap,'Q',0,Ref)
            
            # Saturated conditions at the condenser
            hdew_cond = PropsSI('H','P',psat_cond,'Q',1,Ref)
            Tdew_cond = PropsSI('T','P',psat_cond,'Q',1,Ref)
            hbubble_cond = PropsSI('H','P',psat_cond,'Q',0,Ref)
            Tbubble_cond = PropsSI('T','P',psat_cond,'Q',0,Ref)
    
    
            pcond_range1 = np.linspace(pmin,psat_cond,1000)
            pcond_range2 = np.linspace(psat_cond,pmax,1000)
            hcond_1 = 0.0*pcond_range1
            hcond_2 = 0.0*pcond_range2
            pevap_range1 = np.linspace(pmin,psat_evap,1000)
            pevap_range2 = np.linspace(psat_evap,pmax,1000)
            hevap_1 = 0.0*pevap_range1
            hevap_2 = 0.0*pevap_range2
            
            for i in range(len(pcond_range1)):
                
                try:
                    hcond_1[i] = PropsSI('H','T',Tdew_cond,'P',pcond_range1[i],Ref)/1000
                except:
                    hcond_1[i] = PropsSI('H','P',psat_cond,'Q',1,Ref)/1000
                
                try:
                    hcond_2[i] = PropsSI('H','T',Tbubble_cond,'P',pcond_range2[i],Ref)/1000
                except:
                    hcond_2[i] = PropsSI('H','P',psat_cond,'Q',0,Ref)/1000
                
                try:    
                    hevap_1[i] = PropsSI('H','T',Tdew_evap,'P',pevap_range1[i],Ref)/1000
                except:
                    hevap_1[i] = PropsSI('H','P',psat_evap,'Q',1,Ref)/1000
                    
                try:
                    hevap_2[i] = PropsSI('H','T',Tbubble_evap,'P',pevap_range2[i],Ref)/1000
                except:
                    hevap_2[i] = PropsSI('H','P',psat_evap,'Q',0,Ref)/1000
 
            # Isolines condensing temperature
            self.axes.semilogy(hcond_1,pcond_range1/1000,'k-.',lw=1, alpha = 0.3)
            self.axes.semilogy(hcond_2,pcond_range2/1000,'k-.',lw=1, alpha = 0.3)
            self.axes.semilogy([hdew_cond/1000,hbubble_cond/1000],[psat_cond/1000,psat_cond/1000],'k-.',lw=1, alpha = 0.3)               
    
        # Isolines evaporating temperature
        self.axes.semilogy(hevap_1,pevap_range1/1000,'k-.',lw=1, alpha = 0.3)
        self.axes.semilogy(hevap_2,pevap_range2/1000,'k-.',lw=1, alpha = 0.3)
        self.axes.semilogy([hdew_evap/1000,hbubble_evap/1000],[psat_evap/1000,psat_evap/1000],'k-.',lw=1, alpha = 0.3)
        
        # Critical isothermal
        self.axes.semilogy(h_lim,p_lim/1000,'k-.',lw=1, alpha = 0.3)
        
        # ph saturated line plot
        self.axes.semilogy(hsatL,psatL,'k', lw = 1.5, alpha = 0.5)
        self.axes.semilogy(hsatV,psatV,'k', lw = 1.5, alpha = 0.5)

        # Axes labels
        self.axes.set_xlabel('Enthalpy [kJ/kg]')
        self.axes.set_ylabel('Pressure [kPa]')
        
        # Check if the cycle is transcritical
        if Cycle.Compressor.Tout_r > T_critical  and Cycle.Compressor.pout_r > p_critical: 
            # Axes limits
            self.axes.set_xlim(0.7*Cycle.GasCooler.hout_r/1000,1.1*Cycle.Compressor.hout_r/1000)        
            self.axes.set_ylim(10*pmin/1000,2*Cycle.GasCooler.psat_r/1000)
        
        else:
            # Axes limits
            self.axes.set_xlim(0.7*Cycle.Condenser.hout_r/1000,1.2*Cycle.Compressor.hout_r/1000)
            self.axes.set_ylim(10*pmin/1000,0.5*pmax/1000)

        # Check cycle type to plot the cycle
        if Cycle.CycleType=='1INJ_Econ' or Cycle.CycleType=='1INJ_FlashTank':
            h_comp=np.r_[Cycle.Compressor.hin_r/1000.,Cycle.Compressor.h_31/1000,Cycle.Compressor.hinj_r/1000.,Cycle.Compressor.hout_r/1000.]
            p_comp=np.r_[Cycle.Compressor.pin_r/1000,Cycle.Compressor.pinj_r/1000,Cycle.Compressor.pinj_r/1000,Cycle.Compressor.pout_r/1000]
            self.axes.semilogy(h_comp,p_comp,'k',lw = 2)
            self.axes.semilogy(h_comp[0],p_comp[0],'ko', mfc = 'w')
            self.axes.semilogy(h_comp[1],p_comp[1],'ko', mfc = 'w')
            self.axes.semilogy(h_comp[2],p_comp[2],'ko', mfc = 'w')
            self.axes.semilogy(h_comp[3],p_comp[3],'ko', mfc = 'w')
            self.axes.text(h_comp[0],p_comp[0],' 1',ha='left',va='top')
            #self.axes.text(h_comp[1],p_comp[1],' 11',ha='left',va='bottom')
            self.axes.text(h_comp[2],p_comp[2],' 1inj',ha='left',va='bottom')
            self.axes.text(h_comp[3],p_comp[3],' 2',ha='left',va='bottom')
            
            h_s_comp=np.r_[Cycle.Compressor.hin_r/1000.,Cycle.Compressor.h_31s/1000,Cycle.Compressor.h_41s/1000]
            p_s_comp=np.r_[Cycle.Compressor.pin_r/1000.,Cycle.Compressor.pinj_r/1000,Cycle.Compressor.pout_r/1000.]
            #self.axes.plot(h_s_comp,p_s_comp,'r--')
            
            #self.axes.plot(Cycle.Compressor.h_31/1000,Cycle.Compressor.pinj_r/1000,'bo')
            #self.axes.plot(Cycle.Compressor.h_31s/1000,Cycle.Compressor.pinj_r/1000,'bo')
            #self.axes.plot(Cycle.Compressor.h_22s/1000,Cycle.Compressor.pinj_r/1000,'bo')
            #self.axes.plot(Cycle.Compressor.h_22/1000,Cycle.Compressor.pinj_r/1000,'bo')
            #self.axes.plot(Cycle.Compressor.h_32s/1000,Cycle.Compressor.pout_r/1000,'bo')
            #self.axes.plot(Cycle.Compressor.h_41s/1000,Cycle.Compressor.pout_r/1000,'ro')
            
        else:        
            h_comp=np.r_[Cycle.Compressor.hin_r/1000.,Cycle.Compressor.hout_r/1000.]
            p_comp=np.r_[Cycle.Compressor.pin_r/1000,Cycle.Compressor.pout_r/1000]
            self.axes.semilogy(h_comp,p_comp,'k', lw = 2)
            self.axes.semilogy(h_comp[0],p_comp[0],'ko', mfc = 'w')
            self.axes.semilogy(h_comp[1],p_comp[1],'ko', mfc = 'w')
            self.axes.text(h_comp[0],p_comp[0],' 1',ha='left',va='top')
            self.axes.text(h_comp[1],p_comp[1],' 2',ha='left',va='bottom')

        if Cycle.CycleType=='Secondary' and Cycle.Mode=='HP':
            h_evap=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Evaporator.hout_r/1000]
            p_evap=np.r_[Cycle.Evaporator.psat_r/1000, Cycle.Evaporator.psat_r/1000]
            self.axes.semilogy(h_evap,p_evap,'k', lw = 2)
            self.axes.semilogy(h_evap[0],p_evap[0],'ko', mfc = 'w')
            self.axes.text(h_evap[0],p_evap[0],'4$\quad\quad$',ha='right',va='bottom')
        elif Cycle.CycleType=='1INJ_Econ' or Cycle.CycleType=='1INJ_FlashTank':
            h_cond=np.r_[Cycle.Condenser.hin_r/1000.,Cycle.Condenser.hout_r/1000.]
            p_cond=np.r_[Cycle.Compressor.pout_r/1000, Cycle.Compressor.pout_r/1000]
            self.axes.semilogy(h_cond,p_cond,'k', lw = 2)
            self.axes.semilogy(h_cond[1],p_cond[1],'ko',mfc='w')
            #self.axes.semilogy(h_cond[1],p_cond[1],'3$\quad\quad$',ha='right',va='bottom')
        elif Cycle.CycleType == 'CO2-DX' or Cycle.CycleType == 'CO2-LSHX-DX' or Cycle.CycleType == 'Transcrit-DX':
            h_gascool=np.r_[Cycle.Compressor.hout_r/1000.,Cycle.GasCooler.hout_r/1000.]
            p_gascool=np.r_[Cycle.GasCooler.psat_r/1000, (Cycle.GasCooler.psat_r-Cycle.GasCooler.DP_r)/1000]
            self.axes.semilogy(h_gascool,p_gascool,'k', lw = 2)
            self.axes.semilogy(h_gascool[1],p_gascool[1],'ko', mfc = 'w')
            self.axes.text(h_gascool[1],p_gascool[1],'3$\quad\quad$',ha='right',va='bottom')            
        else:
            h_cond=np.r_[Cycle.Condenser.hin_r/1000.,Cycle.Condenser.hout_r/1000.]
            p_cond=np.r_[Cycle.Condenser.psat_r/1000, (Cycle.Condenser.psat_r-Cycle.Condenser.DP_r)/1000]
            self.axes.semilogy(h_cond,p_cond,'k', lw = 2)
            self.axes.semilogy(h_cond[1],p_cond[1],'ko', mfc = 'w')
            self.axes.text(h_cond[1],p_cond[1],'3$\quad\quad$',ha='right',va='bottom')
        
        if Cycle.CycleType=="Secondary":
            if Cycle.IHXType=='Coaxial':
                IHX=Cycle.CoaxialIHX
                h_IHX=np.r_[IHX.hin_r/1000., IHX.hout_r/1000.]
                p_IHX=np.r_[IHX.pin_r/1000, IHX.pin_r]
                h_XV=np.r_[IHX.hin_r/1000.,Cycle.Condenser.hout_r/1000.]
                p_XV=np.r_[IHX.pin_r/1000,Cycle.Condenser.psat_r/1000]
            elif Cycle.IHXType=='PHE':
                IHX=Cycle.PHEIHX
                
                if Cycle.Mode=='AC':
                    h_IHX=np.r_[IHX.hin_c/1000., IHX.hout_c/1000.]
                    p_IHX=np.r_[IHX.pin_c/1000, IHX.pin_c/1000]
                    h_XV=np.r_[IHX.hin_c/1000.,Cycle.Condenser.hout_r/1000.]
                    p_XV=np.r_[IHX.pin_c/1000,Cycle.Condenser.psat_r/1000]
                    self.axes.semilogy(h_IHX[0],p_IHX[0],'ko', mfc = 'w')
                    self.axes.text(h_IHX[0],p_IHX[0],'4$\quad\quad$',ha='right',va='top')
                else:
                    h_IHX=np.r_[IHX.hin_h/1000., IHX.hout_h/1000.]
                    p_IHX=np.r_[IHX.pin_h/1000, IHX.pin_h/1000]
                    h_XV=np.r_[Cycle.Evaporator.hin_r/1000.,IHX.hout_h/1000.]
                    p_XV=np.r_[Cycle.Evaporator.psat_r/1000,IHX.pin_h/1000]
                    self.axes.semilogy(h_IHX[1],p_IHX[1],'ko', mfc = 'w')
                    self.axes.text(h_IHX[1],p_IHX[1],'3$\quad\quad$',ha='right',va='top')
            else:
                raise ValueError('Secondary loop system must have a coaxial or PHE heat exchanger')
            
            self.axes.semilogy(h_IHX,p_IHX,'k', lw = 2)
            
            self.axes.semilogy(h_XV,p_XV,'b', lw= 2)
            
        elif Cycle.CycleType=="DX":
            h_evap=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Evaporator.hout_r/1000.]
            p_evap=np.r_[Cycle.Evaporator.psat_r/1000,Cycle.Evaporator.psat_r/1000]
            self.axes.semilogy(h_evap,p_evap,'k', lw = 2)
            self.axes.semilogy(h_evap[0],p_evap[0],'ko', mfc = 'w')
            self.axes.text(h_evap[0],p_evap[0],'4$\quad\quad$',ha='right',va='top')
            
            h_XV=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Condenser.hout_r/1000.]
            p_XV=np.r_[Cycle.Evaporator.psat_r/1000,Cycle.Condenser.psat_r/1000]
            self.axes.semilogy(h_XV,p_XV,'k', lw = 2)
        elif Cycle.CycleType == 'CO2-DX' or Cycle.CycleType == 'CO2-LSHX-DX' or Cycle.CycleType == 'Transcrit-DX':
            h_evap=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Evaporator.hout_r/1000.]
            p_evap=np.r_[Cycle.Evaporator.psat_r/1000,Cycle.Evaporator.psat_r/1000]
            self.axes.semilogy(h_evap,p_evap,'k', lw = 2)
            self.axes.semilogy(h_evap[0],p_evap[0],'ko', mfc = 'w')
            self.axes.text(h_evap[0],p_evap[0],'4$\quad\quad$',ha='right',va='top')
            
            h_XV=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.GasCooler.hout_r/1000.]
            p_XV=np.r_[Cycle.Evaporator.psat_r/1000,Cycle.GasCooler.psat_r/1000]
            self.axes.semilogy(h_XV,p_XV,'k', lw = 2)                        
        elif Cycle.CycleType=='1INJ_Econ':  
            h_evap=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Evaporator.hout_r/1000.]
            p_evap=np.r_[Cycle.Evaporator.psat_r/1000,Cycle.Evaporator.psat_r/1000]
            self.axes.semilogy(h_evap,p_evap,'k', lw = 2)
            self.axes.semilogy(h_evap[0],p_evap[0],'ko', mfc = 'w')
            self.axes.text(h_evap[0],p_evap[0],'7$\quad\quad$',ha='right',va='top')

            h_PHEHX = np.r_[Cycle.PHEHX.hout_h/1000.,Cycle.Condenser.hout_r/1000.]
            p_PHEHX = np.r_[Cycle.Compressor.pout_r/1000,Cycle.Compressor.pout_r/1000]
            self.axes.semilogy(h_PHEHX,p_PHEHX,'k', lw = 2)
            self.axes.semilogy(h_PHEHX[0],p_PHEHX[0],'ko', mfc = 'w')
            self.axes.text(h_PHEHX[0],p_PHEHX[0],'4$\quad\quad$',ha='right',va='top')
            h_XV=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.PHEHX.hout_h/1000.]
            p_XV=np.r_[Cycle.Evaporator.psat_r/1000,(Cycle.PHEHX.pin_h-Cycle.PHEHX.DP_h)/1000.]
            self.axes.plot(h_XV,p_XV,'k', lw = 2)   

            h_XV2=np.r_[Cycle.Condenser.hout_r/1000.,Cycle.Condenser.hout_r/1000.]
            p_XV2=np.r_[Cycle.Condenser.psat_r/1000.,Cycle.PHEHX.pin_c/1000.]
            #self.axes.plot(h_XV[1],p_XV[1],'bo')
            #self.axes.plot(h_XV2,p_XV2,'b')                  
            
            h_1INJ=np.r_[Cycle.PHEHX.hout_h/1000,Cycle.PHEHX.hout_c/1000.]
            p_1INJ=np.r_[Cycle.PHEHX.pin_c/1000,(Cycle.PHEHX.pin_c-Cycle.PHEHX.DP_c)/1000.]
            self.axes.semilogy(h_1INJ,p_1INJ,'k', lw = 2)
            self.axes.semilogy(h_1INJ[0],p_1INJ[0],'ko', mfc = 'w')
            self.axes.text(h_1INJ[0],p_1INJ[0],'5$\quad\quad$',ha='right',va='top')

        elif Cycle.CycleType=='1INJ_FlashTank':  

            h_evap=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Evaporator.hout_r/1000.]
            p_evap=np.r_[Cycle.Evaporator.psat_r/1000,Cycle.Evaporator.psat_r/1000]
            self.axes.semilogy(h_evap,p_evap,'k', lw = 2)
            self.axes.semilogy(h_evap[0],p_evap[0],'ko', mfc = 'w')
            self.axes.text(h_evap[0],p_evap[0],'9$\quad\quad$',ha='right',va='top')

            h_XV=np.r_[Cycle.FlashTank.hin/1000.,Cycle.Condenser.hout_r/1000.]
            p_XV=np.r_[Cycle.FlashTank.pin/1000,Cycle.Compressor.pout_r/1000]
            self.axes.plot(h_XV,p_XV,'k', lw = 2) 

            h_FlashTank_liquid = np.r_[Cycle.FlashTank.hin/1000.,Cycle.FlashTank.hout/1000.]
            p_FlashTank_liquid = np.r_[Cycle.FlashTank.pin/1000,Cycle.FlashTank.pin/1000]
            self.axes.semilogy(h_FlashTank_liquid,p_FlashTank_liquid,'k', lw = 2)
            self.axes.semilogy(h_FlashTank_liquid[0],p_FlashTank_liquid[0],'ko', mfc = 'w')
            self.axes.text(h_FlashTank_liquid[0],p_FlashTank_liquid[0],'6$\quad\quad$',ha='right',va='top')
            
            
            h_1INJ=np.r_[Cycle.FlashTank.hin/1000,Cycle.Compressor.hinj_r/1000.]
            p_1INJ=np.r_[Cycle.FlashTank.pin/1000,(Cycle.Compressor.pinj_r)/1000.]
            self.axes.semilogy(h_1INJ,p_1INJ,'k', lw = 2)
            self.axes.semilogy(h_1INJ[0],p_1INJ[0],'ro', mfc = 'w')
            self.axes.text(h_1INJ[0],p_1INJ[0],'8$\quad\quad$',ha='right',va='top')

            h_XV2=np.r_[Cycle.FlashTank.hout/1000.,Cycle.Evaporator.hin_r/1000.]
            p_XV2=np.r_[Cycle.FlashTank.pin/1000.,Cycle.Evaporator.psat_r/1000.]
            self.axes.semilogy(h_XV[1],p_XV[1],'ko', mfc = 'w')
            self.axes.semilogy(h_XV2,p_XV2,'k', lw = 2) 
            self.axes.semilogy(h_XV2[0],p_XV2[0],'ko', mfc = 'w')
            self.axes.text(h_XV2[0],p_XV2[0],'8$\quad\quad$',ha='right',va='top')

        plt.show()

    
    def PHEjector(self,Ejector,**kwargs):
  
        Ref=Ejector.Ref

        Tmin=PropsSI('Tmin',Ref)
        Tmax=PropsSI('Tmax',Ref)  
        pmin=PropsSI('pmin',Ref)
        pmax=PropsSI('pmax',Ref)
        

        fig = plt.figure(figsize=(6,4))
        self.axes=fig.add_subplot(111)
        
        Tsat = np.linspace(Tmin,PropsSI(Ref,"Tcrit")-0.0000000001,1000)
        (hsatL,hsatV)=(0.0*Tsat,0.0*Tsat)
        (psatL,psatV)=(0.0*Tsat,0.0*Tsat)
        for i in np.arange(len(Tsat)):
            psatL[i] = PropsSI('P','T',Tsat[i],'Q',0,Ref)/1000
            psatV[i] = PropsSI('P','T',Tsat[i],'Q',1,Ref)/1000
            hsatL[i] = PropsSI('H','T',Tsat[i],'Q',0,Ref)/1000
            hsatV[i] = PropsSI('H','T',Tsat[i],'Q',1,Ref)/1000

        h_xd = np.zeros_like((Tsat))
        psat = np.linspace(pmin,0.999*PropsSI(Ref,"pcrit"),1000)
        for i in np.arange(len(psat)):
            h_xd[i] = PropsSI('H','P',psat[i],'Q',Ejector.x_d,Ref)/1000

        # ph saturated line plot
        self.axes.semilogy(hsatL,psatL,'k', lw = 1.5, alpha = 0.5)
        self.axes.semilogy(hsatV,psatV,'k', lw = 1.5, alpha = 0.5)
        self.axes.semilogy(h_xd,psat/1000,'k', lw = 1.5, alpha = 0.3)
        
        self.axes.text(Ejector.h_d/1000,0.8*Ejector.psat_ev/1000,'%0.2f$\quad\quad$' % Ejector.x_d,ha='right',va='top',backgroundcolor='w', color='k',alpha=0.3)        
        
        # Axes labels
        self.axes.set_xlabel('h [kJ/kg]')
        self.axes.set_ylabel('p [kPa]')       

        # Axes limits
        self.axes.set_xlim(0.7*Ejector.hout_gc/1000,1.1*Ejector.hout_ev/1000)
        self.axes.set_ylim(0.7*Ejector.psat_ev/1000,1.2*Ejector.pout_gc/1000)        

        # Gas cooler outlet - ejector inlet
        h_gc_mb = np.r_[Ejector.hout_gc/1000,Ejector.h_mb/1000]
        p_gc_mb = np.r_[Ejector.pout_gc/1000,Ejector.p_mb/1000]
        self.axes.semilogy(h_gc_mb,p_gc_mb,'ko-', mfc = 'w')
        self.axes.text(h_gc_mb[0],p_gc_mb[0],'gc$\quad\quad$',ha='left',va='bottom')
        self.axes.text(h_gc_mb[1],p_gc_mb[1],'mb$\quad\quad$',ha='right',va='top')

        # Evaporator outlet
        hout_ev = Ejector.hout_ev/1000
        psat_ev = Ejector.psat_ev/1000
        self.axes.semilogy(hout_ev,psat_ev,'ko', mfc = 'w')
        self.axes.text(hout_ev,psat_ev,'ev$\quad\quad$',ha='left',va='bottom')

        # Ejector suction stream to mix
        h_sb_mix = np.r_[Ejector.h_sb/1000,Ejector.h_mix/1000]
        p_sb_mix = np.r_[Ejector.p_sb/1000,Ejector.p_mix/1000]
        self.axes.semilogy(h_sb_mix,p_sb_mix,'ko-', mfc = 'w')
        self.axes.text(h_sb_mix[0],p_sb_mix[0],'sb$\quad\quad$',ha='right',va='top')
        self.axes.text(h_sb_mix[1],p_sb_mix[1],'mix$\quad\quad$',ha='right',va='top')        

        # Ejector motion stream to mix
        h_mb_mix = np.r_[Ejector.h_mb/1000,Ejector.h_mix/1000]
        p_mb_mix = np.r_[Ejector.p_mb/1000,Ejector.p_mix/1000]
        self.axes.semilogy(h_mb_mix,p_mb_mix,'ko-', mfc = 'w')   
        
        # Ejector discharge
        h_d = Ejector.h_d/1000
        p_d = Ejector.p_d/1000
        self.axes.semilogy(h_d,p_d,'or', mfc = 'w')
        self.axes.text(h_d,p_d,'d$\quad\quad$',ha='left',va='bottom')        

        self.axes.yaxis.set_major_formatter(ScalarFormatter())
        self.axes.yaxis.set_minor_formatter(ScalarFormatter())
        
        fig.savefig('Ejector.png',dpi=600)
        plt.show()
        
        