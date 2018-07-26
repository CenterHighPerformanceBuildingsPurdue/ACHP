import sys,os

from matplotlib.figure import Figure
from matplotlib import pyplot as plt
import numpy as np
from CoolProp.CoolProp import PropsSI

class PlotsClass():
    
    def __init__(self, *args, **kwds):
        pass    
        
    def TSOverlay(self,Cycle,**kwargs):
        Ref=Cycle.Ref
        Tmin=PropsSI('Tmin',Ref)
        Tmax=PropsSI('Tmax',Ref)   
        

        fig = plt.figure(figsize=(6,4))
        self.axes=fig.add_subplot(111)
    
        Tsat = np.linspace(Tmin,PropsSI(Ref,"Tcrit")-0.0000000001,1000)
        (ssatL,ssatV)=(0.0*Tsat,0.0*Tsat)
        hcrit=np.linspace(0,0,2)
        pcrit=np.linspace(0,0,2)
        for i in np.arange(len(Tsat)):
            ssatL[i] = PropsSI('S','T',Tsat[i],'Q',0,Ref)/1000
            ssatV[i] = PropsSI('S','T',Tsat[i],'Q',1,Ref)/1000
    
        self.axes.plot(ssatL,Tsat,'k')
        self.axes.plot(ssatV,Tsat,'k')
        
        self.axes.set_xlabel('Entropy [kJ/kg-K]')
        self.axes.set_ylabel('Temperature [K]')
        
        s_comp=np.r_[Cycle.Compressor.sin_r/1000.,Cycle.Compressor.sout_r/1000.]
        T_comp=np.r_[Cycle.Compressor.Tin_r,Cycle.Compressor.Tout_r]
        self.axes.plot(s_comp,T_comp,'b')
        self.axes.plot(s_comp[0],T_comp[0],'bo')
        self.axes.plot(s_comp[1],T_comp[1],'bo')
        self.axes.text(s_comp[0],T_comp[0],' 1',ha='left',va='top')
        self.axes.text(s_comp[1],T_comp[1],' 2',ha='left',va='bottom')
        
        #DX systems and secondary loop in cooling mode have condenser    
        if (Cycle.CycleType=='Secondary' and Cycle.Mode=='AC') or Cycle.CycleType=='DX':
            sL=PropsSI('S','T',Cycle.Tdew_cond,'Q',0,Ref)/1000
            sV=PropsSI('S','T',Cycle.Tdew_cond,'Q',1,Ref)/1000
            s_cond=np.r_[Cycle.Condenser.sin_r/1000., sV,sL,Cycle.Condenser.sout_r/1000.]
            T_cond=np.r_[Cycle.Condenser.Tin_r, Cycle.Tdew_cond,Cycle.Tdew_cond,Cycle.Condenser.Tout_r]
            self.axes.plot(s_cond,T_cond,'b')
            self.axes.plot(s_cond[3],T_cond[3],'bo')
            self.axes.text(s_cond[3],T_cond[3],'3$\quad\quad$',ha='right',va='bottom')
            
            self.axes.plot([s_cond[3],s_cond[0]],[Cycle.Condenser.Tin_a,Cycle.Condenser.Tout_a],'r')
            self.axes.text(0.5*s_cond[0]+0.5*s_cond[3],Cycle.Condenser.Tin_a,'Outdoor Air',backgroundcolor='w',ha='center',va='center')
        elif Cycle.CycleType=='Secondary' and Cycle.Mode=='HP':
            sV=PropsSI('S','T',Cycle.Tdew_evap,'Q',1,Ref)/1000
            s_evap=np.r_[Cycle.Evaporator.sin_r/1000., sV,Cycle.Evaporator.sout_r/1000.]
            T_evap=np.r_[Cycle.Evaporator.Tin_r, Cycle.Tdew_evap,Cycle.Evaporator.Tout_r]
            self.axes.plot(s_evap,T_evap,'b')
            self.axes.plot(s_evap[0],T_evap[0],'bo')
            self.axes.text(s_evap[0],T_evap[0],'4$\quad\quad$',ha='right',va='bottom')
            
            self.axes.plot([s_evap[2],s_evap[0]],[Cycle.Evaporator.Tin_a,Cycle.Evaporator.Tout_a],'r')
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
                    self.axes.plot(s_IHX[0],T_IHX[0],'bo')
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
                    self.axes.plot(s_IHX[3],T_IHX[3],'bo')
            else:
                raise ValueError('Secondary loop system must have a coaxial or PHE heat exchanger')
        
            self.axes.plot(s_IHX,T_IHX,'b')
            
            
            self.axes.plot(s_XV,T_XV,'b')            
            self.axes.plot([s_IHX[2],s_IHX[0]],[Cycle.CoolingCoil.Tin_a,Cycle.CoolingCoil.Tout_a],'r')
            self.axes.text(0.5*s_IHX[2]+0.5*s_IHX[0],Cycle.CoolingCoil.Tin_a,'Indoor Air',backgroundcolor='w',ha='center',va='center')
            self.axes.text(0.5*s_IHXg[0]+0.5*s_IHXg[1],np.mean(T_IHXg),'Glycol IHX',backgroundcolor='w',ha='center',va='center')            
            self.axes.plot(s_IHXg,T_IHXg,'g-.')
            
        elif Cycle.CycleType=="DX":
            sL=PropsSI('S','T',Cycle.Tdew_evap,'Q',0,Ref)/1000
            sV=PropsSI('S','T',Cycle.Tdew_evap,'Q',1,Ref)/1000
            s_evap=np.r_[Cycle.Evaporator.sin_r/1000., sV,Cycle.Evaporator.sout_r/1000.]
            T_evap=np.r_[Cycle.Evaporator.Tsat_r, Cycle.Tdew_evap,Cycle.Evaporator.Tout_r]
            self.axes.plot(s_evap,T_evap,'b')
            self.axes.plot(s_evap[0],T_evap[0],'bo')
            self.axes.text(s_evap[0],T_evap[0],'4$\quad\quad$',ha='right',va='top')
            
            s_XV=np.r_[Cycle.Evaporator.sin_r/1000.,Cycle.Condenser.sout_r/1000.]
            T_XV=np.r_[Cycle.Evaporator.Tsat_r,Cycle.Condenser.Tout_r]
            self.axes.plot(s_XV,T_XV,'b')
            
            self.axes.plot([s_evap[2],s_evap[0]],[Cycle.Evaporator.Tin_a,Cycle.Evaporator.Tout_a],'r')
            self.axes.text(0.5*s_evap[2]+0.5*s_evap[0],Cycle.Evaporator.Tin_a,'Indoor Air',backgroundcolor='w',ha='center',va='center')
            
        plt.show()
        
    def PHOverlay(self,Cycle,**kwargs):
        Ref=Cycle.Ref
        Tmin=PropsSI('Tmin',Ref)
        Tmax=PropsSI('Tmax',Ref)  

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
    
        self.axes.plot(hsatL,psatL,'k')
        self.axes.plot(hsatV,psatV,'k')
        
        self.axes.set_xlabel('Enthalpy [kJ/kg]')
        self.axes.set_ylabel('Pressure [kPa]')
        
        if Cycle.CycleType=='1INJ_Econ' or Cycle.CycleType=='1INJ_FlashTank':
            h_comp=np.r_[Cycle.Compressor.hin_r/1000.,Cycle.Compressor.h_31/1000,Cycle.Compressor.hinj_r/1000.,Cycle.Compressor.hout_r/1000.]
            p_comp=np.r_[Cycle.Compressor.pin_r/1000,Cycle.Compressor.pinj_r/1000,Cycle.Compressor.pinj_r/1000,Cycle.Compressor.pout_r/1000]
            self.axes.plot(h_comp,p_comp,'b')
            self.axes.plot(h_comp[0],p_comp[0],'bo')
            self.axes.plot(h_comp[1],p_comp[1],'bo')
            self.axes.plot(h_comp[2],p_comp[2],'bo')
            self.axes.plot(h_comp[3],p_comp[3],'bo')
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
            self.axes.plot(h_comp,p_comp,'b')
            self.axes.plot(h_comp[0],p_comp[0],'bo')
            self.axes.plot(h_comp[1],p_comp[1],'bo')
            self.axes.text(h_comp[0],p_comp[0],' 1',ha='left',va='top')
            self.axes.text(h_comp[1],p_comp[1],' 2',ha='left',va='bottom')

        if Cycle.CycleType=='Secondary' and Cycle.Mode=='HP':
            h_evap=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Evaporator.hout_r/1000]
            p_evap=np.r_[Cycle.Evaporator.psat_r/1000, Cycle.Evaporator.psat_r/1000]
            self.axes.plot(h_evap,p_evap,'b')
            self.axes.plot(h_evap[0],p_evap[0],'bo')
            self.axes.text(h_evap[0],p_evap[0],'4$\quad\quad$',ha='right',va='bottom')
        elif Cycle.CycleType=='1INJ_Econ' or Cycle.CycleType=='1INJ_FlashTank':
            h_cond=np.r_[Cycle.Condenser.hin_r/1000.,Cycle.Condenser.hout_r/1000.]
            p_cond=np.r_[Cycle.Compressor.pout_r/1000, Cycle.Compressor.pout_r/1000]
            self.axes.plot(h_cond,p_cond,'b')
            self.axes.plot(h_cond[1],p_cond[1],'bo')
            self.axes.text(h_cond[1],p_cond[1],'3$\quad\quad$',ha='right',va='bottom')
        else:
            h_cond=np.r_[Cycle.Condenser.hin_r/1000.,Cycle.Condenser.hout_r/1000.]
            p_cond=np.r_[Cycle.Condenser.psat_r/1000, Cycle.Condenser.psat_r/1000]
            self.axes.plot(h_cond,p_cond,'b')
            self.axes.plot(h_cond[1],p_cond[1],'bo')
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
                    self.axes.plot(h_IHX[0],p_IHX[0],'bo')
                    self.axes.text(h_IHX[0],p_IHX[0],'4$\quad\quad$',ha='right',va='top')
                else:
                    h_IHX=np.r_[IHX.hin_h/1000., IHX.hout_h/1000.]
                    p_IHX=np.r_[IHX.pin_h/1000, IHX.pin_h/1000]
                    h_XV=np.r_[Cycle.Evaporator.hin_r/1000.,IHX.hout_h/1000.]
                    p_XV=np.r_[Cycle.Evaporator.psat_r/1000,IHX.pin_h/1000]
                    self.axes.plot(h_IHX[1],p_IHX[1],'bo')
                    self.axes.text(h_IHX[1],p_IHX[1],'3$\quad\quad$',ha='right',va='top')
            else:
                raise ValueError('Secondary loop system must have a coaxial or PHE heat exchanger')
            
            self.axes.plot(h_IHX,p_IHX,'b')
            
            self.axes.plot(h_XV,p_XV,'b')
            
        elif Cycle.CycleType=="DX":
            h_evap=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Evaporator.hout_r/1000.]
            p_evap=np.r_[Cycle.Evaporator.psat_r/1000,Cycle.Evaporator.psat_r/1000]
            self.axes.plot(h_evap,p_evap,'b')
            self.axes.plot(h_evap[0],p_evap[0],'bo')
            self.axes.text(h_evap[0],p_evap[0],'4$\quad\quad$',ha='right',va='top')
            
            h_XV=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Condenser.hout_r/1000.]
            p_XV=np.r_[Cycle.Evaporator.psat_r/1000,Cycle.Condenser.psat_r/1000]
            self.axes.plot(h_XV,p_XV,'b')
            
        elif Cycle.CycleType=='1INJ_Econ':  
            h_evap=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Evaporator.hout_r/1000.]
            p_evap=np.r_[Cycle.Evaporator.psat_r/1000,Cycle.Evaporator.psat_r/1000]
            self.axes.plot(h_evap,p_evap,'b')
            self.axes.plot(h_evap[0],p_evap[0],'bo')
            self.axes.text(h_evap[0],p_evap[0],'7$\quad\quad$',ha='right',va='top')

            h_PHEHX = np.r_[Cycle.PHEHX.hout_h/1000.,Cycle.Condenser.hout_r/1000.]
            p_PHEHX = np.r_[Cycle.Compressor.pout_r/1000,Cycle.Compressor.pout_r/1000]
            self.axes.plot(h_PHEHX,p_PHEHX,'b')
            self.axes.plot(h_PHEHX[0],p_PHEHX[0],'bo')
            self.axes.text(h_PHEHX[0],p_PHEHX[0],'4$\quad\quad$',ha='right',va='top')
            h_XV=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.PHEHX.hout_h/1000.]
            p_XV=np.r_[Cycle.Evaporator.psat_r/1000,(Cycle.PHEHX.pin_h-Cycle.PHEHX.DP_h)/1000.]
            self.axes.plot(h_XV,p_XV,'b')   

            h_XV2=np.r_[Cycle.Condenser.hout_r/1000.,Cycle.Condenser.hout_r/1000.]
            p_XV2=np.r_[Cycle.Condenser.psat_r/1000.,Cycle.PHEHX.pin_c/1000.]
            #self.axes.plot(h_XV[1],p_XV[1],'bo')
            #self.axes.plot(h_XV2,p_XV2,'b')                  
            
            h_1INJ=np.r_[Cycle.PHEHX.hout_h/1000,Cycle.PHEHX.hout_c/1000.]
            p_1INJ=np.r_[Cycle.PHEHX.pin_c/1000,(Cycle.PHEHX.pin_c-Cycle.PHEHX.DP_c)/1000.]
            self.axes.plot(h_1INJ,p_1INJ,'b')
            self.axes.plot(h_1INJ[0],p_1INJ[0],'bo')
            self.axes.text(h_1INJ[0],p_1INJ[0],'5$\quad\quad$',ha='right',va='top')

        elif Cycle.CycleType=='1INJ_FlashTank':  

            h_evap=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.Evaporator.hout_r/1000.]
            p_evap=np.r_[Cycle.Evaporator.psat_r/1000,Cycle.Evaporator.psat_r/1000]
            self.axes.plot(h_evap,p_evap,'b')
            self.axes.plot(h_evap[0],p_evap[0],'bo')
            self.axes.text(h_evap[0],p_evap[0],'9$\quad\quad$',ha='right',va='top')

            h_XV=np.r_[Cycle.FlashTank.hin/1000.,Cycle.Condenser.hout_r/1000.]
            p_XV=np.r_[Cycle.FlashTank.pin/1000,Cycle.Compressor.pout_r/1000]
            self.axes.plot(h_XV,p_XV,'b') 

            h_FlashTank_liquid = np.r_[Cycle.FlashTank.hin/1000.,Cycle.FlashTank.hout/1000.]
            p_FlashTank_liquid = np.r_[Cycle.FlashTank.pin/1000,Cycle.FlashTank.pin/1000]
            self.axes.plot(h_FlashTank_liquid,p_FlashTank_liquid,'b')
            self.axes.plot(h_FlashTank_liquid[0],p_FlashTank_liquid[0],'bo')
            self.axes.text(h_FlashTank_liquid[0],p_FlashTank_liquid[0],'6$\quad\quad$',ha='right',va='top')
            
            
            h_1INJ=np.r_[Cycle.FlashTank.hin/1000,Cycle.Compressor.hinj_r/1000.]
            p_1INJ=np.r_[Cycle.FlashTank.pin/1000,(Cycle.Compressor.pinj_r)/1000.]
            self.axes.plot(h_1INJ,p_1INJ,'b')
            self.axes.plot(h_1INJ[0],p_1INJ[0],'ro')
            self.axes.text(h_1INJ[0],p_1INJ[0],'8$\quad\quad$',ha='right',va='top')

            h_XV2=np.r_[Cycle.FlashTank.hout/1000.,Cycle.Evaporator.hin_r/1000.]
            p_XV2=np.r_[Cycle.FlashTank.pin/1000.,Cycle.Evaporator.psat_r/1000.]
            self.axes.plot(h_XV[1],p_XV[1],'bo')
            self.axes.plot(h_XV2,p_XV2,'b') 
            self.axes.plot(h_XV2[0],p_XV2[0],'bo')
            self.axes.text(h_XV2[0],p_XV2[0],'8$\quad\quad$',ha='right',va='top')

        plt.show()


