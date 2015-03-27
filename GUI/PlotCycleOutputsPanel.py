#Developer Note: Detailed mode has been deprecated due to development time constraints
#Can be re-enabled, but some debugging remains

#Imports
import sys,os
sys.path.append(os.path.join('..','Props'))
import wx
from matplotlib.figure import Figure
import numpy as np
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.patches import FancyArrowPatch


class PlotCycleOutputsPanel(wx.Panel):
    def __init__(self, *args, **kwds):
        # begin wxGlade: PlotCycleOutputsPanel.__init__
        kwds["style"] = wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, *args, **kwds)

        self.__set_properties()
        self.__do_layout()
        # end wxGlade
        self.figure=Figure(figsize=(7,5),dpi=100,facecolor='w')
        self.axes=self.figure.add_subplot(111)
        self.axes=[0,0,1,1]
        self.canvas=FigureCanvas(self,wx.ID_ANY,self.figure)

    def __set_properties(self):
        # begin wxGlade: PlotCycleOutputsPanel.__set_properties
        pass
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: PlotCycleOutputsPanel.__do_layout
        pass
        # end wxGlade
        
    def box(self,H,W,x,y):
        X=np.r_[x+W/2.0,x+W/2.0,x-W/2.0,x-W/2.0,x+W/2.0]
        Y=np.r_[y-H/2.0,y+H/2.0,y+H/2.0,y-H/2.0,y-H/2.0]
        return (X,Y)
    
    def TwoTriangles(self,H,W,x,y):
        X=np.r_[x-W/2.0,x+W/2.0,x+W/2.0,x-W/2.0,x-W/2.0]
        Y=np.r_[y-H/2.0,y+H/2.0,y-H/2.0,y+H/2.0,y-H/2.0]
        return (X,Y)
    
    
    def circle(self,x0,y0,r,N):
        t=np.linspace(0,2*np.pi,N)
        x=x0+r*np.cos(t)
        y=y0+r*np.sin(t)
        return (x,y)
    
    def bendyArrow(self,x1,y1,x2,y2,dir1,dir2):
        x3=x2
        y3=y2
        if dir1 == 'up':
            y2=y3
            x2=x1
        elif dir1 =='down':
            y2=y3
            x2=x1
        elif dir1 == 'left' or dir1 == 'right':
            y2=y1
            x2=x3
        xL=np.r_[x1,x2]
        yL=np.r_[y1,y2]
        xA=np.r_[x2,x3]
        yA=np.r_[y2,y3]
        return (xL,yL,xA,yA)
    
    def TwoBendArrow(self,x1,y1,x2,y2,orientation,parallelCoord):
        x4=x2
        y4=y2
        if orientation == 'h':
            y2=parallelCoord
            y3=y2
            x2=x1
            x3=x4
        elif orientation =='v':
            x2=parallelCoord
            x3=x2
            y2=y1
            y3=y4
        xL=np.r_[x1,x2,x3]
        yL=np.r_[y1,y2,y3]
        xA=np.r_[x3,x4]
        yA=np.r_[y3,y4]
        return (xL,yL,xA,yA)
        
    def doPlot(self,Cycle,**kwargs):
        if 'dpi' in kwargs:
            dpi=kwargs['dpi']
        else:
            dpi=120
            
        if 'detailed' in kwargs:
            detailed=kwargs['detailed']
        else:
            detailed='no'
            
        if 'figsize' in kwargs:
            self.figure=Figure(figsize=kwargs['figsize'],dpi=dpi,facecolor='w')
        
        self.axes=self.figure.add_subplot(111)
        self.axes=[0,0,1,1]
        self.canvas=FigureCanvas(self,wx.ID_ANY,self.figure)
            
        self.figure.clf()
        self.axes=self.figure.add_subplot(111)
        self.canvas=FigureCanvas(self,wx.ID_ANY,self.figure)
        
        xIndoor=0.0
        yIndoor=0.0
        HIndoor=4
        WIndoor=1.5
        
        xPump=3
        yPump=4
        rPump=0.75
        
        xIHX=6.0
        yIHX=0.0
        HIHX=4
        WIHX=1.5
        
        xComp=9
        yComp=4
        rComp=0.75
        
        xCondenser=12.0
        yCondenser=0.0
        HCondenser=4
        WCondenser=1.5
        
        xXV=9.0
        yXV=-5.0
        WXV=0.5
        HXV=0.5
        
        #Indoor Coil to Pump/Compressor
        if Cycle.CycleType=='Secondary':
            if Cycle.IHXType=='PHE':
                Cycle.IHX=Cycle.PHEIHX
            
            (xIndoor_,yIndoor_)=self.box(HIndoor,WIndoor,xIndoor,yIndoor)
            self.axes.plot(xIndoor_,yIndoor_,'k',lw=2)
            self.axes.fill(xIndoor_,yIndoor_,'w',hatch='x',lw=0.2)
            if detailed=='no':
                if Cycle.Mode=='HP':
                    label='Heating Coil\nQ$_{net}$: %0.1f W' %(Cycle.CoolingCoil.Q-Cycle.CoolingCoil.Fins.Air.FanPower)
                    self.axes.text(xIndoor+WIndoor/2.0+0.1,yIndoor,label,ha='left',va='center',size=10)
                else:   
                    label='Cooling Coil\nQ$_{net}$: %0.1f W' %(Cycle.CoolingCoil.Q-Cycle.CoolingCoil.Fins.Air.FanPower)
                    self.axes.text(xIndoor+WIndoor/2.0+0.1,yIndoor,label,ha='left',va='center',size=10)
            else:
                label='Cooling Coil\nQ$_{net}$: %0.1f W\n$T_{out,a}$=%0.2f K' %(Cycle.CoolingCoil.Q-Cycle.CoolingCoil.Fins.Air.FanPower,Cycle.CoolingCoil.Tout_a)
                self.axes.text(xIndoor+WIndoor/2.0+0.1,yIndoor,label,ha='left',va='center',size=10)
                label='$\dot V_{a}$=%0.3f $m^3/s$\n$T_{in,a}$ = %0.2f K \n $\phi_{in}$ = %0.2f %%' %(Cycle.CoolingCoil.Fins.Air.Vdot_ha,Cycle.CoolingCoil.Tin_a,Cycle.CoolingCoil.Fins.Air.RH*100.)
                self.axes.text(xIndoor-WIndoor/2.0-0.1,yIndoor,label,ha='right',va='center',size=10)
                
            (xL,yL,xA,yA)=self.bendyArrow(xIndoor,HIndoor/2.0,xPump-rPump,yPump,'up','down')
            self.axes.plot(xL,yL,'k',lw=1)
            self.axes.plot(xA,yA,'k',lw=1)
            if detailed=='yes':
                self.axes.plot(xIndoor,yPump,'ko',ms=8)
                label='T=%0.1f K'%(Cycle.CoolingCoil.Tout_g)
                self.axes.text(xIndoor+0.1,yPump+0.1,label,size=10,ha='left',va='bottom')
    
            #Pump
            (xPump_,yPump_)=self.circle(xPump,yPump,rPump,50)
            self.axes.plot(xPump_,yPump_,'k',lw=2)
            label='P: %0.1f W' %(Cycle.Pump.W)
            self.axes.text(xPump,yPump+rPump+0.1,label,ha='center',va='bottom',size=10)
            self.axes.text(xPump,yPump-rPump-0.1,'Pump',ha='center',va='top',size=10)
            
            (xL,yL,xA,yA)=self.bendyArrow(xPump+rPump,yPump,xIHX-WIHX/6.0,HIHX/2.0,'right','down')
            self.axes.plot(xL,yL,'k',lw=1)
            self.axes.plot(xA,yA,'k',lw=1)
            if detailed=='yes':
                self.axes.plot(xIHX-WIHX/6.0,yPump,'ko',ms=8)
                label='T=%0.1f K'%(Cycle.IHX.Tin_g)
                self.axes.text(xIHX-WIHX/6.0-0.1,yPump+0.1,label,size=10,ha='right',va='bottom')
            
            (xIHX_,yIHX_)=self.box(HIHX,WIHX,xIHX,yIHX)
            self.axes.plot(xIHX_,yIHX_,'k',lw=2)
            label='IHX\nQ: %0.1f W' %(Cycle.IHX.Q)
            self.axes.text(xIHX+WIHX/2.0+0.1,yIHX,label,ha='left',va='center',size=10)        
            
            (xL,yL,xA,yA)=self.bendyArrow(xIHX+WIHX/6.0,HIHX/2.0,xComp-rComp,yComp,'up','down')
            self.axes.plot(xL,yL,'k',lw=1)
            self.axes.plot(xA,yA,'k',lw=1)
            
            (xL,yL,xA,yA)=self.bendyArrow(xIHX+WIHX/6.0,-HIHX/2.0,xXV-WXV/2.0,yXV,'up','down')
            self.axes.plot(xL,yL,'k',lw=1)
            self.axes.plot(xA,yA,'k',lw=1)
                        
            (xL,yL,xA,yA)=self.TwoBendArrow(xIHX-WIHX/6.0,-HIHX/2.0,xIndoor,-HIndoor/2.0,'h',yXV)
            self.axes.plot(xL,yL,'k',lw=1)
            self.axes.plot(xA,yA,'k',lw=1)
            if detailed=='yes':
                self.axes.plot(xIHX-WIHX/6.0,yXV,'ko',ms=8)
                label='T=%0.1f K'%(Cycle.IHX.Tout_g)
                self.axes.text(xIHX-WIHX/6.0-0.1,yXV+0.1,label,size=10,ha='right',va='bottom')
                
                self.axes.plot(xIndoor,yXV,'ko',ms=8,)
                label='T=%0.1f K'%(Cycle.CoolingCoil.Tin_g)
                self.axes.text(xIndoor+0.1,yXV+0.1,label,size=10,ha='left',va='bottom')
                
                label='$\dot m_{fluid}$ = %0.4f kg/s' %(Cycle.Pump.mdot_g)
                self.axes.text((xIndoor+xIHX)/2.0,yXV/2.0,label,size=10,ha='center',va='center')
                
                label='$\dot m_{ref}$ = %0.4f kg/s' %(Cycle.Compressor.mdot_r)
                self.axes.text((xCondenser+xIHX)/2.0,yXV/2.0,label,size=10,ha='center',va='center')
                
                self.axes.plot(xIHX+WIHX/6.0,yXV,'ko',ms=8)
                label='x=%0.3f\nT=%0.1f K'%(Cycle.IHXxin_r,Cycle.IHXTin_r)
                self.axes.text(xIHX+WIHX/6.0+0.1,yXV+0.1,label,size=10,ha='left',va='bottom')
                
                self.axes.plot(xIHX+WIHX/6.0,yComp,'ko',ms=8)
                label='T=%0.1f K'%(Cycle.IHX.Tout_r)
                self.axes.text(xIHX+WIHX/6.0+0.1,yComp+0.1,label,size=10,ha='left',va='bottom')
            
            if Cycle.Mode=='AC':
                self.axes.add_patch(FancyArrowPatch((xIHX+WIHX/6.0,-HIHX/2.0),(xIHX+WIHX/6.0,HIHX/2.0),arrowstyle='-|>',mutation_scale=20,fc='b',ec='b',lw=2))
                self.axes.add_patch(FancyArrowPatch((xIHX-WIHX/6.0,HIHX/2.0),(xIHX-WIHX/6.0,-HIHX/2.0),arrowstyle='-|>',mutation_scale=20,fc='b',ec='b',lw=2))
            else:
                self.axes.add_patch(FancyArrowPatch((xIHX+WIHX/6.0,+HIHX/2.0),(xIHX+WIHX/6.0,-HIHX/2.0),arrowstyle='-|>',mutation_scale=20,fc='b',ec='b',lw=2))
                self.axes.add_patch(FancyArrowPatch((xIHX-WIHX/6.0,-HIHX/2.0),(xIHX-WIHX/6.0,+HIHX/2.0),arrowstyle='-|>',mutation_scale=20,fc='b',ec='b',lw=2))
            
        else:
            (xIndoor_,yIndoor_)=self.box(HIndoor,WIndoor,xIndoor,yIndoor)
            self.axes.plot(xIndoor_,yIndoor_,'k',lw=2)
            self.axes.fill(xIndoor_,yIndoor_,'w',hatch='x',lw=0.2)
            if detailed=='no':
                label='DX Coil\nQ$_{net}$: %0.1f W' %(Cycle.Evaporator.Q-Cycle.Evaporator.Fins.Air.FanPower)
                self.axes.text(xIndoor+WIndoor/2.0+0.1,yIndoor,label,ha='left',va='center',size=10)
            else:
                label='DX Coil\nQ$_{net}$: %0.1f W\n$T_{out,a}$=%0.2f K' %(Cycle.Evaporator.Q-Cycle.Evaporator.Fins.Air.FanPower,Cycle.Evaporator.Tout_a)
                self.axes.text(xIndoor+WIndoor/2.0+0.1,yIndoor,label,ha='left',va='center',size=10)
                label='$\dot V_{a}$=%0.3f $m^3/s$\n$T_{in,a}$ = %0.2f K \n $\phi_{in}$ = %0.2f %%' %(Cycle.Evaporator.Fins.Air.Vdot_ha,Cycle.Evaporator.Tin_a,Cycle.Evaporator.Fins.Air.RH*100.)
                self.axes.text(xIndoor-WIndoor/2.0-0.1,yIndoor,label,ha='right',va='center',size=10)
        
            (xL,yL,xA,yA)=self.bendyArrow(xIndoor,HIndoor/2.0,xComp-rComp,yComp,'up','down')
            self.axes.plot(xL,yL,'k',lw=1)
            self.axes.plot(xA,yA,'k',lw=1)
            
            (xL,yL,xA,yA)=self.bendyArrow(xXV-WXV/2.0,yXV,xIndoor,-HIndoor/2.0,'left','down')
            self.axes.plot(xL,yL,'k',lw=1)
            self.axes.plot(xA,yA,'k',lw=1)
            
            if detailed=='yes':
                label='$\dot m_{ref}$ = %0.4f kg/s' %(Cycle.Compressor.mdot_r)
                self.axes.text((xCondenser+xIndoor)/2.0,yXV/2.0,label,size=10,ha='center',va='center')
                
                self.axes.plot(xIndoor,yXV,'ko',ms=8)
                label='x=%0.3f\nT=%0.1f K'%(Cycle.Evaporator.xin_r,Cycle.Evaporator.Tsat_r)
                self.axes.text(xIndoor+0.1,yXV+0.1,label,size=10,ha='left',va='bottom')
                
                self.axes.plot(xIndoor,yComp,'ko',ms=8)
                label='T=%0.1f K'%(Cycle.Evaporator.Tout_r)
                self.axes.text(xIndoor+0.1,yComp+0.1,label,size=10,ha='left',va='bottom')
            
        (xCompressor,yCompressor)=self.circle(xComp,yComp,rComp,50)
        self.axes.plot(xCompressor,yCompressor,'k',lw=2)
        label='Q: %0.2f W\nP: %0.2f W' %(-Cycle.Compressor.W*Cycle.Compressor.fp,Cycle.Compressor.W)
        self.axes.text(xComp,yComp+rComp+0.1,label,ha='center',va='bottom',size=10)
        self.axes.text(xComp,yComp-rComp-0.1,'Comp',ha='center',va='top',size=10)
        
        (xL,yL,xA,yA)=self.bendyArrow(xComp+rComp,yComp,xCondenser,HCondenser/2.0,'right','down')
        self.axes.plot(xL,yL,'k',lw=1)
        self.axes.plot(xA,yA,'k',lw=1)
        if detailed=='yes':
            self.axes.plot(xCondenser,yComp,'ko',ms=8)
            label='T=%0.1f K'%(Cycle.Compressor.Tout_r)
            self.axes.text(xCondenser-0.1,yComp+0.1,label,size=10,ha='right',va='bottom')
        
        (xCondenser_,yCondenser_)=self.box(HCondenser,WCondenser,xCondenser,yCondenser)
        self.axes.plot(xCondenser_,yCondenser_,'k',lw=2)
        self.axes.fill(xCondenser_,yCondenser_,'w',hatch='+',lw=0.2)
            
        if detailed=='no':
            if Cycle.CycleType=='Secondary' and Cycle.Mode=='HP':
                label='Evaporator\nQ: %0.1f W' %(Cycle.Evaporator.Q)
                self.axes.text(xCondenser+WCondenser/2.0+0.1,yCondenser,label,ha='left',va='center',size=10)
            else:
                label='Condenser\nQ: %0.1f W' %(Cycle.Condenser.Q)
                self.axes.text(xCondenser+WCondenser/2.0+0.1,yCondenser,label,ha='left',va='center',size=10)
        else:
            label='Condenser\nQ: %0.1f W' %(Cycle.Condenser.Q)
            self.axes.text(xCondenser-WCondenser/2.0-0.1,yCondenser,label,ha='right',va='center',size=10)
            
            label='$\dot V_{a}$=%0.3f $m^3/s$\n$T_{in,a}$ = %0.2f K \n $\phi_{in}$ = %0.2f %%' %(Cycle.Condenser.Fins.Air.Vdot_ha,Cycle.Condenser.Tin_a,Cycle.Condenser.Fins.Air.RH*100.)
            self.axes.text(xCondenser+WCondenser/2.0+0.1,yCondenser,label,ha='left',va='center',size=10)
                
            self.axes.plot(xCondenser,yComp,'ko',ms=8)
            label='T=%0.1f K'%(Cycle.Compressor.Tout_r)
            self.axes.text(xCondenser-0.1,yComp+0.1,label,size=10,ha='right',va='bottom')
            
            self.axes.plot(xCondenser,yXV,'ko',ms=8)
            label='T=%0.1f K'%(Cycle.Condenser.Tout_r)
            self.axes.text(xCondenser-0.1,yXV+0.1,label,size=10,ha='right',va='bottom')
        
        #ExpansionValve
        (xXV_,yXV_)=self.TwoTriangles(HXV,WXV,xXV,yXV)
        self.axes.plot(xXV_,yXV_,'k',lw=2)
        self.axes.text(xXV,yXV+HXV/2.0+0.1,'XV',ha='center',va='bottom',size=10)
            
        (xL,yL,xA,yA)=self.bendyArrow(xCondenser,-HCondenser/2.0,xXV+WXV/2.0,yXV,'down','left')
        self.axes.plot(xL,yL,'k',lw=1)
        self.axes.plot(xA,yA,'k',lw=1)
        
        self.axes.axis('off')
        self.axes.axis('equal')

# end of class PlotCycleOutputsPanel


