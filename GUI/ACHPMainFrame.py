
import wx,os,sys
sys.path.append('..')

from LoadGUI import LoadGUI,BindEvents
from PlotCycleOutputsPanel import PlotCycleOutputsPanel
from MPLPanel import MPLPanel
import itertools
from ACHPFileIO import GUI2DXCycleInputs,CycleOutputs2GUI,GUI2SecondaryCycleInputs
import numpy as np
from CoolProp.CoolProp import Props
from CircuitsClass import CircuitsClass
from TubeSelectClass import TubeSelectClass
from ARIInfoClass import ARIInfoClass
from SchematicFrameClass import SchematicFrameClass

# Import Psyco if available
try:
    import psyco
    psyco.profile()
except ImportError:
    pass

global ChangeSystemType

class ACHPMainFrame(wx.Frame):
    
    def __init__(self,*args,**kwargs):
        super(ACHPMainFrame, self).__init__(None, title='ACHP', 
            size=(800, 600))
        # Run external code in LoadGUI module to 
        # form all the elements of GUI
        LoadGUI(self)
        BindEvents(self)
        self.ChangeSolverMethod()
        self.ChangeSystemType()
        
        #If the name of a config file is passed as a command line argument, load it
        # Must be in the format -O CONFIGGILE.cfg or -R CONFIGFILE.cfg
        if len(sys.argv)>1:
            if '-R' in sys.argv:
                #Read in the config file passed in, and also run the file
                self.AutoRun=True
                i=sys.argv.index('-R')
                self.ReadConfigFile(file=sys.argv[i+1])
                self.RunCode()
                self.Close()
            elif '-O' in sys.argv:
                #Read in the config file passed in, but don't run the file
                self.AutoRun=True
                i=sys.argv.index('-O')
                self.ReadConfigFile(file=sys.argv[i+1])
            else:
                self.AutoRun=False
        
    def RunCode(self, event=None): 
        
        
            
        #Check what type of solver is to be used
        if self.optCycleCharge.GetValue()==True:
            ImposedVariable='Charge'
        else:
            ImposedVariable='Subcooling'
            
        # Check whether a Direct-Expansion system (DX) or secondary loop system 
        # (Secondary) and then load up the run parameters for the appropriate 
        # mode
        if self.optCycleDX.GetValue()==True:
            Cycle=GUI2DXCycleInputs(self)
            Cycle.ImposedVariable=ImposedVariable
            Cycle.CycleType='DX'
        else:
            Cycle=GUI2SecondaryCycleInputs(self)
            Cycle.ImposedVariable=ImposedVariable
            Cycle.CycleType='Secondary'
        
        # Check that if you are using a REFPROP refrigerant that 
        # REFPROP is installed in a good location
        # Important!! REFPROP only works on windows
        if not os.path.exists('c:/Program Files/REFPROP/refprop.dll') and Cycle.Ref.startswith('REFPROP'):
            dlgString='Error: refprop.dll not found in folder\nc:/Program Files/REFPROP/refprop.dll'
            dlg = wx.MessageDialog(self,dlgString,'Run Error',wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            
        if self.radCycleMode.GetStringSelection()=='Cooling Mode':
            Cycle.Mode='AC'
        else: 
            Cycle.Mode='HP'
            
##         if Cycle.Mode=='HP' and Cycle.CycleType=='Secondary' and not os.path.exists('unlock.xyz'):
##             dlg = wx.MessageDialog(self,"Sorry but secondary loops in heating mode not yet supported")
##             dlg.ShowModal()
##             return
            
        #Actually do the solving, either with a parametric study or fixed point
        if self.radSolverSelect.GetStringSelection()=='Parametric Study':
            #It is a parametric solve
            self.ParametricStudy(Cycle)
        else:
            # Cycle here could either be a DXSolver Class or a 
            # SecondaryCycleClass. Both have a PreconditionedSolve() method 
            Cycle.PreconditionedSolve()
            
            #Update all the fields in the GUI
            CycleOutputs2GUI(self,Cycle)
            
            #Make T-s and p-h plots
            kwargs={
                'Tmin': 260,
                'Tmax': 350,
                'Ref' : Cycle.Ref
            }
            self.pltTs.TSOverlay(Cycle,**kwargs)
            kwargs={
                'Tmin': 260,
                'Tmax': 350,
                'Ref' : Cycle.Ref
            }
            self.pltPh.PHOverlay(Cycle,**kwargs)
        
            #Plot the overview
            self.pnlSchematic.doPlot(Cycle)
            
            #Keep a copy of Cycle
            self.Cycle=Cycle
        
    def ParametricStudy(self,Cycle,**kwargs):
        
        NVariables=4
        VariableList=[]
        ValueList=[]
        for i in range(1,NVariables+1):
            #get checkbox
            chkbox=getattr(self,'chkParaVariable'+str(i))
            #See if checkbox is checked
            if chkbox.GetValue()==True:
                #get combobox
                cmbbox=getattr(self,'cmbParaVariable'+str(i))
                #Add variable
                VariableList.append(str(cmbbox.GetValue()))
                Min=getattr(self,'txtParaVariable'+str(i)+'Min').GetValue()
                Max=getattr(self,'txtParaVariable'+str(i)+'Max').GetValue()
                N=getattr(self,'txtParaVariable'+str(i)+'N').GetValue()
                if N=='L':
                    #Use comma-separated list of values stored in Min
                    values=Min.split(',')
                    #Convert list elements to floats
                    values=map(float,values)
                    #Append to list of values
                    ValueList.append(values)
                else:
                    values=np.linspace(float(Min),float(Max),int(N))
                    ValueList.append(values)
                    
        #Create a permutation of all the possibilities using the function product
        # An example:
        #  >> import itertools
        #  >> A=['r','g','b']
        #  >> B=['x','y']
        #  >> print list(itertools.product(A,B))
        # should yield 
        # [('r', 'x'), ('r', 'y'), ('r', 'z'), ('g', 'x'), ('g', 'y'), 
        #    ('g', 'z'), ('b', 'x'), ('b', 'y'), ('b', 'z')]
        #
        # Also see http://stackoverflow.com/questions/798854/all-combinations-of-a-list-of-lists
        # Doing *ValueList unpacks the list of numpy-like arrays and converts
        # product([values1, values2, values3]) to product(values1,values2,values3)
        #
        runList=list(itertools.product(*ValueList)) 
        
        # Check if file already exists, don't do anything if the user doesn't
        # want to overwrite it
        if os.path.exists(str(self.txtParaPath.GetValue)):
            dlg = wx.MessageDialog(self, "Parametric Ouput File exists. Append?", "Overwrite confirmation", wx.YES_NO|wx.NO_DEFAULT|wx.ICON_QUESTION)
            if not (dlg.ShowModal() == wx.ID_YES):
                return
        
        if len(runList[0]) == 0:
            dlgString='Parametric mode selected but no parametric variables selected'
            dlg = wx.MessageDialog(self,dlgString,'Run Error',wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            return
        
        #Loop over the runList
        for i in range(len(runList)):
            #For each run, set all the variables that are needed
            for j in range(len(VariableList)):
                #Get a string representation of the variable name
                VariableString=self.ParamVariables[VariableList[j]]
                #Drill down through the variable name, ensuring that it is valid
                # -----------------------
                #Split to list of fields
                fields=VariableString.strip().split('.')
                item=Cycle
                if len(fields)==0:
                    raise ValueError()
                if fields[0]=='Cycle':
                    fields.pop(0)
                while len(fields)>1:
                    item=getattr(item,fields[0])
                    fields.pop(0)
                # The variable item now points to the owner class of variable
                # Update the value
                setattr(item,fields[0],runList[i][j])
                self.ParamValues[VariableList[j]]=runList[i][j]
                
            #Now solve the given point
            Cycle.PreconditionedSolve()
            #Save it to file
            if i==0:
                HeaderOn=True
            else:
                HeaderOn=False
            self.SaveRun(Cycle,path=Cycle.ParaPath,type='para',Header=HeaderOn)

    def ChangeSolverMethod(self, event=None):
        ParamItems=[
            self.chkParaVariable1,self.cmbParaVariable1,
            self.txtParaVariable1Min,self.txtParaVariable1Max,
            self.txtParaVariable1N,
            self.chkParaVariable2,self.cmbParaVariable2,
            self.txtParaVariable2Min,self.txtParaVariable2Max,
            self.txtParaVariable2N,
            self.chkParaVariable3,self.cmbParaVariable3,
            self.txtParaVariable3Min,self.txtParaVariable3Max,
            self.txtParaVariable3N,
            self.chkParaVariable4,self.cmbParaVariable4,
            self.txtParaVariable4Min,self.txtParaVariable4Max,
            self.txtParaVariable4N]
        
        for item in ParamItems:
            if self.radSolverSelect.GetStringSelection()=='Parametric Study':
                item.Enable()
            else:
                item.Disable()

    def ChangeCompressorModel(self, event=None): 
        CompMapItems=[self.lblCompM1, self.lblCompM2, self.lblCompM3,
            self.lblCompM4, self.lblCompM5, self.lblCompM6, self.lblCompM7,
            self.lblCompM8, self.lblCompM9, self.lblCompM10, self.txtCompM1,
            self.txtCompM2, self.txtCompM3, self.txtCompM4, self.txtCompM5,
            self.txtCompM6, self.txtCompM7, self.txtCompM8, self.txtCompM9,
            self.txtCompM10, self.lblCompP1, self.lblCompP2, self.lblCompP3,
            self.lblCompP4, self.lblCompP5, self.lblCompP6, self.lblCompP7,
            self.lblCompP8, self.lblCompP9, self.lblCompP10, self.txtCompP1,
            self.txtCompP2, self.txtCompP3, self.txtCompP4, self.txtCompP5,
            self.txtCompP6, self.txtCompP7, self.txtCompP8, self.txtCompP9,
            self.txtCompP10]
        EfficiencyItems=[self.lblCompIsenEff, self.txtCompIsenEff,
            self.lblCompVolEff, self.txtCompVolEff, self.lblCompVdot,
            self.txtCompVdot]

        for item in CompMapItems:
            if self.optCompMapModel.GetValue() == 1:
                item.Enable()
            else:
                item.Disable()
        for item in EfficiencyItems:
            if self.optCompMapModel.GetValue() == 0:
                item.Enable()
            else:
                item.Disable()
        
            
    def LoadCoeffs(self, event):
        dlg = wx.MessageDialog(self,"Comma separated values with first column the mass flow coeffs, \nsecond column the power coefficients.  Based on saturation temps in degF, \nmass flows in lbm/hr, and power in W")
        if (dlg.ShowModal() == wx.ID_OK):
            dlg = wx.FileDialog(self, "Open CSV file", "comps", "",
                               "CSV files (*.csv)|*.csv|All Files|*.*", wx.OPEN)
            if (dlg.ShowModal() == wx.ID_OK):
                fileName = dlg.GetFilename()
                dirName = dlg.GetDirectory()
                print fileName,dirName
            
                A=np.loadtxt(dirName + os.sep + fileName,delimiter=',')
                
                self.txtCompM1.SetValue("%0.9f " % A[0,0])
                self.txtCompM2.SetValue("%0.9f " % A[1,0])
                self.txtCompM3.SetValue("%0.9f " % A[2,0])
                self.txtCompM4.SetValue("%0.9f " % A[3,0])
                self.txtCompM5.SetValue("%0.9f " % A[4,0])
                self.txtCompM6.SetValue("%0.9f " % A[5,0])
                self.txtCompM7.SetValue("%0.9f " % A[6,0])
                self.txtCompM8.SetValue("%0.9f " % A[7,0])
                self.txtCompM9.SetValue("%0.9f " % A[8,0])
                self.txtCompM10.SetValue("%0.9f " % A[9,0])
                
                self.txtCompP1.SetValue("%0.9f " % A[0,1])
                self.txtCompP2.SetValue("%0.9f " % A[1,1])
                self.txtCompP3.SetValue("%0.9f " % A[2,1])
                self.txtCompP4.SetValue("%0.9f " % A[3,1])
                self.txtCompP5.SetValue("%0.9f " % A[4,1])
                self.txtCompP6.SetValue("%0.9f " % A[5,1])
                self.txtCompP7.SetValue("%0.9f " % A[6,1])
                self.txtCompP8.SetValue("%0.9f " % A[7,1])
                self.txtCompP9.SetValue("%0.9f " % A[8,1])
                self.txtCompP10.SetValue("%0.9f " % A[9,1])
                
                self.txtCompCoeffsFile.SetValue(fileName)
        event.Skip()
        
    def ChangeImposedVariable(self, event=None):
        if self.optCycleSubcooling.GetValue()==True:
            self.txtCycleSubcooling.Enable()
            self.txtCycleCharge.Disable()
        else:
            self.txtCycleSubcooling.Disable()
            self.txtCycleCharge.Enable()
        
    def ChangeSystemType(self, event=None): 
        SecondaryItems=[
            self.PumpIHXInputsPanel,
            self.PumpIHXOutputsPanel,
            self.CoolingCoilOutputsPanel  ]
        
        DXItems=[self.EvaporatorOutputsPanel ]
        
        for item in SecondaryItems:
            if self.optCycleDX.GetValue()==True:
                item.Disable()
            else:
                item.Enable()
        for item in DXItems:
            if self.optCycleDX.GetValue()==True:
                item.Enable()
            else:
                item.Disable()
                
        if self.radCycleMode.GetStringSelection()=='Cooling Mode':
            if self.optCycleDX.GetValue()==True:
                # It is a DX cycle
                self.bmpCycle.SetBitmap(wx.Bitmap(os.path.join('imgs','DX-Cycle-AC.png')))
                self.lblCycleDTsh.SetLabel('Evap outlet superheat [K]')
            else:
                self.bmpCycle.SetBitmap(wx.Bitmap(os.path.join('imgs','Secondary-Cycle-AC.png')))
                self.lblCycleDTsh.SetLabel('IHX outlet superheat [K]')
                self.EvaporatorOutputsPanel.Disable()
                self.CondenserOutputsPanel.Enable()
                 
        elif self.radCycleMode.GetStringSelection()=='Heating Mode':
            if self.optCycleDX.GetValue()==True:
                # It is a DX cycle
                self.bmpCycle.SetBitmap(wx.Bitmap(os.path.join('imgs','DX-Cycle-HP.png')))
                self.lblCycleDTsh.SetLabel('Evap outlet superheat [K]')
                self.EvaporatorOutputsPanel.Enable()
                self.CondenserOutputsPanel.Enable()
            else:
                self.bmpCycle.SetBitmap(wx.Bitmap(os.path.join('imgs','Secondary-Cycle-HP.png')))
                self.lblCycleDTsh.SetLabel('IHX outlet superheat [K]')
                self.EvaporatorOutputsPanel.Enable()
                self.CondenserOutputsPanel.Disable()
                

    def WriteConfigFile(self, event=None, fName=None,OverWrite=False,**kwargs): # wxGlade: MyFrame.<event_handler>
        """
        OverWrite : True/False   [Whether to prompt on overwrite] 
        """
        if fName == None:
            #Ask for a file name
            dlg = wx.FileDialog(self, "Save Config file", "configs", "",
                   "CFG files (*.cfg)|*.cfg|All Files|*.*", wx.SAVE)
            if (dlg.ShowModal() == wx.ID_OK):
                fileName = dlg.GetFilename()
                dirName = dlg.GetDirectory()
                #Join the absolute path to the file together
                fName=os.path.join(dirName,fileName)
        #If they cancelled, path is empty, quit
        if fName==None:
            return
        
        #Check if it exists, confirm overwrite,otherwise quit
        if os.path.exists(fName) and OverWrite==False:
            dlg = wx.MessageDialog(self, "File exists. Overwrite?", "Overwrite confirmation", wx.YES_NO|wx.NO_DEFAULT|wx.ICON_QUESTION)
            
            if (dlg.ShowModal() == wx.ID_NO):
                #Don't do anything and quit function
                return
        
        #Actually write the file
        f=open(fName,'w')
        list=dir(self) #get a listing of all things in this class
        for i in range(len(list)):
            name=str(list[i])
            ## Look at the beginning of the item name
            if name.startswith(('txt','cmb','opt','chk')):
                f.write(' %s : %s\n' %(name,getattr(self,name).GetValue()))
            elif name.startswith('rad'):
                f.write(' %s : %s\n' %(name,getattr(self,name).GetSelection()))
        f.close()
        
    def ReadConfigFile(self, event=None, **kwargs): 
        
        if 'file' not in kwargs:
            dlg = wx.FileDialog(self, "Load Config file", "configs", "",
                                   "CFG files (*.cfg)|*.cfg|All Files|*.*", wx.OPEN)
            if (dlg.ShowModal() == wx.ID_OK):
                fileName = dlg.GetFilename()
                dirName = dlg.GetDirectory()
            else:
                return
            file=dirName+os.sep+fileName
        else:
            file=kwargs['file']
            
        notFound=[]
        f=open(file,'r')
        for line in f:
            (name,value)=str.split(line,':')
            name=str.strip(name,'\n')
            value=str.strip(value,'\n')
            name=str.lstrip(name,' ')
            name=str.rstrip(name,' ')
            value=str.lstrip(value,' ')
            value=str.rstrip(value,' ')
            value=str.rstrip(value)
            try: 
                if name.startswith(('cmb','txt')):
                    getattr(self,name).SetValue(value)
                elif name.startswith('rad'):
                    getattr(self,name).SetSelection(int(value))                    
                elif name.startswith(('opt','chk')):
                    if value.startswith('True'):
                        getattr(self,name).SetValue(True)
                    else:
                        getattr(self,name).SetValue(False)
            except:
                notFound.append(name)
            
        f.close()
        if len(notFound)>0:
            print 'These GUI items were not loaded: ',notFound
        if len(notFound)>0:
            nfString=",".join(["%s" % k for k in notFound])
            dlg = wx.MessageDialog(self,nfString+" not found")
            
        self.ChangeSystemType()
        self.ChangeImposedVariable()
        self.ChangeSolverMethod()
        #event.Skip()

    def openHelp(self, event): # wxGlade: MyFrame.<event_handler>
        os.startfile("help/help.pdf")

    def CoolingCoilFinSelect(self, event): 
        if self.cmbCoolingCoilFink.GetValue()=='Aluminum':
            self.txtCoolingCoilFink.SetValue('237')
        elif self.cmbCoolingCoilFink.GetValue()=='Copper':
            self.txtCoolingCoilFink.SetValue('401')
        
    def ShowCoolingCoilCircuits(self, event): 
        Circuits = CircuitsClass(self, -1, "")
        print float(self.txtCoolingCoilTubesNcircuit.GetValue())
        kwargs={'Ncircuits' : float(self.txtCoolingCoilTubesNcircuit.GetValue()),
                'pl' : float(self.txtCoolingCoilTubesPl.GetValue()),
                'pt' : float(self.txtCoolingCoilTubesPt.GetValue()),
                'Nbank' : float(self.txtCoolingCoilTubesNbank.GetValue()),
                'Ntubes_bank' : float(self.txtCoolingCoilTubesNtubes.GetValue())
                }
        Circuits.pnlCircuits.doPlot(**kwargs)
        wx.GetApp().SetTopWindow(Circuits)
        Circuits.Show()
        
    def ShowCondenserCircuits(self, event): 
        Circuits = CircuitsClass(self, -1, "")
        print float(self.txtCoolingCoilTubesNcircuit.GetValue())
        kwargs={'Ncircuits' : float(self.txtCondenserTubesNcircuit.GetValue()),
                'pl' : float(self.txtCondenserTubesPl.GetValue()),
                'pt' : float(self.txtCondenserTubesPt.GetValue()),
                'Nbank' : float(self.txtCondenserTubesNbank.GetValue()),
                'Ntubes_bank' : float(self.txtCondenserTubesNtubes.GetValue())
                }
        Circuits.pnlCircuits.doPlot(**kwargs)
        wx.GetApp().SetTopWindow(Circuits)
        Circuits.Show()

    def SelectLineSetReturnTube(self, event): 
        TubeSelect = TubeSelectClass(self, -1, "")
        TubeSelect.setTarget(self.txtLineSetOD_return,self.txtLineSetID_return)
        wx.GetApp().SetTopWindow(TubeSelect)
        TubeSelect.Show()

    def SelectLineSetSupplyTube(self, event): 
        TubeSelect = TubeSelectClass(self, -1, "")
        TubeSelect.setTarget(self.txtLineSetOD_supply,self.txtLineSetID_supply)
        wx.GetApp().SetTopWindow(TubeSelect)
        TubeSelect.Show()

    def SelectEvaporatorTube(self, event): 
        TubeSelect = TubeSelectClass(self, -1, "")
        TubeSelect.setTarget(self.txtEvaporatorTubesOD,self.txtEvaporatorTubesID)
        wx.GetApp().SetTopWindow(TubeSelect)
        TubeSelect.Show()

    def SelectCoolingCoilTube(self, event): 
        TubeSelect = TubeSelectClass(self, -1, "")
        TubeSelect.setTarget(self.txtCoolingCoilTubesOD,self.txtCoolingCoilTubesID)
        wx.GetApp().SetTopWindow(TubeSelect)
        TubeSelect.Show()
        
    def SaveRun(self,Cycle,type=None,Header=False,path='test.csv'):
        
        if Cycle.CycleType=='Secondary':
            if Cycle.IHXType=='Coaxial':
                IHX=Cycle.CoaxialIHX
            else:
                IHX=Cycle.PHEIHX
                
            if Cycle.Mode=='AC':
                Lists=[(Cycle,'Cycle'),
                   (Cycle.Compressor,'Compressor'),
                   (Cycle.Condenser,'Condenser')]
            else:
                Lists=[(Cycle,'Cycle'),
                   (Cycle.Compressor,'Compressor'),
                   (Cycle.Evaporator,'Evaporator')]
                
            Lists+=[(IHX,'Internal Heat Exchanger'),
                   (Cycle.CoolingCoil,'Cooling Coil'),
                   (Cycle.LineSetSupply,'Line Set Supply'),
                   (Cycle.LineSetReturn,'Line Set Return'),
                   (Cycle.Pump,'Pump')]
        elif Cycle.CycleType=='DX':
            Lists=[(Cycle,'Cycle'),
                   (Cycle.Compressor,'Compressor'),
                   (Cycle.Condenser,'Condenser'),
                   (Cycle.Evaporator,'Evaporator'),
                   (Cycle.LineSetSupply,'Line Set Supply'),
                   (Cycle.LineSetReturn,'Line Set Return'), ]
        else:
            raise ValueError
        
        parahead=''
        paraunits=''
        paravals=''
        paratitle=''
        if type=='para':
            f=open(path,'a')
            #Collect the terms for the parametric study
            
            for key in self.ParamUnits.keys():
                if key in self.ParamValues.keys():
                    parahead+=','+key
                    paraunits+=','+self.ParamUnits[key]
                    paravals+=','+str(self.ParamValues[key])
                    paratitle+=','+'Parametric Study'
            parahead=parahead.lstrip(',')
            paraunits=paraunits.lstrip(',')
            paravals=paravals.lstrip(',')
            paratitle=paratitle.lstrip(',')
        else:
            f=open(path,'w')
        
        head=''
        units=''
        vals=''
        title=''
        
        for component,desc in Lists:
            if not hasattr(component,'OutputList'):
                raise AttributeError('Component '+component.__name__+'must have the function OutputList')
            for headVal,unitVal,value in component.OutputList():
                title+=','+desc
                head+=','+headVal
                units+=','+unitVal
                vals+=','+str(value)
        head=parahead+head.strip()
        units=paraunits+units.strip()
        vals=paravals+vals.strip()
        title=paratitle+title.strip()
            
        if Header==True:
            f.write('%s\n' % (title))
            f.write('%s\n' % (head))
            f.write('%s\n' % (units))
            f.write('%s\n' % (vals))
        else:
            f.write('%s\n' % (vals))
        f.close()
        
        return
        
    def writeRun(self, event):
        dlg = wx.FileDialog(self, "Save Output file", ".", "",
                               "Output files (*.achpo)|*.achpo|All Files|*.*", wx.SAVE)
        if (dlg.ShowModal() == wx.ID_OK):
            fileName = dlg.GetFilename()
            dirName = dlg.GetDirectory()
            if os.path.exists(dirName+'\\'+fileName):
                dlg = wx.MessageDialog(self, "File exists. Overwrite?", "Overwrite confirmation", wx.YES_NO|wx.NO_DEFAULT|wx.ICON_QUESTION)
                if (dlg.ShowModal() == wx.ID_YES):
                    self.SaveRun(self.Cycle,path=dirName+'\\'+fileName,Header=True)
            else:
                self.SaveRun(self.Cycle,path=dirName+'\\'+fileName,Header=True)

    def ARIMapInfo(self, event): # wxGlade: MyFrame.<event_handler>
        ARIInfo = ARIInfoClass(self, -1, "")
        wx.GetApp().SetTopWindow(ARIInfo)
        ARIInfo.Show()

    def SelectLineSetInsulation(self, event): # wxGlade: MyFrame.<event_handler>
        if self.cmbLineSetInsulk.GetStringSelection()=='Armaflex':
            self.txtLineSetInsulk.SetValue('0.036')
        elif self.cmbLineSetInsulk.GetStringSelection()=='Polystyrene':
            self.txtLineSetInsulk.SetValue('0.027')

    def SelectParaPath(self, event): # wxGlade: MyFrame.<event_handler>
        print "Event handler `SelectParaPath' not implemented"
        event.Skip()

    def FileQuit(self, event): # wxGlade: MyFrame.<event_handler>
        self.Close()

    def modifyEvaporatorFink(self, event): 
        self.cmbEvaporatorFink.SetValue('User Defined')

    def modifyCoolingCoilFink(self, event): 
        self.cmbCoolingCoilFink.SetValue('User Defined')

    def modifyLineSetInsulation(self, event): 
        self.cmbLineSetInsulk.SetValue('User Defined')

    def openSchematic(self, event): 
        #This is temporarily deprecated
        dpi=100
        Schematic=SchematicFrameClass(self,-1,"")
        Schematic.Maximize()
        (w,h)=Schematic.GetClientSize()
        Schematic.pnlSchematicDetailed.doPlot(self.Cycle,figsize=(w/dpi,h/dpi),dpi=dpi,detailed='yes')
        
        wx.GetApp().SetTopWindow(Schematic)
        Schematic.Show()

if __name__=='__main__':
    #Running code for debugging purposes
    
    ACHPApp = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    
    MainFrame = ACHPMainFrame(None, -1)
    ACHPApp.SetTopWindow(MainFrame)
    MainFrame.Show()

    ACHPApp.MainLoop()
    