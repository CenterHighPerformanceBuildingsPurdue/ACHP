
import wx,os

# This file takes the place of wxGlade and is used to form the USER interface
# using the wxPython wrappers of wxWidgets objects
#

def LoadParametricParams(path):
    paramsList=[]
    lines=open(path,'r').readlines()
    for line in lines:
        #Strip whitespace
        line=line.strip()
        #Skip lines starting with #
        if not line.lstrip().startswith('#') and len(line)>0:
            (desc,units,val)=line.strip().split('::')
            desc=desc.strip()
            units=units.strip()
            val=val.strip()
            paramsList.append((desc,units,val))
    return paramsList
    
def LoadGUI(self):
    
    ##Imports
    from PlotCycleOutputsPanel import PlotCycleOutputsPanel
    from MPLPanel import MPLPanel
    
    _icon = wx.EmptyIcon()
    _icon.CopyFromBitmap(wx.Bitmap(os.path.join('imgs','sc.ico'), wx.BITMAP_TYPE_ANY))
    self.SetIcon(_icon)
    
    def MenuBar(self):
        #################################
        ####       Menu Bar         #####
        #################################
        
        # Menu Bar
        self.MenuBar = wx.MenuBar()
        
        self.File = wx.Menu()
        self.FileSave = wx.MenuItem(self.File, wx.NewId(), "Save Config File...\tCtrl+S", "", wx.ITEM_NORMAL)
        self.File.AppendItem(self.FileSave)
        self.FileOpen = wx.MenuItem(self.File, wx.NewId(), "Open Config File...\tCtrl+O", "", wx.ITEM_NORMAL)
        self.File.AppendItem(self.FileOpen)
        self.menuFileQuit = wx.MenuItem(self.File, wx.NewId(), "Quit\tCtrl+Q", "", wx.ITEM_NORMAL)
        self.File.AppendItem(self.menuFileQuit)
        self.MenuBar.Append(self.File, "File")
        
        self.Solve = wx.Menu()
        self.SolveSolve = wx.MenuItem(self.Solve, wx.NewId(), "Solve\tF5", "", wx.ITEM_NORMAL)
        self.Solve.AppendItem(self.SolveSolve)
        self.MenuBar.Append(self.Solve, "Solve")
        
        self.Help = wx.Menu()
        self.HelpHelp = wx.MenuItem(self.Help, wx.NewId(), "Help...\tCtrl+H", "", wx.ITEM_NORMAL)
        self.Help.AppendItem(self.HelpHelp)
        self.MenuBar.Append(self.Help, "Help")
        
        #Actually set it
        self.SetMenuBar(self.MenuBar)
    def Panels(self):
        #################################
        ####  Define all the panels #####
        #################################
    
        #Main notebook:
        self.nbIO = wx.Notebook(self, -1, style=0)
        #Create panel for input notebook, and attach it to notebook
        InputsPanel = wx.Panel(self.nbIO, -1)
        self.nbIO.AddPage(InputsPanel, "Input")
        #Create panel for input notebook, and attach it to notebook
        self.SolversPanel = wx.Panel(self.nbIO, -1)
        self.nbIO.AddPage(self.SolversPanel, "Solvers")
        #Create panel for output notebook, and attach it to notebook
        OutputsPanel = wx.Panel(self.nbIO, -1)
        self.nbIO.AddPage(OutputsPanel, "Output")
        
        ##Inputs tab panels:
        self.nbInputs = wx.Notebook(InputsPanel, -1, style=0)
        hboxCycleInputs = wx.BoxSizer(wx.VERTICAL)
        hboxCycleInputs.Add(self.nbInputs)
        InputsPanel.SetSizer(hboxCycleInputs)
        
        self.CycleInputsPanel=wx.Panel(self.nbInputs,-1)
        self.nbInputs.AddPage(self.CycleInputsPanel, "Cycle")
        
        self.CoolingCoilInputsPanel=wx.Panel(self.nbInputs,-1)
        self.nbInputs.AddPage(self.CoolingCoilInputsPanel, "Indoor Coil")
        
        self.CondenserInputsPanel=wx.Panel(self.nbInputs,-1)
        self.nbInputs.AddPage(self.CondenserInputsPanel, "Outdoor Coil")
        
        self.PumpIHXInputsPanel=wx.Panel(self.nbInputs,-1)
        self.nbInputs.AddPage(self.PumpIHXInputsPanel, "Pump && IHX")
        
        self.CompressorInputsPanel=wx.Panel(self.nbInputs,-1)
        self.nbInputs.AddPage(self.CompressorInputsPanel, "Compressor")
        
        self.LineSetInputsPanel=wx.Panel(self.nbInputs,-1)
        self.nbInputs.AddPage(self.LineSetInputsPanel, "Line Set")        
        
        ##Solver notebook:
        self.nbSolvers = wx.Notebook(self.SolversPanel, -1, style=0)
        hboxSolvers = wx.BoxSizer(wx.VERTICAL)
        hboxSolvers.Add(self.nbSolvers)
        self.SolversPanel.SetSizer(hboxSolvers)
        
        self.SolverMethodPanel=wx.Panel(self.nbSolvers,-1)
        self.nbSolvers.AddPage(self.SolverMethodPanel, "Main")
        
        ##Output notebook:
        self.nbOutputs = wx.Notebook(OutputsPanel, -1, style=0)
        hboxCycleOutputs = wx.BoxSizer(wx.VERTICAL)
        hboxCycleOutputs.Add(self.nbOutputs)
        OutputsPanel.SetSizer(hboxCycleOutputs)
        
        self.MainCycleOutputsPanel=wx.Panel(self.nbOutputs,-1)
        self.nbOutputs.AddPage(self.MainCycleOutputsPanel, "Main")
        
        self.TsPlotPanel=wx.Panel(self.nbOutputs,-1)
        self.nbOutputs.AddPage(self.TsPlotPanel, "T-s Plot")
        
        self.phPlotPanel=wx.Panel(self.nbOutputs,-1)
        self.nbOutputs.AddPage(self.phPlotPanel, "p-h Plot")
        
        self.CompressorOutputsPanel=wx.Panel(self.nbOutputs,-1)
        self.nbOutputs.AddPage(self.CompressorOutputsPanel, "Compressor")
        
        self.LineSetOutputsPanel=wx.Panel(self.nbOutputs,-1)
        self.nbOutputs.AddPage(self.LineSetOutputsPanel, "Line Set")
        
        self.CondenserOutputsPanel=wx.Panel(self.nbOutputs,-1)
        self.nbOutputs.AddPage(self.CondenserOutputsPanel, "Condenser")
        
        self.PumpIHXOutputsPanel=wx.Panel(self.nbOutputs,-1)
        self.nbOutputs.AddPage(self.PumpIHXOutputsPanel, "IHX && Pump")
        
        self.CoolingCoilOutputsPanel=wx.Panel(self.nbOutputs,-1)
        self.nbOutputs.AddPage(self.CoolingCoilOutputsPanel, "Cooling Coil")
        
        self.EvaporatorOutputsPanel=wx.Panel(self.nbOutputs,-1)
        self.nbOutputs.AddPage(self.EvaporatorOutputsPanel, "Evaporator")
    def CycleInputs(self):
        ##--------------------------------------------
        ##        Cycle Inputs
        ##--------------------------------------------
        hboxMainInputs = wx.GridBagSizer()
        
        self.lblPrimaryCycleInputs = wx.StaticText(self.CycleInputsPanel, label="Primary Cycle")
        self.lblPrimaryCycleInputs.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        hboxMainInputs.Add(self.lblPrimaryCycleInputs,pos=(0,0))
        
        self.lblSecondaryCycleInputs = wx.StaticText(self.CycleInputsPanel, label="Secondary Cycle")
        self.lblSecondaryCycleInputs.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        hboxMainInputs.Add(self.lblSecondaryCycleInputs,pos=(0,3))
        
        self.lblCycleMode = wx.StaticText(self.CycleInputsPanel, label="Cycle Mode")
        self.lblCycleMode.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        hboxMainInputs.Add(self.lblCycleMode,pos=(0,5))
        
        self.radCycleMode = wx.RadioBox(self.CycleInputsPanel, -1, "", choices=["Cooling Mode","Heating Mode"], majorDimension=0, style=wx.RA_SPECIFY_ROWS)        
        
        self.lineCycleInputs1=wx.StaticLine(self.CycleInputsPanel, -1)
        
        hboxMainInputs.Add(self.radCycleMode, pos=(2,5))
        
        hboxMainInputs.Add(self.lineCycleInputs1,pos=(1,0),span=(1,6),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblRefrigerant = wx.StaticText(self.CycleInputsPanel, -1, "Primary Refrigerant")
        self.cmbRefrigerant = wx.ComboBox(self.CycleInputsPanel, -1, choices=["R290", "R410A", "R404A", "R134a"], style=wx.CB_DROPDOWN)
        self.optCycleCharge = wx.RadioButton(self.CycleInputsPanel, -1, "Primary Charge [kg]", style=wx.RB_GROUP)
        self.txtCycleCharge = wx.TextCtrl(self.CycleInputsPanel, -1, "0.8")
        self.optCycleSubcooling = wx.RadioButton(self.CycleInputsPanel, -1, "Subcooling [K]")
        self.txtCycleSubcooling = wx.TextCtrl(self.CycleInputsPanel, -1, "7")
        self.lblCycleDTsh = wx.StaticText(self.CycleInputsPanel, -1, "\"Evaporator\" outlet superheat [K]")
        self.txtCycleDTsh = wx.TextCtrl(self.CycleInputsPanel, -1, "5")
        
        self.optSecFluid = wx.RadioButton(self.CycleInputsPanel, -1, "Yes/ Fluid", style=wx.RB_GROUP)
        self.cmbSecFluid = wx.ComboBox(self.CycleInputsPanel, -1, choices=["Water", "EG-20%", "EG-30%", "HC-10", "PG-20%", "PG-30%", "Methanol-10%", "Methanol-20%", "Methanol-30%", "Methanol-40%", "NH3/H2O-10%", "NH3/H2O-20%"], style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.optCycleDX = wx.RadioButton(self.CycleInputsPanel, -1, "No/ DX Primary")
        
        
        fgsPrimaryInputs = wx.FlexGridSizer(4, 2, 3, 5)
        fgsPrimaryInputs.AddMany([
                     (self.lblRefrigerant),(self.cmbRefrigerant), 
                     (self.lblCycleDTsh),(self.txtCycleDTsh),
                     (self.optCycleCharge),(self.txtCycleCharge),
                     (self.optCycleSubcooling),(self.txtCycleSubcooling),
                     ])
        hboxMainInputs.Add(fgsPrimaryInputs, pos=(2,0))
        
        fgsSecondaryInputs = wx.FlexGridSizer(2, 2, 3, 5)
        fgsSecondaryInputs.AddMany([
                     (self.optSecFluid),(self.cmbSecFluid), 
                     (self.optCycleDX),((0,0)),
                     ])
        hboxMainInputs.Add(fgsSecondaryInputs, pos=(2,3))
        
        self.bmpCycle = wx.StaticBitmap(self.CycleInputsPanel, -1, wx.Bitmap(os.path.join('imgs','Secondary-Cycle-AC.png'), wx.BITMAP_TYPE_ANY))
        hboxMainInputs.Add(self.bmpCycle, pos=(5,0),span=(5,20))
        
        listREFPROP=["REFPROP-1BUTENE","REFPROP-ACETONE",
            "REFPROP-AMMONIA","REFPROP-ARGON","REFPROP-BENZENE","REFPROP-BUTANE",
            "REFPROP-C12","REFPROP-C2BUTENE","REFPROP-C4F10","REFPROP-C5F12",
            "REFPROP-CF3I","REFPROP-CO","REFPROP-CO2","REFPROP-COS",
            "REFPROP-CYCLOHEX","REFPROP-CYCLOPRO","REFPROP-D2","REFPROP-D2O",
            "REFPROP-DECANE","REFPROP-DME","REFPROP-ETHANE","REFPROP-ETHANOL",
            "REFPROP-ETHYLENE","REFPROP-FLUORINE","REFPROP-H2S","REFPROP-HELIUM",
            "REFPROP-HEPTANE","REFPROP-HEXANE","REFPROP-HMX","REFPROP-HYDROGEN",
            "REFPROP-IBUTENE","REFPROP-IHEXANE","REFPROP-IPENTANE","REFPROP-ISOBUTAN",
            "REFPROP-KRYPTON","REFPROP-METHANE","REFPROP-METHANOL","REFPROP-N2O",
            "REFPROP-NEON","REFPROP-NEOPENTN","REFPROP-NF3","REFPROP-NITROGEN",
            "REFPROP-NONANE","REFPROP-OCTANE","REFPROP-OXYGEN","REFPROP-PARAHYD",
            "REFPROP-PENTANE","REFPROP-PROPANE","REFPROP-PROPYLEN","REFPROP-PROPYNE",
            "REFPROP-R11","REFPROP-R113","REFPROP-R114","REFPROP-R115","REFPROP-R116",
            "REFPROP-R12","REFPROP-R123","REFPROP-R1234YF","REFPROP-R124","REFPROP-R125",
            "REFPROP-R13","REFPROP-R134A","REFPROP-R14","REFPROP-R141B","REFPROP-R142B",
            "REFPROP-R143A","REFPROP-R152A","REFPROP-R21","REFPROP-R218","REFPROP-R22",
            "REFPROP-R227EA","REFPROP-R23","REFPROP-R236EA","REFPROP-R236FA",
            "REFPROP-R245CA","REFPROP-R245FA","REFPROP-R32","REFPROP-R365MFC",
            "REFPROP-R41","REFPROP-RC318","REFPROP-SF6","REFPROP-SO2",
            "REFPROP-T2BUTENE","REFPROP-TOLUENE","REFPROP-WATER","REFPROP-XENON","REFPROP-MIX:R32[0.697615]&R125[0.302385]"]
            
            #Can add the Pseudo-pure back if you can get them to work properly...
            
        self.cmbRefrigerant.AppendItems(listREFPROP)
        
        self.CycleInputsPanel.SetSizer(hboxMainInputs)
    def PumpIHXInputs(self):
        """
        --------------------------------------------
                         Pump Inputs
        --------------------------------------------
        """
        
        hboxPumpIHXInputs = wx.GridBagSizer()
        
        
        self.lblPumpInputs = wx.StaticText(self.PumpIHXInputsPanel, label="Pump Inputs")
        self.lblPumpInputs.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        hboxPumpIHXInputs.Add(self.lblPumpInputs,pos=(0,0))
        
        self.lineCoolingCoil1=wx.StaticLine(self.PumpIHXInputsPanel, -1)
        hboxPumpIHXInputs.Add(self.lineCoolingCoil1,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
                
        
        self.lblPumpp = wx.StaticText(self.PumpIHXInputsPanel, -1, "Nominal Loop Pressure [kPa]")
        self.txtPumpp = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "200")
        self.lblPumpp.SetToolTipString("Pressure of the secondary loop.  Not currently used since properties are only a function of temperature.")
        self.lblPumpmdot = wx.StaticText(self.PumpIHXInputsPanel, -1, "Mass Flow Rate [kg/s]")
        self.txtPumpmdot = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.38")
        self.lblPumpmdot.SetToolTipString("Mass flow rate of secondary fluid")
        self.lblPumpEfficiency = wx.StaticText(self.PumpIHXInputsPanel, -1, "Overall Efficiency (Motor+Pump) [-]")
        self.txtPumpEfficiency = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.5")
        self.lblPumpEfficiency.SetToolTipString("Combined efficiency of motor and pump set")
        self.lblPumpEfficiency.SetMinSize((180, -1))

        fgsPumpInputs = wx.FlexGridSizer(3, 2, 3, 5)
        fgsPumpInputs.AddMany([
                     (self.lblPumpp),(self.txtPumpp), 
                     (self.lblPumpmdot),(self.txtPumpmdot),
                     (self.lblPumpEfficiency),(self.txtPumpEfficiency)
                     ])
        hboxPumpIHXInputs.Add(fgsPumpInputs, pos=(2,0))
        
        self.optIHXUseCoaxial = wx.RadioButton(self.PumpIHXInputsPanel, -1, "Use Coaxial IHX")
        self.lblIHXLength = wx.StaticText(self.PumpIHXInputsPanel, -1, "Length of IHX [m]")
        self.txtIHXLength = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "15")
        self.lblIHXAnnOD = wx.StaticText(self.PumpIHXInputsPanel, -1, "Annulus OD [m]")
        self.txtIHXAnnOD = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.08")
        self.lblIHXAnnID = wx.StaticText(self.PumpIHXInputsPanel, -1, "Annulus ID [m]")
        self.txtIHXAnnID = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.03415")
        self.lblIHXTubeID = wx.StaticText(self.PumpIHXInputsPanel, -1, "Tube ID [m]")
        self.txtIHXTubeID = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.0278")
        self.lblIHXLength.SetToolTipString("Length of the internal heat exchanger")
        self.lblIHXAnnOD.SetToolTipString("Outer diameter of the wetted portion of the annulus")
        self.lblIHXAnnID.SetToolTipString("Inner diameter of the wetted portion of the annulus")
        self.lblIHXTubeID.SetToolTipString("Inner diameter of the tube at the center of the IHX")
        
        self.optIHXUsePHE = wx.RadioButton(self.PumpIHXInputsPanel, -1, "Use Plate IHX")
        self.lblIHXPlateBp = wx.StaticText(self.PumpIHXInputsPanel, -1, "Width [-]")
        self.txtIHXPlateBp = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.117")
        self.lblIHXPlateLp = wx.StaticText(self.PumpIHXInputsPanel, -1, "Port-Port centerline distance [-]")
        self.txtIHXPlateLp = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.300")
        self.lblIHXPlateN = wx.StaticText(self.PumpIHXInputsPanel, -1, "Number of plates [-]")
        self.txtIHXPlateN = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "10")
        self.lblIHXPlateAmplitude = wx.StaticText(self.PumpIHXInputsPanel, -1, "Plate amplitude [m]")
        self.txtIHXPlateAmplitude = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.001")
        self.lblIHXPlateThickness = wx.StaticText(self.PumpIHXInputsPanel, -1, "Plate thickness [m]")
        self.txtIHXPlateThickness = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.0005")
        self.lblIHXPlateConductivity = wx.StaticText(self.PumpIHXInputsPanel, -1, "Plate conductivity [W/m-K]")
        self.txtIHXPlateConductivity = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.03415")
        self.lblIHXPlateWavelength = wx.StaticText(self.PumpIHXInputsPanel, -1, "Plate wavelength [m]")
        self.txtIHXPlateWavelength = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "0.0278")
        self.lblIHXInclinationAngle = wx.StaticText(self.PumpIHXInputsPanel, -1, "Inclination Angle [deg]")
        self.txtIHXInclinationAngle = wx.TextCtrl(self.PumpIHXInputsPanel, -1, "45")
        self.radIHXChannelSelect = wx.RadioBox(self.PumpIHXInputsPanel, -1, "Additional Channel", choices=["Hot","Cold"], majorDimension=0, style=wx.RA_SPECIFY_ROWS)
        
        self.lblPumpInputs = wx.StaticText(self.PumpIHXInputsPanel, label="IHX Inputs")
        self.lblPumpInputs.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        hboxPumpIHXInputs.Add(self.lblPumpInputs,pos=(3,0))
        
        self.lineCoolingCoil1=wx.StaticLine(self.PumpIHXInputsPanel, -1)
        hboxPumpIHXInputs.Add(self.lineCoolingCoil1,pos=(4,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        fgsIHXInputs = wx.FlexGridSizer(cols = 2, vgap = 3, hgap = 5)
        fgsIHXInputs.AddMany([
                     (self.optIHXUseCoaxial),(0,0),
                     (self.lblIHXLength),(self.txtIHXLength),
                     (self.lblIHXAnnOD),(self.txtIHXAnnOD),
                     (self.lblIHXAnnID),(self.txtIHXAnnID),
                     (self.lblIHXTubeID),(self.txtIHXTubeID),
                     ])
        hboxPumpIHXInputs.Add(fgsIHXInputs, pos=(5,0))
        
        fgsPlateIHXInputs = wx.FlexGridSizer(cols = 2, vgap = 3, hgap = 5)
        fgsPlateIHXInputs.AddMany([
                    (self.optIHXUsePHE),(0,0),
                    (self.lblIHXPlateBp),(self.txtIHXPlateBp),(self.lblIHXPlateLp),(self.txtIHXPlateLp),
                    (self.lblIHXPlateN),(self.txtIHXPlateN),
                    (self.lblIHXPlateAmplitude),(self.txtIHXPlateAmplitude),
                    (self.lblIHXPlateThickness),(self.txtIHXPlateThickness),
                    (self.lblIHXPlateConductivity),(self.txtIHXPlateConductivity),
                    (self.lblIHXPlateWavelength),(self.txtIHXPlateWavelength),
                    (self.lblIHXInclinationAngle),(self.txtIHXInclinationAngle),
                    (self.radIHXChannelSelect),(0,0)
                     ])
        hboxPumpIHXInputs.Add(fgsPlateIHXInputs, pos=(5,3))
        
        self.bmpIHXSchem = wx.StaticBitmap(self.PumpIHXInputsPanel, -1, wx.Bitmap(os.path.join('imgs','IHX.png'), wx.BITMAP_TYPE_ANY))
        hboxPumpIHXInputs.Add(self.bmpIHXSchem,pos=(6,0))
        
        self.PumpIHXInputsPanel.SetSizer(hboxPumpIHXInputs)      
    def CoolingCoilInputs(self):
        ### ----
        ### Indoor Coil
        ### ----

        sizerCoolingCoilInputs = wx.GridBagSizer()
        
        self.lblCoolingCoilFins = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Fins")
        self.lblCoolingCoilFins.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        sizerCoolingCoilInputs.Add(self.lblCoolingCoilFins,pos=(0,0))
        
        self.lineCoolingCoil1=wx.StaticLine(self.CoolingCoilInputsPanel, -1)
        sizerCoolingCoilInputs.Add(self.lineCoolingCoil1,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        #Fins
        self.lblCoolingCoilFinType = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Fin Type")
        self.cmbCoolingCoilFinType = wx.ComboBox(self.CoolingCoilInputsPanel, -1, choices=["Wavy Lanced"], style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.lblCoolingCoilFinFPI = wx.StaticText(self.CoolingCoilInputsPanel, -1, "FPI [1/in]")
        self.txtCoolingCoilFinFPI = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "14.5")
        self.lblCoolingCoilFinpd = wx.StaticText(self.CoolingCoilInputsPanel, -1, "pd [m]")
        self.txtCoolingCoilFinpd = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "0.001")
        self.lblCoolingCoilFinxf = wx.StaticText(self.CoolingCoilInputsPanel, -1, "xf [m]")
        self.txtCoolingCoilFinxf = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "0.001")
        self.lblCoolingCoilFint = wx.StaticText(self.CoolingCoilInputsPanel, -1, "t [m]")
        self.txtCoolingCoilFint = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "0.00011")
        self.lblCoolingCoilFink = wx.StaticText(self.CoolingCoilInputsPanel, -1, "k [W/m-K]")
        self.txtCoolingCoilFink = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "237")
        self.cmbCoolingCoilFink = wx.ComboBox(self.CoolingCoilInputsPanel, -1, choices=["Aluminum", "Copper", "User Specified"], style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.lblCoolingCoilTubesNtubes = wx.StaticText(self.CoolingCoilInputsPanel, -1, "# Tubes / bank [-]")
        self.txtCoolingCoilTubesNtubes = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "32")
        self.lblCoolingCoilTubesNbank = wx.StaticText(self.CoolingCoilInputsPanel, -1, "# bank [-]")
        self.txtCoolingCoilTubesNbank = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "4")
        self.lblCoolingCoilTubesL = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Tube Length [m]")
        self.txtCoolingCoilTubesL = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "0.452")
        self.lblCoolingCoilTubesOD = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Tube OD [m]")
        self.txtCoolingCoilTubesOD = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "0.009525")
        self.cmdSelectCoolingCoilTube = wx.Button(self.CoolingCoilInputsPanel, -1, "Select")
        self.lblCoolingCoilTubesID = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Tube ID [m]")
        self.txtCoolingCoilTubesID = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "0.0089154")
        self.lblCoolingCoilTubesPl = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Longitud. Pitch [m]")
        self.txtCoolingCoilTubesPl = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "0.0254")
        self.lblCoolingCoilTubesPt = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Transverse Pitch [m]")
        self.txtCoolingCoilTubesPt = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "0.0219964")
        self.lblCoolingCoilTubesNcircuit = wx.StaticText(self.CoolingCoilInputsPanel, -1, "# circuit [-]")
        self.txtCoolingCoilTubesNcircuit = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "4")
        self.cmdShowCoolingCoilCircuits = wx.Button(self.CoolingCoilInputsPanel, -1, "Show Circuits ...")
        
        #Tubes
        self.lblCoolingCoilFinType.SetMinSize((120, -1))
        self.lblCoolingCoilFinType.SetToolTipString("The type of fins employed")
        self.cmbCoolingCoilFinType.SetSelection(0)
        self.lblCoolingCoilFinFPI.SetToolTipString("Number of fins per inch")
        self.lblCoolingCoilFinpd.SetToolTipString("Twice the amplitude of the wave of the fin")
        self.lblCoolingCoilFinxf.SetToolTipString("Half the wavelength of the wave of the fin")
        self.lblCoolingCoilFint.SetToolTipString("Thickness of the fin")
        self.lblCoolingCoilFink.SetToolTipString("Fin Conductivity")
        self.cmbCoolingCoilFink.SetSelection(0)
        self.lblCoolingCoilTubesNtubes.SetMinSize((120, -1))
        self.lblCoolingCoilTubesNtubes.SetToolTipString("Number of tubes per bank (AKA row) of the heat exchanger.  See circuiting figure for clarification")
        self.lblCoolingCoilTubesNbank.SetToolTipString("Number of banks (AKA rows) of tubes in the heat exchanger (see coil circuiting diagram for clarification)")
        self.lblCoolingCoilTubesL.SetToolTipString("Length of one tube.  If looking face-on at the coil installed in duct with the tubes horizontally, the width of the duct")
        self.lblCoolingCoilTubesOD.SetToolTipString("Outer diameter of the tubes")
        self.lblCoolingCoilTubesID.SetToolTipString("Inner diameter of the tubes")
        self.lblCoolingCoilTubesPl.SetToolTipString("Longitudinal pitch of the tubes.  Also called horizontal pitch.  If you look edge-on at the coil with flow of air from left to right, the horizontal spacing between banks")
        self.lblCoolingCoilTubesPt.SetToolTipString("Transverse pitch of tubes.  Also called vertical pitch.  Vertical spacing of tubes if looked at end-on with flow from left to right")
        self.lblCoolingCoilTubesNcircuit.SetToolTipString("Number of circuits that the flow is distributed betwen for the working fluid.  See circuiting diagram for clarity")
        
        fgsCoolingCoilInputs1 = wx.FlexGridSizer(6, 3, 2, 2)
        fgsCoolingCoilInputs1.AddMany([(self.lblCoolingCoilFinType), (self.cmbCoolingCoilFinType), ((0,0)),
                     (self.lblCoolingCoilFinFPI), (self.txtCoolingCoilFinFPI), ((0,0)),
                     (self.lblCoolingCoilFinpd), (self.txtCoolingCoilFinpd), ((0,0)),
                     (self.lblCoolingCoilFinxf),(self.txtCoolingCoilFinxf),  ((0,0)),
                     (self.lblCoolingCoilFint),(self.txtCoolingCoilFint),  ((0,0)),
                     (self.lblCoolingCoilFink),(self.txtCoolingCoilFink),(self.cmbCoolingCoilFink)
                     ])
        fgsCoolingCoilInputs1.AddGrowableRow(2, 1)
        fgsCoolingCoilInputs1.AddGrowableCol(1, 1)
        sizerCoolingCoilInputs.Add(fgsCoolingCoilInputs1,pos=(2,0))
        self.bmpCoolingCoilFins = wx.StaticBitmap(self.CoolingCoilInputsPanel, -1, wx.Bitmap("imgs/Fins.png", wx.BITMAP_TYPE_ANY))
        sizerCoolingCoilInputs.Add(self.bmpCoolingCoilFins,pos=(2,4))
        
        self.lblCoolingCoilTubes = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Tubes")
        self.lblCoolingCoilTubes.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        sizerCoolingCoilInputs.Add(self.lblCoolingCoilTubes,pos=(3,0))
        
        self.lineCoolingCoil2=wx.StaticLine(self.CoolingCoilInputsPanel, -1)
        sizerCoolingCoilInputs.Add(self.lineCoolingCoil2,pos=(4,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        fgsCoolingCoilInputs2 = wx.FlexGridSizer(8, 3, 2, 2)
        fgsCoolingCoilInputs2.AddMany([
                (self.lblCoolingCoilTubesNtubes), (self.txtCoolingCoilTubesNtubes), ((0,0)),
                (self.lblCoolingCoilTubesNbank), (self.txtCoolingCoilTubesNbank), ((0,0)),
                (self.lblCoolingCoilTubesL), (self.txtCoolingCoilTubesL), ((0,0)),
                (self.lblCoolingCoilTubesOD),(self.txtCoolingCoilTubesOD),  (self.cmdSelectCoolingCoilTube),
                (self.lblCoolingCoilTubesID),(self.txtCoolingCoilTubesID),  ((0,0)),
                (self.lblCoolingCoilTubesPl),(self.txtCoolingCoilTubesPl),  ((0,0)),
                (self.lblCoolingCoilTubesPt),(self.txtCoolingCoilTubesPt),  ((0,0)),
                (self.lblCoolingCoilTubesNcircuit),(self.txtCoolingCoilTubesNcircuit),  (self.cmdShowCoolingCoilCircuits),
            ])
        fgsCoolingCoilInputs2.AddGrowableRow(2, 1)
        fgsCoolingCoilInputs2.AddGrowableCol(1, 1)
        sizerCoolingCoilInputs.Add(fgsCoolingCoilInputs2,pos=(5,0))
        
        self.bmpCoolingCoilTubes = wx.StaticBitmap(self.CoolingCoilInputsPanel, -1, wx.Bitmap("imgs/Tubes.png", wx.BITMAP_TYPE_ANY))
        sizerCoolingCoilInputs.Add(self.bmpCoolingCoilTubes,pos=(5,4))
        
        #Air         
        self.lblCoolingCoilAirVdot = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Vdot [m^3/s]")
        self.txtCoolingCoilAirVdot = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "0.56319")
        self.lblCoolingCoilAirTdb = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Tdb [K]")
        self.txtCoolingCoilAirTdb = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "297.039")
        self.lblCoolingCoilAirp = wx.StaticText(self.CoolingCoilInputsPanel, -1, "p [kPa]")
        self.txtCoolingCoilAirp = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "101.325")
        self.lblCoolingCoilAirRH = wx.StaticText(self.CoolingCoilInputsPanel, -1, "RH [-]")
        self.txtCoolingCoilAirRH = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "0.5")
        self.lblCoolingCoilPower = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Fan Power [W]")
        self.txtCoolingCoilPower = wx.TextCtrl(self.CoolingCoilInputsPanel, -1, "260")
        self.lblCoolingCoilAirVdot.SetMinSize((120, -1))
        self.lblCoolingCoilAirVdot.SetToolTipString("Volumetric flow of humid air into the coil in cubic meters of air per second.")
        self.lblCoolingCoilAirTdb.SetToolTipString("Dry bulb temperature")
        self.lblCoolingCoilAirp.SetToolTipString("Ambient pressure (absolute)")
        self.lblCoolingCoilAirRH.SetToolTipString("Relative humidity of the air from 0 to 1  (0% to 100% RH)")
        self.lblCoolingCoilPower.SetToolTipString("Fan power of the blower")

        self.lblCoolingCoilAir = wx.StaticText(self.CoolingCoilInputsPanel, -1, "Air")
        self.lblCoolingCoilAir.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        sizerCoolingCoilInputs.Add(self.lblCoolingCoilAir,pos=(6,0))
        
        self.lineCoolingCoil3=wx.StaticLine(self.CoolingCoilInputsPanel, -1)
        sizerCoolingCoilInputs.Add(self.lineCoolingCoil3,pos=(7,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        fgsCoolingCoilInputs3 = wx.FlexGridSizer(5, 2, 2, 2)
        fgsCoolingCoilInputs3.AddMany([
                (self.lblCoolingCoilAirVdot), (self.txtCoolingCoilAirVdot),
                (self.lblCoolingCoilAirTdb), (self.txtCoolingCoilAirTdb),
                (self.lblCoolingCoilAirp),(self.txtCoolingCoilAirp),  
                (self.lblCoolingCoilAirRH),(self.txtCoolingCoilAirRH),
                (self.lblCoolingCoilPower),(self.txtCoolingCoilPower),
            ])
        fgsCoolingCoilInputs3.AddGrowableRow(2, 1)
        fgsCoolingCoilInputs3.AddGrowableCol(1, 1)
        sizerCoolingCoilInputs.Add(fgsCoolingCoilInputs3,pos=(8,0))
        
        self.CoolingCoilInputsPanel.SetSizer(sizerCoolingCoilInputs)
    def CondenserInputs(self):
        ##--------------------------------------------
        ##        Condenser Inputs
        ##--------------------------------------------
        
        sizerCondenserInputs = wx.GridBagSizer()
        
        self.lblCondenserFins = wx.StaticText(self.CondenserInputsPanel, -1, "Fins")
        self.lblCondenserFins.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        sizerCondenserInputs.Add(self.lblCondenserFins,pos=(0,0))
        
        self.lineCondenser1=wx.StaticLine(self.CondenserInputsPanel, -1)
        sizerCondenserInputs.Add(self.lineCondenser1,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblCondenserFinType = wx.StaticText(self.CondenserInputsPanel, -1, "Fin Type")
        self.cmbCondenserFinType = wx.ComboBox(self.CondenserInputsPanel, -1, choices=["Wavy Lanced"], style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.lblCondenserFinFPI = wx.StaticText(self.CondenserInputsPanel, -1, "FPI [1/in]")
        self.txtCondenserFinFPI = wx.TextCtrl(self.CondenserInputsPanel, -1, "25")
        self.lblCondenserFinpd = wx.StaticText(self.CondenserInputsPanel, -1, "pd [m]")
        self.txtCondenserFinpd = wx.TextCtrl(self.CondenserInputsPanel, -1, "0.001")
        self.lblCondenserFinxf = wx.StaticText(self.CondenserInputsPanel, -1, "xf [m]")
        self.txtCondenserFinxf = wx.TextCtrl(self.CondenserInputsPanel, -1, "0.001")
        self.lblCondenserFint = wx.StaticText(self.CondenserInputsPanel, -1, "t [m]")
        self.txtCondenserFint = wx.TextCtrl(self.CondenserInputsPanel, -1, "0.00011")
        self.lblCondenserFink = wx.StaticText(self.CondenserInputsPanel, -1, "k [W/m-K]")
        self.txtCondenserFink = wx.TextCtrl(self.CondenserInputsPanel, -1, "237")
        self.cmbCondenserFink = wx.ComboBox(self.CondenserInputsPanel, -1, choices=["Aluminum", "Copper", "User Specified"], style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.lblCondenserFinType.SetMinSize((120, -1))
        self.lblCondenserFinType.SetToolTipString("The type of fins employed")
        self.cmbCondenserFinType.SetSelection(0)
        self.cmbCondenserFink.SetSelection(0)
        self.lblCondenserFinFPI.SetToolTipString("Number of fins per inch")
        self.lblCondenserFinpd.SetToolTipString("Twice the amplitude of the wave of the fin")
        self.lblCondenserFinxf.SetToolTipString("Half the wavelength of the wave of the fin")
        self.lblCondenserFint.SetToolTipString("Thickness of the fin")
        self.lblCondenserFink.SetToolTipString("Fin Conductivity")
        
        fgsCondenserInputs1 = wx.FlexGridSizer(6, 3, 2, 2)
        fgsCondenserInputs1.AddMany([(self.lblCondenserFinType), (self.cmbCondenserFinType), ((0,0)),
                     (self.lblCondenserFinFPI), (self.txtCondenserFinFPI), ((0,0)),
                     (self.lblCondenserFinpd), (self.txtCondenserFinpd), ((0,0)),
                     (self.lblCondenserFinxf),(self.txtCondenserFinxf),  ((0,0)),
                     (self.lblCondenserFint),(self.txtCondenserFint),  ((0,0)),
                     (self.lblCondenserFink),(self.txtCondenserFink),(self.cmbCondenserFink)
                     ])
        fgsCondenserInputs1.AddGrowableRow(2, 1)
        fgsCondenserInputs1.AddGrowableCol(1, 1)
        sizerCondenserInputs.Add(fgsCondenserInputs1,pos=(2,0))
        self.bmpCondenserFins = wx.StaticBitmap(self.CondenserInputsPanel, -1, wx.Bitmap("imgs/Fins.png", wx.BITMAP_TYPE_ANY))
        sizerCondenserInputs.Add(self.bmpCondenserFins,pos=(2,4))
        
        self.lblCondenserTubes = wx.StaticText(self.CondenserInputsPanel, -1, "Tubes")
        self.lblCondenserTubes.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        sizerCondenserInputs.Add(self.lblCondenserTubes,pos=(3,0))
        
        self.lineCondenser2=wx.StaticLine(self.CondenserInputsPanel, -1)
        sizerCondenserInputs.Add(self.lineCondenser2,pos=(4,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)

        self.lblCondenserTubesNtubes = wx.StaticText(self.CondenserInputsPanel, -1, "# Tubes/bank [-]")
        self.txtCondenserTubesNtubes = wx.TextCtrl(self.CondenserInputsPanel, -1, "41")
        self.lblCondenserTubesNbank = wx.StaticText(self.CondenserInputsPanel, -1, "# bank [-]")
        self.txtCondenserTubesNbank = wx.TextCtrl(self.CondenserInputsPanel, -1, "1")
        self.lblCondenserTubesL = wx.StaticText(self.CondenserInputsPanel, -1, "Tube Length [m]")
        self.txtCondenserTubesL = wx.TextCtrl(self.CondenserInputsPanel, -1, "2.286")
        self.lblCondenserTubesOD = wx.StaticText(self.CondenserInputsPanel, -1, "Tube OD [m]")
        self.txtCondenserTubesOD = wx.TextCtrl(self.CondenserInputsPanel, -1, "0.007")
        self.lblCondenserTubesID = wx.StaticText(self.CondenserInputsPanel, -1, "Tube ID [m]")
        self.txtCondenserTubesID = wx.TextCtrl(self.CondenserInputsPanel, -1, "0.0063904")
        self.lblCondenserTubesPl = wx.StaticText(self.CondenserInputsPanel, -1, "Longit. Pitch [m]")
        self.txtCondenserTubesPl = wx.TextCtrl(self.CondenserInputsPanel, -1, "0.01905")
        self.lblCondenserTubesPt = wx.StaticText(self.CondenserInputsPanel, -1, "Transverse Pitch [m]")
        self.txtCondenserTubesPt = wx.TextCtrl(self.CondenserInputsPanel, -1, "0.022225")
        self.lblCondenserTubesNcircuit = wx.StaticText(self.CondenserInputsPanel, -1, "# circuit [-]")
        self.txtCondenserTubesNcircuit = wx.TextCtrl(self.CondenserInputsPanel, -1, "5")
        self.cmdShowCondenserCircuits = wx.Button(self.CondenserInputsPanel, -1, "Show Circuits ...")       
        self.lblCondenserTubesNtubes.SetMinSize((120, -1))
        self.lblCondenserTubesNtubes.SetToolTipString("Number of tubes per bank (AKA row) of the heat exchanger.  See circuiting figure for clarification")
        self.lblCondenserTubesNbank.SetToolTipString("Number of banks (AKA rows) of tubes in the heat exchanger (see coil circuiting diagram for clarification)")
        self.lblCondenserTubesL.SetToolTipString("Length of one tube.  If looking face-on at the coil installed in duct with the tubes horizontally, the width of the duct")
        self.lblCondenserTubesOD.SetToolTipString("Outer diameter of the tubes")
        self.lblCondenserTubesID.SetToolTipString("Inner diameter of the tubes")
        self.lblCondenserTubesPl.SetToolTipString("Longitudinal pitch of the tubes.  Also called horizontal pitch.  If you look edge-on at the coil with flow of air from left to right, the horizontal spacing between banks")
        self.lblCondenserTubesPt.SetToolTipString("Transverse pitch of tubes.  Also called vertical pitch.  Vertical spacing of tubes if looked at end-on with flow from left to right")
        self.lblCondenserTubesNcircuit.SetToolTipString("Number of circuits that the flow is distributed betwen for the working fluid.  See circuiting diagram for clarity") 
        
        fgsCondenserInputs2 = wx.FlexGridSizer(8, 3, 2, 2)
        fgsCondenserInputs2.AddMany([
                (self.lblCondenserTubesNtubes), (self.txtCondenserTubesNtubes), ((0,0)),
                (self.lblCondenserTubesNbank), (self.txtCondenserTubesNbank), ((0,0)),
                (self.lblCondenserTubesL), (self.txtCondenserTubesL), ((0,0)),
                (self.lblCondenserTubesOD),(self.txtCondenserTubesOD),  ((0,0)),
                (self.lblCondenserTubesID),(self.txtCondenserTubesID),  ((0,0)),
                (self.lblCondenserTubesPl),(self.txtCondenserTubesPl),  ((0,0)),
                (self.lblCondenserTubesPt),(self.txtCondenserTubesPt),  ((0,0)),
                (self.lblCondenserTubesNcircuit),(self.txtCondenserTubesNcircuit),  (self.cmdShowCondenserCircuits),
            ])
        fgsCondenserInputs2.AddGrowableRow(2, 1)
        fgsCondenserInputs2.AddGrowableCol(1, 1)
        sizerCondenserInputs.Add(fgsCondenserInputs2,pos=(5,0))
        
        self.bmpCondenserTubes = wx.StaticBitmap(self.CondenserInputsPanel, -1, wx.Bitmap("imgs/Tubes.png", wx.BITMAP_TYPE_ANY))
        sizerCondenserInputs.Add(self.bmpCondenserTubes,pos=(5,4))

        self.lblCondenserAirVdot = wx.StaticText(self.CondenserInputsPanel, -1, "Vdot air [m^3/s]")
        self.txtCondenserAirVdot = wx.TextCtrl(self.CondenserInputsPanel, -1, "1.7934")
        self.lblCondenserAirTdb = wx.StaticText(self.CondenserInputsPanel, -1, "Dry Bulb Temp [K]")
        self.txtCondenserAirTdb = wx.TextCtrl(self.CondenserInputsPanel, -1, "308.15")
        self.lblCondenserAirp = wx.StaticText(self.CondenserInputsPanel, -1, "p_atm [kPa]")
        self.txtCondenserAirp = wx.TextCtrl(self.CondenserInputsPanel, -1, "101.325")
        self.lblCondenserAirRH = wx.StaticText(self.CondenserInputsPanel, -1, "Rel. Hum. [-]")
        self.txtCondenserAirRH = wx.TextCtrl(self.CondenserInputsPanel, -1, "0.5")
        self.lblCondenserPower = wx.StaticText(self.CondenserInputsPanel, -1, "Rated Fan Power [W]")
        self.txtCondenserPower = wx.TextCtrl(self.CondenserInputsPanel, -1, "438")
        self.lblCondenserAirVdot.SetMinSize((120, -1))
        self.lblCondenserAirVdot.SetToolTipString("Volumetric flow of humid air into the coil in cubic meters of air per second.")
        self.lblCondenserAirTdb.SetToolTipString("Dry bulb temperature")
        self.lblCondenserAirp.SetToolTipString("Ambient pressure (absolute)")
        self.lblCondenserAirRH.SetToolTipString("Relative humidity of the air from 0 to 1  (0% to 100% RH)")
        self.lblCondenserPower.SetToolTipString("Fan power of the blower")

        self.lblCondenserAir = wx.StaticText(self.CondenserInputsPanel, -1, "Air")
        self.lblCondenserAir.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        sizerCondenserInputs.Add(self.lblCondenserAir,pos=(6,0))
        
        self.lineCondenser3=wx.StaticLine(self.CondenserInputsPanel, -1)
        sizerCondenserInputs.Add(self.lineCondenser3,pos=(7,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        fgsCondenserInputs3 = wx.FlexGridSizer(5, 2, 2, 2)
        fgsCondenserInputs3.AddMany([
                (self.lblCondenserAirVdot), (self.txtCondenserAirVdot),
                (self.lblCondenserAirTdb), (self.txtCondenserAirTdb),
                (self.lblCondenserAirp),(self.txtCondenserAirp),  
                (self.lblCondenserAirRH),(self.txtCondenserAirRH),
                (self.lblCondenserPower),(self.txtCondenserPower),
            ])
        fgsCondenserInputs3.AddGrowableRow(2, 1)
        fgsCondenserInputs3.AddGrowableCol(1, 1)
        sizerCondenserInputs.Add(fgsCondenserInputs3,pos=(8,0))
        
        self.CondenserInputsPanel.SetSizer(sizerCondenserInputs)
    
    def EvaporatorInputs(self):
        self.sizerEvaporatorInputs = wx.GridBagSizer()
        
        self.lblEvaporatorFins = wx.StaticText(self.EvaporatorInputsPanel, -1, "Fins")
        self.lblEvaporatorFins.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerEvaporatorInputs.Add(self.lblEvaporatorFins,pos=(0,0))
        
        self.lineEvaporator1=wx.StaticLine(self.EvaporatorInputsPanel, -1)
        self.sizerEvaporatorInputs.Add(self.lineEvaporator1,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblEvaporatorFinType = wx.StaticText(self.EvaporatorInputsPanel, -1, "Fin Type")
        self.cmbEvaporatorFinType = wx.ComboBox(self.EvaporatorInputsPanel, -1, choices=["Wavy Lanced"], style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.lblEvaporatorFinFPI = wx.StaticText(self.EvaporatorInputsPanel, -1, "FPI [1/in]")
        self.txtEvaporatorFinFPI = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "14.5")
        self.lblEvaporatorFinpd = wx.StaticText(self.EvaporatorInputsPanel, -1, "pd [m]")
        self.txtEvaporatorFinpd = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "0.001")
        self.lblEvaporatorFinxf = wx.StaticText(self.EvaporatorInputsPanel, -1, "xf [m]")
        self.txtEvaporatorFinxf = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "0.001")
        self.lblEvaporatorFint = wx.StaticText(self.EvaporatorInputsPanel, -1, "t [m]")
        self.txtEvaporatorFint = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "0.00011")
        self.lblEvaporatorFink = wx.StaticText(self.EvaporatorInputsPanel, -1, "k [W/m-K]")
        self.txtEvaporatorFink = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "237")
        self.cmbEvaporatorFink = wx.ComboBox(self.EvaporatorInputsPanel, -1, choices=["Aluminum", "Copper", "User Specified"], style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.lblEvaporatorFinType.SetToolTipString("The type of fins employed")
        self.cmbEvaporatorFinType.SetSelection(0)
        self.lblEvaporatorFinFPI.SetToolTipString("Number of fins per inch")
        self.lblEvaporatorFinpd.SetToolTipString("Twice the amplitude of the wave of the fin")
        self.lblEvaporatorFinxf.SetToolTipString("Half the wavelength of the wave of the fin")
        self.lblEvaporatorFint.SetToolTipString("Thickness of the fin")
        self.lblEvaporatorFink.SetMinSize((120, -1))
        self.lblEvaporatorFink.SetToolTipString("Fin Conductivity")
        self.cmbEvaporatorFink.SetSelection(0)

        fgsEvaporatorInputs1 = wx.FlexGridSizer(6, 3, 2, 2)
        fgsEvaporatorInputs1.AddMany([(self.lblEvaporatorFinType), (self.cmbEvaporatorFinType), ((0,0)),
                     (self.lblEvaporatorFinFPI), (self.txtEvaporatorFinFPI), ((0,0)),
                     (self.lblEvaporatorFinpd), (self.txtEvaporatorFinpd), ((0,0)),
                     (self.lblEvaporatorFinxf),(self.txtEvaporatorFinxf),  ((0,0)),
                     (self.lblEvaporatorFint),(self.txtEvaporatorFint),  ((0,0)),
                     (self.lblEvaporatorFink),(self.txtEvaporatorFink),(self.cmbEvaporatorFink)
                     ])
        fgsEvaporatorInputs1.AddGrowableRow(2, 1)
        fgsEvaporatorInputs1.AddGrowableCol(1, 1)
        self.sizerEvaporatorInputs.Add(fgsEvaporatorInputs1,pos=(2,0))
        self.bmpEvaporatorFins = wx.StaticBitmap(self.EvaporatorInputsPanel, -1, wx.Bitmap("imgs/Fins.png", wx.BITMAP_TYPE_ANY))
        self.sizerEvaporatorInputs.Add(self.bmpEvaporatorFins,pos=(2,4))
        
        self.lblEvaporatorTubes = wx.StaticText(self.EvaporatorInputsPanel, -1, "Tubes")
        self.lblEvaporatorTubes.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerEvaporatorInputs.Add(self.lblEvaporatorTubes,pos=(3,0))
        
        self.lineEvaporator2=wx.StaticLine(self.EvaporatorInputsPanel, -1)
        self.sizerEvaporatorInputs.Add(self.lineEvaporator2,pos=(4,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        

        self.lblEvaporatorTubesNtubes = wx.StaticText(self.EvaporatorInputsPanel, -1, "# Tubes / bank [-]")
        self.txtEvaporatorTubesNtubes = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "32")
        self.lblEvaporatorTubesNbank = wx.StaticText(self.EvaporatorInputsPanel, -1, "# bank [-]")
        self.txtEvaporatorTubesNbank = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "3")
        self.lblEvaporatorTubesL = wx.StaticText(self.EvaporatorInputsPanel, -1, "Tube Length [m]")
        self.txtEvaporatorTubesL = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "0.452")
        self.lblEvaporatorTubesOD = wx.StaticText(self.EvaporatorInputsPanel, -1, "Tube OD [m]")
        self.txtEvaporatorTubesOD = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "0.009525")
        self.cmdSelectEvaporatorTube = wx.Button(self.EvaporatorInputsPanel, -1, "Select")
        self.lblEvaporatorTubesID = wx.StaticText(self.EvaporatorInputsPanel, -1, "Tube ID [m]")
        self.txtEvaporatorTubesID = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "0.0089154")
        self.lblEvaporatorTubesPl = wx.StaticText(self.EvaporatorInputsPanel, -1, "Longit. Pitch [m]")
        self.txtEvaporatorTubesPl = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "0.0254")
        self.lblEvaporatorTubesPt = wx.StaticText(self.EvaporatorInputsPanel, -1, "Transverse Pitch [m]")
        self.txtEvaporatorTubesPt = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "0.0219964")
        self.lblEvaporatorTubesNcircuit = wx.StaticText(self.EvaporatorInputsPanel, -1, "# circuit [-]")
        self.txtEvaporatorTubesNcircuit = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "6")
        self.cmdShowEvaporatorCircuits = wx.Button(self.EvaporatorInputsPanel, -1, "Show Circuits ...")
        self.lblEvaporatorTubesNtubes.SetToolTipString("Number of tubes per bank (AKA row) of the heat exchanger.  See circuiting figure for clarification")
        self.lblEvaporatorTubesNbank.SetToolTipString("Number of banks (AKA rows) of tubes in the heat exchanger (see coil circuiting diagram for clarification)")
        self.lblEvaporatorTubesL.SetToolTipString("Length of one tube.  If looking face-on at the coil installed in duct with the tubes horizontally, the width of the duct")
        self.lblEvaporatorTubesOD.SetToolTipString("Outer diameter of the tubes")
        self.lblEvaporatorTubesID.SetToolTipString("Inner diameter of the tubes")
        self.lblEvaporatorTubesPl.SetToolTipString("Longitudinal pitch of the tubes.  Also called horizontal pitch.  If you look edge-on at the coil with flow of air from left to right, the horizontal spacing between banks")
        self.lblEvaporatorTubesPt.SetMinSize((120, -1))
        self.lblEvaporatorTubesPt.SetToolTipString("Transverse pitch of tubes.  Also called vertical pitch.  Vertical spacing of tubes if looked at end-on with flow from left to right")
        self.lblEvaporatorTubesNcircuit.SetToolTipString("Number of circuits that the flow is distributed betwen for the working fluid.  See circuiting diagram for clarity")

        fgsEvaporatorInputs2 = wx.FlexGridSizer(8, 3, 2, 2)
        fgsEvaporatorInputs2.AddMany([
                (self.lblEvaporatorTubesNtubes), (self.txtEvaporatorTubesNtubes), ((0,0)),
                (self.lblEvaporatorTubesNbank), (self.txtEvaporatorTubesNbank), ((0,0)),
                (self.lblEvaporatorTubesL), (self.txtEvaporatorTubesL), ((0,0)),
                (self.lblEvaporatorTubesOD),(self.txtEvaporatorTubesOD),  (self.cmdSelectEvaporatorTube),
                (self.lblEvaporatorTubesID),(self.txtEvaporatorTubesID),  ((0,0)),
                (self.lblEvaporatorTubesPl),(self.txtEvaporatorTubesPl),  ((0,0)),
                (self.lblEvaporatorTubesPt),(self.txtEvaporatorTubesPt),  ((0,0)),
                (self.lblEvaporatorTubesNcircuit),(self.txtEvaporatorTubesNcircuit),  (self.cmdShowEvaporatorCircuits),
            ])
        fgsEvaporatorInputs2.AddGrowableRow(2, 1)
        fgsEvaporatorInputs2.AddGrowableCol(1, 1)
        self.sizerEvaporatorInputs.Add(fgsEvaporatorInputs2,pos=(5,0))
        
        self.bmpEvaporatorTubes = wx.StaticBitmap(self.EvaporatorInputsPanel, -1, wx.Bitmap("imgs/Tubes.png", wx.BITMAP_TYPE_ANY))
        self.sizerEvaporatorInputs.Add(self.bmpEvaporatorTubes,pos=(5,4))
                
        self.lblEvaporatorAirVdot = wx.StaticText(self.EvaporatorInputsPanel, -1, "Vdot [m^3/s]")
        self.txtEvaporatorAirVdot = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "0.56319")
        self.lblEvaporatorAirTdb = wx.StaticText(self.EvaporatorInputsPanel, -1, "Tdb [K]")
        self.txtEvaporatorAirTdb = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "297.039")
        self.lblEvaporatorAirp = wx.StaticText(self.EvaporatorInputsPanel, -1, "p [kPa]")
        self.txtEvaporatorAirp = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "101.325")
        self.lblEvaporatorAirRH = wx.StaticText(self.EvaporatorInputsPanel, -1, "RH [-]")
        self.txtEvaporatorAirRH = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "0.5")
        self.lblEvaporatorPower = wx.StaticText(self.EvaporatorInputsPanel, -1, "Fan Power [W]")
        self.txtEvaporatorPower = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "260")
        self.lblEvaporatorAirVdot.SetMinSize((120, -1))
        self.lblEvaporatorAirVdot.SetToolTipString("Volumetric flow rate of air through indoor heat exchanger")
        self.txtEvaporatorAirVdot.SetToolTipString("Volumetric flow of humid air into the coil in cubic meters of air per second.")
        self.lblEvaporatorAirTdb.SetToolTipString("Dry bulb temperature")
        self.lblEvaporatorAirp.SetToolTipString("Ambient pressure (absolute)")
        self.lblEvaporatorAirRH.SetToolTipString("Relative humidity of the air from 0 to 1  (0% to 100% RH)")
        self.lblEvaporatorPower.SetToolTipString("Fan power of the blower")

        self.lblEvaporatorAir = wx.StaticText(self.EvaporatorInputsPanel, -1, "Air")
        self.lblEvaporatorAir.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerEvaporatorInputs.Add(self.lblEvaporatorAir,pos=(6,0))
        
        self.lineEvaporator3=wx.StaticLine(self.EvaporatorInputsPanel, -1)
        self.sizerEvaporatorInputs.Add(self.lineEvaporator3,pos=(7,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        fgsEvaporatorInputs3 = wx.FlexGridSizer(5, 2, 2, 2)
        fgsEvaporatorInputs3.AddMany([
                (self.lblEvaporatorAirVdot), (self.txtEvaporatorAirVdot),
                (self.lblEvaporatorAirTdb), (self.txtEvaporatorAirTdb),
                (self.lblEvaporatorAirp),(self.txtEvaporatorAirp),  
                (self.lblEvaporatorAirRH),(self.txtEvaporatorAirRH),
                (self.lblEvaporatorPower),(self.txtEvaporatorPower),
            ])
        fgsEvaporatorInputs3.AddGrowableRow(2, 1)
        fgsEvaporatorInputs3.AddGrowableCol(1, 1)
        self.sizerEvaporatorInputs.Add(fgsEvaporatorInputs3,pos=(8,0))
        

        self.lblEvaporatorRef = wx.StaticText(self.EvaporatorInputsPanel, -1, "Refrigerant")
        self.lblEvaporatorRef.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerEvaporatorInputs.Add(self.lblEvaporatorRef,pos=(6,4))

        self.lblEvaporatorDTsh = wx.StaticText(self.EvaporatorInputsPanel, -1, "Superheat [K]")
        self.txtEvaporatorDTsh = wx.TextCtrl(self.EvaporatorInputsPanel, -1, "3.5")
        
        fgsEvaporatorInputs4 = wx.FlexGridSizer(1,2, 2, 2)
        fgsEvaporatorInputs4.AddMany([(self.lblEvaporatorDTsh),( self.txtEvaporatorDTsh)])
        self.sizerEvaporatorInputs.Add(fgsEvaporatorInputs4,pos=(8,4))

        self.EvaporatorInputsPanel.SetSizer(self.sizerEvaporatorInputs)
    
    ###  Compressor
    def CompressorInputs(self):
        self.sizerCompressorInputs = wx.GridBagSizer()
        
        self.lblCompressor = wx.StaticText(self.CompressorInputsPanel, -1, "Compressor")
        self.lblCompressor.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerCompressorInputs.Add(self.lblCompressor,pos=(0,0))
        
        self.lineCompressor1=wx.StaticLine(self.CompressorInputsPanel, -1)
        self.sizerCompressorInputs.Add(self.lineCompressor1,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblCompfp = wx.StaticText(self.CompressorInputsPanel, -1, "Heat Loss Fraction")
        self.txtCompfp = wx.TextCtrl(self.CompressorInputsPanel, -1, "0.15")
        self.lblCompVdot_ratio = wx.StaticText(self.CompressorInputsPanel, -1, "Compressor Displacement \nscale factor")
        self.txtCompVdot_ratio = wx.TextCtrl(self.CompressorInputsPanel, -1, "1.2")
        self.lblCompfp.SetToolTipString("Fraction of the electrical power input which is lost to the ambient")
        self.lblCompVdot_ratio.SetToolTipString("Scale factor to the compressor displacement")
        
        fgsCompressorInputs1 = wx.FlexGridSizer(2,2, 2, 2)
        fgsCompressorInputs1.AddMany([
                (self.lblCompfp), (self.txtCompfp),
                (self.lblCompVdot_ratio), (self.txtCompVdot_ratio) 
            ])
        self.sizerCompressorInputs.Add(fgsCompressorInputs1,pos=(2,0))
        
        self.lblCompressor2 = wx.StaticText(self.CompressorInputsPanel, -1, "Compressor Model")
        self.lblCompressor2.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerCompressorInputs.Add(self.lblCompressor2,pos=(3,0))
        
        self.lineCompressor2=wx.StaticLine(self.CompressorInputsPanel, -1)
        self.sizerCompressorInputs.Add(self.lineCompressor2,pos=(4,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.optCompMapModel = wx.RadioButton(self.CompressorInputsPanel, -1, "Compressor Map")
        self.sizerCompressorInputs.Add(self.optCompMapModel,pos=(5,0),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblCompM1 = wx.StaticText(self.CompressorInputsPanel, -1, "M1")
        self.txtCompM1 = wx.TextCtrl(self.CompressorInputsPanel, -1, "30.25660922")
        self.lblCompP1 = wx.StaticText(self.CompressorInputsPanel, -1, "P1")
        self.txtCompP1 = wx.TextCtrl(self.CompressorInputsPanel, -1, "280.6815683")
        self.lblCompM2 = wx.StaticText(self.CompressorInputsPanel, -1, "M2")
        self.txtCompM2 = wx.TextCtrl(self.CompressorInputsPanel, -1, "-0.568382876")
        self.lblCompP2 = wx.StaticText(self.CompressorInputsPanel, -1, "P2")
        self.txtCompP2 = wx.TextCtrl(self.CompressorInputsPanel, -1, "-81.89665598")
        self.lblCompM3 = wx.StaticText(self.CompressorInputsPanel, -1, "M3")
        self.txtCompM3 = wx.TextCtrl(self.CompressorInputsPanel, -1, "2.831710808")
        self.lblCompP3 = wx.StaticText(self.CompressorInputsPanel, -1, "P3")
        self.txtCompP3 = wx.TextCtrl(self.CompressorInputsPanel, -1, "44.78234097")
        self.lblCompM4 = wx.StaticText(self.CompressorInputsPanel, -1, "M4")
        self.txtCompM4 = wx.TextCtrl(self.CompressorInputsPanel, -1, "0.103625357")
        self.lblCompP4 = wx.StaticText(self.CompressorInputsPanel, -1, "P4")
        self.txtCompP4 = wx.TextCtrl(self.CompressorInputsPanel, -1, "2.076283422")
        self.lblCompM5 = wx.StaticText(self.CompressorInputsPanel, -1, "M5")
        self.txtCompM5 = wx.TextCtrl(self.CompressorInputsPanel, -1, "0.015111817")
        self.lblCompP5 = wx.StaticText(self.CompressorInputsPanel, -1, "P5")
        self.txtCompP5 = wx.TextCtrl(self.CompressorInputsPanel, -1, "0.206102738")
        self.lblCompM6 = wx.StaticText(self.CompressorInputsPanel, -1, "M6")
        self.txtCompM6 = wx.TextCtrl(self.CompressorInputsPanel, -1, "-0.022000643")
        self.lblCompP6 = wx.StaticText(self.CompressorInputsPanel, -1, "P6")
        self.txtCompP6 = wx.TextCtrl(self.CompressorInputsPanel, -1, "-0.36128051")
        self.lblCompM7 = wx.StaticText(self.CompressorInputsPanel, -1, "M7")
        self.txtCompM7 = wx.TextCtrl(self.CompressorInputsPanel, -1, "-0.000851757")
        self.lblCompP7 = wx.StaticText(self.CompressorInputsPanel, -1, "P7")
        self.txtCompP7 = wx.TextCtrl(self.CompressorInputsPanel, -1, "-0.015962071")
        self.lblCompM8 = wx.StaticText(self.CompressorInputsPanel, -1, "M8")
        self.txtCompM8 = wx.TextCtrl(self.CompressorInputsPanel, -1, "1.79561E-05")
        self.lblCompP8 = wx.StaticText(self.CompressorInputsPanel, -1, "P8")
        self.txtCompP8 = wx.TextCtrl(self.CompressorInputsPanel, -1, "-0.002293473")
        self.lblCompM9 = wx.StaticText(self.CompressorInputsPanel, -1, "M9")
        self.txtCompM9 = wx.TextCtrl(self.CompressorInputsPanel, -1, "-6.2467E-05")
        self.lblCompP9 = wx.StaticText(self.CompressorInputsPanel, -1, "P9")
        self.txtCompP9 = wx.TextCtrl(self.CompressorInputsPanel, -1, "-0.000758934")
        self.lblCompM10 = wx.StaticText(self.CompressorInputsPanel, -1, "M10")
        self.txtCompM10 = wx.TextCtrl(self.CompressorInputsPanel, -1, "4.07877E-05")
        self.lblCompP10 = wx.StaticText(self.CompressorInputsPanel, -1, "P10")
        self.txtCompP10 = wx.TextCtrl(self.CompressorInputsPanel, -1, "0.001701369")
        self.cmdLoadCoeffs = wx.Button(self.CompressorInputsPanel, -1, "Load Coefficients")
        self.cmdARIInfo = wx.Button(self.CompressorInputsPanel, -1, "Info...")
        self.lblCompCoeffsFile = wx.StaticText(self.CompressorInputsPanel, -1, "File:")
        self.txtCompCoeffsFile = wx.TextCtrl(self.CompressorInputsPanel, -1, "", style=wx.TE_READONLY)
        
        self.optCompMapModel.SetValue(1)
        self.cmdLoadCoeffs.SetToolTipString("Load ARI coefficients from file")
        self.cmdARIInfo.SetToolTipString("More information about the ARI correlation")
        self.txtCompCoeffsFile.SetToolTipString("Currently loaded file")
        
        fgsCompressorInputs2 = wx.FlexGridSizer(13,4, 2, 2)
        fgsCompressorInputs2.AddMany([
                (self.lblCompM1), (self.txtCompM1), (self.lblCompP1), (self.txtCompP1),
                (self.lblCompM2), (self.txtCompM2), (self.lblCompP2), (self.txtCompP2),
                (self.lblCompM3), (self.txtCompM3), (self.lblCompP3), (self.txtCompP3),
                (self.lblCompM4), (self.txtCompM4), (self.lblCompP4), (self.txtCompP4),
                (self.lblCompM5), (self.txtCompM5), (self.lblCompP5), (self.txtCompP5),
                (self.lblCompM6), (self.txtCompM6), (self.lblCompP6), (self.txtCompP6),
                (self.lblCompM7), (self.txtCompM7), (self.lblCompP7), (self.txtCompP7),
                (self.lblCompM8), (self.txtCompM8), (self.lblCompP8), (self.txtCompP8),
                (self.lblCompM9), (self.txtCompM9), (self.lblCompP9), (self.txtCompP9),
                (self.lblCompM10), (self.txtCompM10), (self.lblCompP10), (self.txtCompP10),
                ((0,0)),(self.cmdLoadCoeffs) , ((0,0)), (self.cmdARIInfo),
                (self.lblCompCoeffsFile),(self.txtCompCoeffsFile,1,wx.EXPAND)
            ])
        fgsCompressorInputs2.AddGrowableRow(2, 1)
        fgsCompressorInputs2.AddGrowableCol(1, 1)
        self.sizerCompressorInputs.Add(fgsCompressorInputs2,pos=(5,1))
        
        self.lineCompressor3=wx.StaticLine(self.CompressorInputsPanel, -1)
        self.sizerCompressorInputs.Add(self.lineCompressor3,pos=(6,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.optCompEffModel = wx.RadioButton(self.CompressorInputsPanel, -1, "Efficiency")
        self.sizerCompressorInputs.Add(self.optCompEffModel,pos=(7,0),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblCompIsenEff = wx.StaticText(self.CompressorInputsPanel, -1, "Overall Isentropic Efficiency [-]")
        self.txtCompIsenEff = wx.TextCtrl(self.CompressorInputsPanel, -1, "0.6")
        self.lblCompVolEff = wx.StaticText(self.CompressorInputsPanel, -1, "Volumetric Efficiency [-]")
        self.txtCompVolEff = wx.TextCtrl(self.CompressorInputsPanel, -1, "0.9")
        self.lblCompVdot = wx.StaticText(self.CompressorInputsPanel, -1, "Displacement Rate [m^3/h]")
        self.txtCompVdot = wx.TextCtrl(self.CompressorInputsPanel, -1, "400")
        
        self.optCompEffModel.Enable(False)
        self.lblCompIsenEff.Enable(False)
        self.txtCompIsenEff.Enable(False)
        self.lblCompVolEff.Enable(False)
        self.txtCompVolEff.Enable(False)
        self.lblCompVdot.Enable(False)
        self.txtCompVdot.Enable(False)
        
        fgsCompressorInputs3 = wx.FlexGridSizer(3,2, 2, 2)
        fgsCompressorInputs3.AddMany([
            (self.lblCompIsenEff),(self.txtCompIsenEff),
            (self.lblCompVolEff),(self.txtCompVolEff),
            (self.lblCompVdot),(self.txtCompVdot)
            ])
        self.sizerCompressorInputs.Add(fgsCompressorInputs3,pos=(7,1))
        
        self.lineCompressor4=wx.StaticLine(self.CompressorInputsPanel, -1)
        self.sizerCompressorInputs.Add(self.lineCompressor4,pos=(8,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.CompressorInputsPanel.SetSizer(self.sizerCompressorInputs)
    def LineSetInputs(self):
        self.sizerLineSetInputs=wx.GridBagSizer()
        
        self.lblLineSet = wx.StaticText(self.LineSetInputsPanel, -1, "Line Set")
        self.lblLineSet.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerLineSetInputs.Add(self.lblLineSet,pos=(0,0))
        
        self.lineLineSet1=wx.StaticLine(self.LineSetInputsPanel, -1)
        self.sizerLineSetInputs.Add(self.lineLineSet1,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblLineSetL = wx.StaticText(self.LineSetInputsPanel, -1, "Line Length [m]")
        self.txtLineSetL = wx.TextCtrl(self.LineSetInputsPanel, -1, "5")
        self.lblLineSetOD_supply = wx.StaticText(self.LineSetInputsPanel, -1, "Supply[Liquid] Tube OD [m]")
        self.txtLineSetOD_supply = wx.TextCtrl(self.LineSetInputsPanel, -1, "0.04826")
        self.cmdLineSetSupplySelect = wx.Button(self.LineSetInputsPanel, -1, "Select")
        self.lblLineSetID_supply = wx.StaticText(self.LineSetInputsPanel, -1, "Supply[Liquid] Tube ID [m]")
        self.txtLineSetID_supply = wx.TextCtrl(self.LineSetInputsPanel, -1, "0.04089")
        self.lblLineSetOD_return = wx.StaticText(self.LineSetInputsPanel, -1, "Return[Suction] Tube OD [m]")
        self.txtLineSetOD_return = wx.TextCtrl(self.LineSetInputsPanel, -1, "0.04826")
        self.cmdLineSetReturnSelect = wx.Button(self.LineSetInputsPanel, -1, "Select")
        self.lblLineSetID_return = wx.StaticText(self.LineSetInputsPanel, -1, "Return[Suction] Tube ID [m]")
        self.txtLineSetID_return = wx.TextCtrl(self.LineSetInputsPanel, -1, "0.04089")
        self.lblLineSetTubek = wx.StaticText(self.LineSetInputsPanel, -1, "Tube Conductivity [W/m-K]")
        self.txtLineSetTubek = wx.TextCtrl(self.LineSetInputsPanel, -1, "0.19")
        self.lblLineSetInsult = wx.StaticText(self.LineSetInputsPanel, -1, "Insulation Thickness [m]")
        self.txtLineSetInsult = wx.TextCtrl(self.LineSetInputsPanel, -1, "0.02")
        self.lblLineSetInsulk = wx.StaticText(self.LineSetInputsPanel, -1, "Insulation Conductivity [W/m-K]")
        self.txtLineSetInsulk = wx.TextCtrl(self.LineSetInputsPanel, -1, "0.036")
        self.cmbLineSetInsulk = wx.ComboBox(self.LineSetInputsPanel, -1, choices=["Armaflex", "Polystyrene", "User Defined"], style=wx.CB_DROPDOWN)
        self.lblLineSeth_air = wx.StaticText(self.LineSetInputsPanel, -1, "Air Effective HTC [W/m^2-K]")
        self.txtLineSeth_air = wx.TextCtrl(self.LineSetInputsPanel, -1, "6")
        self.lblLineSetT_air = wx.StaticText(self.LineSetInputsPanel, -1, "Air Effective Temp [K]")
        self.txtLineSetT_air = wx.TextCtrl(self.LineSetInputsPanel, -1, "297")
        self.lblLineSetL.SetToolTipString("Length of the tube from the condensing unit to the coil installed indoors in the duct work.  ")
        self.lblLineSetOD_supply.SetToolTipString("Outer diameter of the liquid line for DX systems, supply line for secondary systems")
        self.lblLineSetID_supply.SetToolTipString("Inner diameter of the suction line for DX systems, return line for secondary systems")
        self.lblLineSetOD_return.SetToolTipString("Outer diameter of the suction line for DX systems, return line for secondary systems")
        self.lblLineSetID_return.SetToolTipString("Inner diameter of the suction line for DX systems, return line for secondary systems")
        self.lblLineSetTubek.SetToolTipString("Tube Conductivity")
        self.lblLineSetInsult.SetToolTipString("Thickness of the insulation")
        self.lblLineSetInsulk.SetToolTipString("Thermal conductivity of the insulation")
        self.cmbLineSetInsulk.SetSelection(0)
        self.lblLineSeth_air.SetToolTipString("Combined heat transfer coefficient for radiation and convection")
        self.lblLineSetT_air.SetToolTipString("Effective temperature of the environment that the line set interacts with thermally")
        
        fgsLineSetInputs1 = wx.FlexGridSizer(10,3, 2, 2)
        fgsLineSetInputs1.AddMany([
                                   
       (self.lblLineSetL),(self.txtLineSetL), ((0,0)),
        (self.lblLineSetOD_supply),(self.txtLineSetOD_supply),(self.cmdLineSetSupplySelect),
        (self.lblLineSetID_supply),(self.txtLineSetID_supply), ((0,0)),
        (self.lblLineSetOD_return),(self.txtLineSetOD_return),(self.cmdLineSetReturnSelect),
        (self.lblLineSetID_return),(self.txtLineSetID_return), ((0,0)),
        (self.lblLineSetTubek),(self.txtLineSetTubek), ((0,0)),
        (self.lblLineSetInsult),(self.txtLineSetInsult), ((0,0)),
        (self.lblLineSetInsulk),(self.txtLineSetInsulk),(self.cmbLineSetInsulk),
        (self.lblLineSeth_air),(self.txtLineSeth_air), ((0,0)),
        (self.lblLineSetT_air),(self.txtLineSetT_air), ((0,0)),
            ])
        self.sizerLineSetInputs.Add(fgsLineSetInputs1,pos=(2,0))
        
        self.bmpLineSet = wx.StaticBitmap(self.LineSetInputsPanel, -1, wx.Bitmap(os.path.join('imgs','LineSet.png'), wx.BITMAP_TYPE_ANY))
        self.sizerLineSetInputs.Add(self.bmpLineSet,pos=(2,1))
        
        self.LineSetInputsPanel.SetSizer(self.sizerLineSetInputs)
    def SolversInputs(self):
        self.sizerParaInputs = wx.GridBagSizer()
        
        self.radSolverSelect = wx.RadioBox(self.SolverMethodPanel, -1, "Solver Method", choices=["Specified Design Point ", "Parametric Study"], majorDimension=0, style=wx.RA_SPECIFY_ROWS)
        self.sizerParaInputs.Add(self.radSolverSelect,pos=(0,0),border=3)
        
        self.lblParaVariable = wx.StaticText(self.SolverMethodPanel, -1, "Variable")
        
        self.lblParaMinValue = wx.StaticText(self.SolverMethodPanel, -1, "Min Value / List")
        self.lblParaMaxValue = wx.StaticText(self.SolverMethodPanel, -1, "Max Value", style=wx.ALIGN_CENTRE)
        self.lblParaNstep = wx.StaticText(self.SolverMethodPanel, -1, " Number Steps", style=wx.ALIGN_CENTRE)
        self.chkParaVariable1 = wx.CheckBox(self.SolverMethodPanel, -1, "")
        self.cmbParaVariable1 = wx.ComboBox(self.SolverMethodPanel, -1, choices=["R290", "R410A", "R404A", "R134a"], style=wx.CB_DROPDOWN|wx.CB_READONLY)        
        self.txtParaVariable1Min = wx.TextCtrl(self.SolverMethodPanel, -1, "292.6")
        self.txtParaVariable1Max = wx.TextCtrl(self.SolverMethodPanel, -1, "324.8")
        self.txtParaVariable1N = wx.TextCtrl(self.SolverMethodPanel, -1, "6")
        self.chkParaVariable2 = wx.CheckBox(self.SolverMethodPanel, -1, "")
        self.cmbParaVariable2 = wx.ComboBox(self.SolverMethodPanel, -1, choices=["R290", "R410A", "R404A", "R134a"], style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.txtParaVariable2Min = wx.TextCtrl(self.SolverMethodPanel, -1, "290")
        self.txtParaVariable2Max = wx.TextCtrl(self.SolverMethodPanel, -1, "295")
        self.txtParaVariable2N = wx.TextCtrl(self.SolverMethodPanel, -1, "2")
        self.chkParaVariable3 = wx.CheckBox(self.SolverMethodPanel, -1, "")
        self.cmbParaVariable3 = wx.ComboBox(self.SolverMethodPanel, -1, choices=["R290", "R410A", "R404A", "R134a"], style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.txtParaVariable3Min = wx.TextCtrl(self.SolverMethodPanel, -1, "0.2146,0.3552,0.3851,0.5108,0.6832")
        self.txtParaVariable3Max = wx.TextCtrl(self.SolverMethodPanel, -1, "")
        self.txtParaVariable3N = wx.TextCtrl(self.SolverMethodPanel, -1, "L")
        self.chkParaVariable4 = wx.CheckBox(self.SolverMethodPanel, -1, "")
        self.cmbParaVariable4 = wx.ComboBox(self.SolverMethodPanel, -1, choices=["R290", "R410A", "R404A", "R134a"], style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.txtParaVariable4Min = wx.TextCtrl(self.SolverMethodPanel, -1, "0.2146,0.3552,0.3851,0.5108,0.6832")
        self.txtParaVariable4Max = wx.TextCtrl(self.SolverMethodPanel, -1, "")
        self.txtParaVariable4N = wx.TextCtrl(self.SolverMethodPanel, -1, "L")
        
        self.lblParaPath = wx.StaticText(self.SolverMethodPanel, -1, "Parametric study output filename:")
        self.txtParaPath = wx.TextCtrl(self.SolverMethodPanel, -1, "Para1.csv")
        self.cmdSelectParaPath = wx.Button(self.SolverMethodPanel, -1, "Select ...")
        
        #Load up the list of possible parametric variables
        paramsList=LoadParametricParams(os.path.join('parametric','params.txt'))
        
        self.ParamUnits={}
        self.ParamVariables={}
        self.ParamValues={}
        for combo in [self.cmbParaVariable1,self.cmbParaVariable2,self.cmbParaVariable3,self.cmbParaVariable4]:
            combo.Clear()
            for (desc,units,variable) in paramsList:
                combo.Append(desc)
                #Use description as the key in dictionaries
                self.ParamUnits[desc]=units
                self.ParamVariables[desc]=variable
        
        fgsParaInputs2 = wx.FlexGridSizer(7,5, 2, 2)
        fgsParaInputs2.AddMany([
            ((0,0)),    (self.lblParaVariable), (self.lblParaMinValue), (self.lblParaMaxValue), (self.lblParaNstep),
        (self.chkParaVariable1), (self.cmbParaVariable1), (self.txtParaVariable1Min),(self.txtParaVariable1Max),(self.txtParaVariable1N),
        (self.chkParaVariable2), (self.cmbParaVariable2), (self.txtParaVariable2Min),(self.txtParaVariable2Max),(self.txtParaVariable2N),
        (self.chkParaVariable3), (self.cmbParaVariable3), (self.txtParaVariable3Min),(self.txtParaVariable3Max),(self.txtParaVariable3N),
        (self.chkParaVariable4), (self.cmbParaVariable4), (self.txtParaVariable4Min),(self.txtParaVariable4Max),(self.txtParaVariable4N),
            ])
        self.sizerParaInputs.Add(fgsParaInputs2,pos=(2,0))
        
        fgsParaInputs3 = wx.FlexGridSizer(1,3, 2, 2)
        fgsParaInputs3.AddMany([
            (self.lblParaPath),(self.txtParaPath),(self.cmdSelectParaPath)
            ])
        self.sizerParaInputs.Add(fgsParaInputs3,pos=(4,0))

        self.SolverMethodPanel.SetSizer(self.sizerParaInputs)
    def MainCycleOutputs(self):
        self.sizerMainCycleOutputs=wx.GridBagSizer()
        
        self.lblMainCycleOutputs = wx.StaticText(self.MainCycleOutputsPanel, -1, "Main Outputs")
        self.lblMainCycleOutputs.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerMainCycleOutputs.Add(self.lblMainCycleOutputs,pos=(0,0))
        
        self.lineMainCycleOutputs=wx.StaticLine(self.MainCycleOutputsPanel, -1)
        self.sizerMainCycleOutputs.Add(self.lineMainCycleOutputs,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblCycleCOP = wx.StaticText(self.MainCycleOutputsPanel, -1, "COP (Q/W)[-]")
        self.txtCycleCOP = wx.TextCtrl(self.MainCycleOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCycleCOPeff = wx.StaticText(self.MainCycleOutputsPanel, -1, "COP w/ fans,pump [-]")
        self.txtCycleCOPeff = wx.TextCtrl(self.MainCycleOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCycleSHR = wx.StaticText(self.MainCycleOutputsPanel, -1, "Sensible Heat Ratio [-]")
        self.txtCycleSHR = wx.TextCtrl(self.MainCycleOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCycleChargeOutput = wx.StaticText(self.MainCycleOutputsPanel, -1, "Charge [kg]")
        self.txtCycleChargeOutput = wx.TextCtrl(self.MainCycleOutputsPanel, -1, "")
        self.lblCycleTsatCond = wx.StaticText(self.MainCycleOutputsPanel, -1, "Condensing Temp [K/C]")
        self.txtCycleTsatCond = wx.TextCtrl(self.MainCycleOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCycleTsatEvap = wx.StaticText(self.MainCycleOutputsPanel, -1, "Evaporation Temp [K/C]")
        self.txtCycleTsatEvap = wx.TextCtrl(self.MainCycleOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCycleMainFreezeProtection = wx.StaticText(self.MainCycleOutputsPanel, -1, "Freeze protection [K/C]")
        self.txtCycleMainFreezeProtection = wx.TextCtrl(self.MainCycleOutputsPanel, -1, "", style=wx.TE_READONLY)
#        self.cmdOpenSchematic = wx.Button(self.MainCycleOutputsPanel, -1, "Detailed Schematic")
        self.cmdWriteRun = wx.Button(self.MainCycleOutputsPanel, -1, "Write To File...")
        
        fgsMainCycleOutputs = wx.FlexGridSizer(8,2, 2, 2)
        fgsMainCycleOutputs.AddMany([
            (self.lblCycleCOP),(self.txtCycleCOP),
            (self.lblCycleCOPeff),(self.txtCycleCOPeff),
            (self.lblCycleSHR),(self.txtCycleSHR),
            (self.lblCycleChargeOutput),(self.txtCycleChargeOutput),
            (self.lblCycleTsatCond),(self.txtCycleTsatCond),
            (self.lblCycleTsatEvap),(self.txtCycleTsatEvap),
            (self.lblCycleMainFreezeProtection),(self.txtCycleMainFreezeProtection),
            ((0,0)),(self.cmdWriteRun),                
        ])
        self.sizerMainCycleOutputs.Add(fgsMainCycleOutputs,pos=(2,0))
        
        self.pnlSchematic = PlotCycleOutputsPanel(self.MainCycleOutputsPanel, -1)
        self.sizerMainCycleOutputs.Add(self.pnlSchematic,pos=(3,0))
    
        self.MainCycleOutputsPanel.SetSizer(self.sizerMainCycleOutputs)
    
    def PlotsOutputs(self):
        fgsTsPlot = wx.BoxSizer(wx.VERTICAL)
        self.pltTs = MPLPanel(self.TsPlotPanel, -1)
        fgsTsPlot.Add((self.pltTs))
        self.TsPlotPanel.SetSizer(fgsTsPlot)
        
        fgsPhPlot = wx.BoxSizer(wx.VERTICAL)
        self.pltPh = MPLPanel(self.phPlotPanel, -1)
        fgsPhPlot.Add((self.pltPh))
        self.phPlotPanel.SetSizer(fgsPhPlot)
    def CompressorOutputs(self):
        self.sizerCompressorOutputs=wx.GridBagSizer()
        
        self.lblCompressor = wx.StaticText(self.CompressorOutputsPanel, -1, "Compressor")
        self.lblCompressor.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerCompressorOutputs.Add(self.lblLineSet,pos=(0,0))
        
        self.lineLineSet1=wx.StaticLine(self.CompressorOutputsPanel, -1)
        self.sizerCompressorOutputs.Add(self.lineLineSet1,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblCompressormdot = wx.StaticText(self.CompressorOutputsPanel, -1, "Mass flow rate [kg/s]")
        self.txtCompressormdot = wx.TextCtrl(self.CompressorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCompressorPower = wx.StaticText(self.CompressorOutputsPanel, -1, "Power [W]")
        self.txtCompressorPower = wx.TextCtrl(self.CompressorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCompressorTin = wx.StaticText(self.CompressorOutputsPanel, -1, "Tin [K / C]")
        self.txtCompressorTin = wx.TextCtrl(self.CompressorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCompressorTout = wx.StaticText(self.CompressorOutputsPanel, -1, "Tout [K / C]")
        self.txtCompressorTout = wx.TextCtrl(self.CompressorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCompressorhin = wx.StaticText(self.CompressorOutputsPanel, -1, "Enthalpy in [J/kg]")
        self.txtCompressorhin = wx.TextCtrl(self.CompressorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCompressorhout = wx.StaticText(self.CompressorOutputsPanel, -1, "Enthalpy out [J/kg]")
        self.txtCompressorhout = wx.TextCtrl(self.CompressorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCompressoretaoi = wx.StaticText(self.CompressorOutputsPanel, -1, "Overall Isentropic Efficiency [-]")
        self.txtCompressoretaoi = wx.TextCtrl(self.CompressorOutputsPanel, -1, "", style=wx.TE_READONLY)
        
        fgsParaInputs3 = wx.FlexGridSizer(cols = 2, vgap = 2, hgap = 2)
        fgsParaInputs3.AddMany([
            (self.lblCompressormdot),(self.txtCompressormdot),
        (self.lblCompressorPower),(self.txtCompressorPower),
        (self.lblCompressorTin),(self.txtCompressorTin),
        (self.lblCompressorTout),(self.txtCompressorTout),
        (self.lblCompressorhin),(self.txtCompressorhin),
        (self.lblCompressorhout),(self.txtCompressorhout),
        (self.lblCompressoretaoi),(self.txtCompressoretaoi)
            ])
        self.sizerCompressorOutputs.Add(fgsParaInputs3,pos=(2,0))
        
        self.CompressorOutputsPanel.SetSizer(self.sizerCompressorOutputs)
    def LineSetOutputs(self):
        self.sizerLineSetOutputs=wx.GridBagSizer()
        
        self.lblSupplyLineSet = wx.StaticText(self.LineSetOutputsPanel, -1, "Supply")
        self.lblSupplyLineSet.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerLineSetOutputs.Add(self.lblSupplyLineSet,pos=(0,0))
        
        self.lineLineSetOutputs1=wx.StaticLine(self.LineSetOutputsPanel, -1)
        self.sizerLineSetOutputs.Add(self.lineLineSetOutputs1,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblLineSetSupplyOutputsQ = wx.StaticText(self.LineSetOutputsPanel, -1, "Heat Transfer [W]")
        self.txtLineSetSupplyOutputsQ = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetSupplyOutputsDP = wx.StaticText(self.LineSetOutputsPanel, -1, "Pressure Drop [Pa]")
        self.txtLineSetSupplyOutputsDP = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetSupplyOutputsRe = wx.StaticText(self.LineSetOutputsPanel, -1, "Reynolds Number [-]")
        self.txtLineSetSupplyOutputsRe = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetSupplyOutputsh = wx.StaticText(self.LineSetOutputsPanel, -1, "Mean HTC [W/m^2-K]")
        self.txtLineSetSupplyOutputsh = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetSupplyOutputsTin = wx.StaticText(self.LineSetOutputsPanel, -1, "Inlet Temp [K/C]")
        self.txtLineSetSupplyOutputsTin = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetSupplyOutputsTout = wx.StaticText(self.LineSetOutputsPanel, -1, "OutletTemp [K/C]")
        self.txtLineSetSupplyOutputsTout = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetSupplyOutputsCharge = wx.StaticText(self.LineSetOutputsPanel, -1, "Charge [kg]")
        self.txtLineSetSupplyOutputsCharge = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        
        fgsParaInputs1 = wx.FlexGridSizer(7,2, 2, 2)
        fgsParaInputs1.AddMany([
            (self.lblLineSetSupplyOutputsQ),(self.txtLineSetSupplyOutputsQ),
        (self.lblLineSetSupplyOutputsDP),(self.txtLineSetSupplyOutputsDP),
        (self.lblLineSetSupplyOutputsRe),(self.txtLineSetSupplyOutputsRe),
        (self.lblLineSetSupplyOutputsh),(self.txtLineSetSupplyOutputsh),
        (self.lblLineSetSupplyOutputsTin),(self.txtLineSetSupplyOutputsTin),
        (self.lblLineSetSupplyOutputsTout),(self.txtLineSetSupplyOutputsTout),
        (self.lblLineSetSupplyOutputsCharge),(self.txtLineSetSupplyOutputsCharge),
            ])
        self.sizerLineSetOutputs.Add(fgsParaInputs1,pos=(3,0))
        
        self.lblReturnLineSet = wx.StaticText(self.LineSetOutputsPanel, -1, "Return")
        self.lblReturnLineSet.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerLineSetOutputs.Add(self.lblReturnLineSet,pos=(0,1))
        
        self.lblLineSetReturnOutputsQ = wx.StaticText(self.LineSetOutputsPanel, -1, "Heat Transfer [W]")
        self.txtLineSetReturnOutputsQ = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetReturnOutputsDP = wx.StaticText(self.LineSetOutputsPanel, -1, "Pressure Drop [Pa]")
        self.txtLineSetReturnOutputsDP = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetReturnOutputsRe = wx.StaticText(self.LineSetOutputsPanel, -1, "Reynolds Number [-]")
        self.txtLineSetReturnOutputsRe = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetReturnOutputsh = wx.StaticText(self.LineSetOutputsPanel, -1, "Mean HTC [W/m^2-K]")
        self.txtLineSetReturnOutputsh = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetReturnOutputsTin = wx.StaticText(self.LineSetOutputsPanel, -1, "Inlet Temp [K/C]")
        self.txtLineSetReturnOutputsTin = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetReturnOutputsTout = wx.StaticText(self.LineSetOutputsPanel, -1, "Outlet Temp [K/C]")
        self.txtLineSetReturnOutputsTout = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblLineSetReturnOutputsCharge = wx.StaticText(self.LineSetOutputsPanel, -1, "Charge [kg]")
        self.txtLineSetReturnOutputsCharge = wx.TextCtrl(self.LineSetOutputsPanel, -1, "", style=wx.TE_READONLY)
        
        fgsParaInputs2 = wx.FlexGridSizer(7,2, 2, 2)
        fgsParaInputs2.AddMany([
            (self.lblLineSetReturnOutputsQ),(self.txtLineSetReturnOutputsQ),
            (self.lblLineSetReturnOutputsDP),(self.txtLineSetReturnOutputsDP),
            (self.lblLineSetReturnOutputsRe),(self.txtLineSetReturnOutputsRe),
            (self.lblLineSetReturnOutputsh),(self.txtLineSetReturnOutputsh),
            (self.lblLineSetReturnOutputsTin),(self.txtLineSetReturnOutputsTin),
            (self.lblLineSetReturnOutputsTout),(self.txtLineSetReturnOutputsTout),
            (self.lblLineSetReturnOutputsCharge),(self.txtLineSetReturnOutputsCharge),
            ])
        self.sizerLineSetOutputs.Add(fgsParaInputs2,pos=(3,1))
        
        self.LineSetOutputsPanel.SetSizer(self.sizerLineSetOutputs)
    
    def CondenserOutputs(self):
        self.sizerCondenserOutputs=wx.GridBagSizer()
        
        self.lblCondenserOutputs = wx.StaticText(self.CondenserOutputsPanel, -1, "Refrigerant-Side")
        self.lblCondenserOutputs.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerCondenserOutputs.Add(self.lblCondenserOutputs,pos=(0,0))
        
        self.lineCondenserOutputs=wx.StaticLine(self.CondenserOutputsPanel, -1)
        self.sizerCondenserOutputs.Add(self.lineCondenserOutputs,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblCondenserQ = wx.StaticText(self.CondenserOutputsPanel, -1, "Heat Transfer [W]")
        self.txtCondenserQ = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserQ_superheat = wx.StaticText(self.CondenserOutputsPanel, -1, "--> Heat Transfer Superheated [W]")
        self.txtCondenserQ_superheat = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserQ_2phase = wx.StaticText(self.CondenserOutputsPanel, -1, "--> Heat Transfer Two-Phase [W]")
        self.txtCondenserQ_2phase = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserQ_subcool = wx.StaticText(self.CondenserOutputsPanel, -1, "--> Heat Transfer Subcooled [W]")
        self.txtCondenserQ_subcool = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserCharge = wx.StaticText(self.CondenserOutputsPanel, -1, "Charge [kg]")
        self.txtCondenserCharge = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserCharge_superheat = wx.StaticText(self.CondenserOutputsPanel, -1, "--> Charge Superheated [kg]")
        self.txtCondenserCharge_superheat = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserCharge_2phase = wx.StaticText(self.CondenserOutputsPanel, -1, "--> Charge Two-Phase [kg]")
        self.txtCondenserCharge_2phase = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserCharge_subcool = wx.StaticText(self.CondenserOutputsPanel, -1, "--> Charge Subcooled [kg]")
        self.txtCondenserCharge_subcool = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserh_superheat = wx.StaticText(self.CondenserOutputsPanel, -1, "Ref. Mean HTC Superheat [W/m^2-K]")
        self.txtCondenserh_superheat = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserh_2phase = wx.StaticText(self.CondenserOutputsPanel, -1, "Ref. Mean HTC Two-Phase [W/m^2-K]")
        self.txtCondenserh_2phase = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserh_subcool = wx.StaticText(self.CondenserOutputsPanel, -1, "Ref. Mean HTC Subcooled [W/m^2-K]")
        self.txtCondenserh_subcool = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserw_superheat = wx.StaticText(self.CondenserOutputsPanel, -1, "Volume Fraction Superheat [-]")
        self.txtCondenserw_superheat = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserw_2phase = wx.StaticText(self.CondenserOutputsPanel, -1, "Volume Fraction Two-Phase[-]")
        self.txtCondenserw_2phase = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserw_subcool = wx.StaticText(self.CondenserOutputsPanel, -1, "Volume Fraction Subcooled [-]")
        self.txtCondenserw_subcool = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserTin_r = wx.StaticText(self.CondenserOutputsPanel, -1, "Inlet Temperature [K/C]")
        self.txtCondenserTin_r = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserTout_r = wx.StaticText(self.CondenserOutputsPanel, -1, "Outlet Temperature [K/C]")
        self.txtCondenserTout_r = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserDP_r = wx.StaticText(self.CondenserOutputsPanel, -1, "Pressure Drop Total [Pa]")
        self.txtCondenserDP_r = wx.TextCtrl(self.CondenserOutputsPanel, -1, "")
        self.lblCondenserDP_r_superheat = wx.StaticText(self.CondenserOutputsPanel, -1, "Pressure Drop Superheat [Pa]")
        self.txtCondenserDP_r_superheat = wx.TextCtrl(self.CondenserOutputsPanel, -1, "")
        self.lblCondenserDP_r_2phase = wx.StaticText(self.CondenserOutputsPanel, -1, "Pressure Drop Two-Phase [Pa]")
        self.txtCondenserDP_r_2phase = wx.TextCtrl(self.CondenserOutputsPanel, -1, "")
        self.lblCondenserDP_r_subcool = wx.StaticText(self.CondenserOutputsPanel, -1, "Pressure Drop Subcooled [Pa]")
        self.txtCondenserDP_r_subcool = wx.TextCtrl(self.CondenserOutputsPanel, -1, "")
        
        fgsCondenserOutputs1 = wx.FlexGridSizer(20,2, 2, 2)
        fgsCondenserOutputs1.AddMany([
            (self.lblCondenserQ),(self.txtCondenserQ),
            (self.lblCondenserQ_superheat),(self.txtCondenserQ_superheat),
            (self.lblCondenserQ_2phase),(self.txtCondenserQ_2phase),
            (self.lblCondenserQ_subcool),(self.txtCondenserQ_subcool),
            (self.lblCondenserCharge),(self.txtCondenserCharge),
            (self.lblCondenserCharge_superheat),(self.txtCondenserCharge_superheat),
            (self.lblCondenserCharge_2phase),(self.txtCondenserCharge_2phase),
            (self.lblCondenserCharge_subcool),(self.txtCondenserCharge_subcool),
            (self.lblCondenserh_superheat),(self.txtCondenserh_superheat),
            (self.lblCondenserh_2phase),(self.txtCondenserh_2phase),
            (self.lblCondenserh_subcool),(self.txtCondenserh_subcool),
            (self.lblCondenserw_superheat),(self.txtCondenserw_superheat),
            (self.lblCondenserw_2phase),(self.txtCondenserw_2phase),
            (self.lblCondenserw_subcool),(self.txtCondenserw_subcool),
            (self.lblCondenserTin_r),(self.txtCondenserTin_r),
            (self.lblCondenserTout_r),(self.txtCondenserTout_r),
            (self.lblCondenserDP_r),(self.txtCondenserDP_r),
            (self.lblCondenserDP_r_superheat),(self.txtCondenserDP_r_superheat),
            (self.lblCondenserDP_r_2phase),(self.txtCondenserDP_r_2phase),
            (self.lblCondenserDP_r_subcool),(self.txtCondenserDP_r_subcool),
            ])
        self.sizerCondenserOutputs.Add(fgsCondenserOutputs1,pos=(2,0))
        
        self.lblCondenserOutputs2 = wx.StaticText(self.CondenserOutputsPanel, -1, "Air-Side")
        self.lblCondenserOutputs2.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerCondenserOutputs.Add(self.lblCondenserOutputs2,pos=(0,1))
        
        self.lblCondenserAirh_a = wx.StaticText(self.CondenserOutputsPanel, -1, "Air Mean HTC [W/m^2]")
        self.txtCondenserAirh_a = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserAirA_a = wx.StaticText(self.CondenserOutputsPanel, -1, "Total air-side surface area [m^2]")
        self.txtCondenserAirA_a = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserAirmdot_a = wx.StaticText(self.CondenserOutputsPanel, -1, "Air mass flow rate [kg/s]")
        self.txtCondenserAirmdot_a = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserAireta_a = wx.StaticText(self.CondenserOutputsPanel, -1, "External surface efficiency [-]")
        self.txtCondenserAireta_a = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCondenserAirdP_a = wx.StaticText(self.CondenserOutputsPanel, -1, "Air-side Pressure drop [Pa]")
        self.txtCondenserAirdP_a = wx.TextCtrl(self.CondenserOutputsPanel, -1, "", style=wx.TE_READONLY)            
        
        fgsCondenserOutputs2 = wx.FlexGridSizer(5,2, 2, 2)
        fgsCondenserOutputs2.AddMany([
            (self.lblCondenserAirh_a),(self.txtCondenserAirh_a),
            (self.lblCondenserAirA_a),(self.txtCondenserAirA_a),
            (self.lblCondenserAirmdot_a),(self.txtCondenserAirmdot_a),
            (self.lblCondenserAireta_a),(self.txtCondenserAireta_a),
            (self.lblCondenserAirdP_a),(self.txtCondenserAirdP_a),
            ])
        self.sizerCondenserOutputs.Add(fgsCondenserOutputs2,pos=(2,1))
        
        self.CondenserOutputsPanel.SetSizer(self.sizerCondenserOutputs)
        
    def CoolingCoilOutputs(self):
        self.sizerCoolingCoilOutputs=wx.GridBagSizer()
        
        self.lblCoolingCoilOutputsGlycol = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Glycol-Side")
        self.lblCoolingCoilOutputsGlycol.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerCoolingCoilOutputs.Add(self.lblCoolingCoilOutputsGlycol,pos=(0,0))
        
        self.lineCoolingCoilOutputs=wx.StaticLine(self.CoolingCoilOutputsPanel, -1)
        self.sizerCoolingCoilOutputs.Add(self.lineCoolingCoilOutputs,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblCoolingCoilQ = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Heat Transfer [W]")
        self.txtCoolingCoilQ = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilh_g = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Mean HTC [W/m^2-K]")
        self.txtCoolingCoilh_g = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilRe_g = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Mean Re [-]")
        self.txtCoolingCoilRe_g = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilDP_g = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Pressure drop [Pa]")
        self.txtCoolingCoilDP_g = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilTin_g = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Inlet Temperature [K/C]")
        self.txtCoolingCoilTin_g = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilTout_g = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Outlet Temperature [K/C]")
        self.txtCoolingCoilTout_g = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
    
        
        fgsCoolingCoil1 = wx.FlexGridSizer(20,2, 2, 2)
        fgsCoolingCoil1.AddMany([
            (self.lblCoolingCoilQ),(self.txtCoolingCoilQ),
            (self.lblCoolingCoilh_g),(self.txtCoolingCoilh_g),
            (self.lblCoolingCoilRe_g),(self.txtCoolingCoilRe_g),
            (self.lblCoolingCoilDP_g),(self.txtCoolingCoilDP_g),
            (self.lblCoolingCoilTin_g),(self.txtCoolingCoilTin_g),
            (self.lblCoolingCoilTout_g),(self.txtCoolingCoilTout_g),
            ])
        self.sizerCoolingCoilOutputs.Add(fgsCoolingCoil1,pos=(2,0))
        
        self.lblCoolingCoilOutputsAir = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Air-Side")
        self.lblCoolingCoilOutputsAir.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerCoolingCoilOutputs.Add(self.lblCoolingCoilOutputsAir,pos=(0,1))
        
        self.lblCoolingCoilAirh_a = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Mean HTC [W/m^2]")
        self.txtCoolingCoilAirh_a = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilAirA_a = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Total surface area [m^2]")
        self.txtCoolingCoilAirA_a = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilAirmdot_a = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Air mass flow rate [kg/s]")
        self.txtCoolingCoilAirmdot_a = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilAireta_a = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Surface efficiency [-]")
        self.txtCoolingCoilAireta_a = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilAirdP_a = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Air-side Pressure drop [Pa]")
        self.txtCoolingCoilAirdP_a = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilAirSHR = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "SHR [-]")
        self.txtCoolingCoilAirSHR = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilAirf_dry = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Dry Fraction [-]")
        self.txtCoolingCoilAirf_dry = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblCoolingCoilAirTout = wx.StaticText(self.CoolingCoilOutputsPanel, -1, "Air Outlet Temperature [K/C]")
        self.txtCoolingCoilAirTout = wx.TextCtrl(self.CoolingCoilOutputsPanel, -1, "", style=wx.TE_READONLY)
        
        fgsCoolingCoil2 = wx.FlexGridSizer(8,2, 2, 2)
        fgsCoolingCoil2.AddMany([
            (self.lblCoolingCoilAirh_a),(self.txtCoolingCoilAirh_a),
            (self.lblCoolingCoilAirA_a),(self.txtCoolingCoilAirA_a),
            (self.lblCoolingCoilAirmdot_a),(self.txtCoolingCoilAirmdot_a),
            (self.lblCoolingCoilAireta_a),(self.txtCoolingCoilAireta_a),
            (self.lblCoolingCoilAirdP_a),(self.txtCoolingCoilAirdP_a),
            (self.lblCoolingCoilAirSHR),(self.txtCoolingCoilAirSHR),
            (self.lblCoolingCoilAirf_dry),(self.txtCoolingCoilAirf_dry),
            (self.lblCoolingCoilAirTout),(self.txtCoolingCoilAirTout),
            ])
        self.sizerCoolingCoilOutputs.Add(fgsCoolingCoil2,pos=(2,1))
        
        self.CoolingCoilOutputsPanel.SetSizer(self.sizerCoolingCoilOutputs)
        
    def EvaporatorOutputs(self):
        self.sizerEvaporatorOutputs=wx.GridBagSizer()
        
        self.lblEvaporatorOutputsRef = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Refrigerant-Side")
        self.lblEvaporatorOutputsRef.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerEvaporatorOutputs.Add(self.lblEvaporatorOutputsRef,pos=(0,0))
        
        self.lineEvaporatorOutputs=wx.StaticLine(self.EvaporatorOutputsPanel, -1)
        self.sizerEvaporatorOutputs.Add(self.lineEvaporatorOutputs,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblEvaporatorQ = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Heat Transfer [W]")
        self.txtEvaporatorQ = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorQ_superheat = wx.StaticText(self.EvaporatorOutputsPanel, -1, "--> Heat Transfer Superheated [W]")
        self.txtEvaporatorQ_superheat = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorQ_2phase = wx.StaticText(self.EvaporatorOutputsPanel, -1, "--> Heat Transfer Two-Phase [W]")
        self.txtEvaporatorQ_2phase = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorQ_subcool = wx.StaticText(self.EvaporatorOutputsPanel, -1, "--> Heat Transfer Subcooled [W]")
        self.txtEvaporatorQ_subcool = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorCharge = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Charge [kg]")
        self.txtEvaporatorCharge = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorCharge_superheat = wx.StaticText(self.EvaporatorOutputsPanel, -1, "--> Charge Superheated [kg]")
        self.txtEvaporatorCharge_superheat = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorCharge_2phase = wx.StaticText(self.EvaporatorOutputsPanel, -1, "--> Charge Two-Phase [kg]")
        self.txtEvaporatorCharge_2phase = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorCharge_subcool = wx.StaticText(self.EvaporatorOutputsPanel, -1, "--> Charge Subcooled [kg]")
        self.txtEvaporatorCharge_subcool = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorh_r_superheat = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Mean HTC Superheat [W/m^2-K]")
        self.txtEvaporatorh_r_superheat = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorh_r_2phase = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Mean HTC Two-Phase [W/m^2-K]")
        self.txtEvaporatorh_r_2phase = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorh_r_subcool = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Mean HTC Subcooled [W/m^2-K]")
        self.txtEvaporatorh_r_subcool = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorw_superheat = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Volume Fraction Superheat [-]")
        self.txtEvaporatorw_superheat = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorw_2phase = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Volume Fraction Two-Phase[-]")
        self.txtEvaporatorw_2phase = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorw_subcool = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Volume Fraction Subcooled [-]")
        self.txtEvaporatorw_subcool = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorDP_r = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Pressure Drop Total [Pa]")
        self.txtEvaporatorDP_r = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "")
        self.lblEvaporatorDP_r_superheat = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Pressure Drop Superheat [Pa]")
        self.txtEvaporatorDP_r_superheat = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "")
        self.lblEvaporatorDP_r_2phase = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Pressure Drop Two-Phase [Pa]")
        self.txtEvaporatorDP_r_2phase = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "")
        self.lblEvaporatorDP_r_subcool = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Pressure Drop Subcooled [Pa]")
        self.txtEvaporatorDP_r_subcool = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "")
        
        fgsEvaporator1 = wx.FlexGridSizer(20,2, 2, 2)
        fgsEvaporator1.AddMany([
            (self.lblEvaporatorQ),(self.txtEvaporatorQ),
            (self.lblEvaporatorQ_superheat),(self.txtEvaporatorQ_superheat),
            (self.lblEvaporatorQ_2phase),(self.txtEvaporatorQ_2phase),
            (self.lblEvaporatorQ_subcool),(self.txtEvaporatorQ_subcool),
            (self.lblEvaporatorCharge),(self.txtEvaporatorCharge),
            (self.lblEvaporatorCharge_superheat),(self.txtEvaporatorCharge_superheat),
            (self.lblEvaporatorCharge_2phase),(self.txtEvaporatorCharge_2phase),
            (self.lblEvaporatorCharge_subcool),(self.txtEvaporatorCharge_subcool),
            (self.lblEvaporatorh_r_superheat),(self.txtEvaporatorh_r_superheat),
            (self.lblEvaporatorh_r_2phase),(self.txtEvaporatorh_r_2phase),
            (self.lblEvaporatorh_r_subcool),(self.txtEvaporatorh_r_subcool),
            (self.lblEvaporatorw_superheat),(self.txtEvaporatorw_superheat),
            (self.lblEvaporatorw_2phase),(self.txtEvaporatorw_2phase),
            (self.lblEvaporatorw_subcool),(self.txtEvaporatorw_subcool),
            (self.lblEvaporatorDP_r),(self.txtEvaporatorDP_r),
            (self.lblEvaporatorDP_r_superheat),(self.txtEvaporatorDP_r_superheat),
            (self.lblEvaporatorDP_r_2phase),(self.txtEvaporatorDP_r_2phase),
            (self.lblEvaporatorDP_r_subcool),(self.txtEvaporatorDP_r_subcool),
            ])
        self.sizerEvaporatorOutputs.Add(fgsEvaporator1,pos=(2,0))
        
        self.lblEvaporatorOutputsAir = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Air-Side")
        self.lblEvaporatorOutputsAir.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerEvaporatorOutputs.Add(self.lblEvaporatorOutputsAir,pos=(0,1))
        
        self.lblEvaporatorh_a = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Mean HTC [W/m^2-K]")
        self.txtEvaporatorh_a = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorA_a = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Total Area [m^2]")
        self.txtEvaporatorA_a = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatormdot_a = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Air mass flow rate [kg/s]")
        self.txtEvaporatormdot_a = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatoreta_a = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Surface Efficiency [-]")
        self.txtEvaporatoreta_a = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatordP_a = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Air-side pressure drop [Pa]")
        self.txtEvaporatordP_a = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorSHR = wx.StaticText(self.EvaporatorOutputsPanel, -1, "SHR [-]")
        self.txtEvaporatorSHR = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorfdry_2phase = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Dry fraction (2phase) [-]")
        self.txtEvaporatorfdry_2phase = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorfdry_superheat = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Dry fraction (superheat) [-]")
        self.txtEvaporatorfdry_superheat = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblEvaporatorTout_a = wx.StaticText(self.EvaporatorOutputsPanel, -1, "Outlet Air Temperature [K/C]")
        self.txtEvaporatorTout_a = wx.TextCtrl(self.EvaporatorOutputsPanel, -1, "", style=wx.TE_READONLY)
        
        fgsEvaporator2 = wx.FlexGridSizer(cols = 2, hgap = 2, vgap = 2)
        fgsEvaporator2.AddMany([
            (self.lblEvaporatorh_a),(self.txtEvaporatorh_a),
            (self.lblEvaporatorA_a),(self.txtEvaporatorA_a),
            (self.lblEvaporatormdot_a),(self.txtEvaporatormdot_a),
            (self.lblEvaporatoreta_a),(self.txtEvaporatoreta_a),
            (self.lblEvaporatordP_a),(self.txtEvaporatordP_a),
            (self.lblEvaporatorSHR),(self.txtEvaporatorSHR),
            (self.lblEvaporatorfdry_2phase),(self.txtEvaporatorfdry_2phase),
            (self.lblEvaporatorfdry_superheat),(self.txtEvaporatorfdry_superheat),
            (self.lblEvaporatorTout_a),(self.txtEvaporatorTout_a),
            ])
        self.sizerEvaporatorOutputs.Add(fgsEvaporator2,pos=(2,1))
        
    
        self.EvaporatorOutputsPanel.SetSizer(self.sizerEvaporatorOutputs)
    def PumpIHXOutputs(self):    
        
        self.sizerPumpIHXOutputs=wx.GridBagSizer()
        
        self.lblPumpIHXRef = wx.StaticText(self.PumpIHXOutputsPanel, -1, "IHX Refrigerant-Side")
        self.lblPumpIHXRef.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerPumpIHXOutputs.Add(self.lblPumpIHXRef,pos=(0,0))
        
        self.linePumpIHXOutputs=wx.StaticLine(self.PumpIHXOutputsPanel, -1)
        self.sizerPumpIHXOutputs.Add(self.linePumpIHXOutputs,pos=(1,0),span=(1,5),border=3,flag=wx.EXPAND|wx.BOTTOM)
        
        self.lblIHXQ = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Heat Transfer [W]")
        self.txtIHXQ = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXQ_superheat = wx.StaticText(self.PumpIHXOutputsPanel, -1, "--> Heat Transfer Superheated [W]")
        self.txtIHXQ_superheat = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXQ_2phase = wx.StaticText(self.PumpIHXOutputsPanel, -1, "--> Heat Transfer Two-Phase [W]")
        self.txtIHXQ_2phase = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXQ_subcool = wx.StaticText(self.PumpIHXOutputsPanel, -1, "--> Heat Transfer Subcooled [W]")
        self.txtIHXQ_subcool = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXCharge = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Charge [kg]")
        self.txtIHXCharge = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXCharge_superheat = wx.StaticText(self.PumpIHXOutputsPanel, -1, "--> Charge Superheated [kg]")
        self.txtIHXCharge_superheat = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXCharge_2phase = wx.StaticText(self.PumpIHXOutputsPanel, -1, "--> Charge Two-Phase [kg]")
        self.txtIHXCharge_2phase = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXCharge_subcool = wx.StaticText(self.PumpIHXOutputsPanel, -1, "--> Charge Subcooled [kg]")
        self.txtIHXCharge_subcool = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXh_r_superheat = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Mean HTC Superheat [W/m^2-K]")
        self.txtIHXh_r_superheat = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXh_r_2phase = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Mean HTC Two-Phase [W/m^2-K]")
        self.txtIHXh_r_2phase = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXh_r_subcool = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Mean HTC Subcooled [W/m^2-K]")
        self.txtIHXh_r_subcool = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXw_superheat = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Volume Fraction Superheat [-]")
        self.txtIHXw_superheat = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXw_2phase = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Volume Fraction Two-Phase [-]")
        self.txtIHXw_2phase = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXw_subcool = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Volume Fraction Subcooled [-]")
        self.txtIHXw_subcool = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXTin_r = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Inlet Temperature [K/C]")
        self.txtIHXTin_r = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXTout_r = wx.StaticText(self.PumpIHXOutputsPanel, -1, "OutletTemperature [K/C]")
        self.txtIHXTout_r = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXDP_r = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Pressure Drop Total [Pa]")
        self.txtIHXDP_r = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "")
        self.lblIHXDP_r_superheat = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Pressure Drop Superheat [Pa]")
        self.txtIHXDP_r_superheat = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "")
        self.lblIHXDP_r_2phase = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Pressure Drop Two-Phase [Pa]")
        self.txtIHXDP_r_2phase = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "")
        self.lblIHXDP_r_subcool = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Pressure Drop Subcooled [Pa]")
        self.txtIHXDP_r_subcool = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "")
        
        fgsPumpIHX1 = wx.FlexGridSizer(20,2, 2, 2)
        fgsPumpIHX1.AddMany([
            (self.lblIHXQ),(self.txtIHXQ),
            (self.lblIHXQ_superheat),(self.txtIHXQ_superheat),
            (self.lblIHXQ_2phase),(self.txtIHXQ_2phase),
            (self.lblIHXQ_subcool),(self.txtIHXQ_subcool),
            (self.lblIHXCharge),(self.txtIHXCharge),
            (self.lblIHXCharge_superheat),(self.txtIHXCharge_superheat),
            (self.lblIHXCharge_2phase),(self.txtIHXCharge_2phase),
            (self.lblIHXCharge_subcool),(self.txtIHXCharge_subcool),
            (self.lblIHXh_r_superheat),(self.txtIHXh_r_superheat),
            (self.lblIHXh_r_2phase),(self.txtIHXh_r_2phase),
            (self.lblIHXh_r_subcool),(self.txtIHXh_r_subcool),
            (self.lblIHXw_superheat),(self.txtIHXw_superheat),
            (self.lblIHXw_2phase),(self.txtIHXw_2phase),
            (self.lblIHXw_subcool),(self.txtIHXw_subcool),
            (self.lblIHXTin_r),(self.txtIHXTin_r),
            (self.lblIHXTout_r),(self.txtIHXTout_r),
            (self.lblIHXDP_r),(self.txtIHXDP_r),
            (self.lblIHXDP_r_superheat),(self.txtIHXDP_r_superheat),
            (self.lblIHXDP_r_2phase),(self.txtIHXDP_r_2phase),
            (self.lblIHXDP_r_subcool),(self.txtIHXDP_r_subcool),
            ])
        self.sizerPumpIHXOutputs.Add(fgsPumpIHX1,pos=(2,0))
        
        self.lblPumpIHXGly = wx.StaticText(self.PumpIHXOutputsPanel, -1, "IHX Glycol-Side")
        self.lblPumpIHXGly.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerPumpIHXOutputs.Add(self.lblPumpIHXGly,pos=(0,1))
        
        self.lblIHXh_g = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Mean HTC [W/m^2-K]")
        self.txtIHXh_g = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXRe_g = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Mean Reynolds # [-]")
        self.txtIHXRe_g = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXDP_g = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Pressure Drop [Pa]")
        self.txtIHXDP_g = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXTin_g = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Inlet Temperature [K/C]")
        self.txtIHXTin_g = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblIHXTout_g = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Outlet Temperature [K/C]")
        self.txtIHXTout_g = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        
        
        fgsPumpIHX2 = wx.FlexGridSizer(8,2, 2, 2)
        fgsPumpIHX2.AddMany([
                (self.lblIHXh_g),(self.txtIHXh_g),
                (self.lblIHXRe_g),(self.txtIHXRe_g),
                (self.lblIHXDP_g),(self.txtIHXDP_g),
                (self.lblIHXTin_g),(self.txtIHXTin_g),
                (self.lblIHXTout_g),(self.txtIHXTout_g),
            ])
        self.sizerPumpIHXOutputs.Add(fgsPumpIHX2,pos=(2,1))
        
        self.lblPumpIHXPump = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Pump")
        self.lblPumpIHXPump.SetFont(wx.Font(-1, wx.DEFAULT, wx.NORMAL, wx.BOLD, 0, ""))
        self.sizerPumpIHXOutputs.Add(self.lblPumpIHXPump,pos=(0,2))
    
        self.lblPumpDP = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Pressure Lift [Pa]")
        self.txtPumpDP = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        self.lblPumpPower = wx.StaticText(self.PumpIHXOutputsPanel, -1, "Power [W]")
        self.txtPumpPower = wx.TextCtrl(self.PumpIHXOutputsPanel, -1, "", style=wx.TE_READONLY)
        
        fgsPumpIHX3 = wx.FlexGridSizer(2,2, 2, 2)
        fgsPumpIHX3.AddMany([
                (self.lblPumpDP),(self.txtPumpDP),
                (self.lblPumpPower),(self.txtPumpPower),
            ])
        self.sizerPumpIHXOutputs.Add(fgsPumpIHX3,pos=(2,2))
        
        self.PumpIHXOutputsPanel.SetSizer(self.sizerPumpIHXOutputs)

    MenuBar(self)
    Panels(self)
    CycleInputs(self)
    PumpIHXInputs(self)
    CoolingCoilInputs(self)
    CondenserInputs(self)
#    EvaporatorInputs(self)
    CompressorInputs(self)
    LineSetInputs(self)
    SolversInputs(self)
    MainCycleOutputs(self)
    PlotsOutputs(self)
    CompressorOutputs(self)
    LineSetOutputs(self)
    CondenserOutputs(self)
    CoolingCoilOutputs(self)
    EvaporatorOutputs(self)
    PumpIHXOutputs(self)
    
def BindEvents(self):
    
    self.Bind(wx.EVT_MENU, self.WriteConfigFile, self.FileSave)
    self.Bind(wx.EVT_MENU, self.ReadConfigFile, self.FileOpen)
    self.Bind(wx.EVT_MENU, self.FileQuit, self.menuFileQuit)
    self.Bind(wx.EVT_MENU, self.openHelp, self.HelpHelp)
    self.Bind(wx.EVT_MENU, self.RunCode, self.SolveSolve)
    self.Bind(wx.EVT_RADIOBUTTON, self.ChangeImposedVariable, self.optCycleCharge)
    self.Bind(wx.EVT_RADIOBUTTON, self.ChangeImposedVariable, self.optCycleSubcooling)
    self.Bind(wx.EVT_RADIOBUTTON, self.ChangeSystemType, self.optSecFluid)
    self.Bind(wx.EVT_RADIOBUTTON, self.ChangeSystemType, self.optCycleDX)
    self.Bind(wx.EVT_BUTTON, self.SelectLineSetSupplyTube, self.cmdLineSetSupplySelect)
    self.Bind(wx.EVT_BUTTON, self.SelectLineSetReturnTube, self.cmdLineSetReturnSelect)
    self.Bind(wx.EVT_TEXT, self.modifyLineSetInsulation, self.txtLineSetInsulk)
    self.Bind(wx.EVT_COMBOBOX, self.SelectLineSetInsulation, self.cmbLineSetInsulk)
    self.Bind(wx.EVT_TEXT, self.modifyCoolingCoilFink, self.txtCoolingCoilFink)
    self.Bind(wx.EVT_COMBOBOX, self.CoolingCoilFinSelect, self.cmbCoolingCoilFink)
    self.Bind(wx.EVT_BUTTON, self.SelectCoolingCoilTube, self.cmdSelectCoolingCoilTube)
    self.Bind(wx.EVT_BUTTON, self.ShowCoolingCoilCircuits, self.cmdShowCoolingCoilCircuits)
    self.Bind(wx.EVT_BUTTON, self.ShowCondenserCircuits, self.cmdShowCondenserCircuits)
    self.Bind(wx.EVT_RADIOBUTTON, self.ChangeCompressorModel, self.optCompMapModel)
    self.Bind(wx.EVT_BUTTON, self.LoadCoeffs, self.cmdLoadCoeffs)
    self.Bind(wx.EVT_BUTTON, self.ARIMapInfo, self.cmdARIInfo)
    self.Bind(wx.EVT_RADIOBUTTON, self.ChangeCompressorModel, self.optCompEffModel)
    self.Bind(wx.EVT_RADIOBOX, self.ChangeSolverMethod, self.radSolverSelect)
#    self.Bind(wx.EVT_RADIOBOX, self.ChangeMode, self.radCycleMode)
    self.Bind(wx.EVT_RADIOBOX, self.ChangeSystemType, self.radCycleMode)
    self.Bind(wx.EVT_BUTTON, self.SelectParaPath, self.cmdSelectParaPath)
#    self.Bind(wx.EVT_BUTTON, self.openSchematic, self.cmdOpenSchematic)
    self.Bind(wx.EVT_BUTTON, self.writeRun, self.cmdWriteRun)
    
    # end wxGlade
    
    self.ChangeSystemType()
    self.ReadConfigFile(file=os.path.join('configs','Default.cfg'))