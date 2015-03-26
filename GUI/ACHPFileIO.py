from CoolProp.CoolProp import Props
from ACHP.Cycle import DXCycleClass,SecondaryCycleClass
from math import pi

def setField(Cycle,fieldName,value):
    """
    
    """
    S=Cycle
    fields=fieldName.split('.')
    while len(fields)>1:
        #Update current item
        S=getattr(S,fields[0])
        #Update field
        fieldName=fieldName.partition('.')[2] #Keep only that to the right of first '.'
        #Update the fields list
        fields=fieldName.split('.')
        print type(S),fieldName,len(fieldName.split('.'))
    
    print type(S),fields[0],value
    setattr(S,fields[0],value)
                
def GUI2DXCycleInputs(GUI):
    
    #Create the basic structure of the DXCycle
    Cycle=DXCycleClass()
            
    M=[0,0,0,0,0,0,0,0,0,0] #Empty list
    M[0]=float(GUI.txtCompM1.GetValue())
    M[1]=float(GUI.txtCompM2.GetValue())
    M[2]=float(GUI.txtCompM3.GetValue())
    M[3]=float(GUI.txtCompM4.GetValue())
    M[4]=float(GUI.txtCompM5.GetValue())
    M[5]=float(GUI.txtCompM6.GetValue())
    M[6]=float(GUI.txtCompM7.GetValue())
    M[7]=float(GUI.txtCompM8.GetValue())
    M[8]=float(GUI.txtCompM9.GetValue())
    M[9]=float(GUI.txtCompM10.GetValue())
    
    P=[0,0,0,0,0,0,0,0,0,0] #Empty list
    P[0]=float(GUI.txtCompP1.GetValue())
    P[1]=float(GUI.txtCompP2.GetValue())
    P[2]=float(GUI.txtCompP3.GetValue())
    P[3]=float(GUI.txtCompP4.GetValue())
    P[4]=float(GUI.txtCompP5.GetValue())
    P[5]=float(GUI.txtCompP6.GetValue())
    P[6]=float(GUI.txtCompP7.GetValue())
    P[7]=float(GUI.txtCompP8.GetValue())
    P[8]=float(GUI.txtCompP9.GetValue())
    P[9]=float(GUI.txtCompP10.GetValue())
    
    Cycle.Compressor.M=M
    Cycle.Compressor.P=P

    Cycle.Compressor.fp                 =float(GUI.txtCompfp.GetValue())
    Cycle.Compressor.Vdot_ratio         =float(GUI.txtCompVdot_ratio.GetValue())
    Cycle.Evaporator.DT_sh              =float(GUI.txtCycleDTsh.GetValue())
    
    if GUI.radCycleMode.GetStringSelection()=='Cooling Mode':
        Cycle.Condenser.Fins.Air.Vdot_ha    =float(GUI.txtCondenserAirVdot.GetValue())
        Cycle.Condenser.Fins.Air.Tdb        =float(GUI.txtCondenserAirTdb.GetValue())
        Cycle.Condenser.Fins.Air.p          =float(GUI.txtCondenserAirp.GetValue())
        Cycle.Condenser.Fins.Air.RH         =float(GUI.txtCondenserAirRH.GetValue())
        Cycle.Condenser.Fins.Air.Tmean      =float(GUI.txtCondenserAirTdb.GetValue())
        Cycle.Condenser.Fins.Air.RHmean     =float(GUI.txtCondenserAirRH.GetValue())
        
        Cycle.Condenser.Fins.Tubes.NTubes_per_bank   =int(GUI.txtCondenserTubesNtubes.GetValue())
        Cycle.Condenser.Fins.Tubes.Nbank    =int(GUI.txtCondenserTubesNbank.GetValue())
        Cycle.Condenser.Fins.Tubes.Ltube    =float(GUI.txtCondenserTubesL.GetValue())
        Cycle.Condenser.Fins.Tubes.OD       =float(GUI.txtCondenserTubesOD.GetValue())
        Cycle.Condenser.Fins.Tubes.ID       =float(GUI.txtCondenserTubesID.GetValue())
        Cycle.Condenser.Fins.Tubes.Pl       =float(GUI.txtCondenserTubesPl.GetValue())
        Cycle.Condenser.Fins.Tubes.Pt       =float(GUI.txtCondenserTubesPt.GetValue())
        Cycle.Condenser.Fins.Tubes.Ncircuits =int(GUI.txtCondenserTubesNcircuit.GetValue())
        Cycle.Condenser.Fins.Fins.FPI       =float(GUI.txtCondenserFinFPI.GetValue())
        Cycle.Condenser.Fins.Fins.Pd        =float(GUI.txtCondenserFinpd.GetValue())
        Cycle.Condenser.Fins.Fins.xf        =float(GUI.txtCondenserFinxf.GetValue())
        Cycle.Condenser.Fins.Fins.t         =float(GUI.txtCondenserFint.GetValue())
        Cycle.Condenser.Fins.Fins.k_fin     =float(GUI.txtCondenserFink.GetValue())
        Cycle.Condenser.Fins.Air.FanPower   =float(GUI.txtCondenserPower.GetValue())
        Cycle.Condenser.Verbosity           =0
        
        Cycle.Evaporator.Fins.Air.Vdot_ha   =float(GUI.txtCoolingCoilAirVdot.GetValue())
        Cycle.Evaporator.Fins.Air.Tdb       =float(GUI.txtCoolingCoilAirTdb.GetValue())
        Cycle.Evaporator.Fins.Air.p         =float(GUI.txtCoolingCoilAirp.GetValue())
        Cycle.Evaporator.Fins.Air.RH        =float(GUI.txtCoolingCoilAirRH.GetValue())
        Cycle.Evaporator.Fins.Air.Tmean      =float(GUI.txtCoolingCoilAirTdb.GetValue())
        Cycle.Evaporator.Fins.Air.RHmean     =float(GUI.txtCoolingCoilAirRH.GetValue())
        Cycle.Evaporator.Fins.Tubes.NTubes_per_bank  =float(GUI.txtCoolingCoilTubesNtubes.GetValue())
        Cycle.Evaporator.Fins.Tubes.Nbank   =float(GUI.txtCoolingCoilTubesNbank.GetValue())
        Cycle.Evaporator.Fins.Tubes.Ltube   =float(GUI.txtCoolingCoilTubesL.GetValue())
        Cycle.Evaporator.Fins.Tubes.OD      =float(GUI.txtCoolingCoilTubesOD.GetValue())
        Cycle.Evaporator.Fins.Tubes.ID      =float(GUI.txtCoolingCoilTubesID.GetValue())
        Cycle.Evaporator.Fins.Tubes.Pl      =float(GUI.txtCoolingCoilTubesPl.GetValue())
        Cycle.Evaporator.Fins.Tubes.Pt      =float(GUI.txtCoolingCoilTubesPt.GetValue())
        Cycle.Evaporator.Fins.Tubes.Ncircuits=float(GUI.txtCoolingCoilTubesNcircuit.GetValue())
        Cycle.Evaporator.Fins.Fins.FPI      =float(GUI.txtCoolingCoilFinFPI.GetValue())
        Cycle.Evaporator.Fins.Fins.Pd       =float(GUI.txtCoolingCoilFinpd.GetValue())
        Cycle.Evaporator.Fins.Fins.xf       =float(GUI.txtCoolingCoilFinxf.GetValue())
        Cycle.Evaporator.Fins.Fins.t        =float(GUI.txtCoolingCoilFint.GetValue())
        Cycle.Evaporator.Fins.Fins.k_fin    =float(GUI.txtCoolingCoilFink.GetValue())
        Cycle.Evaporator.Fins.Air.FanPower  =float(GUI.txtCoolingCoilPower.GetValue())
        Cycle.Evaporator.DT_sh              =float(GUI.txtCycleDTsh.GetValue())
        Cycle.Evaporator.Verbosity           =0
        
    elif GUI.radCycleMode.GetStringSelection()=='Heating Mode':
        Cycle.Evaporator.Fins.Air.Vdot_ha    =float(GUI.txtCondenserAirVdot.GetValue())
        Cycle.Evaporator.Fins.Air.Tdb        =float(GUI.txtCondenserAirTdb.GetValue())
        Cycle.Evaporator.Fins.Air.p          =float(GUI.txtCondenserAirp.GetValue())
        Cycle.Evaporator.Fins.Air.RH         =float(GUI.txtCondenserAirRH.GetValue())
        Cycle.Evaporator.Fins.Air.Tmean      =float(GUI.txtCondenserAirTdb.GetValue())
        Cycle.Evaporator.Fins.Air.RHmean     =float(GUI.txtCondenserAirRH.GetValue())
        
        Cycle.Evaporator.Fins.Tubes.NTubes_per_bank   =int(GUI.txtCondenserTubesNtubes.GetValue())
        Cycle.Evaporator.Fins.Tubes.Nbank    =int(GUI.txtCondenserTubesNbank.GetValue())
        Cycle.Evaporator.Fins.Tubes.Ltube    =float(GUI.txtCondenserTubesL.GetValue())
        Cycle.Evaporator.Fins.Tubes.OD       =float(GUI.txtCondenserTubesOD.GetValue())
        Cycle.Evaporator.Fins.Tubes.ID       =float(GUI.txtCondenserTubesID.GetValue())
        Cycle.Evaporator.Fins.Tubes.Pl       =float(GUI.txtCondenserTubesPl.GetValue())
        Cycle.Evaporator.Fins.Tubes.Pt       =float(GUI.txtCondenserTubesPt.GetValue())
        Cycle.Evaporator.Fins.Tubes.Ncircuits =int(GUI.txtCondenserTubesNcircuit.GetValue())
        Cycle.Evaporator.Fins.Fins.FPI       =float(GUI.txtCondenserFinFPI.GetValue())
        Cycle.Evaporator.Fins.Fins.Pd        =float(GUI.txtCondenserFinpd.GetValue())
        Cycle.Evaporator.Fins.Fins.xf        =float(GUI.txtCondenserFinxf.GetValue())
        Cycle.Evaporator.Fins.Fins.t         =float(GUI.txtCondenserFint.GetValue())
        Cycle.Evaporator.Fins.Fins.k_fin     =float(GUI.txtCondenserFink.GetValue())
        Cycle.Evaporator.Fins.Air.FanPower            =float(GUI.txtCondenserPower.GetValue())
        Cycle.Evaporator.Verbosity           =0
        
        Cycle.Condenser.Fins.Air.Vdot_ha   =float(GUI.txtCoolingCoilAirVdot.GetValue())
        Cycle.Condenser.Fins.Air.Tdb       =float(GUI.txtCoolingCoilAirTdb.GetValue())
        Cycle.Condenser.Fins.Air.p         =float(GUI.txtCoolingCoilAirp.GetValue())
        Cycle.Condenser.Fins.Air.RH        =float(GUI.txtCoolingCoilAirRH.GetValue())
        Cycle.Condenser.Fins.Air.Tmean      =float(GUI.txtCoolingCoilAirTdb.GetValue())
        Cycle.Condenser.Fins.Air.RHmean     =float(GUI.txtCoolingCoilAirRH.GetValue())
        Cycle.Condenser.Fins.Tubes.NTubes_per_bank  =float(GUI.txtCoolingCoilTubesNtubes.GetValue())
        Cycle.Condenser.Fins.Tubes.Nbank   =float(GUI.txtCoolingCoilTubesNbank.GetValue())
        Cycle.Condenser.Fins.Tubes.Ltube   =float(GUI.txtCoolingCoilTubesL.GetValue())
        Cycle.Condenser.Fins.Tubes.OD      =float(GUI.txtCoolingCoilTubesOD.GetValue())
        Cycle.Condenser.Fins.Tubes.ID      =float(GUI.txtCoolingCoilTubesID.GetValue())
        Cycle.Condenser.Fins.Tubes.Pl      =float(GUI.txtCoolingCoilTubesPl.GetValue())
        Cycle.Condenser.Fins.Tubes.Pt      =float(GUI.txtCoolingCoilTubesPt.GetValue())
        Cycle.Condenser.Fins.Tubes.Ncircuits=float(GUI.txtCoolingCoilTubesNcircuit.GetValue())
        Cycle.Condenser.Fins.Fins.FPI      =float(GUI.txtCoolingCoilFinFPI.GetValue())
        Cycle.Condenser.Fins.Fins.Pd       =float(GUI.txtCoolingCoilFinpd.GetValue())
        Cycle.Condenser.Fins.Fins.xf       =float(GUI.txtCoolingCoilFinxf.GetValue())
        Cycle.Condenser.Fins.Fins.t        =float(GUI.txtCoolingCoilFint.GetValue())
        Cycle.Condenser.Fins.Fins.k_fin    =float(GUI.txtCoolingCoilFink.GetValue())
        Cycle.Condenser.Fins.Air.FanPower  =float(GUI.txtCoolingCoilPower.GetValue())
        Cycle.Condenser.Verbosity           =0
        
    
    Cycle.LineSetSupply.L                 =float(GUI.txtLineSetL.GetValue())
    Cycle.LineSetSupply.OD                =float(GUI.txtLineSetOD_supply.GetValue())
    Cycle.LineSetSupply.ID                =float(GUI.txtLineSetID_supply.GetValue())
    Cycle.LineSetSupply.t_insul           =float(GUI.txtLineSetInsult.GetValue())
    Cycle.LineSetSupply.k_tube            =float(GUI.txtLineSetTubek.GetValue())
    Cycle.LineSetSupply.k_insul           =float(GUI.txtLineSetInsulk.GetValue())
    Cycle.LineSetSupply.h_air             =float(GUI.txtLineSeth_air.GetValue())
    Cycle.LineSetSupply.T_air             =float(GUI.txtLineSetT_air.GetValue())
    
    Cycle.LineSetReturn.L                 =float(GUI.txtLineSetL.GetValue())
    Cycle.LineSetReturn.OD                =float(GUI.txtLineSetOD_return.GetValue())
    Cycle.LineSetReturn.ID                =float(GUI.txtLineSetID_return.GetValue())
    Cycle.LineSetReturn.t_insul           =float(GUI.txtLineSetInsult.GetValue())
    Cycle.LineSetReturn.k_tube            =float(GUI.txtLineSetTubek.GetValue())
    Cycle.LineSetReturn.k_insul           =float(GUI.txtLineSetInsulk.GetValue())
    Cycle.LineSetReturn.h_air             =float(GUI.txtLineSeth_air.GetValue())
    Cycle.LineSetReturn.T_air             =float(GUI.txtLineSetT_air.GetValue())
    
    Cycle.Charge_target=float(GUI.txtCycleCharge.GetValue())
    Cycle.DT_sc_target=float(GUI.txtCycleSubcooling.GetValue())
    Cycle.Mode='AC'
    
    Cycle.Ref=str(GUI.cmbRefrigerant.GetValue())
    Cycle.Verbosity = 10
    Cycle.ParaPath = GUI.txtParaPath.GetValue()
        
    return Cycle

def CycleOutputs2GUI(GUI,Cycle):
    
    ## These outputs are common to both DX and Secondary Loop systems
    GUI.txtCycleCOP.SetValue("%0.3f " % Cycle.COP )
    GUI.txtCycleCOPeff.SetValue("%0.3f " % Cycle.COSP )
    GUI.txtCycleChargeOutput.SetValue("%0.3f " % Cycle.Charge )
    GUI.txtCycleTsatCond.SetValue("%0.2f / %0.2f" % (Cycle.Tdew_cond,Cycle.Tdew_cond-273.15) )
    GUI.txtCycleTsatEvap.SetValue("%0.2f / %0.2f" % (Cycle.Tdew_evap,Cycle.Tdew_evap-273.15) )
    
    GUI.txtCompressorhin.SetValue("%0.2f" % (Cycle.Compressor.hin_r ))
    GUI.txtCompressorhout.SetValue("%0.2f" % (Cycle.Compressor.hout_r ))
    GUI.txtCompressoretaoi.SetValue("%0.4f" % (Cycle.Compressor.eta_oi ))
    GUI.txtCompressorPower.SetValue("%0.2f" % (Cycle.Compressor.W ))
    GUI.txtCompressormdot.SetValue("%0.6f" % (Cycle.Compressor.mdot_r ))
    GUI.txtCompressorTin.SetValue("%0.2f / %0.2f" % (Cycle.Compressor.Tin_r,Cycle.Compressor.Tin_r-273.15) )
    GUI.txtCompressorTout.SetValue("%0.2f / %0.2f" % (Cycle.Compressor.Tout_r,Cycle.Compressor.Tout_r-273.15) )
        
    if Cycle.Mode=='AC':
        GUI.txtCondenserQ.SetValue("%0.2f" % (Cycle.Condenser.Q ))
        GUI.txtCondenserQ_superheat.SetValue("%0.2f" % (Cycle.Condenser.Q_superheat ))
        GUI.txtCondenserQ_2phase.SetValue("%0.2f" % (Cycle.Condenser.Q_2phase ))
        GUI.txtCondenserQ_subcool.SetValue("%0.2f" % (Cycle.Condenser.Q_subcool ))
        GUI.txtCondenserCharge.SetValue("%0.5f" % (Cycle.Condenser.Charge ))
        GUI.txtCondenserCharge_superheat.SetValue("%0.5f" % (Cycle.Condenser.Charge_superheat ))
        GUI.txtCondenserCharge_2phase.SetValue("%0.5f" % (Cycle.Condenser.Charge_2phase ))
        GUI.txtCondenserCharge_subcool.SetValue("%0.5f" % (Cycle.Condenser.Charge_subcool ))
        GUI.txtCondenserAirh_a.SetValue("%0.2f" % (Cycle.Condenser.Fins.h_a ))
        GUI.txtCondenserAireta_a.SetValue("%0.4f" % (Cycle.Condenser.Fins.eta_a ))
        GUI.txtCondenserAirA_a.SetValue("%0.2f" % (Cycle.Condenser.Fins.A_a ))
        GUI.txtCondenserAirmdot_a.SetValue("%0.4f" % (Cycle.Condenser.Fins.mdot_da ))
        GUI.txtCondenserAirdP_a.SetValue("%0.4f" % (Cycle.Condenser.Fins.dP_a ))
        GUI.txtCondenserh_superheat.SetValue("%0.2f" % (Cycle.Condenser.h_r_superheat ))
        GUI.txtCondenserh_2phase.SetValue("%0.2f" % (Cycle.Condenser.h_r_2phase ))
        GUI.txtCondenserh_subcool.SetValue("%0.2f" % (Cycle.Condenser.h_r_subcool ))
        GUI.txtCondenserw_superheat.SetValue("%0.3f" % (Cycle.Condenser.w_superheat ))
        GUI.txtCondenserw_2phase.SetValue("%0.3f" % (Cycle.Condenser.w_2phase ))
        GUI.txtCondenserw_subcool.SetValue("%0.3f" % (Cycle.Condenser.w_subcool ))
        GUI.txtCondenserTin_r.SetValue("%0.2f / %0.2f" % (Cycle.Condenser.Tin_r,Cycle.Condenser.Tin_r-273.15 ))
        GUI.txtCondenserTout_r.SetValue("%0.2f / %0.2f" % (Cycle.Condenser.Tout_r,Cycle.Condenser.Tout_r-273.15 ))
        GUI.txtCondenserDP_r.SetValue("%0.3f" % (Cycle.Condenser.DP_r ))
        GUI.txtCondenserDP_r_superheat.SetValue("%0.3f" % (Cycle.Condenser.DP_r_superheat ))
        GUI.txtCondenserDP_r_2phase.SetValue("%0.3f" % (Cycle.Condenser.DP_r_2phase ))
        GUI.txtCondenserDP_r_subcool.SetValue("%0.3f" % (Cycle.Condenser.DP_r_subcool ))
    elif Cycle.Mode=='HP':
        GUI.txtEvaporatorQ.SetValue("%0.2f" % (Cycle.Evaporator.Q ))
        GUI.txtEvaporatorQ_superheat.SetValue("%0.2f" % (Cycle.Evaporator.Q_superheat ))
        GUI.txtEvaporatorQ_2phase.SetValue("%0.2f" % (Cycle.Evaporator.Q_2phase ))
        GUI.txtEvaporatorQ_subcool.SetValue("N/A")
        GUI.txtEvaporatorCharge.SetValue("%0.5f" % (Cycle.Evaporator.Charge ))
        GUI.txtEvaporatorCharge_superheat.SetValue("%0.5f" % (Cycle.Evaporator.Charge_superheat ))
        GUI.txtEvaporatorCharge_2phase.SetValue("%0.5f" % (Cycle.Evaporator.Charge_2phase ))
        GUI.txtEvaporatorCharge_subcool.SetValue("N/A" )
        GUI.txtEvaporatorh_r_superheat.SetValue("%0.2f" % (Cycle.Evaporator.h_r_superheat ))
        GUI.txtEvaporatorh_r_2phase.SetValue("%0.2f" % (Cycle.Evaporator.h_r_2phase ))
        GUI.txtEvaporatorh_r_subcool.SetValue("N/A")
        GUI.txtEvaporatorw_superheat.SetValue("%0.3f" % (Cycle.Evaporator.w_superheat ))
        GUI.txtEvaporatorw_2phase.SetValue("%0.3f" % (Cycle.Evaporator.w_2phase ))
        GUI.txtEvaporatorw_subcool.SetValue("N/A" )
        GUI.txtEvaporatorh_a.SetValue("%0.3f" % (Cycle.Evaporator.Fins.h_a ))
        GUI.txtEvaporatoreta_a.SetValue("%0.3f" % (Cycle.Evaporator.Fins.eta_a ))
        GUI.txtEvaporatordP_a.SetValue("%0.3f" % (Cycle.Evaporator.Fins.dP_a ))
        GUI.txtEvaporatormdot_a.SetValue("%0.3f" % (Cycle.Evaporator.Fins.mdot_da ))
        GUI.txtEvaporatorA_a.SetValue("%0.3f" % (Cycle.Evaporator.Fins.A_a ))
        GUI.txtEvaporatorSHR.SetValue("%0.3f" % (Cycle.Evaporator.SHR ))
        GUI.txtEvaporatorfdry_2phase.SetValue("%0.3f" % (Cycle.Evaporator.fdry_2phase ))
        GUI.txtEvaporatorfdry_superheat.SetValue("%0.3f" % (Cycle.Evaporator.fdry_superheat ))
        GUI.txtEvaporatorTout_a.SetValue("%0.2f / %0.2f" % (Cycle.Evaporator.Tout_a,Cycle.Evaporator.Tout_a-273.15 ))
        GUI.txtEvaporatorDP_r.SetValue("%0.3f" % (Cycle.Evaporator.DP_r ))
        GUI.txtEvaporatorDP_r_superheat.SetValue("%0.3f" % (Cycle.Evaporator.DP_r_superheat ))
        GUI.txtEvaporatorDP_r_2phase.SetValue("%0.3f" % (Cycle.Evaporator.DP_r_2phase ))
        GUI.txtEvaporatorDP_r_subcool.SetValue("N/A" )
        
    
    GUI.txtLineSetSupplyOutputsQ.SetValue("%0.3f" % (Cycle.LineSetSupply.Q ))
    GUI.txtLineSetSupplyOutputsDP.SetValue("%0.3f" % (Cycle.LineSetSupply.DP ))
    GUI.txtLineSetSupplyOutputsRe.SetValue("%0.3f" % (Cycle.LineSetSupply.Re_fluid ))
    GUI.txtLineSetSupplyOutputsh.SetValue("%0.3f" % (Cycle.LineSetSupply.h_fluid ))
    GUI.txtLineSetSupplyOutputsCharge.SetValue("%0.3f" % (Cycle.LineSetSupply.Charge ))
    GUI.txtLineSetSupplyOutputsTin.SetValue("%0.2f / %0.2f" % (Cycle.LineSetSupply.Tin,Cycle.LineSetSupply.Tin-273.15 ))
    GUI.txtLineSetSupplyOutputsTout.SetValue("%0.2f / %0.2f" % (Cycle.LineSetSupply.Tout,Cycle.LineSetSupply.Tout-273.15 ))
    
    GUI.txtLineSetReturnOutputsQ.SetValue("%0.3f" % (Cycle.LineSetReturn.Q ))
    GUI.txtLineSetReturnOutputsDP.SetValue("%0.3f" % (Cycle.LineSetReturn.DP ))
    GUI.txtLineSetReturnOutputsRe.SetValue("%0.3f" % (Cycle.LineSetReturn.Re_fluid ))
    GUI.txtLineSetReturnOutputsh.SetValue("%0.3f" % (Cycle.LineSetReturn.h_fluid ))
    GUI.txtLineSetReturnOutputsCharge.SetValue("%0.3f" % (Cycle.LineSetReturn.Charge ))
    GUI.txtLineSetReturnOutputsTin.SetValue("%0.2f / %0.2f" % (Cycle.LineSetReturn.Tin,Cycle.LineSetReturn.Tin-273.15 ))
    GUI.txtLineSetReturnOutputsTout.SetValue("%0.2f / %0.2f" % (Cycle.LineSetReturn.Tout,Cycle.LineSetReturn.Tout-273.15 ))
        
    if Cycle.CycleType=='DX':
        GUI.txtCycleSHR.SetValue("%0.4f" % (Cycle.Evaporator.SHR ))
        GUI.txtCycleMainFreezeProtection.SetValue("N/A" )
        
        GUI.txtEvaporatorQ.SetValue("%0.2f" % (Cycle.Evaporator.Q ))
        GUI.txtEvaporatorQ_superheat.SetValue("%0.2f" % (Cycle.Evaporator.Q_superheat ))
        GUI.txtEvaporatorQ_2phase.SetValue("%0.2f" % (Cycle.Evaporator.Q_2phase ))
        GUI.txtEvaporatorQ_subcool.SetValue("N/A")
        GUI.txtEvaporatorCharge.SetValue("%0.5f" % (Cycle.Evaporator.Charge ))
        GUI.txtEvaporatorCharge_superheat.SetValue("%0.5f" % (Cycle.Evaporator.Charge_superheat ))
        GUI.txtEvaporatorCharge_2phase.SetValue("%0.5f" % (Cycle.Evaporator.Charge_2phase ))
        GUI.txtEvaporatorCharge_subcool.SetValue("N/A" )
        GUI.txtEvaporatorh_r_superheat.SetValue("%0.2f" % (Cycle.Evaporator.h_r_superheat ))
        GUI.txtEvaporatorh_r_2phase.SetValue("%0.2f" % (Cycle.Evaporator.h_r_2phase ))
        GUI.txtEvaporatorh_r_subcool.SetValue("N/A")
        GUI.txtEvaporatorw_superheat.SetValue("%0.3f" % (Cycle.Evaporator.w_superheat ))
        GUI.txtEvaporatorw_2phase.SetValue("%0.3f" % (Cycle.Evaporator.w_2phase ))
        GUI.txtEvaporatorw_subcool.SetValue("N/A" )
        GUI.txtEvaporatorh_a.SetValue("%0.3f" % (Cycle.Evaporator.Fins.h_a ))
        GUI.txtEvaporatoreta_a.SetValue("%0.3f" % (Cycle.Evaporator.Fins.eta_a ))
        GUI.txtEvaporatordP_a.SetValue("%0.3f" % (Cycle.Evaporator.Fins.dP_a ))
        GUI.txtEvaporatormdot_a.SetValue("%0.3f" % (Cycle.Evaporator.Fins.mdot_da ))
        GUI.txtEvaporatorA_a.SetValue("%0.3f" % (Cycle.Evaporator.Fins.A_a ))
        GUI.txtEvaporatorSHR.SetValue("%0.3f" % (Cycle.Evaporator.SHR ))
        GUI.txtEvaporatorfdry_2phase.SetValue("%0.3f" % (Cycle.Evaporator.fdry_2phase ))
        GUI.txtEvaporatorfdry_superheat.SetValue("%0.3f" % (Cycle.Evaporator.fdry_superheat ))
        GUI.txtEvaporatorTout_a.SetValue("%0.2f / %0.2f" % (Cycle.Evaporator.Tout_a,Cycle.Evaporator.Tout_a-273.15 ))
        GUI.txtEvaporatorDP_r.SetValue("%0.3f" % (Cycle.Evaporator.DP_r ))
        GUI.txtEvaporatorDP_r_superheat.SetValue("%0.3f" % (Cycle.Evaporator.DP_r_superheat ))
        GUI.txtEvaporatorDP_r_2phase.SetValue("%0.3f" % (Cycle.Evaporator.DP_r_2phase ))
        GUI.txtEvaporatorDP_r_subcool.SetValue("N/A" )
    
    elif Cycle.CycleType=='Secondary':
        
        GUI.txtCycleSHR.SetValue("%0.4f" % (Cycle.CoolingCoil.SHR ))
        #Freeze protection with secondary working fluid
        Tf=Props('F','T',275,'P',101,Cycle.CoolingCoil.Ref_g)
        GUI.txtCycleMainFreezeProtection.SetValue("%0.2f / %0.2f" % (Tf,Tf-273.15 ))
        
        GUI.txtPumpPower.SetValue("%0.2f" % (Cycle.Pump.W ))
        GUI.txtPumpDP.SetValue("%0.2f" % (Cycle.Pump.DP_g ))
        
        GUI.txtCoolingCoilQ.SetValue("%0.2f" % (Cycle.CoolingCoil.Q ))
        GUI.txtCoolingCoilh_g.SetValue("%0.2f" % (Cycle.CoolingCoil.h_g ))
        GUI.txtCoolingCoilRe_g.SetValue("%0.2f" % (Cycle.CoolingCoil.Re_g ))
        GUI.txtCoolingCoilDP_g.SetValue("%0.2f" % (Cycle.CoolingCoil.DP_g ))
        GUI.txtCoolingCoilTin_g.SetValue("%0.2f / %0.2f" % (Cycle.CoolingCoil.Tin_g,Cycle.CoolingCoil.Tin_g-273.15) )
        GUI.txtCoolingCoilTout_g.SetValue("%0.2f / %0.2f" % (Cycle.CoolingCoil.Tout_g,Cycle.CoolingCoil.Tout_g-273.15) )
        GUI.txtCoolingCoilAirdP_a.SetValue("%0.4f" % (Cycle.CoolingCoil.Fins.dP_a ))
        GUI.txtCoolingCoilAirh_a.SetValue("%0.4f" % (Cycle.CoolingCoil.Fins.h_a ))
        GUI.txtCoolingCoilAireta_a.SetValue("%0.4f" % (Cycle.CoolingCoil.Fins.eta_a ))
        GUI.txtCoolingCoilAirmdot_a.SetValue("%0.4f" % (Cycle.CoolingCoil.Fins.mdot_da ))
        GUI.txtCoolingCoilAirA_a.SetValue("%0.4f" % (Cycle.CoolingCoil.Fins.A_a ))
        GUI.txtCoolingCoilAirf_dry.SetValue("%0.4f" % (Cycle.CoolingCoil.f_dry ))
        GUI.txtCoolingCoilAirSHR.SetValue("%0.4f" % (Cycle.CoolingCoil.SHR ))
        GUI.txtCoolingCoilAirTout.SetValue("%0.2f / %0.2f" % (Cycle.CoolingCoil.Tout_a, Cycle.CoolingCoil.Tout_a-273.15 ))
        
        if Cycle.IHXType=='Coaxial':
            IHX=Cycle.CoaxialIHX
            GUI.txtIHXQ.SetValue("%0.2f" % (IHX.Q ))
            GUI.txtIHXQ_superheat.SetValue("%0.2f" % (IHX.Q_superheat ))
            GUI.txtIHXQ_2phase.SetValue("%0.2f" % (IHX.Q_2phase ))
            GUI.txtIHXQ_subcool.SetValue("%0.2f" % (IHX.Q_subcool ))
            GUI.txtIHXCharge.SetValue("%0.5f" % (IHX.Charge_r))
            GUI.txtIHXCharge_superheat.SetValue("%0.5f" % (IHX.Charge_r_superheat ))
            GUI.txtIHXCharge_2phase.SetValue("%0.5f" % (IHX.Charge_r_2phase ))
            GUI.txtIHXCharge_subcool.SetValue("%0.5f" % (IHX.Charge_r_subcool ))
            GUI.txtIHXh_r_superheat.SetValue("%0.2f" % (IHX.h_r_superheat ))
            GUI.txtIHXh_r_2phase.SetValue("%0.2f" % (IHX.h_r_2phase ))
            GUI.txtIHXh_r_subcool.SetValue("%0.2f" % (IHX.h_r_subcool ))
            GUI.txtIHXw_superheat.SetValue("%0.3f" % (IHX.w_superheat ))
            GUI.txtIHXw_2phase.SetValue("%0.3f" % (IHX.w_2phase ))
            GUI.txtIHXw_subcool.SetValue("%0.3f" % (IHX.w_subcool ))
            GUI.txtIHXTout_r.SetValue("%0.2f / %0.2f" % (IHX.Tout_r,IHX.Tout_r-273.15 ))
            GUI.txtIHXTin_r.SetValue("%0.2f / %0.2f" % (IHX.Tin_r,IHX.Tin_r-273.15 ))
            GUI.txtIHXTout_g.SetValue("%0.2f / %0.2f" % (IHX.Tout_g,IHX.Tout_g-273.15 ))
            GUI.txtIHXTin_g.SetValue("%0.2f / %0.2f" % (IHX.Tin_g,IHX.Tin_g-273.15 ))
            GUI.txtIHXDP_r.SetValue("%0.3f" % (IHX.DP_r ))
            GUI.txtIHXDP_r_superheat.SetValue("%0.3f" % (IHX.DP_r_superheat ))
            GUI.txtIHXDP_r_2phase.SetValue("%0.3f" % (IHX.DP_r_2phase ))
            GUI.txtIHXDP_r_subcool.SetValue("%0.3f" % (IHX.DP_r_subcool ))
            
            GUI.txtIHXh_g.SetValue("%0.2f" % (IHX.h_g))
            GUI.txtIHXRe_g.SetValue("%0.2f" % (IHX.Re_g))
            GUI.txtIHXDP_g.SetValue("%0.2f" % (IHX.DP_g))
        elif Cycle.IHXType=='PHE':
            IHX=Cycle.PHEIHX
            GUI.txtIHXQ.SetValue("%0.2f" % (IHX.Q ))
            GUI.txtIHXQ_superheat.SetValue("%0.2f" % (IHX.Q_superheated_c ))
            GUI.txtIHXQ_2phase.SetValue("%0.2f" % (IHX.Q_2phase_c ))
            GUI.txtIHXQ_subcool.SetValue("%0.2f" % (IHX.Q_subcooled_c ))
            GUI.txtIHXCharge.SetValue("%0.5f" % (IHX.Charge_c))
            GUI.txtIHXCharge_superheat.SetValue("%0.5f" % (IHX.Charge_superheated_c ))
            GUI.txtIHXCharge_2phase.SetValue("%0.5f" % (IHX.Charge_2phase_c ))
            GUI.txtIHXCharge_subcool.SetValue("%0.5f" % (IHX.Charge_subcooled_c ))
#            GUI.txtIHXh_r_superheat.SetValue("%0.2f" % (IHX.h_r_superheat ))
#            GUI.txtIHXh_r_2phase.SetValue("%0.2f" % (IHX.h_r_2phase ))
#            GUI.txtIHXh_r_subcool.SetValue("%0.2f" % (IHX.h_r_subcool ))
            GUI.txtIHXw_superheat.SetValue("%0.3f" % (IHX.w_superheated_c ))
            GUI.txtIHXw_2phase.SetValue("%0.3f" % (IHX.w_2phase_c ))
            GUI.txtIHXw_subcool.SetValue("%0.3f" % (IHX.w_subcooled_c ))
            GUI.txtIHXTout_r.SetValue("%0.2f / %0.2f" % (IHX.Tout_c,IHX.Tout_c-273.15 ))
            GUI.txtIHXTin_r.SetValue("%0.2f / %0.2f" % (IHX.Tin_c,IHX.Tin_c-273.15 ))
            GUI.txtIHXTout_g.SetValue("%0.2f / %0.2f" % (IHX.Tout_h,IHX.Tout_h-273.15 ))
            GUI.txtIHXTin_g.SetValue("%0.2f / %0.2f" % (IHX.Tin_h,IHX.Tin_h-273.15 ))
            GUI.txtIHXDP_r.SetValue("%0.3f" % (IHX.DP_c ))
            GUI.txtIHXDP_r_superheat.SetValue("%0.3f" % (IHX.DP_superheated_c ))
            GUI.txtIHXDP_r_2phase.SetValue("%0.3f" % (IHX.DP_2phase_c ))
            GUI.txtIHXDP_r_subcool.SetValue("%0.3f" % (IHX.DP_subcooled_c ))
            GUI.txtIHXDP_g.SetValue("%0.2f" % (IHX.DP_h))
        else:
            raise ValueError('Secondary loop system must have a coaxial or PHE heat exchanger')

def GUI2SecondaryCycleInputs(GUI):
    #Create the basic structure of the DXCycle
    Cycle=SecondaryCycleClass()
    
    Cycle.Mode='AC'
    Cycle.Ref=str(GUI.cmbRefrigerant.GetValue())
    Cycle.SecLoopFluid=str(GUI.cmbSecFluid.GetValue())
    
    
            
    M=[0,0,0,0,0,0,0,0,0,0] #Empty list
    M[0]=float(GUI.txtCompM1.GetValue())
    M[1]=float(GUI.txtCompM2.GetValue())
    M[2]=float(GUI.txtCompM3.GetValue())
    M[3]=float(GUI.txtCompM4.GetValue())
    M[4]=float(GUI.txtCompM5.GetValue())
    M[5]=float(GUI.txtCompM6.GetValue())
    M[6]=float(GUI.txtCompM7.GetValue())
    M[7]=float(GUI.txtCompM8.GetValue())
    M[8]=float(GUI.txtCompM9.GetValue())
    M[9]=float(GUI.txtCompM10.GetValue())
    
    P=[0,0,0,0,0,0,0,0,0,0] #Empty list
    P[0]=float(GUI.txtCompP1.GetValue())
    P[1]=float(GUI.txtCompP2.GetValue())
    P[2]=float(GUI.txtCompP3.GetValue())
    P[3]=float(GUI.txtCompP4.GetValue())
    P[4]=float(GUI.txtCompP5.GetValue())
    P[5]=float(GUI.txtCompP6.GetValue())
    P[6]=float(GUI.txtCompP7.GetValue())
    P[7]=float(GUI.txtCompP8.GetValue())
    P[8]=float(GUI.txtCompP9.GetValue())
    P[9]=float(GUI.txtCompP10.GetValue())
    
    Cycle.Compressor.M=M
    Cycle.Compressor.P=P

    Cycle.Compressor.fp                 =float(GUI.txtCompfp.GetValue())
    Cycle.Compressor.Vdot_ratio         =float(GUI.txtCompVdot_ratio.GetValue())
    
    
    Cycle.CoolingCoil.Fins.Air.Vdot_ha   =float(GUI.txtCoolingCoilAirVdot.GetValue())
    Cycle.CoolingCoil.Fins.Air.Tdb       =float(GUI.txtCoolingCoilAirTdb.GetValue())
    Cycle.CoolingCoil.Fins.Air.p         =float(GUI.txtCoolingCoilAirp.GetValue())
    Cycle.CoolingCoil.Fins.Air.RH        =float(GUI.txtCoolingCoilAirRH.GetValue())
    Cycle.CoolingCoil.Fins.Air.Tmean      =float(GUI.txtCoolingCoilAirTdb.GetValue())
    Cycle.CoolingCoil.Fins.Air.RHmean     =float(GUI.txtCoolingCoilAirRH.GetValue())
    Cycle.CoolingCoil.Fins.Tubes.NTubes_per_bank  =float(GUI.txtCoolingCoilTubesNtubes.GetValue())
    Cycle.CoolingCoil.Fins.Tubes.Nbank   =float(GUI.txtCoolingCoilTubesNbank.GetValue())
    Cycle.CoolingCoil.Fins.Tubes.Ltube   =float(GUI.txtCoolingCoilTubesL.GetValue())
    Cycle.CoolingCoil.Fins.Tubes.OD      =float(GUI.txtCoolingCoilTubesOD.GetValue())
    Cycle.CoolingCoil.Fins.Tubes.ID      =float(GUI.txtCoolingCoilTubesID.GetValue())
    Cycle.CoolingCoil.Fins.Tubes.Pl      =float(GUI.txtCoolingCoilTubesPl.GetValue())
    Cycle.CoolingCoil.Fins.Tubes.Pt      =float(GUI.txtCoolingCoilTubesPt.GetValue())
    Cycle.CoolingCoil.Fins.Tubes.Ncircuits=float(GUI.txtCoolingCoilTubesNcircuit.GetValue())
    Cycle.CoolingCoil.Fins.Fins.FPI      =float(GUI.txtCoolingCoilFinFPI.GetValue())
    Cycle.CoolingCoil.Fins.Fins.Pd       =float(GUI.txtCoolingCoilFinpd.GetValue())
    Cycle.CoolingCoil.Fins.Fins.xf       =float(GUI.txtCoolingCoilFinxf.GetValue())
    Cycle.CoolingCoil.Fins.Fins.t        =float(GUI.txtCoolingCoilFint.GetValue())
    Cycle.CoolingCoil.Fins.Fins.k_fin    =float(GUI.txtCoolingCoilFink.GetValue())
    Cycle.CoolingCoil.Fins.Air.FanPower  =float(GUI.txtCoolingCoilPower.GetValue())
    Cycle.CoolingCoil.Verbosity           =0
    Cycle.CoolingCoil.pin_g               =float(GUI.txtPumpp.GetValue())
    Cycle.CoolingCoil.Ref_g               =str(GUI.cmbSecFluid.GetValue())
    
    Cycle.Pump.mdot_g                    =float(GUI.txtPumpmdot.GetValue())
    Cycle.Pump.eta                       =float(GUI.txtPumpEfficiency.GetValue())
    Cycle.Pump.pin_g                     =float(GUI.txtPumpp.GetValue())
    Cycle.Pump.Ref_g                     =str(GUI.cmbSecFluid.GetValue())
    
    Cycle.LineSetSupply.L                 =float(GUI.txtLineSetL.GetValue())
    Cycle.LineSetSupply.OD                =float(GUI.txtLineSetOD_supply.GetValue())
    Cycle.LineSetSupply.ID                =float(GUI.txtLineSetID_supply.GetValue())
    Cycle.LineSetSupply.t_insul           =float(GUI.txtLineSetInsult.GetValue())
    Cycle.LineSetSupply.k_tube            =float(GUI.txtLineSetTubek.GetValue())
    Cycle.LineSetSupply.k_insul           =float(GUI.txtLineSetInsulk.GetValue())
    Cycle.LineSetSupply.h_air             =float(GUI.txtLineSeth_air.GetValue())
    Cycle.LineSetSupply.T_air             =float(GUI.txtLineSetT_air.GetValue())
    Cycle.LineSetSupply.Ref               =str(GUI.cmbSecFluid.GetValue())
    Cycle.LineSetSupply.pin               =float(GUI.txtPumpp.GetValue())
    
    Cycle.LineSetReturn.L                 =float(GUI.txtLineSetL.GetValue())
    Cycle.LineSetReturn.OD                =float(GUI.txtLineSetOD_supply.GetValue())
    Cycle.LineSetReturn.ID                =float(GUI.txtLineSetID_supply.GetValue())
    Cycle.LineSetReturn.t_insul           =float(GUI.txtLineSetInsult.GetValue())
    Cycle.LineSetReturn.k_tube            =float(GUI.txtLineSetTubek.GetValue())
    Cycle.LineSetReturn.k_insul           =float(GUI.txtLineSetInsulk.GetValue())
    Cycle.LineSetReturn.h_air             =float(GUI.txtLineSeth_air.GetValue())
    Cycle.LineSetReturn.T_air             =float(GUI.txtLineSetT_air.GetValue())
    Cycle.LineSetReturn.Ref               =str(GUI.cmbSecFluid.GetValue())
    Cycle.LineSetReturn.pin               =float(GUI.txtPumpp.GetValue())
    
    if GUI.optIHXUseCoaxial.GetValue()==1:
        Cycle.CoaxialIHX.ID_i                =float(GUI.txtIHXTubeID.GetValue())
        Cycle.CoaxialIHX.OD_i                =float(GUI.txtIHXAnnID.GetValue())
        Cycle.CoaxialIHX.ID_o                =float(GUI.txtIHXAnnOD.GetValue())
        Cycle.CoaxialIHX.L                   =float(GUI.txtIHXLength.GetValue())
        Cycle.CoaxialIHX.pin_g               =float(GUI.txtPumpp.GetValue())
        Cycle.CoaxialIHX.Ref_g               =str(GUI.cmbSecFluid.GetValue())
        Cycle.CoaxialIHX.Ref_r               =str(GUI.cmbRefrigerant.GetValue())
        Cycle.CoaxialIHX.Verbosity           =0
        Cycle.IHXType                        ='Coaxial'
    else: 
        Cycle.PHEIHX.Bp                      =float(GUI.txtIHXPlateBp.GetValue())
        Cycle.PHEIHX.Lp                      =float(GUI.txtIHXPlateLp.GetValue())
        Cycle.PHEIHX.Nplates                 =int(GUI.txtIHXPlateN.GetValue())
        Cycle.PHEIHX.PlateAmplitude          =float(GUI.txtIHXPlateAmplitude.GetValue())
        Cycle.PHEIHX.PlateThickness          =float(GUI.txtIHXPlateThickness.GetValue())
        Cycle.PHEIHX.PlateConductivity       =float(GUI.txtIHXPlateConductivity.GetValue())
        Cycle.PHEIHX.PlateWavelength         =float(GUI.txtIHXPlateWavelength.GetValue())
        Cycle.PHEIHX.InclinationAngle        =float(GUI.txtIHXInclinationAngle.GetValue())/180.0*pi
        Cycle.IHXType                        ='PHE'
        Cycle.PHEIHX.Verbosity               =0
        Cycle.PHEIHX.DT_sh               =float(GUI.txtCycleDTsh.GetValue())

        if GUI.radIHXChannelSelect.GetSelection()==0:
                Cycle.PHEIHX.MoreChannels='Hot'
        elif GUI.radIHXChannelSelect.GetSelection()==1:
            Cycle.PHEIHX.MoreChannels='Cold'
        
        if GUI.radCycleMode.GetStringSelection()=='Cooling Mode':
            Cycle.PHEIHX.pin_h                   =Cycle.Pump.pin_g
            Cycle.PHEIHX.Ref_c                   =str(GUI.cmbRefrigerant.GetValue())
            Cycle.PHEIHX.Ref_h                   =str(GUI.cmbSecFluid.GetValue())
        elif GUI.radCycleMode.GetStringSelection()=='Heating Mode':
            Cycle.PHEIHX.pin_c                   =Cycle.Pump.pin_g
            Cycle.PHEIHX.Ref_c                   =str(GUI.cmbSecFluid.GetValue())
            Cycle.PHEIHX.Ref_h                   =str(GUI.cmbRefrigerant.GetValue())
            
    if GUI.radCycleMode.GetStringSelection()=='Cooling Mode':
        Cycle.Condenser.Fins.Air.Vdot_ha    =float(GUI.txtCondenserAirVdot.GetValue())
        Cycle.Condenser.Fins.Air.Tdb        =float(GUI.txtCondenserAirTdb.GetValue())
        Cycle.Condenser.Fins.Air.p          =float(GUI.txtCondenserAirp.GetValue())
        Cycle.Condenser.Fins.Air.RH         =float(GUI.txtCondenserAirRH.GetValue())
        Cycle.Condenser.Fins.Air.Tmean      =float(GUI.txtCondenserAirTdb.GetValue())
        Cycle.Condenser.Fins.Air.RHmean     =float(GUI.txtCondenserAirRH.GetValue())
        
        Cycle.Condenser.Fins.Tubes.NTubes_per_bank   =int(GUI.txtCondenserTubesNtubes.GetValue())
        Cycle.Condenser.Fins.Tubes.Nbank    =int(GUI.txtCondenserTubesNbank.GetValue())
        Cycle.Condenser.Fins.Tubes.Ltube    =float(GUI.txtCondenserTubesL.GetValue())
        Cycle.Condenser.Fins.Tubes.OD       =float(GUI.txtCondenserTubesOD.GetValue())
        Cycle.Condenser.Fins.Tubes.ID       =float(GUI.txtCondenserTubesID.GetValue())
        Cycle.Condenser.Fins.Tubes.Pl       =float(GUI.txtCondenserTubesPl.GetValue())
        Cycle.Condenser.Fins.Tubes.Pt       =float(GUI.txtCondenserTubesPt.GetValue())
        Cycle.Condenser.Fins.Tubes.Ncircuits =int(GUI.txtCondenserTubesNcircuit.GetValue())
        Cycle.Condenser.Fins.Fins.FPI       =float(GUI.txtCondenserFinFPI.GetValue())
        Cycle.Condenser.Fins.Fins.Pd        =float(GUI.txtCondenserFinpd.GetValue())
        Cycle.Condenser.Fins.Fins.xf        =float(GUI.txtCondenserFinxf.GetValue())
        Cycle.Condenser.Fins.Fins.t         =float(GUI.txtCondenserFint.GetValue())
        Cycle.Condenser.Fins.Fins.k_fin     =float(GUI.txtCondenserFink.GetValue())
        Cycle.Condenser.Fins.Air.FanPower   =float(GUI.txtCondenserPower.GetValue())
        Cycle.Condenser.Verbosity           =0
    elif GUI.radCycleMode.GetStringSelection()=='Heating Mode':
        #In heating mode, the "condenser" from cooling mode is the evaporator
        Cycle.Evaporator.Fins.Air.Vdot_ha    =float(GUI.txtCondenserAirVdot.GetValue())
        Cycle.Evaporator.Fins.Air.Tdb        =float(GUI.txtCondenserAirTdb.GetValue())
        Cycle.Evaporator.Fins.Air.p          =float(GUI.txtCondenserAirp.GetValue())
        Cycle.Evaporator.Fins.Air.RH         =float(GUI.txtCondenserAirRH.GetValue())
        Cycle.Evaporator.Fins.Air.Tmean      =float(GUI.txtCondenserAirTdb.GetValue())
        Cycle.Evaporator.Fins.Air.RHmean     =float(GUI.txtCondenserAirRH.GetValue())
        
        Cycle.Evaporator.Fins.Tubes.NTubes_per_bank   =int(GUI.txtCondenserTubesNtubes.GetValue())
        Cycle.Evaporator.Fins.Tubes.Nbank    =int(GUI.txtCondenserTubesNbank.GetValue())
        Cycle.Evaporator.Fins.Tubes.Ltube    =float(GUI.txtCondenserTubesL.GetValue())
        Cycle.Evaporator.Fins.Tubes.OD       =float(GUI.txtCondenserTubesOD.GetValue())
        Cycle.Evaporator.Fins.Tubes.ID       =float(GUI.txtCondenserTubesID.GetValue())
        Cycle.Evaporator.Fins.Tubes.Pl       =float(GUI.txtCondenserTubesPl.GetValue())
        Cycle.Evaporator.Fins.Tubes.Pt       =float(GUI.txtCondenserTubesPt.GetValue())
        Cycle.Evaporator.Fins.Tubes.Ncircuits =int(GUI.txtCondenserTubesNcircuit.GetValue())
        Cycle.Evaporator.Fins.Fins.FPI       =float(GUI.txtCondenserFinFPI.GetValue())
        Cycle.Evaporator.Fins.Fins.Pd        =float(GUI.txtCondenserFinpd.GetValue())
        Cycle.Evaporator.Fins.Fins.xf        =float(GUI.txtCondenserFinxf.GetValue())
        Cycle.Evaporator.Fins.Fins.t         =float(GUI.txtCondenserFint.GetValue())
        Cycle.Evaporator.Fins.Fins.k_fin     =float(GUI.txtCondenserFink.GetValue())
        Cycle.Evaporator.Fins.Air.FanPower   =float(GUI.txtCondenserPower.GetValue())
        Cycle.Evaporator.Verbosity           =0
        
        Cycle.Evaporator.DT_sh              =float(GUI.txtCycleDTsh.GetValue())
    else:
        raise ValueError()
    
    
    Cycle.Charge_target=float(GUI.txtCycleCharge.GetValue())
    Cycle.DT_sc_target=float(GUI.txtCycleSubcooling.GetValue())
    
    Cycle.Verbosity = 10
    Cycle.ParaPath = GUI.txtParaPath.GetValue()
    
    return Cycle
    
if __name__=='__main__':
    import sys,os
    sys.path.append('..')
    import SecondaryCycle
    InputStructure=SecondaryCycle.CycleInputVals()
    InputFile=os.path.join('configs','Default.cfg')
    ConfigFile2InputFields(InputFile, InputStructure)
