from __future__ import division, print_function, absolute_import
from CoolProp.CoolProp import HAPropsSI 
from math import sqrt,pi,tanh,exp,cos,log,sin
from ACHP.ACHPTools import ValidateFields
from scipy.optimize import fsolve

def IsFinsClass(MicroFins):
    '''
    Returns the Fins class if Fins is an instance of the FinInputs class, False otherwise
    
    Convenience function for the Validator 
    '''
    if isinstance(MicroFins,MicroFinInputs):
        return MicroFins
    else:
        return False
    
class FinsVals:
    pass    
class TubesVals:
    pass
class AirVals:
    pass
class LouversVals:
    pass
        
class MicroFinInputs():
    """ 
    Empty Class for fin data
    """
    def __init__(self):
        #Don't do anything
        self.Tubes=TubesVals()
        self.Air=AirVals()
        self.Fins=FinsVals()
        self.Louvers=LouversVals()
        
    def __repr__(self):
        string="Tubes::\n"
        for field in self.Tubes.__dict__.keys():
            string+=field+": "+repr(self.Tubes.__dict__[field])+"\n"
        string+="Fins::\n"
        for field in self.Fins.__dict__.keys():
            string+=field+": "+repr(self.Fins.__dict__[field])+"\n"
        string+="Air::\n"
        for field in self.Air.__dict__.keys():
            string+=field+": "+repr(self.Air.__dict__[field])+"\n"
        string+="Louvers::\n"
        for field in self.Louvers.__dict__.keys():
            string+=field+": "+repr(self.Louvers.__dict__[field])+"\n"
        return string
    
    def Validate(self):
        """
        Check that all required fields are included, and no extra fields not listed here are added
        
        This is quite strict, but for the best to avoid typos
        """       
        reqFields=[
            ('RH',float,0,1),
            ('Tdb',float,-80+273.15,200+273.15),
            ('FanPower',float,0,4000),
            ('p',float,10,10000000),
            ('Vdot_ha',float,0.001,10)
        ]
        optFields=['RHmean','Tmean']
        d=dict(self.Air.__dict__) #The current values
        ValidateFields(d,reqFields,optFields)
        
        reqFields=[
            ('FPI',float,0.1,100),
            ('Lf',float,0.0001,1),
            ('t',float,0.00001,0.01),
            ('k_fin',float,0.01,10000)
        ]
        optFields=None
        d=dict(self.Fins.__dict__) #The current values
        ValidateFields(d,reqFields,optFields)
        
        reqFields=[
            ('NTubes',float,0.1,100),
            ('Nbank',float,1,50),
            ('Npass',float,1,50),
            ('Nports',float,1,50),
            ('Ltube',float,0.001,10),
            ('Td',float,0.0001,1),
            ('Ht',float,0.0001,1),
            ('b',float,0.0001,1),
            ('tw',float,0.00001,0.01),
            ('twp',float,0.00001,0.01),
            ('beta',float,0.00001,100),
            ('kw',float,0.01,10000)
        ]
        optFields=None
        d=dict(self.Tubes.__dict__) #The current values
        ValidateFields(d,reqFields,optFields)

        reqFields=[
            ('Lalpha',float,1,89),
            ('lp',float,0.0001,1),
        ]
        optFields=['Llouv']
        d=dict(self.Louvers.__dict__) #The current values
        ValidateFields(d,reqFields,optFields)

def MultiLouveredMicroFins(Inputs):
    """
    # The analysis is based on Lee book 2010 "Thermal design:Themoelectric, .." 
    # In addition, correlation from
    # Man-Hoe Kim and Clark W. Bullard, 2002, "Air-side thermal hydraulic
    # performance of multi-louvered fin aluminum heat exchanger", International 
    # journal of refrigeration
    
    
    |<------------------- Lf ------------------->| 
     
     ____________________________________________
    |                                            |     =
    |    ||     ||    ||        ||    ||   ||    |     |
    |    ||     ||    ||        ||    ||   ||    |   Llouv
    |    ||     ||    ||        ||    ||   ||    |     |
    |    ||     ||    ||        ||    ||   ||    |     =
    |____________________________________________|
                                      ||
                                   -->||<--
                                      lp
       
          /     /     /---------\     \    \
    _____/     /     /           \     \    \____
        
      
     Lf: Fin length (flow depth)
     Llouv: Louver cut length
     Lalpha: Louver angle
     
     
    
             |-  Td  -|
              ________     
             |________|
       =
       |      
       |
       b
       |
       |
       =      ________   
             |________|  Ht
     
     
     Td: tube outside width (depth)
     Ht: tube outside height (major diameter)
     b: tube spacing                     
    
    """
    Lalpha =      Inputs.Louvers.Lalpha         #Louver angle, in degree
    lp =          Inputs.Louvers.lp             #Louver pitch
    
    delta =       Inputs.Fins.t                 #Fin thickness
    Lf =          Inputs.Fins.Lf                #Fin length (flow depth)
    k_fin =       Inputs.Fins.k_fin             #Thermal conductivity
    FPI =         Inputs.Fins.FPI               #Fin per inch (fin density)
    
    Ntubes =      Inputs.Tubes.NTubes           #Number of tubes
    Nbank =       Inputs.Tubes.Nbank            #Number of banks
    L3 =          Inputs.Tubes.Ltube            #length of a single tube    
    L2 =          Inputs.Tubes.Td               #Tube outside width (depth)
    Ht =          Inputs.Tubes.Ht               #Tube outside height (major diameter)
    b =           Inputs.Tubes.b                #Tube spacing       
    
    Vdot_ha =     Inputs.Air.Vdot_ha            #Air volume flow rate, m^3/s
    p =           Inputs.Air.p                  #Air pressure, Pa
    
    if hasattr(Inputs.Air,'Tmean'):
        Temp = Inputs.Air.Tmean
    else:
        Temp = Inputs.Air.Tdb
    
    if hasattr(Inputs.Air,'RHmean'):
        RHin = Inputs.Air.RHmean
    else:
        RHin = Inputs.Air.RH

        
    # Check that cs_cp is defined, if so, set it to the value passed in
    if (hasattr(Inputs,'cs_cp') and Inputs.cs_cp>0) or (hasattr(Inputs,'WetDry') and Inputs.WetDry=='Wet'):
        isWet=True
        cs_cp=Inputs.Air.cs_cp
    else:
        isWet=False
        cs_cp=1.0

    #Fins per meter (fin density) [1/m]
    FPM = FPI / 0.0254
    #Fin pitch (distance between centerlines of fins)
    pf = 1 / FPM
    #Fin height
    sf = sqrt(b**2 + pf**2) 
    #Louver cut length
    if hasattr(Inputs.Louvers,'Llouv'):
        Llouv = Inputs.Louvers.Llouv
    else:
        Llouv = 0.85 *sf
    
    #Fin pitch
    pt = Ht + b
    
    #Louver height
    lh = lp * sin(Lalpha*pi/180)
    
    #Air passages
    Npg = Ntubes + 1 #sometime there are no tubes on the edges of HX, so >>> Npg = Ntubes_bank + 1, otherwise Npg = Ntubes_bank - 1
    #Height of heat exchanger (core width)
    L1 = Npg * b + Ntubes * Ht
    #Total number of fins (per bank)
    nf = L3/pf * Npg
    #Primary area =  Tube outside surface area - Fin base area (per bank)
    Ap = (2*(L2 - Ht) + pi * Ht)*L3 *Ntubes - 2*delta*L2*nf
    #Total number of louvers (per bank)
    nlouv = (Lf/lp - 1)*nf
    #Total fin area = Fin area + Louver edge area (per bank)
    Af = 2 * (sf*Lf + sf*delta)*nf + 2*Llouv*delta*nlouv
    #Total surface area on air-side
    At = (Af + Ap) * Nbank
    
    #Minimum free-flow area on air-side = area spacing between tubes - fin and louver edge area
    Ac = b*L3*Npg - (delta*(sf-Llouv) +Llouv*lh)*nf
    #Frontal area on air-side
    Afr = L1*L3
    #Hydraulic diameter on air-side
    Dh = 4*Ac*L2*Nbank / At
    #Porosity on air-side
    sigma = Ac/Afr
    #Volume of HX on air-side
    Vhx = L2*L3*b*Npg*Nbank
    #Surface area density on air-side
    beta = At/Vhx
    
    
    #Evaluate the mass flow rate based on inlet conditions
    # To convert a parameter from per kg_{humid air} to per kg_{dry air}, divide by (1+W)
    W=HAPropsSI('W','T',Inputs.Air.Tdb,'P',p,'R',Inputs.Air.RH)
    v_da=HAPropsSI('V','T',Inputs.Air.Tdb,'P',p,'W',W)
    h_da=HAPropsSI('H','T',Inputs.Air.Tdb,'P',p,'W',W)
    rho_ha = 1 / v_da*(1+W) #[m^3/kg_ha]
    rho_da = 1 / v_da #[m^3/kg_da]
    mdot_ha = Vdot_ha * rho_ha #[kg_ha/s]
    mdot_da = Vdot_ha * rho_da #[kg_da/s]
    
    #mass flux on air-side
    G = mdot_ha / Ac #[kg/m^2-s]
    #maximum velocity on air-side
    umax = G / rho_ha * Afr/Ac #[m/s]
    
    #Use a forward difference to calculate cp from cp=dh/dT
    dT=0.0001 #[K]
    cp_da=(HAPropsSI('H','T',Inputs.Air.Tdb+dT,'P', p, 'W',W)-h_da)/dT #*1000 #[J/kg_da/K]
    cp_ha=cp_da/(1+W) #[J/kg_ha/K]
    
    #Transport properties of humid air from CoolProp
    mu_ha=HAPropsSI('M','T',Inputs.Air.Tdb,'P',p,'W',W)
    k_ha=HAPropsSI('K','T',Inputs.Air.Tdb,'P',p,'W',W)
    
    #Dimensionless Groups
    Pr = cp_ha * mu_ha / k_ha #Prandlt's number
    Re_Dh = G * Dh / mu_ha #Reynold's number on air-side
    Re_lp = rho_ha * umax * lp / mu_ha
    
    #Heat transfer
    
    #Colburn j-Factor based on Chang & Wang, "A Generalized Heat Transfer 
    #Correlation for Louver Fin Geometry." Int. J. Heat Mass Transfer, 40 (1997): 553-554
    #j = pow(Re_lp,-0.49) * pow((Lalpha/90.0),0.27) * pow(pf/lp,-0.14) * pow(b/lp,-0.29) * pow(Lf/lp,-0.23) * pow(Llouv/lp,0.68) * pow(pt/lp,-0.28) * pow(delta/lp,-0.05)
    #h_a = j * G * cp_ha / pow(Pr,2.0/3.0)
    
    #Colburn j-Factor based on Kim & Bullard, "Air-side thermal hydraulic performance 
    #of multi-louvered fin aluminum heat exchanger" Int. J. Refrigeration, 25 (2002): 390-400
    j = pow(Re_lp,-0.487) * pow((Lalpha/90.0),0.257) * pow(pf/lp,-0.13) * pow(b/lp,-0.29) * pow(Lf/lp,-0.235) * pow(Llouv/lp,0.68) * pow(pt/lp,-0.279) * pow(delta/lp,-0.05)
    h_a = j * rho_ha * umax * cp_ha / pow(Pr,2.0/3.0)
    
    #Air-side pressure drop correlations based on Chang and el., "A Generalized Friction
    #Correlation for Louver Fin Geometry." (2000) Int. J. Heat Mass Transfer, 43, 2237-2243
    if (Re_lp<150):
        fa = 14.39 * Re_lp**(-0.805*pf/sf) * pow(log(1.0 + pf/lp),3.04)
        fb = pow(log((delta/pf)**0.48 + 0.9),-1.453) * pow(Dh/lp,-3.01) * pow(log(0.5*Re_lp),-3.01)
        fc = pow(pf/Llouv,-0.308) * pow(Lf/Llouv,-0.308) *exp(-0.1167*pt/Ht) * pow(Lalpha,0.35)
    else:
        fa = 4.97 * pow(Re_lp,(0.6049 - 1.064/Lalpha**0.2)) * pow(log((delta/pf)**0.5 + 0.9),-0.527)
        fb = pow((Dh/lp)*log(0.3*Re_lp),-2.966) * pow(pf/Llouv,-0.7931*pt/b)
        fc = pow(pt/Ht, -0.0446) * pow(log(1.2 + (lp/pf)**1.4),-3.553) * pow(Lalpha,-0.477)    
    #Fanning friction factor    
    f = fa*fb*fc
    #Air-side pressure drop
    #neglecting contraction effect, momentum effect, expansion effect (Kc, Ke ,..etc)
    #this assumption is valid for compact HX based on book of Shah and 
    #Sekulic 2003,"Fundamentals of Heat Exchanger Design"
    DeltaP_air= f * At/Ac * G**2 / (2.0*rho_ha)
    
    
    #Air-side pressure drop including momentum, expansion and contraction effects
    #so the coefficients are send back for give structure
    #Assume Reynold's number is turbulent due to louver fin
    Re_d = 10**7
    #Contraction and expansion coefficient
    C_tube = 4.374e-4 * exp(6.737*sqrt(sigma)) + 0.621
    #Calculate Friction factor, velocity distribution coefficient (triangular tube), 
    if (Re_d >= 2300):
        f_d = 0.049*Re_d**(-0.2) #Friction factor
        Kd_tube = 1.09068*(4*f_d) + 0.05884*sqrt(4*f_d) + 1
        Kd_tri = 1 + 1.29*(Kd_tube-1)    
    else:
        f_d = 64/Re_d
        Kd_tube = 1.33
        Kd_tri = 1.43
        
    #Expansion coefficient
    Ke_tri = 1 - 2*Kd_tri*sigma + sigma**2
    #Contraction coefficient   
    Kc_tri = (1 - 2*C_tube + C_tube**2 *(2*Kd_tri - 1))/C_tube**2
    
    #calcs needed for specific fin types (Based on Kim & Bullard 2002 paper)
    mf = sqrt(2 * h_a * cs_cp / (k_fin * delta) *(1 + delta/Lf) ) #cs_cp is the correction for heat/mass transfer
    #characteristic length (Based on Kim & Bullard 2002 paper)
    Ls = sf/2 - delta
    #Finned surface efficiency
    eta_f = tanh(mf * Ls) / (mf * Ls) #Can be included for wet surface >>> *cos(0.1 * m * Ls)
    #overall surface efficiency
    eta_o = 1 - Af / At * (1 - eta_f)

    
    #write necessary values back into the given structure
    Inputs.Llouv=Llouv
    Inputs.A_a=At
    Inputs.A_a_c=Ac
    Inputs.cp_da=cp_da
    Inputs.cp_ha=cp_ha
    if isWet==True:
        Inputs.eta_a_wet=eta_o
    else:
        Inputs.eta_a=eta_o
    Inputs.h_a=h_a
    Inputs.mdot_ha=mdot_ha
    Inputs.mdot_da=mdot_da
    Inputs.f_a=f
    Inputs.dP_a=DeltaP_air
    Inputs.Re=Re_Dh
    """ADD NEW"""
    Inputs.G_air = G
    Inputs.rho_i_air = rho_ha
    Inputs.sigma = sigma
    Inputs.Kc_tri = Kc_tri
    Inputs.Ke_tri = Ke_tri
    
        
if __name__=='__main__':
    
    LouversFinsTubes=MicroFinInputs()
    
    LouversFinsTubes.Tubes.NTubes=61.354           #Number of tubes
    LouversFinsTubes.Tubes.Nbank=1                 #Number of banks (set to 1 for now!)
    LouversFinsTubes.Tubes.Npass=3                 #Number of passes (per bank)
    LouversFinsTubes.Tubes.Nports=1                #Number of rectangular ports
    LouversFinsTubes.Tubes.Ltube=0.30213           #length of a single tube
    LouversFinsTubes.Tubes.Td=0.0333               #Tube outside width (depth)
    LouversFinsTubes.Tubes.Ht= 0.002               #Tube outside height (major diameter)
    LouversFinsTubes.Tubes.b=0.00635               #Tube spacing     
    LouversFinsTubes.Tubes.tw=0.0003               #Tube wall thickness     
    LouversFinsTubes.Tubes.twp=0.0004              #Port wall thickness     
    LouversFinsTubes.Tubes.beta=1                  #Port (channel) aspect ratio (=width/height)
    LouversFinsTubes.Tubes.kw=117                  #wall thermal conductivity
    
    LouversFinsTubes.Fins.FPI=11.0998              #Fin per inch
    LouversFinsTubes.Fins.Lf=0.0333                #Fin length
    LouversFinsTubes.Fins.t=0.000152               #Fin thickness
    LouversFinsTubes.Fins.k_fin=117                #Fin thermal conductivity
    
    LouversFinsTubes.Air.Vdot_ha=1.05              #Air volume flow rate in m^3/s
    LouversFinsTubes.Air.Tdb=298                   #Air inlet temperature, K
    LouversFinsTubes.Air.p=100000                  #Air pressure in Pa
    LouversFinsTubes.Air.RH=0.5                    #Air inlet relative humidity
    LouversFinsTubes.Air.FanPower=327.36           #Fan power, Watts
    
    LouversFinsTubes.Louvers.Lalpha=20             #Louver angle, in degree
    LouversFinsTubes.Louvers.lp=0.001              #Louver pitch
    LouversFinsTubes.Louvers.Llouv=0.005737        #Louver cut length
    
    LouversFinsTubes.Validate()
    
    print (LouversFinsTubes)  #just print our inputs
    MultiLouveredMicroFins(LouversFinsTubes)  #calculate
    print ("Multi-Louvered Micro fins:","eta_a is:"+str(LouversFinsTubes.eta_a)+", dP_a is:"+str(LouversFinsTubes.dP_a)+" Pa")