from __future__ import division, print_function, absolute_import
from math import log,exp
from scipy import optimize
from numpy import float64,isnan,isinf,fabs
from math import floor


def s_oil(Liq,T):
    
    T0 = 273.15
    P0 = 101.325

    #Liquid Properties
    if Liq == "PAO":
        s_L = 1.940 * log(T/T0)
    elif Liq == "PAG":
        # PAG 0-OB-1020 from Tribology Data Handbook
        # T in K, cp in kJ/kg-K
        try:
            s_L=2.74374E-03*(T-T0)+1.08646*log(T/T0)
        except:
            a=4
    elif Liq == "POE150":
        s_L = 2.30 * log(T/T0)
    elif Liq == "POE32":
        s_L = 2.30 * log(T/T0)
    elif Liq == "Duratherm_LT":
        "specific entropy [kJ/kg-K] of Duratherm LT given T in K"  
        "T0=298.15 K"
        s_L = (3.4014*(T-298)+1094.3*log(T/298.0))/1000  #
    elif Liq == "ACD100FY":
        cl_A = 1.304
        cl_B = 1.035e-3
        cl_C = 2.801e-6
        S_RSV = 1 #[kJ/kg-K]
        s_L = cl_A*log(T/T0) + cl_B*(T-T0) + cl_C/2*(T**2-T0**2) #+ S_RSV
    else:
        print ("Invalid fluid")

    return s_L*1000  #J/kgK


def u_oil(Liq,T):
    
    if Liq == "Duratherm_LT":
        #internal energy [kJ/kg] of Duratherm LT given T in K
        u = (3.4014/2*pow(T,2)+1094.3*T)/1000  #LT
    elif Liq == "Zerol60":
        u = (5.186/2*pow(T,2)+337.116*T)/1000  #Zerol60
    elif Liq == "POE":
        u = (2.0935/2*T**2 + 1186.7*T)/1000 #POE 150 SUS
    else:
        print ("Invalid Fluid")
    
    return u*1000  #J/kg


def h_oil(Liq,T,P):
    """
    Enthalpy of the mixture as a function of temperature [K]
    and pressure [kPa].  Output in J/kg
    """
    T0 = 273.15
    P0 = 101.325
    h = 0
    h_L = 0

    if Liq =='PAO':
        h_L = 1.940*(T-T0)+(P-P0)/849
    elif Liq=='PAG':
        # PAG 0-OB-1020 from Tribology Data Handbook
        rho_L=-0.726923*float64(T)+1200.22;
        h_L=2.74374E-03*(float64(T)**2-T0**2)/2.0+1.08646*(float64(T)-T0)+(float64(P)-P0)/rho_L;
    elif Liq == 'POE':
        # From Totten, p 261, cp=0.55 cal/g-C --> 2.30 kJ/kg-K
        h_L = 2.30*(T-T0)+(P-P0)/930
    elif Liq == 'Duratherm_LT':
        #the specific enthalpy of Duratherm LT [kJ/kg-k]"
        h_L = u_oil(Liq ,T)/1000 + (P-P0)/rho_oil(Liq ,T)
    elif Liq == "ACD100FY":
        cl_A = 1.304
        cl_B = 1.035e-3
        cl_C = 2.801e-6
        H_RSV = 200 #[kJ/kg]
        h_L = cl_A*(T-T0) + cl_B/2*(T**2-T0**2) + cl_C/3*(T**3-T0**3) + H_RSV

    else:
        print ("Invalid fluid")

    return h_L*1000


def rho_oil(Liq,T):
    
    """
    Input:
    T : temperature [K]
    """

    #Liquid Properties
    if Liq == 'PAO':
        rho_L=849
    elif Liq == 'PAG':
        # PAG 0-OB-1020 from Tribology Data Handbook
        rho_L=-0.726923*T+1200.22
    elif Liq == 'POE150':
        rho_L= -0.7*T+1186 ##POE 150 SUS
    elif Liq == 'POE32':
        # Bell 2011 PhD
        rho_L = (-0.00074351165*(T -273.15) + 0.9924395)*1000
    elif Liq == 'Duratherm_LT':
        #density [kg/m^3] of Duratherm LT given T in K"
        rho_L = -0.6793*T + 1012.4 
    elif Liq =="Water":
        # Water props from Yaws
        rhol_A=0.3471     
        rhol_B=0.274      
        rhol_n=0.28571    
        rhol_Tc=647.13
        rho_L=rhol_A/pow(rhol_B,pow(1-T/rhol_Tc,rhol_n))*1000
    elif Liq == 'ACD100FY':
        rho_L = 1.18631056e+03 -7.31369048e-01*T  #T [K]
    else:
        print ("Invalid fluid")

    return rho_L


def cp_oil(Liq,T):

    if Liq ==  'PAO':
        cp_L = 1.940
    elif Liq == 'POE150':
        cp_L = 2.30
    elif Liq == 'POE32':
        cp_L = 2.30
    elif Liq == 'PAG':
        # PAG 0-OB-1020 from Tribology Data Handbook
        # T in K, cp in kJ/kg-K
        cp_L=2.74374E-03*T+1.08646;
    elif Liq == 'Duratherm_LT':
        #specific heat [kJ/kg-K] of Duratherm LT given T in K
        c_L = (3.4014*T + 1094.3)/1000
    elif Liq == "ACD100FY":
        # 273 < T [K] < 387
        cp_L = 1.304 + 1.035e-3*T+2.801e-6*T**2   #[kJ/kg-K]
    return cp_L*1000
    
    
    
def mu_oil(Liq,T):

    "returns the viscosity given temp in K, [Pa-s]"
    mu_l=0
    
    if Liq == 'Duratherm_LT': 
    
        mu_l = 8e12*pow(T,-6.001)  #LT

    elif Liq == 'Zerol60': 
        mu_l = 1.0*(-0.0001235*T + 0.04808) #Zerol60 Vincent
    
    elif Liq == 'POE150':
        #POE equation only valid from 60 C to 120 C
        #POE 150 SUS
        mu_l = 0.000000000517*T**4 - 0.000000795840*T**3 + 0.000460766590*T**2 - 0.118976538068*T + 11.571730524692   
    elif Liq == 'POE32':
        # 32-3MAF POE 
        rho_l = rho_oil(Liq,T)
        mu_l = 0.0002389593*(log(T)**2 - 0.1927238779*log(T) + 40.3718884485)*rho_l*1e-6
            
    elif Liq == "ACD100FY":
        # 313 < T [K] < 387
        #mu_l = 1.603e+46*pow(T,-19.69) #T [K]   
        mu_l = 2.0022e+31*pow(T,-1.2961e+01)
    else:
        print ("Invalid fluid")
        
    return mu_l


def k_oil(Liq,T):

    "Thermal conductivity [W/m-K] for T in Kelvin"

    if Liq == 'Duratherm_LT': 
        k_l = -9e-5*T + 0.1223

    elif Liq == 'Zerol60': 
        k_l =  0.1153 #Zerol60 Vincent
        #k_l = 0.170 #[W/m-K] Ian appendix C4
    
    elif Liq == 'POE':
        #POE 150 SUS
        k_l = 0.138                   
    
    elif Liq == "ACD100FY":
        k_l = 0.13    #TODO: check value of conductivity
    
    else:
        print ("Invalid fluid")
        
    return k_l

    
def Solubility_Ref_in_Liq(Ref,Liq,T,p):
    x_Ref=0.0
    Tmax=273.15
    Tmin=273.15
    error=False

    if Ref=='R744' and Liq=='Water':
        Tmin=273.15+40
        Tmax=273.15+100
        x_CO2=(8.47565584E-01-6.36371603E-03*T
        +1.84478093E-05*T**2-2.21498670E-08*T**3
        +8.76225186E-12*T**4-5.10876533E-05*p
        +3.57823191E-09*p**2-1.37995983E-13*p**3
        +2.05230067E-18*p**4+2.52713006E-07*T*p
        -7.95947422E-12*T*p**2+8.30553293E-17*T*p**3
        -6.90483618E-10*T**2*p+2.36780380E-14*T**2*p**2
        -2.84627591E-19*T**2*p**3+6.06021155E-13*T**3*p
        -2.13606178E-17*T**3*p**2+2.75290092E-22*T**3*p**3)
        x_Ref=x_CO2
    elif Ref=='R744' and Liq =='PAG':
        Tmin=273.15+40
        Tmax=273.15+100
        #kPa to bar
        p/=100
        x_40c_pag=+9.099439359-1.777194511588e+01*p**(0.96261)+15.8417564431*p-0.00785617429085*p**2+9.29355174098e-06*p**3
        x_100c_pag=+326.6783782-3.230924121929e+02*p**(0.0067777)+0.284830938229*p-0.00120728637563*p**2+2.90768469832e-06*p**3
        x_Ref=((x_100c_pag-x_40c_pag)/60*(T-313.15)+x_40c_pag)/100
    
    elif Ref=='R744' and Liq=='PAO':
        Tmin=273.15+40
        Tmax=273.15+100
        p/=100
        x_40c_pao=+0.08197431533-1.627640711357e+01*p**(0.99157)+16.0542422174*p-0.00239828193696*p**2+1.88834705177e-06*p**3
        x_100c_pao=+139.8515113-1.401455321692e+02*p**(0.001384)+0.117510290184*p-0.000176084086229*p**2+4.04345850972e-07*p**3
        x_Ref=((x_100c_pao-x_40c_pao)/60*(T-313.15)+x_40c_pao)/100

    elif Ref=='R744' and Liq=='POE':
        Tmin=273.15+40
        Tmax=273.15+100
        p/=100
        x_40c_poe=-0.07347897661-9.364025183662e+00*p**(0.99171)+9.34087009398*p+0.00230160520305*p**2-1.73258741724e-05*p**3
        x_100c_poe=-0.161283727+1.134483722323e+01*p**(0.99097)-10.9522884211*p+0.00316840572421*p**2-9.72225995465e-06*p**3
        x_Ref=((x_100c_poe-x_40c_poe)/60*(T-313.15)+x_40c_poe)/100
    elif Ref=='R407C' and Liq=='POE32':
        """
        Harms 2002 "Charge Inventory system modeling and validation for unitary air conditioners"
        xi = solubility
        xi = m_diss/(m_diss + m_oil)
        sigmoid curves p = -(c1/2)+c1/(1-exp(xi/c2))
        
        for 0 < xi < 0.25
        p [kPa]
        T [C]
        """        
        xi = p/(990 + 91.9*T + 0.633*T**2)
        if xi > 0.25:
            raise Exception('Check solubility range R407C in POE')
        x_Ref = xi
        
    elif Ref=='R410A' and Liq=='POE32':
        x_Ref = 0.0
        
    else:
        print ("Ref/Liquid [%s/%s] not implemented" %(Ref,Liq))
        x_Ref=0.0
    if T<Tmin or T>Tmax:
#        print "Error: T[%0.2f] out of range [%0.2f,%0.2f]" %(T,Tmin,Tmax)
        error=True
    return x_Ref,error

if __name__== "__main__":

    print ('cp_l:',cp_oil("PAO",373.))
    print ('ACD1000FY s_l:', s_oil("ACD100FY",370))
    print ('POE s_l:', s_oil("POE",370))
    print ('Duratherm_LT s_l', s_oil("Duratherm_LT",370))
    print ('POE h_l:', h_oil("POE",370,200))
    print ('ACD1000FY h_l:', h_oil("ACD100FY",340,200))
    print ('Duratherm_LT h_l', h_oil("Duratherm_LT",370,200))
    print ('POE u_l:', u_oil("POE",370))
    print ('Duratherm_LT u_l', u_oil("Duratherm_LT",370))
    print ('ACD1000FY mu_l:', mu_oil("ACD100FY",370))
    print ('POE mu_l:', mu_oil("POE",370))

    
    
    