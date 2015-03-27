from __future__ import division #Make integer 3/2 give 1.5 in python 2.x
from CoolProp.CoolProp import PropsSI, IsFluidType #,UseSaturationLUT ,T_hp, 
from math import pi,log,sqrt,exp,cos,sin,tan,log10
from scipy.integrate import quad,quadrature,trapz,simps,fixed_quad
from scipy.optimize import brentq,fsolve
import numpy as np

try:
    import psyco
    psyco.full()
except ImportError:
    pass

#Machine precision
machine_eps=np.finfo(np.float).eps

def Phase_ph(Ref,p,h,Tbubble,Tdew,rhosatL,rhosatV):
    """
    Convenience function to return just the Phase rather than densities and phase
    """
    (rhoL,rhoV,Phase)=TrhoPhase_ph(Ref,p,h,Tbubble,Tdew,rhosatL,rhosatV)
    return Phase

def TrhoPhase_ph(Ref,p,h,Tbubble,Tdew,rhosatL=None,rhosatV=None):
    """
    Convenience function to find temperature, density, and phase of fluid as a function of pressure and enthalpy
    """
    #UseSaturationLUT(1)
    #h/=1000  #convert J/kg to kJ/kg since CoolProp uses kJ/kg >>> modification: CoolProp 5.x uses J/kg
    
    if IsFluidType(Ref,'Brine')==1:
        #It is subcooled
        # Use a guess of 10 degrees below max temp
        #Tguess=PropsSI('M','T',0,'P',0, Ref)-10
        T=PropsSI('T','H',h,'P',p,Ref)  #T_hp(Ref,h,p,Tguess)
        rho=PropsSI('D','T',T,'P',p,Ref)
        return T,rho,'Subcooled'
        
    #Check if it is supercritical
    if p>PropsSI(Ref,'pcrit'):                                                  #Updated from PropsSI('E','T',0,'P',0, Ref):
        return None,None,'Supercritical'
    #It is not supercritical
    else:
        if rhosatL==None:
            rhosatL=PropsSI('D','T',Tbubble,'Q',0.0,Ref)
            rhosatV=PropsSI('D','T',Tdew,'Q',1.0,Ref)
        vsatL=1/rhosatL
        vsatV=1/rhosatV
        hsatL=PropsSI('H','T',Tbubble,'D',rhosatL,Ref)
        hsatV=PropsSI('H','T',Tdew,'D',rhosatV,Ref)
        
        if h>hsatV:
            #It's superheated
            cp=PropsSI('C','T',Tdew,'D',rhosatV,Ref)
            Tguess=Tdew+(h-hsatV)/cp
            T=PropsSI('T','H',h,'P',p,Ref)   #T_hp(Ref,h,p,Tguess)
            rho=PropsSI('D','T',T,'P',p,Ref)
            return T,rho,'Superheated'
        elif h<hsatL:
            # It's subcooled
            cp=PropsSI('C','T',Tbubble,'D',rhosatL,Ref)
            Tguess=Tbubble-(hsatL-h)/cp
            T=PropsSI('T','H',h,'P',p, Ref)  #T_hp(Ref,h,p,Tguess)
            rho=PropsSI('D','T',T,'P',p,Ref)
            return T,rho,'Subcooled'
        else:
            #It's two-phase
            x=(h-hsatL)/(hsatV-hsatL)
            v=x*vsatV+(1-x)*vsatL
            T=x*Tdew+(1-x)*Tbubble
            rho=1/v
            return T,rho,'TwoPhase'

def TwoPhaseDensity(Ref,xmin,xmax,Tdew,Tbubble,slipModel='Zivi'):
    rhog=PropsSI('D', 'T', Tdew, 'Q', 1, Ref)
    rhof=PropsSI('D', 'T', Tbubble, 'Q', 0, Ref)

    if slipModel=='Zivi':
        S=pow(rhof/rhog,0.3333)
    elif slipModel=='Homogeneous':
        S=1
    else:
        raise ValueError("slipModel must be either 'Zivi' or 'Homogeneous'")
    C=S*rhog/rhof

    if xmin+5*machine_eps<0 or xmax-5*machine_eps>1.0:
        raise ValueError('Quality must be between 0 and 1')
    #Avoid the zero and one qualities (undefined integral)
    if xmin==xmax:
        alpha_average=1/(1+C*(1-xmin)/xmin)
    else:    
        if xmin>=1.0:
            alpha_average=1.0
        elif xmax<=0.0:
            alpha_average=0.0
        else:
            alpha_average=-(C*(log( ((xmax-1.0)*C-xmax)/((xmin-1.0)*C-xmin) )+xmax-xmin)-xmax+xmin)/(C**2-2*C+1)/(xmax-xmin)
    return alpha_average*rhog + (1-alpha_average)*rhof

def AccelPressureDrop(x_min,x_max,Ref,G,Tbubble,Tdew,rhosatL=None,rhosatV=None,slipModel='Zivi'):
    """
    Accelerational pressure drop
    
    From -dpdz|A=G^2*d[x^2v_g/alpha+(1-x)^2*v_f/(1-alpha)^2]/dz
    
    Integrating over z from 0 to L where x=x_1 at z=0 and x=x_2 at z=L
    
    Maxima code:
        alpha:1/(1+S*rho_g/rho_f*(1-x)/x)$
        num1:x^2/rho_g$
        num2:(1-x)^2/rho_f$
        subst(num1/alpha+num2/(1-alpha),x,1);
        subst(num1/alpha+num2/(1-alpha),x,0);
    """
    if rhosatL==None or rhosatV==None:
        rhosatL=PropsSI('D','T',Tbubble,'Q',0.0,Ref)
        rhosatV=PropsSI('D','T',Tdew,'Q',1.0,Ref)
        
    def f(x,rhoL,rhoV):
        if abs(x)<1e-12:
            return 1/rhosatL
        elif abs(1-x)<1e-12:
            return 1/rhosatV
        else:
            if slipModel=='Zivi':
                S=pow(rhoL/rhoV,1/3)
            elif slipModel=='Homogeneous':
                S=1
            else:
                raise ValueError("slipModel must be either 'Zivi' or 'Homogeneous'")
            alpha=1/(1+S*rhoV/rhoL*(1-x)/x)
            return x**2/rhoV/alpha+(1-x)**2/rhoL/(1-alpha)
    return G**2*(f(x_min,rhosatL,rhosatV)-f(x_max,rhosatL,rhosatV))
        
def LMPressureGradientAvg(x_min,x_max,Ref,G,D,Tbubble,Tdew,C=None,satTransport=None):
    """
    Returns the average pressure gradient between qualities of x_min and x_max.
    
    To obtain the pressure gradient for a given value of x, pass it in as x_min and x_max
    
    Required parameters:
    * x_min : The minimum quality for the range [-]
    * x_max : The maximum quality for the range [-]
    * Ref : String with the refrigerant name
    * G : Mass flux [kg/m^2/s]
    * D : Diameter of tube [m]
    * Tbubble : Bubblepoint temperature of refrigerant [K]
    * Tdew : Dewpoint temperature of refrigerant [K]
    
    Optional parameters:
    * C : The coefficient in the pressure drop
    * satTransport : A dictionary with the keys 'mu_f','mu_g,'v_f','v_g' for the saturation properties.  So they can be calculated once and passed in for a slight improvement in efficiency 
    """
    def LMFunc(x):
        dpdz,alpha=LockhartMartinelli(Ref,G,D,x,Tbubble,Tdew,C,satTransport)
        return dpdz
    
    ## Use Simpson's Rule to calculate the average pressure gradient
    ## Can't use adapative quadrature since function is not sufficiently smooth
    ## Not clear why not sufficiently smooth at x>0.9
    if x_min==x_max:
        return LMFunc(x_min)
    else:
        #Calculate the tranport properties once
        satTransport={}
        satTransport['v_f']=1/PropsSI('D','T',Tbubble,'Q',0.0,Ref)
        satTransport['v_g']=1/PropsSI('D','T',Tdew,'Q',1.0,Ref)
        satTransport['mu_f']=PropsSI('V','T',Tbubble,'Q',0.0,Ref)
        satTransport['mu_g']=PropsSI('V','T',Tdew,'Q',1.0,Ref)
        
        xx=np.linspace(x_min,x_max,30)
        DP=np.zeros_like(xx)
        for i in range(len(xx)):
            DP[i]=LMFunc(xx[i])
        return -simps(DP,xx)/(x_max-x_min)

def LockhartMartinelli(Ref, G, D, x, Tbubble,Tdew,C=None,satTransport=None):
    # Following the method laid out in ME506 notes on 
    # Separated Flow pressure drop calculations

    #Convert the quality, which might come in as a single numpy float value, to a float
    #With the conversion, >20x speedup in the LockhartMartinelli function, not clear why
    x=float(x)
    
    #UseSaturationLUT(1)
    if satTransport==None:
        # Calculate Necessary saturation properties
        v_f=1/PropsSI('D','T',Tbubble,'Q',0.0,Ref)
        v_g=1/PropsSI('D','T',Tdew,'Q',1.0,Ref)
        mu_f=PropsSI('V','T',Tbubble,'Q',0.0,Ref)
        mu_g=PropsSI('V','T',Tdew,'Q',1.0,Ref)
    else:
        #Pull out of the dictionary
        v_f=satTransport['v_f']
        v_g=satTransport['v_g']
        mu_f=satTransport['mu_f']
        mu_g=satTransport['mu_g']
        
    # 1. Find the Reynolds Number for each phase based on the actual flow rate of the individual phase
    Re_g=G*x*D/mu_g
    Re_f=G*(1-x)*D/mu_f

    # 2. Friction factor for each phase
    if x==1: #No liquid
        f_f=0 #Just to be ok until next step
    elif Re_f<1000: #Laminar
        f_f=16/Re_f
    elif Re_f>2000: #Fully-Turbulent
        f_f=0.046/(Re_f**0.2)
    else: # Mixed
        # Weighting factor
        w=(Re_f-1000)/(2000-1000)
        # Linear interpolation between laminar and turbulent
        f_f=w*16.0/Re_f+(1-w)*0.046/(Re_f**0.2)

    if x==0: #No gas
        f_g=0 #Just to be ok until next step
    elif Re_g<1000: #Laminar
        f_g=16.0/Re_g
    elif Re_g>2000: # Fully-Turbulent
        f_g=0.046/(Re_g**0.2)
    else: # Mixed
        # Weighting factor
        w=(Re_g-1000)/(2000-1000)
        # Linear interpolation between laminar and turbulent
        f_g=w*16.0/Re_g+(1-w)*0.046/(Re_g**0.2)

    # 3. Frictional pressure drop based on actual flow rate of each phase
    dpdz_f=2*f_f*G**2*(1-x)**2*v_f/D
    dpdz_g=2*f_g*G**2*x**2*v_g/D

    if x<=0:
        # Entirely liquid
        alpha=0.0
        dpdz=dpdz_f
        return dpdz,alpha
    if x>=1:
        #Entirely vapor
        alpha=1.0
        dpdz=dpdz_g
        return dpdz,alpha

    # 4. Lockhart-Martinelli parameter
    X=sqrt(dpdz_f/dpdz_g)

    # 5. Find the Constant based on the flow Re of each phase
    #    (using 1500 as the transitional Re to ensure continuity)
    
    #Calculate C if not passed in:
    if C==None:
        if Re_f>1500 and Re_g > 1500:
            C=20.0
        elif Re_f<1500 and Re_g>1500:
            C=12.0
        elif Re_f>1500 and Re_g<1500:
            C=10.0
        else:
            C=5.0

    # 6. Two-phase multipliers for each phase
    phi_g2=1+C*X+X**2
    phi_f2=1+C/X+1/X**2
        
    # 7. Find gradient
    if dpdz_g*phi_g2>dpdz_f*phi_f2:
        dpdz=dpdz_g*phi_g2
    else:
        dpdz=dpdz_f*phi_f2

    # 8. Void Fraction
    alpha=1-X/sqrt(X*X+20*X+1)
    
    return dpdz,alpha

def ShahEvaporation_Average(x_min,x_max,Ref,G,D,p,q_flux,Tbubble,Tdew):
    """
    Returns the average pressure gradient between qualities of x_min and x_max.
    
    To obtain the pressure gradient for a given value of x, pass it in as x_min and x_max
    
    Required parameters:
    * x_min : The minimum quality for the range [-]
    * x_max : The maximum quality for the range [-]
    * Ref : String with the refrigerant name
    * G : Mass flux [kg/m^2/s]
    * D : Diameter of tube [m]
    * p : Pressure [Pa]
    * q_flux : Heat transfer flux [W/m^2]
    * Tbubble : Bubblepoint temperature of refrigerant [K]
    * Tdew : Dewpoint temperature of refrigerant [K]
    """
    # ********************************
    #        Necessary Properties
    # ********************************
    v_g = 1 / PropsSI('D', 'T', Tdew, 'Q', 1, Ref)# [m^3/kg]
    v_f = 1 / PropsSI('D', 'T', Tbubble, 'Q', 0, Ref)# [m^3/kg]
    rho_f = 1 / v_f# [kg/m^3]
    rho_g = 1 / v_g# [kg/m^3]
    mu_f = PropsSI('V', 'T', Tbubble, 'Q', 0, Ref)# [kg/m-s] 
    mu_g = PropsSI('V', 'T', Tdew, 'Q', 1, Ref)# [kg/m-s] 
    h_fg = (PropsSI('H', 'T', Tdew, 'Q', 1, Ref) - PropsSI('H', 'T', Tbubble, 'Q', 0, Ref))#*1000 #[J/kg]
    cp_f = PropsSI('C', 'T', Tbubble, 'Q', 0, Ref)#*1000 # [J/kg-K]
    cp_g = PropsSI('C', 'T', Tdew, 'Q', 1, Ref)#*1000 # [J/kg-K]
    k_f = PropsSI('L', 'T', Tbubble, 'Q', 0, Ref)#*1000 # [W/m-K]
    k_g = PropsSI('L', 'T', Tdew, 'Q', 1, Ref)#*1000 # [W/m-K]
    Pr_f = cp_f * mu_f / k_f #[-]
    Pr_g = cp_g * mu_g / k_g #[-]

    g_grav = 9.81 #[m/s^2]

    # Shah evaporation correlation
    Fr_L = G**2 / (rho_f*rho_f * g_grav * D) #[-]
    Bo = q_flux / (G * h_fg) #[-]
    
    if Bo < 0:
        raise ValueError('Heat flux for Shah Evaporation must be positive')
    
    if Bo > 0.0011:
        F = 14.7
    else:
        F = 15.43
    #Pure vapor single-phase heat transfer coefficient
    h_g = 0.023 * (G*D/mu_g)**(0.8) * Pr_g**(0.4) * k_g / D #[W/m^2-K]    
    def ShahEvaporation(x):
        if abs(1-x)<5*machine_eps:
            return h_g
            
        #If the quality is above 0.999, linearly interpolate to avoid division by zero
        if x>0.999:
            h_1=ShahEvaporation(1.0) #Fully fry
            h_999=ShahEvaporation(0.999) #At a quality of 0.999
            return (h_1-h_999)/(1.0-0.999)*(x-0.999)+h_999 #Linear interpolation
        if abs(x)<5*machine_eps:
            h_L = 0.023 * (G*(1 - x)*D/mu_f)**(0.8) * Pr_f**(0.4) * k_f / D #[W/m^2-K]
            return h_L
        else:
            h_L = 0.023 * (G*(1 - x)*D/mu_f)**(0.8) * Pr_f**(0.4) * k_f / D #[W/m^2-K]
        Co = (1 / x - 1)**(0.8) * (rho_g / rho_f)**(0.5) #[-]

        if Fr_L >= 0.04:
            N = Co
        else:
            N = 0.38 * Fr_L**(-0.3) * Co

        psi_cb = 1.8 / N**(0.8)
        if (0.1 < N and N <= 1.0):
            psi_bs = F * (Bo)**(0.5) * exp(2.74 * N**(-0.1))
            psi = max([psi_bs, psi_cb])
        else:
            if (N > 1.0):
                if (Bo > 0.00003):
                    psi_nb = 230 * (Bo)**(0.5)
                else:
                    psi_nb = 1.0 + 46.0 * (Bo)**(0.5)
                psi = max([psi_nb,psi_cb])
            else:
                psi_bs = F * (Bo)**(0.5) * exp(2.47 * N**(-0.15))
                psi = max([psi_bs, psi_cb])
        return psi * h_L #[W/m^2-K]
    
    #Calculate h over the range of x
    x=np.linspace(x_min,x_max,10)
    h=np.zeros_like(x)
    for i in range(len(x)):
        h[i]=ShahEvaporation(x[i])
    
    #if x_min == x_max, or they are really really close to being the same
    if abs(x_max-x_min)<5*machine_eps:
        #return just one of the edge values
        return h[0]
    else:
        #Use Simpson's rule to carry out numerical integration to get average
        return simps(h,x)/(x_max-x_min)
def LongoCondensation(x_avg,G,dh,Ref,TsatL,TsatV):
    rho_L = PropsSI('D', 'T', TsatL, 'Q', 0, Ref) #kg/m^3
    rho_V = PropsSI('D', 'T', TsatV, 'Q', 1, Ref) #kg/m^3
    mu_L = PropsSI('V', 'T', TsatL, 'Q', 0, Ref) #kg/m-s 
    cp_L = PropsSI('C', 'T', TsatL, 'Q', 0, Ref)#*1000 #J/kg-K
    k_L = PropsSI('L', 'T', TsatV, 'Q', 0, Ref)#*1000 #W/m-K
    Pr_L = cp_L * mu_L / k_L #[-]
    
    Re_eq=G*((1-x_avg)+x_avg*sqrt(rho_L/rho_V))*dh/mu_L
    
    if Re_eq<1750:
        Nu=60*Pr_L**(1/3)
    else:
        Nu=((75-60)/(3000-1750)*(Re_eq-1750)+60)*Pr_L**(1/3)
    h=Nu*k_L/dh
    return h
    
def ShahCondensation_Average(x_min,x_max,Ref,G,D,p,TsatL,TsatV):
    # ********************************
    #        Necessary Properties
    #    Calculated outside the quadrature integration for speed
    # ********************************
    mu_f = PropsSI('V', 'T', TsatL, 'Q', 0, Ref) #kg/m-s 
    cp_f = PropsSI('C', 'T', TsatL, 'Q', 0, Ref)#*1000 #J/kg-K
    k_f = PropsSI('L', 'T', TsatV, 'Q', 0, Ref)#*1000 #W/m-K
    Pr_f = cp_f * mu_f / k_f #[-]
    Pstar = p / PropsSI(Ref,'pcrit')
    h_L = 0.023 * (G*D/mu_f)**(0.8) * Pr_f**(0.4) * k_f / D #[W/m^2-K]
    def ShahCondensation(x,Ref,G,D,p):
        return h_L * ((1 - x)**(0.8) + (3.8 * x**(0.76) * (1 - x)**(0.04)) / (Pstar**(0.38)) )
        
    if not x_min==x_max:
        #A proper range is given
        return quad(ShahCondensation,x_min,x_max,args=(Ref,G,D,p))[0]/(x_max-x_min)
    else:
        #A single value is given
        return ShahCondensation(x_min,Ref,G,D,p)
    
def f_h_1phase_Tube(mdot,ID,T, p,Fluid,Phase='Single'):
    """ 
    Convenience function to run annular model for tube.  Tube is a degenerate case of annulus with inner diameter of 0
    
    """
    return f_h_1phase_Annulus(mdot, ID, 0.0, T, p, Fluid, Phase='Single' )

def f_h_1phase_Annulus(mdot, OD, ID, T, p, Fluid, Phase='Single'):
    """
    
    """
    if Phase =="SatVap":
        mu = PropsSI('V', 'T', T, 'Q', 1, Fluid) #kg/m-s
        cp = PropsSI('C', 'T', T, 'Q', 1, Fluid)*1. #*1000. #J/kg-K
        k = PropsSI('L', 'T', T, 'Q', 1, Fluid)*1. #*1000. #W/m-K
        rho = PropsSI('D', 'T', T, 'Q', 1, Fluid) #kg/m^3
    else:
        mu = PropsSI('V', 'T', T, 'P', p, Fluid)  #kg/m-s
        cp = PropsSI('C', 'T', T, 'P', p, Fluid)*1. #*1000. #J/kg-K
        k = PropsSI('L', 'T', T, 'P', p, Fluid)*1. #*1000. #W/m-K
        rho = PropsSI('D', 'T', T, 'P', p, Fluid) #kg/m^3

    Pr = cp * mu / k #[-]

    Dh = OD - ID
    Area=pi*(OD**2-ID**2)/4.0
    u=mdot/(Area*rho)
    Re=rho*u*Dh/mu

    # Friction factor of Churchill (Darcy Friction factor where f_laminar=64/Re)
    e_D = 0
    A = ((-2.457 * log( (7.0 / Re)**(0.9) + 0.27 * e_D)))**16
    B = (37530.0 / Re)**16
    f = 8 * ((8/Re)**12.0 + 1 / (A + B)**(1.5))**(1/12)

    # Heat Transfer coefficient of Gnielinski
    Nu = (f/8)*(Re-1000)*Pr/(1+12.7*sqrt(f/8)*(Pr**(0.66666)-1)) #[-]
    h = k*Nu/Dh #W/m^2-K
    return (f, h, Re)

def f_h_1phase_Channel(mdot,W,H,T,p,Fluid,Phase):

    if Phase=="SatVap":
        mu = PropsSI('V', 'T', T, 'Q', 1, Fluid) #kg/m-s
        cp = PropsSI('C', 'T', T, 'Q', 1, Fluid)*1. #*1000. #J/kg-K
        k = PropsSI('L', 'T', T, 'Q', 1, Fluid)*1. #*1000. #W/m-K
        rho = PropsSI('D', 'T', T, 'Q', 1, Fluid) #kg/m^3
    else:
        mu = PropsSI('V', 'T', T, 'P', p, Fluid) ##kg/m-s
        cp = PropsSI('C', 'T', T, 'P', p, Fluid)*1. #*1000. #J/kg-K
        k = PropsSI('L', 'T', T, 'P', p, Fluid)*1. #*1000. #W/m-K
        rho = PropsSI('D', 'T', T, 'P', p, Fluid) #kg/m^3

    Pr = cp * mu / k #[-]

    Dh = 2*H*W/(H+W)
    Area=W*H
    u=mdot/(Area*rho)
    Re=rho*u*(Dh)/mu

    # Friction factor of Churchill (Darcy Friction factor where f_laminar=64/Re)
    e_D = 0
    A = ((-2.457 * log( (7.0 / Re)**(0.9) + 0.27 * e_D)))**16
    B = (37530.0 / Re)**16
    f = 8 * ((8/Re)**12.0 + 1 / (A + B)**(1.5))**(1/12)

    # Heat Transfer coefficient of Gnielinski for high Re, 
    if Re>1000:
        Nu = (f / 8) * (Re - 1000.) * Pr / (1 + 12.7 * sqrt(f / 8) * (Pr**(0.66666) - 1)) #[-]
    else:
        Nu = 3.66
    h = k * Nu / Dh #W/m^2-K
    return (f, h, Re)

def PHE_1phase_hdP(Inputs,JustGeo=False):
    """   
    Based on the single-phase pressure drop and heat transfer correlations
    in VDI Heat Atlas Chapter N6: Pressure Drop and Heat Transfer in Plate Heat 
    Exchangers by Holger Martin DOI: 10.1007/978-3-540-77877-6_66 Springer Verlag
    ::
        =============================        
        ||   __               __    ||
        ||  /  \             /  \   ||
        || |    |           |    |  ||  ===
        ||  \__/             \__/   ||   |
        ||                          ||   |
        ||             | <-  B   -> ||   |
        ||                          ||   |
        ||                          ||   |
        ||                          ||
        ||                          ||
        ||             |\           ||
        ||             | \          ||   Lp
        ||             |  \         ||  
        ||             |   \        ||
        ||             |phi \       ||
        ||             |     \      ||   |
        ||                          ||   |
        ||   __               __    ||   |
        ||  /  \             /  \   ||   |
        || |    |           |    |  ||  ===
        ||  \__/             \__/   ||
        ||                          ||
        =============================
         | -----      Bp  --------- |
         
         phi is the inclination angle
    """
        
        
    #Plate parameters
    PlateAmplitude = Inputs['PlateAmplitude']
    PlateWavelength = Inputs['PlateWavelength']
    InclinationAngle = Inputs['InclinationAngle']
    Bp = Inputs['Bp']
    Lp = Inputs['Lp']
    if JustGeo==False:        
        Ref = Inputs['Ref']
        T = Inputs['T']
        p = Inputs['p']
        mdot_gap = Inputs['mdot_gap']
    
    
    X=2*pi*PlateAmplitude/PlateWavelength
    PHI=1/6*(1+sqrt(1+X**2)+4*sqrt(1+X**2/2))
    
    #The plane surface between the ports
    A0=Bp*Lp
    
    #The plane surface of one plate
    Ap=PHI*A0
    
    #The volume of one channel
    Vchannel=Bp*Lp*2*PlateAmplitude
    
    #Hydraulic diameter
    dh=4*PlateAmplitude/PHI
    
    if JustGeo==True:
        return {'Ap':Ap,'Vchannel':Vchannel,'Aflow':2*PlateAmplitude*Bp,'Dh':dh,'PHI':PHI}
    else:
        #Also calculate the thermodynamics and pressure drop
        
        #"Glycol" properties
        rho_g=PropsSI('D','T',T,'P',p,Ref)
        eta_g=PropsSI('V','T',T,'P',p,Ref)
        cp_g=PropsSI('C','T',T,'P',p,Ref)#*1000
        k_g=PropsSI('L','T',T,'P',p,Ref)#*1000
        
        Pr_g=cp_g*eta_g/k_g
        
        eta_g_w=eta_g #TODO: allow for temperature dependence?
        w_g=mdot_gap/rho_g/(2*PlateAmplitude*Bp)
        Re_g=rho_g*w_g*dh/eta_g
        
        #Calculate the friction factor zeta
        phi=InclinationAngle
        
        if Re_g<2000:
            zeta0=64/Re_g
            zeta1_0=597/Re_g+3.85
        else:
            zeta0=(1.8*log(Re_g)-1.5)**(-2)
            zeta1_0=39/Re_g**0.289
        
        a=3.8
        b=0.18
        c=0.36
        
        zeta1=a*zeta1_0
        #RHS from Equation 18
        RHS=cos(phi)/sqrt(b*tan(phi)+c*sin(phi)+zeta0/cos(phi))+(1-cos(phi))/sqrt(zeta1)
        zeta=1/RHS**2
        #Hagen number
        Hg=zeta*Re_g**2/2
        
        #Constants for Nu correlation
        c_q=0.122
        q=0.39#q=0.374
        #Nusselt number [-]
        Nu=c_q*Pr_g**(1/3)*(eta_g/eta_g_w)**(1/6)*(2*Hg*sin(2*phi))**(q)
        
        #Heat transfer coefficient [W/m^2-K]
        h=Nu*k_g/dh
        
        #Pressure drop 
        DELTAP=Hg*eta_g**2*Lp/(rho_g*dh**3)
        
        # There are quite a lot of things that might be useful to have access to
        # in outer functions, so pack up parameters into a dictionary
        Outputs={
             'Dh':dh,
             'h':h,
             'Ap':Ap,
             'DELTAP':DELTAP,
             'Re': Re_g,
             'U': w_g,
             'k': k_g,
             'cp': cp_g,
             'Vchannel':Vchannel,
             'Aflow':2*PlateAmplitude*Bp
        }
        return Outputs

def Cooper_PoolBoiling(pstar,Rp,q,M):
    return 55*pstar**(0.12-0.2*log10(Rp))*(-log10(pstar))**(-0.55)*q**(0.67)*M**(-0.5)

def KandlikarPHE(Ref,xmean,G,D,q,Tbubble,Tdew):
    """
    From http://www.rit.edu/kgcoe/mechanical/taleme/Papers/Conference%20Papers/C041.pdf
    
    Not recommended for fluids other than R134a
    """
    rhoG=PropsSI('D','T',Tdew,'Q',1,Ref)
    rhoL=PropsSI('D','T',Tbubble,'Q',0,Ref)
    mu_f = PropsSI('V', 'T', Tbubble, 'Q', 0, Ref) #kg/m-s 
    cp_f = PropsSI('C', 'T', Tbubble, 'Q', 0, Ref)#*1000 #J/kg-K
    k_f = PropsSI('L', 'T', Tdew, 'Q', 0, Ref)#*1000 #W/m-K
    Pr_f = cp_f * mu_f / k_f #[-]
    h_LG = (PropsSI('H','T',Tdew,'Q',1,Ref)-PropsSI('H','T',Tbubble,'Q',0,Ref))#*1000
    alpha_L = 0.023 * (G*D/mu_f)**(0.8) * Pr_f**(0.4) * k_f / D #[W/m^2-K]
    Co=(rhoG/rhoL)**(0.5)*((1-xmean)/xmean)**(0.8)
    Bo=q/(G*h_LG)
#    E_CB=0.512
#    E_NB=0.338
#    F_fl=1.0
    #alpha_r=(2.312*Co**(-0.3)*E_CB+667.3*Bo**(2.8)*F_fl*E_NB)*(1-xmean)**(0.003)*alpha_L
    alpha_r=1.055*(1.056*Co**(-0.4)+1.02*Bo**(0.9))*xmean**(-0.12)*alpha_L**(0.98)
    return alpha_r

def Bertsch_MC(x,Ref,G,Dh,q,L,Tbubble,Tdew):
    p=PropsSI('P','T',Tdew,'Q',1.0,Ref)
    pc=PropsSI(Ref,'pcrit')                                                     #Updated from PropsSI('E','T',0,'P',0,Ref)
    pr=p/pc
    M=PropsSI('M','T',0,'P',0,Ref)
    g=9.81

    if Ref=='R290':
        sig=55.28*(1-Tdew/369.818)**(1.258)/1000.
    elif Ref=='R410A':
        ## From Okada 1999 "Surface Tension of HFC Refrigerant Mixtures"
        sig=62.38*(1-Tdew/344.56)**(1.246)/1000.
    
    k_L=PropsSI('L','T',Tbubble,'Q',0,Ref)#*1000
    k_G=PropsSI('L','T',Tdew,'Q',1,Ref)#*1000
    cp_L=PropsSI('C','T',Tbubble,'Q',0,Ref)#*1000
    cp_G=PropsSI('C','T',Tdew,'Q',1,Ref)#*1000
    rho_L=PropsSI('D','T',Tbubble,'Q',0,Ref)
    rho_G=PropsSI('D','T',Tdew,'Q',1,Ref)
    mu_L=PropsSI('V','T',Tbubble,'Q',0,Ref)
    mu_G=PropsSI('L','T',Tdew,'Q',1,Ref)
    Re_L=G*Dh/mu_L
    Re_G=G*Dh/mu_G
    Pr_L=cp_L*mu_G/k_L
    Pr_G=cp_G*mu_G/k_G

    h_nb=55*(pr)**(0.12)*(-log10(pr))**(-0.55)*M**(-0.5)*q**(0.67)
    h_conv_l=(3.66+(0.0668*Dh/L*Re_L*Pr_L)/(1+0.04*(Dh/L*Re_L*Pr_L)**(2.0/3.0)))*k_L/Dh
    h_conv_g=(3.66+(0.0668*Dh/L*Re_G*Pr_G)/(1+0.04*(Dh/L*Re_G*Pr_G)**(2.0/3.0)))*k_G/Dh
    h_conv_tp=h_conv_l*(1-x)+h_conv_g*x
    Co=sqrt(sig/(g*(rho_L-rho_G)*Dh**2))
    h_TP=h_nb*(1-x)+h_conv_tp*(1.0+80.0*(x**2-x**6)*exp(-0.6*Co))
    return h_TP
if __name__=='__main__':
    DP_vals_acc=[]
    DP_vals_fric=[]
    x_vals=[]
    import numpy as np,pylab
    for x in np.linspace(0.1,1.0,10):  
        DP_vals_acc.append(AccelPressureDrop(x-0.1,x,'R410A',2,250,250))
        DP_vals_fric.append(LMPressureGradientAvg(x-0.1,x,'R410A',0.1,0.01,250,250)*1*1)
        x_vals.append(x)
        
    print "plot shows accelerational pressure drop as f(x) for 0.1 x segments"
    pylab.plot(x_vals, DP_vals_acc)
    pylab.show()
    print "plot shows frictional pressure drop as f(x) for 0.1 x segments of a fictional tube with unit length"
    pylab.plot(x_vals, DP_vals_fric)
    pylab.show()
    
    