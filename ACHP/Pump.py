from __future__ import division, print_function, absolute_import
import CoolProp as CP

class PumpClass():
    def __init__(self,**kwargs):
        #Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    def Update(self,**kwargs):
        #Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        return [
            ('Loop pressure','kPa',self.pin_g),
            ('Overall Efficiency','-',self.eta),
            ('Pump power','W',self.W),
            ('Pressure lift','Pa',self.DP_g),
            ('Mass flow rate','kg/s',self.mdot_g)
         ]
        
    def Calculate(self):
        #AbstractState
        AS_g = self.AS_g
        if hasattr(self,'MassFrac_g'):
            AS_g.set_mass_fractions([self.MassFrac_g])
        elif hasattr(self, 'VoluFrac_SLF'):
            AS_g.set_volu_fractions([self.VoluFrac_g])
        
        AS_g.update(CP.PT_INPUTS, self.pin_g, self.Tin_g)
        rho=AS_g.rhomass() #[kg/m^3]
        self.W=abs(self.DP_g)*(self.mdot_g/rho)/self.eta