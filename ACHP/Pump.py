from __future__ import division
from CoolProp.CoolProp import PropsSI

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
        rho=PropsSI('D','T',self.Tin_g,'P',self.pin_g,self.Ref_g)
        self.W=abs(self.DP_g)*(self.mdot_g/rho)/self.eta