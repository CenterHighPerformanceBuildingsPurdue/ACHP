import sys
sys.path.append('../')

import SecondaryCycle as sc
In=sc.CycleInputVals()
In.CondAirRH=51.0
In.CondAirVdot=1.4
Out=sc.CycleOutputVals()


sc.SecondaryCycle(In,Out)
print Out.COP