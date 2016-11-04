
import cProfile,subprocess
from PHEHX import SamplePHEHX
from SampleCycles import SampleSecondaryLoopSystem,SampleDXACSystem 
from MCECycles import SampleDXMCEACSystem
#command ="""SamplePHEHX()"""
command ="""SampleDXACSystem()"""
#command ="""SampleDXMCEACSystem()"""
#command ="""SampleSecondaryLoopSystem()"""
#command ="""import Condenser; Condenser.SampleCondenser(40)"""
#command="""from TestSatLUT import d; d()"""
cProfile.runctx(command,globals(),locals(),filename="profile.txt")
subprocess.call(['runsnake','profile.txt'])

from time import time

#t1=time()
#SampleDXACSystem()
#t2=time()
#print t2-t1