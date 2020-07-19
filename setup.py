from __future__ import print_function

from distutils.core import setup
import subprocess,shutil,os,sys
import setuptools

#you ran the file, you probably wanted to run install, append that command
if len(sys.argv)==1:
    print("I think you wanted to install")
    sys.argv.append('develop')

setup (name = 'ACHP',
       version = '1.6',
       author      = "Ian Bell,Ammar Bahman,Davide Ziviani",
       author_email='ian.h.bell@gmail.com;bahman.ammar@gmail.com;dziviani@purdue.edu',
       url='http://achp.sourceforge.net',
       description = """ACHP - Air Conditioning Heat Pump Model""",
       packages=['ACHP'],
       zip_safe=False
       )
