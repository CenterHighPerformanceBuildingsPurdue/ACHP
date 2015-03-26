# add the mpl mpl-data folder and rc file
import matplotlib as mpl

from distutils.core import setup
import py2exe,os,sys,subprocess,glob
sys.path.append(os.path.join('..','Props'))
sys.path.append('..')

# Remove the build folder, a bit slower but ensures that build contains the latest
import shutil
## shutil.rmtree("build", ignore_errors=True)
## shutil.rmtree("ACHP", ignore_errors=True)

# my setup.py is based on one generated with gui2exe, so data_files is done a bit differently
data_files = []
includes = []
excludes = ['_gtkagg', '_tkagg', 'bsddb', 'curses', 'pywin.debugger',
            'pywin.debugger.dbgcon', 'pywin.dialogs', 'tcl',
            'Tkconstants', 'Tkinter', 'pydoc', 'doctest', 'test', 'sqlite3',
            'projections','Tcl','PyQT4','qt','backend_qt','backend_qt4','backend_qt4agg',
            'backend_qtagg','backend_cairo','backend_cocoaagg'
            ]
packages = ['pytz']
dll_excludes = ['libgdk-win32-2.0-0.dll', 'libgobject-2.0-0.dll', 'tcl84.dll',
                'tk84.dll','QtGui4.dll','QtCore4.dll']
icon_resources = [(1,'SC.ico')]
bitmap_resources = []
other_resources = []

data_files += mpl.get_py2exe_datafiles()

#Get all the image files
data_files += [('imgs',glob.glob('imgs\\*.png'))]
data_files += [('imgs',glob.glob('imgs\\*.ico'))]
data_files += [('help',glob.glob('help\\help.pdf'))]

#Get all the config files
data_files += [('configs',glob.glob('configs\\*.cfg'))]

#Get all the compressor map files
data_files += [('comps',glob.glob('comps\\*.*'))]

#Get all the compressor map files
data_files += [('parametric',glob.glob('parametric\\*.txt'))]

setup(
    console=[
            {
                "script":'ACHP.py',
                "icon_resources":[(1,"sc.ico")]
            }
        ],
                          # compressed and optimize reduce the size
    options = {"py2exe": {"compressed": 2, 
                          "optimize": 2,
                          "includes": includes,
                          "excludes": excludes,
                          "packages": packages,
                          "dll_excludes": dll_excludes,
                          # using 2 to reduce number of files in dist folder
                          # using 1 is not recommended as it often does not work
                          "bundle_files": 2,
                          "dist_dir": 'dist',
                          "xref": False,
                          "skip_archive": False,
                          "ascii": False,
                          "custom_boot_script": '',
                         }
              },

    # using zipfile to reduce number of files in dist
    zipfile = r'lib\library.zip',
    
    data_files=data_files
)

exe_dir=os.path.join(os.path.abspath(os.path.curdir),'dist')
os.remove(os.path.join(exe_dir,'w9xpopen.exe'))
os.rename('dist','ACHP')
shutil.rmtree('build')
subprocess.call(['7z','a','ACHP1.4.zip','ACHP'])
shutil.rmtree('ACHP')