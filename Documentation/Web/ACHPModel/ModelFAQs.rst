
.. _FAQs:

Frequently Asked Questions
**************************

.. _FAQS-PythonPath:

Working with the Python Path (aka *Argh!! Why can't Python find my files?*)
===========================================================================

One of the biggest hurdles for anyone getting started with Python is dealing with the python path.  (Author's note: I *still* struggle with this!)  Python searches a number of locations when it is looking for a module from an import command.  In a script, if you write::

    from __future__ import print_function
    import sys
    print(sys.path)
    
it will tell you all the places that Python is looking for your files.  From IPython, you might get an output something like this for the first three entries:

.. ipython::

    In [1]: import sys; print(sys.path[0:3])

If you have multiple copies of a module, it will use the first one it finds as it searches from left to right in the list of folders.

Any modules that you have installed using either the::

    python setup.py install
    
command or with an installer executable should be sitting in the Python folder tree, and should ideally not require any further work from you.  `CoolProp <http://coolprop.sourceforge.net>`_, a required package for ACHP for instance, is installed into the::

    C:\\Python26\\Lib\\site-packages\\CoolProp

folder with Python 2.6.x.  And if you open an ipython prompt and type

.. ipython::

    In [1]: import CoolProp; print(CoolProp.__file__)
    
no errors should be generated, and it will spit out the path to the module files.
    
Updating Python Path
--------------------
If you have code in a file that Python cannot find right now, you have a few options:

#. Move the files to somewhere that Python CAN find (i.e. on the python path)
#. Change the PYTHONPATH environmental variable
#. Add the path to sys.path

In general, for most code my preferred method is to update the PYTHONPATH environmental variable for code that I am likely to reference from a number of other scripts (i.e. the ACHP code).  To do this in Windows, click on the windows symbol in the taskbar (the old start button), right click on Computer and go to properties, then advanced system settings, the Advanced tab, and then Environment Variables.  In the sytem variables window, there are two lists of environmental variables, user variables and system variables.  I tend to add the path to the folder that the files are in to both lists just to be sure.  It is not clear how python decides which list to use.  The list of folders is semicolon delimited.  If the environmental variable PYTHONPATH does not exist, just create a new environmental variable.

The nice part of using the PYTHONPATH solution is that you only have to do it once, and then any subsequent times you want the code, python can find it.

If on the other hand you want to just use the sys.path method, to add a path to the python path, just do something like::

    import sys
    sys.path.append('c:\path\to\folder')
    
which will add the file path to the system path for the current execution of code.