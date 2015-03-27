Installation and Setup
============================

Binaries for Windows Users (Easy way):
--------------------------------------
32-bit binaries including all other required packages are available from https://sourceforge.net/projects/achp/files/.  Download the most recent zip file, and expand it somewhere convenient.  If you have other configuration files, place them in the *configs* folder which is next to the ACHP.exe file

This is the easiest way to get started!

Using the source for Developers (Not for the faint of heart):
-------------------------------------------------------------

**Step 1.**  
Download a distribution of Python.  For Windows users I strongly recommend the `Python(x,y) <http://www.pythonxy.com>`_ distribution. Go to downloads and then pick the full install.  I currently use the 2.6.6 version.  When you are installing, I recommend going into the *other* section at the bottom and adding SWIG and mingw.  This will allow you to recompile CoolProp if necessary.  Otherwise the defaults should be fine.

For Mac OSX users, the Enthought Python Distribution `(EPD) <http://www.enthought.com/products/getepd.php>`_ is recommended.  It is free to use in 32-bit mode for academic users.
    
For masochists, you could download Python, and the other packages that you would need are scipy, numpy, matplotlib, and wx.  It is highly recommended to use one of the above distributions.

For linux users, use your package manager to install python and the required packages. 

**Step 2.**
Download a copy of CoolProp.  For Windows users, you can download an installer directly by going to `CoolProp files <http://sourceforge.net/projects/coolprop/files/CoolProp/>`_ and downloading the most recent file for your architecture.  Pretty much the only architecture that is directly supported is 32bit Windows.  But on other architectures, the requisite files are included.

**Step 3.** 
Download a subversion client.  This allows you to interface with the subversion repository that the code is maintained in.  You can use a command line subversion client (http://www.collab.net/downloads/subversion/ and download the command-line client), or use a graphical user interface.  The GUI interface is a lot easier to use.  I recommend `TortoiseSVN <http://tortoisesvn.net/downloads.html>`_.  If you successfully installed CoolProp in the previous step, skip to Step 5.

**Step 4.**
Download a copy of the CoolProp source (everyone except 32-bit Windows).  If you are using a command line subversion client, type::
	
	svn co http://coolprop.svn.sourceforge.net/svnroot/coolprop/trunk coolprop

This will pull all the CoolProp source code in the main development trunk from the servers and put it the folder coolprop.  If you are using TortoiseSVN, in a Windows explorer window, right click and select *SVN Checkout...* then put https://achp.svn.sourceforge.net/svnroot/achp/trunk in the URL of repository and change the path if you like.  Leave it to fully recursive and HEAD revision.  Open a command prompt in the folder that the source resides and type::

	python setup.py install
	
This will install CoolProp and put it in a location that Python will be able to find.

**Step 5.**
Download a copy of the ACHP source code.  If you are using TortoiseSVN, open an explorer window and right-click where you want the code to go.  Select *SVN Checkout...*, set the URL of repository to https://achp.svn.sourceforge.net/svnroot/achp/trunk and change the path if you like.  Leave it to fully recursive and HEAD revision.  In a few moments you will then have the most recent set of code on your computer.

**Step 6.**
You can now peruse the code.  The components can all be found in the PyACHP subfolder and the GUI is in the GUI folder.  If you want to begin development, I recommend the use of Eclipse with the plugin Pydev.  The source tree includes a Pydev project file.  Eclipse is included with the Python(x,y) distribution.  Eclipse can be found in the pythonxy folder in the Windows start menu.  Once Eclipse is opened, to get the Eclipse project opened, right click in the Pydev Package Explorer and select *Import...*.  Then select *Existing Projects into Workspace*, and in the next window, click on the browse button next to the *Select root directory* option button.  Browse to the folder that contains your ACHP source, and it should load up the project.

To run a given component, open the component in the Pydev Package explorer, click on the down-arrow on the button that looks like a play button in the toolbar, and go to Run as... and then select Python Run.  This should run the component in stand-alone mode.  If you want to run the cycle model, open Cycle.py and run it.

To start the graphical user interface, run GUI/ACHPMainFrame.py using the same Python run method as for a component.

**Step 7.** 
To update your code, either right click on the source folder and select SVN Update... if you are using TortoiseSVN, or open a command prompt in the main source folder and type::

	svn update