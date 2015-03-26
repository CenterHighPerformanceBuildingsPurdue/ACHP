Additional Functionality
========================

Using Configuration Files
-------------------------
If you have a configuration that you like and would like to use at a later time in ACHP, you can save a configuration file that can be loaded back into ACHP.  Once you have all the parameters as you would like them, go to *File->Save Config File...* or press Ctrl+s (Apple+s for OSX users).  This will pop up a file selection dialog, and you can save the file somewhere.  

At a later time, you can load the configuration file back into ACHP by going to *File->Load Conf File...* or pressing Ctrl+o (Apple+o for OSX users).  Again you will get a file selection dialog and you can open your configuration file.

The configuration files are simple text files that are delimited with the delimiter ":::".  The configuration files list each of the items that are in the user interface and their current value (Also see note for developers below).

*For Developers*: When defining the names of items in the user interface, the first few letters of the object name determine whether it will be written to the configuration file or not, and whether it will be loaded back into the interface.  The naming conventions from Visual Basic are used, so names starting with *opt* are option boxes, *cmb* are combo boxes, *lbl* are labels, *txt* are textboxes, *rad* are radio boxes, and *chk* are checkboxes.  Refer to :download:`LoadGUI <../../../GUI/LoadGUI.py>` to see how the user interface is constructed, and the functions ReadConfigFile and WriteConfigFile in :download:`ACHPMainFrame.py <../../../GUI/ACHPMainFrame.py>`

Parametric Studies
------------------
One of the features built into versions 1.3 and onwards of ACHP is the flexible multi-dimensional parametric study solver.  With this solver you can run multi-dimensional parametric studies.  In order to activate this mode, go to the solvers Tab and switch the *Solver Method* to Parametric Study.

In each of the rows of the parametric study, you can select a variable that you would like to vary by selecting the variable from the combobox, and turn it on by checking the checkbox.  You must then provide the values that will be used for the variable.  There are two ways of doing this:

	* Provide the minimum value in the *Min Value/List* column, the maximum value in the *Max Value* column, and the integer number of linearly spaced steps in the *Number Steps* column.  This will yield a linearly spaced set of inputs that will be used.
	
	* Provide a comma-separated list of inputs in the *Min Value/List* column, and put the letter "L" in the *Number Steps* column.  This allows you to put your own values in, particularly useful if you are trying to duplicate non-linearly-spaced rating data, which was the motivating factor for this functionality
	
Not all variables make sense for all configurations.  For instance, with a direct expansion air conditioning system, there is no secondary loop, so altering the secondary loop mass flow rate will not change anything, or may raise an error.

*For Developers*: If you want to add other things to the parametric study, you can open the file *parametric/params.txt* (path relative to ACHP.exe) and add or delete entries.  You must know the name of the variable you are trying to modify, so some insight into the naming conventions of ACHP is required.