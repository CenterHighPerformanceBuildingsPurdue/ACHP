import sys,os,subprocess,shutil
    
def BuildLaTeX():
    '''
    Build the HTML documentation
    '''
    options=''
    subprocess.call('sphinx-build %s -P -b latex -d _build/doctrees . _build/latex' % options,shell=True,stderr=subprocess.STDOUT)
    subprocess.call('make all-pdf',cwd='_build/latex',shell=True,stderr=subprocess.STDOUT) 
    
def BuildHTML():
    '''
    Build the HTML documentation
    '''
    options=''
    subprocess.call('sphinx-build %s -P -b html -d _build/doctrees . _build/html' % options,shell=True,stderr=subprocess.STDOUT)
    
if __name__=="__main__":
    if len(sys.argv)==1:
##         BuildLaTeX() #Build LaTeX first because you need the PDF to add to the _static folder
##         #Make the _static folder if it doesn't already exist
##         try:
##             os.mkdir('_static')
##         except OSError:
##             pass
##         shutil.copy2('_build/latex/ACHP.pdf','_static/ACHP.pdf')
        BuildHTML()
    