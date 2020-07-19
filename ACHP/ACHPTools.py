from __future__ import division, print_function, absolute_import
'''
This file holds some functions that don't have any obvious other home
'''
import os
import csv
import sys
import numpy as np

from io import StringIO

# see http://stackoverflow.com/a/14707227/1360263
from contextlib import contextmanager
@contextmanager
def stdout_redirected(new_stdout):
    save_stdout = sys.stdout
    sys.stdout = new_stdout
    try:
        yield None
    finally:
        sys.stdout = save_stdout

def redirected_exec(pyfile, outfile):
    with open(outfile, "w") as f:
        with stdout_redirected(f):
            with open(pyfile, 'r') as fp:    
                exec(fp.read())

def Write2CSV(Class,file,append=False):
    """
    This function takes in a class and a file pointer
    """
    def OL2Strings(OL):
        head=str(OL[0][0])
        units=str(OL[0][1])
        vals=str(OL[0][2])
        for i in range(1,len(OL)):
            head+=','+str(OL[i][0])
            units+=','+str(OL[i][1])
            vals+=','+str(OL[i][2])
        return head,units,vals
    
    def BuildComponentList(ShapeString,string):
        OutString=string
        for i in range(len(ShapeString.split(','))-1):
            OutString+=','+string
        return OutString
        
    #from Cycle import SecondaryCycleClass,DXCycleClass
    # Check if it is an instance of one of the cycle classes - more work required
    # to collect all the component outputs
    if True: #if isinstance(Class,(SecondaryCycleClass,DXCycleClass)):
        #Pull the cycle outputs
        head,units,vals=OL2Strings(Class.OutputList())
        headList=[head]
        unitsList=[units]
        valsList=[vals]
        componentList=[BuildComponentList(units,'Cycle')]
        
        #Loop over the other things that are there
        for item in dir(Class):            
            #If the item has an outputList, collect it
            if hasattr(getattr(Class,item),'OutputList'):
                head,units,vals=OL2Strings(getattr(Class,item).OutputList())
                componentList+=[BuildComponentList(units,item)]
                headList+=[head]
                unitsList+=[units]
                valsList+=[vals]
        components=','.join(componentList)
        head=','.join(headList)     #print head
        units=','.join(unitsList)   #print units
        vals=','.join(valsList)     #print vals
        IsCycle=True
    else:
        head,units,vals=OL2Strings(Class.OutputList())
        IsCycle=False
        
    if type(file)!=type('some string'):
        #A file object was passed in, use it
        fP=file
        firstRow=False
    else:
        if os.path.exists(file):
            firstRow=False
        else:
            firstRow=True
        if append==True:
            fP=open(file,'a')
        else:
            fP=open(file,'w')
            
    if append==True and firstRow==False:
        fP.write(vals+'\n')
    else:           
        if IsCycle==True:
            fP.write(components+'\n')
        fP.write(head+'\n')
        fP.write(units+'\n')
        fP.write(vals+'\n')
    fP.close()
    
def simple_write_to_file_(head,data,file,append='True',newline=True):
    #function to simply write some data into a file
    #quick and dirty, nonspecific, therefore no formatting...
    #append is unused
    fP=open(file,'a')
    if newline:
        fP.write(str(head)+str(data)+'\n')
    else:
        fP.write(str(head)+str(data))
    fP.close

def simple_write_to_file(head,data,file,append='True',newline=True):
    #function to simply write some data comma seperated into a file
    #quick and dirty, nonspecific, therefore no formatting...
    #append is unused, dictionaries not supported, yet
    def writertool(info_bit):
        fP.write(","+str(info_bit))
    
    def callw(input):
        test=str(input)
        if test==input:
            writertool(input)
        elif isinstance(input, dict):
            for key in input:
                callw(key+":")
                callw(input[key])
        else:
            try:
                for num in range(len(input)):
                    callw(input[num])
                callw('<>')
            except:
                writertool(input)
           
    fP=open(file,'a')
    fP.write(str(head))
    callw(data)
    if newline:
        writertool('\n')
    else:
        writertool(',')
    fP.close
    
def Get_dict_write(filename_old,enforce_float=False,write_file=True):

    #change a column-wise file to a row wise file (easier to read from), then return dictionary of data
    #reads in standard python output file and saves it rotated as csv
    #can enforce only numerical data and switch of write to file functionality

    f = csv.reader(open(filename_old))
    fw = open(str(filename_old[:-4])+'_rows.csv','wb')
    w = csv.writer(fw)
    columns = zip(*f)
    dict={}
    for column in columns:
        if not dict.has_key(str(column[0])):
            dict[str(column[0])]={}
        if enforce_float:
            dict[str(column[0])].update({str(column[1]):np.array(column[3:],dtype='float')})
        else:
            try:
                dict[str(column[0])].update({str(column[1]):np.array(column[3:],dtype='float')}) #this is used for cells that contain numerical data
            except:
                dict[str(column[0])].update({str(column[1]):np.array(column[3:])}) #this uis used for cells that contain strings
        dict[str(column[0])].update({str(column[1])+' units':str(column[2])})
        if write_file:
            w.writerow(column)
    return dict

def print_dict(dict_inp):
    #print fict as given by function
    # Get_dict_write
    for key in dict_inp:
        print (key,":",end='')
        for keyword in dict_inp[key]:
            print (" > "+keyword+" < ",end='')
        print ("<<<")


def ValidateFields(d,reqFields,optFields=None):
        '''
        The function checkFields takes in inputs of:
        
        =========   =============================================================
        Variable    Type & Description
        =========   =============================================================
        d           dict of values that are part of structure
        reqFields   list of tuples in the form (fieldname,typepointer,min,max)
        optFields   list of other fieldnames
        =========   =============================================================
        
        required parameters are checked that they 
        * exist
        * can be cast using the typepointer function pointer
        * is within the range (min,max)
        
        if a parameter is on the optional parameters list, it is ok-ed, but not value checked
        
        Additional parameters raise AttributeError
        
        '''
        #make a copy of d
        d=dict(d)
        
        #Required parameters
        for field,typepointer,min,max in reqFields:
            if field in d:
                #See if you can do a type cast using the conversion function pointer
                # You should get the same value back
                assert typepointer(d[field])==d[field],field+': failed type conversion, should be '+str(typepointer)
                #check the bounds if numeric input
                if typepointer in (float,int):
                    assert d[field]>=min and d[field]<=max,field+' (value: %g) not in the range [%g,%g]'%(d[field],min,max)
                #remove field from dictionary of terms left to check if no errors
                del d[field]
            else:
                raise AttributeError('Required field '+field+' not included')
        #Optional parameters (not strictly checked, just checked their existence)
        if optFields!=None:
            for field in optFields:
                if field in d:
                    del d[field]
        assert len(d)==0,'Unmatched fields found: '+str(d.keys())
        
def get_svn_revision(path=None):
    import re
    rev = None
    if path is None:
        path = "C:/Users/BACHC/Desktop/achp/trunk/PyACHP"
    entries_path = '%s/.svn/entries' % path
    print (entries_path)

    if os.path.exists(entries_path):
        entries = open(entries_path, 'r').read()
        # Versions >= 7 of the entries file are flat text.  The first line is
        # the version number. The next set of digits after 'dir' is the revision.
        if re.match('(\d+)', entries):
            rev_match = re.search('\d+\s+dir\s+(\d+)', entries)
            if rev_match:
                rev = rev_match.groups()[0]
        # Older XML versions of the file specify revision as an attribute of
        # the first entries node.
        else:
            from xml.dom import minidom
            dom = minidom.parse(entries_path)
            rev = dom.getElementsByTagName('entry')[0].getAttribute('revision')
    print ("Warning! This ACHP version tool gives only main revision number - update working copy before usage!")
    if rev:
        return u'SVN-%s' % rev
    return u'SVN-unknown'


def subsample(data, sample_size):
    #subsample data
    sample_size=int(sample_size)
    samples = list(zip(*[iter(data)]*sample_size))   # use 3 for triplets, etc.
    return map(lambda x:sum(x)/float(len(x)), samples)


def smooth_curve(curve_data,N_smooth,exp_max=-1,shift_0=0,fix_first_nonzero=False,plotit=False,t='x for plot'):
    """
    smoothens the curve data for plotting as good as possible while maintaining last and first value
    curve data => np.array, 1D that should be smoothened out
    N_smooth => number of points to smooth over (float)
    exp_max => adjust exponential behavior for average (0='normal' moving average)
    shift_0 => manually fix from where on smoothing is active, e.g. up till where no smootthing is applied
    fix_first_nonzero => if set to true, then automatically determines shift_0 to be where the first nonzero entry is
    plotit => plot results
    t => x-cooordinate for plot
    """
    a=curve_data
    N=N_smooth
    v=np.exp(np.linspace(exp_max, 0., N))
    v=v/v.sum()
    a_v=np.convolve(a,v,'same')
    if fix_first_nonzero==True:
        shift_0=np.nonzero(a != 0)[0][0]
    for n in range(0,len(v)):
        if n!=0:
            v=np.exp(np.linspace(exp_max, 0., n))
            v=v/v.sum()
            a_v[n+shift_0]=np.convolve(a,v,'same')[n+shift_0]
            a_v[len(a)-n-1]=np.convolve(a,v,'same')[len(a)-n-1]
        else:
            a_v[n+shift_0]=a[n+shift_0]
            for i in range(0,n+shift_0):
                a_v[i]=a[i]
            a_v[len(a)-n-1]=a[len(a)-n-1]
    if plotit:
        try:
            np.sin(t)
        except:
            t=np.linspace(0,len(curve_data),len(curve_data))
        import pylab as plt
        plt.plot(t,a,label='original data')
        plt.plot(t,a_v,label='smoothened')
        plt.legend(loc='best',fancybox=True)
        plt.title('curve smoothing')
        plt.show()
    return a_v
    
if __name__=='__main__':
    print (get_svn_revision(None))  
    print (' ')
    
if __name__=='__main__':
    #show how to use curve smoothing
    temps=np.array([-0.5,-0.034782609,0.430434783,0.895652174,1.360869565,1.826086957,2.291304348,2.756521739,3.22173913,3.686956522,4.152173913,4.617391304,5.082608696,5.547826087,6.013043478,6.47826087,6.943478261,7.408695652,7.873913043,8.339130435,8.804347826,9.269565217,9.734782609,10.2,10.66521739,11.13043478,11.59565217,12.06086957,12.52608696,12.99130435,13.45652174,13.92173913,14.38695652,14.85217391,15.3173913,15.7826087,16.24782609,16.71304348,17.17826087,17.64347826,18.10869565,18.57391304,19.03913043,19.50434783,19.96956522,20.43478261,20.9,21.36521739,21.83043478,22.29565217,22.76086957,23.22608696,23.69130435,24.15652174,24.62173913,25.08695652,25.55217391,26.0173913,26.4826087,26.94782609,27.41304348,27.87826087,28.34347826,28.80869565,29.27391304,29.73913043,30.20434783,30.66956522,31.13478261,31.6,32.06521739,32.53043478,32.99565217,33.46086957,33.92608696,34.39130435,34.85652174,35.32173913,35.78695652,36.25217391,36.7173913,37.1826087,37.64782609,38.11304348,38.57826087,39.04347826,39.50869565,39.97391304,40.43913043,40.90434783,41.36956522,41.83478261,42.3])
    data=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.007322468,0.107196952,0.269918828,0.371320001,0.471285827,0.145456344,1.40207551,1.748807824,1.794470974,1.339313933,1.513166286,1.501733385,0.438144108,1.557055465,2.02774747,2.173010628,1.901990502,2.511263394,0.213424172,2.71911854,2.685546322,3.19412648,3.450965543,3.721265054,0.28610417,3.417866297,5.106969924,4.868446692,3.071483957,3.766539738,0.103899522,3.78272759,1.896961631,3.069846495,3.388021634,4.760657037,0,3.73626734,2.295885187,0.762500277,0.281666101,0.312510595])
    smooth_curve(data,10,plotit='True',fix_first_nonzero=True,exp_max=-.2,t=temps)
    
    #show how to save stuff with simple write to file
    data={'some numbers':np.array([1,2,3,4,5,6])}
    head='if someone can save'
    file='blub.csv'
    simple_write_to_file(head,data,file)
    head='more  numbers not only'
    data=['ba',np.array([1,2,3]),[1,2,3]]
    simple_write_to_file(head,data,file)
    
    #subsample to get smaller array
    l =np.array([1, 2, 3, 4,5,6])
    print (subsample(l, 2))
    print (subsample(l, 3))
    print (subsample(l, 4))