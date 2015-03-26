'''
This file holds some functions that don't have any obvious other home
'''
import os
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
        
    from Cycle import SecondaryCycleClass,DXCycleClass
    # Check if it is an instance of one of the cycle classes - more work required
    # to collect all the component outputs
    if isinstance(Class,(SecondaryCycleClass,DXCycleClass)):
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
        head=','.join(headList)
        units=','.join(unitsList)
        vals=','.join(valsList)
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
    print entries_path

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
    print "Warning! This ACHP version tool gives only main revision number - update working copy before usage!"
    if rev:
        return u'SVN-%s' % rev
    return u'SVN-unknown'

if __name__=='__main__':
    print get_svn_revision(None)



