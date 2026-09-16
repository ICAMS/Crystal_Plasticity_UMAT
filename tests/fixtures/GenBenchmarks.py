import os
import numpy as np
import numpy.linalg as npl
import argparse

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
from odbAccess import *  # Access to Symbolic Constants defined in Abaqus

#
#==============================================================================
#
executeOnCaeStartup()
Mdb()

sourcePath = os.getcwd() + r'/source'
UmatPath = sourcePath + r'/umat.f'
PathSave = os.getcwd()

# Test cases
Cases = []

Cases.append('Ux')
Cases.append('Uy')
Cases.append('Uz') 
Cases.append('ShearXY') 
Cases.append('ShearXZ')  
Cases.append('ShearYZ')

## Superalloy: Please activate the option for SuperAlloy in mod_wkcoup.f (i.e. Iwkcoup_int = 1, Iwkcoup_sup = 1)
#Cases.append('SuperAlloy')

## Kinematic hardening: Please active the option for kinematic hardening in mod_wkcoup.f (i.e. Iwkcoup_bk = 1)
Cases.append('Kinematic') 



class OdbData(object):
#""" 
# Python script to analyse Abaqus OBD files

# Homogenzied values are obtained via averaging of element values
# Different element sizes are not considered! 
# Use only for voxelized structures with structured mesh

#Authors:
# Alexander Hartmaier
# ICAMS, Ruhr-Universitat Bochum

#Version: 1.0.0 (2023-03-09)
    
    def __init__(self):
        self.name=[]

        #Times
        self.steptime=[]
        self.totaltime=[]

        #homogenized stress
        self.sigmaH=[]
        self.sigmaHE=[]

        #homogenized strain
        self.etotH=[]
        self.etotHE=[]

        # Temperature
        self.temperature = []

        #tensors
        self.etot1=[]
        self.etot2=[]
        self.etot3=[]
        self.etot4=[]
        self.etot5=[]
        self.etot6=[]

        self.epl1=[]
        self.epl2=[]
        self.epl3=[]
        self.epl4=[]
        self.epl5=[]
        self.epl6=[]
        self.peeq=[]

        self.sig1=[]
        self.sig2=[]
        self.sig3=[]
        self.sig4=[]
        self.sig5=[] 
        self.sig6=[]


    def F_ODB(self, odbName, firstStep, lastStep):
        # Open the ODB
        odb=openOdb(odbName,readOnly=True)      
        
        # Show all setp stored in odb
        allSteps=odb.steps.keys()
        nodeObjectDict=dict()

        # Dictionary definition
        #nodeObjectDict['upper']=Class_NodeObject()
            
        #Access to nodes
        #nodeObjectDict['upper'].abqRegion = odb.rootAssembly.nodeSets['Set-4']
        
        # Loop over all steps
        for stepKey in range(firstStep, lastStep+1):
            # Get current step name
            currentStepName=odb.steps.keys()[stepKey-1]
            #print 'Step: ', currentStepName
            
            # Loop over all frames in current step
            for currentFrame in odb.steps[currentStepName].frames:
                # Time in step for the current frame
                stepTime=currentFrame.frameValue
                self.steptime.append(stepTime)
                print('Analysing step: ',stepTime)
                # Read plastic strain tensor
                epl1, epl2, epl3, epl4, epl5, epl6, eeq = self.average_elmts(
                    currentFrame.fieldOutputs['SDV156'].values,
                    currentFrame.fieldOutputs['SDV157'].values,
                    currentFrame.fieldOutputs['SDV158'].values,
                    currentFrame.fieldOutputs['SDV159'].values,
                    currentFrame.fieldOutputs['SDV160'].values,
                    currentFrame.fieldOutputs['SDV161'].values)
                
                self.epl1.append(epl1)
                self.epl2.append(epl2)
                self.epl3.append(epl3)
                self.epl4.append(epl4)
                self.epl5.append(epl5)
                self.epl6.append(epl6)
                self.peeq.append(eeq)

                # read total strain tensor
#                tens = currentFrame.fieldOutputs['LE']
#                et1, et2, et3, et4, et5, et6, eeq = self.average_elmts(
#                    tens.getScalarField(componentLabel='LE11').values,
#                    tens.getScalarField(componentLabel='LE22').values,
#                    tens.getScalarField(componentLabel='LE33').values,
#                    tens.getScalarField(componentLabel='LE12').values,
#                    tens.getScalarField(componentLabel='LE13').values,
#                    tens.getScalarField(componentLabel='LE23').values,
#                    val_eq = tens.getScalarField(invariant=MAX_PRINCIPAL).values)
#                self.etot1.append(et1)
#                self.etot2.append(et2)
#                self.etot3.append(et3)
#                self.etot4.append(et4)
#                self.etot5.append(et5)
#                self.etot6.append(et6)
#                self.etotHE.append(eeq)

                # read stress tensor
                tens = currentFrame.fieldOutputs['S']
                sig1, sig2, sig3, sig4, sig5, sig6, seq = self.average_elmts(
                    tens.getScalarField(componentLabel='S11').values,
                    tens.getScalarField(componentLabel='S22').values,
                    tens.getScalarField(componentLabel='S33').values,
                    tens.getScalarField(componentLabel='S12').values,
                    tens.getScalarField(componentLabel='S13').values,
                    tens.getScalarField(componentLabel='S23').values,
                    val_eq = tens.getScalarField(invariant=MISES).values)
                self.sig1.append(sig1)
                self.sig2.append(sig2)
                self.sig3.append(sig3)
                self.sig4.append(sig4)
                self.sig5.append(sig5)
                self.sig6.append(sig6)
                self.sigmaHE.append(seq)

                # Read temperature if it exists in field output
                try:
                    temp = currentFrame.fieldOutputs['TEMP'].values
                    self.temperature.append(temp[0].data)
                except:
                    pass
        odb.close()
        return nodeObjectDict

    def average_elmts(self, val1, val2, val3, val4, val5, val6, val_eq=None):
        nd = len(val1)
        eg1 = 0.; eg2 = 0.; eg3 = 0.
        eg4 = 0.; eg5 = 0.; eg6 = 0.
        egv = 0.
        for i in range(nd):
            eg1 += val1[i].data
            eg2 += val2[i].data
            eg3 += val3[i].data
            eg4 += val4[i].data
            eg5 += val5[i].data
            eg6 += val6[i].data
            if val_eq is not None:
                # average equivalent stress or strain given in input
                egv += val_eq[i].data
            else:
                # calculate equivalent strain
                egv += np.sqrt(2.*(val1[i].data**2 + val2[i].data**2 + val3[i].data**2 
                + 0.5*(val4[i].data**2 + val5[i].data**2 + val6[i].data**2))/3.)
        eg1 /= nd
        eg2 /= nd
        eg3 /= nd
        eg4 /= nd
        eg5 /= nd
        eg6 /= nd
        egv /= nd 
        return eg1, eg2, eg3, eg4, eg5, eg6, egv

def writeFiles(Result, odb_trunk):
    if odb_trunk == 'Ux' or odb_trunk == 'Kinematic':
        header_str = '# time (s), plastic strain [11, 22, 33] (.), stress [11] (MPa) \n'
        out_val = zip(Result.steptime, Result.epl1, Result.epl2, Result.epl3, Result.sig1)
            
    elif odb_trunk == 'Uy' or odb_trunk == 'SuperAlloy':
        header_str = '# time (s), plastic strain [11, 22, 33] (.), stress [22] (MPa)\n'
        out_val = zip(Result.steptime, Result.epl1, Result.epl2, Result.epl3, Result.sig2)

    elif odb_trunk == 'Uz':
        header_str = '# time (s), plastic strain [11, 22, 33] (.), stress [33] (MPa)\n'
        out_val = zip(Result.steptime, Result.epl1, Result.epl2, Result.epl3, Result.sig3)
    
    elif odb_trunk == 'ShearXY':
        header_str = '# time (s), plastic strain [12] (.), stress [12] (MPa)\n'
        out_val = zip(Result.steptime, Result.epl4, Result.sig4)

    elif odb_trunk == 'ShearXZ':
        header_str = '# time (s), plastic strain [13] (.), stress [13] (MPa)\n'
        out_val = zip(Result.steptime, Result.epl5, Result.sig5)

    elif odb_trunk == 'ShearYZ':
        header_str = '# time (s), plastic strain [23] (.), stress [23] (MPa)\n'
        out_val = zip(Result.steptime, Result.epl6, Result.sig6)

    with open(odb_trunk + '_Ref.csv','w') as f:
          f.write(header_str)
          np.savetxt(f, out_val, delimiter=',', fmt='%f')
            
def deleteFiles(odb_trunk):
    extensions = ['.odb', '.msg', '.log', '.dat', '.com', '.sta', '.stt','.simlog', '.sim', '.prt', '.2.SMABulk', '.1.SMABulk', '.ipm']
    for extension in extensions:
        name = odb_trunk + extension
        if os.path.exists(name):
            os.remove(name)
#==============================================================================
#==============================================================================
print( 'Starting Script' )


for case in Cases:
    JobName = case
    path = PathSave + r'/' + JobName + '.inp'
    odb_trunk = case
    print(path)
    if not os.path.isfile(path):
        print(path)
        raise FileNotFoundError('File does not exist.')
    mdb.JobFromInputFile(name= JobName, inputFileName= path, userSubroutine = UmatPath, scratch=sourcePath) # userSubroutine
    mdb.jobs[JobName].submit()
    mdb.jobs[JobName].waitForCompletion()

    # Write output
    #Number of steps
    beginStep=1
    endStep=1
    # Definition of Result Files
    odb_name = odb_trunk+'.odb'
    Result=OdbData()
    Result.F_ODB(odb_name, beginStep, endStep)
    
    writeFiles(Result, odb_trunk)

    # Delete abaqus files
    deleteFiles(odb_trunk)

print('Finish \n')
