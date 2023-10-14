#!/usr/env/python
# -*- coding: utf-8 -*-
""" 
 Python script to analyse Abaqus OBD files

 Homogenzied values are obtained via averaging of element values
 Different element sizes are not considered! 
 Use only for voxelized structures with structured mesh

Authors:
 Alexander Hartmaier
 ICAMS, Ruhr-Universität Bochum

Version: 1.0.0 (2023-03-09)

"""
import os
import numpy as np
import numpy.linalg as npl
import argparse

from odbAccess import *  # Access to Symbolic Constants defined in Abaqus
from abaqusConstants import *
#
#==============================================================================
#

# parse command-line argument to get filename
parser = argparse.ArgumentParser(description='Name trunk of Abaqus ODB file')
parser.add_argument('string', metavar='odb_trunk', type=str,
                   help='A string with the ODB filename w/o .odb')
odb_trunk = parser.parse_args().string
odb_name = odb_trunk+'.odb'
if not os.path.isfile(odb_name):
    print(odb_name)
    raise FileNotFoundError('File does not exist.')


class OdbData(object):
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
                tens = currentFrame.fieldOutputs['LE']
                et1, et2, et3, et4, et5, et6, eeq = self.average_elmts(
                    tens.getScalarField(componentLabel='LE11').values,
                    tens.getScalarField(componentLabel='LE22').values,
                    tens.getScalarField(componentLabel='LE33').values,
                    tens.getScalarField(componentLabel='LE12').values,
                    tens.getScalarField(componentLabel='LE13').values,
                    tens.getScalarField(componentLabel='LE23').values,
                    val_eq = tens.getScalarField(invariant=MAX_PRINCIPAL).values)
                self.etot1.append(et1)
                self.etot2.append(et2)
                self.etot3.append(et3)
                self.etot4.append(et4)
                self.etot5.append(et5)
                self.etot6.append(et6)
                self.etotHE.append(eeq)

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



#==============================================================================
#==============================================================================
print( 'Starting Script' )

#Number of steps
beginStep=1
endStep=1
# Definition of Result Files

#Result1
Result1=OdbData()
Result1.F_ODB(odb_name, beginStep, endStep)

header_str = '# time (s), total strain [11, 22, 33, 12, 13, 23] (.), plastic strain [11, 22, 33, 12, 13, 23] (.), stress [11, 22, 33, 12, 13, 23] (MPa), max. principal total strain (.), equiv. stress (MPa)'
lt = len(Result1.temperature)
if lt > 0:
    header_str += ', Temperature (K)\n'
    out_val = zip(Result1.steptime,
              Result1.etot1, Result1.etot2, Result1.etot3,
              Result1.etot4, Result1.etot5, Result1.etot6,
              Result1.epl1, Result1.epl2, Result1.epl3,
              Result1.epl4, Result1.epl5, Result1.epl6,
              Result1.sig1, Result1.sig2, Result1.sig3,
              Result1.sig4, Result1.sig5, Result1.sig6,
              Result1.etotHE, Result1.sigmaHE, Result1.temperature)
else:
    header_str += ', PEEQ\n'
    out_val = zip(Result1.steptime,
              Result1.etot1, Result1.etot2, Result1.etot3,
              Result1.etot4, Result1.etot5, Result1.etot6,
              Result1.epl1, Result1.epl2, Result1.epl3,
              Result1.epl4, Result1.epl5, Result1.epl6,
              Result1.sig1, Result1.sig2, Result1.sig3,
              Result1.sig4, Result1.sig5, Result1.sig6,
              Result1.etotHE, Result1.sigmaHE, Result1.peeq)

with open(odb_trunk+'_sig_eps.csv','w') as f:
    f.write(header_str)
    np.savetxt(f, out_val, delimiter=',', fmt='%f')

print(Result1.sigmaH)

