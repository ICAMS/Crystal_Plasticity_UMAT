#!/usr/env/python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
#
# Martin Boeff
#
#
# -------------------------------------------------------------------

import numpy as np
import os
import time
import argparse

# parse command-line argument to get filename
parser = argparse.ArgumentParser(description='Name Abaqus geometry input file')
parser.add_argument('string', metavar='geom_inp', type=str,
                   help='A string with the input filename')
geom_name = parser.parse_args().string
if not os.path.isfile(geom_name):
    print(geom_name)
    raise FileNotFoundError('File does not exist.')


class InputFile(object):
       
       def __init__(self):
              
              self.nodeNumber=[]
              self.xCoord=[]
              self.yCoord=[]
              self.zCoord=[]


       def F_AssignNodesToCoordinates(self,filename):
              daten=open(filename,'r')
              jumpParameter=0
              content=[]
              CutOff=6
              for zeile in daten:
                     zeile=zeile.strip('\r\n').strip(' ').split('\n')
                     #print zeile
                     content.extend(zeile)
              #print content
              nodeStart=content.index('*Node')+1
        
              #Seachres line where there are no integers are readen any more to get breakout crit
              for item in content[nodeStart:-1]:
                     try:
                            int(item[0])
                     except ValueError:
                            #nodeStop=content.index(str(item))
                            nodeStop = content[nodeStart:-1].index(str(item)) + nodeStart
                            break

              for item in content[nodeStart:nodeStop]:
                     item=item.strip('\r\n').split(',')
                     self.nodeNumber.append(int(item[0]))
                     self.xCoord.append(round(float(item[1]),CutOff))
                     self.yCoord.append(round(float(item[2]),CutOff))
                     self.zCoord.append(round(float(item[3]),CutOff))


class Set(object):
       
       def __init__(self):
              
              self.nodeNumber=[]
              self.xCoord=[]
              self.yCoord=[]
              self.zCoord=[]


       def F_SurfaceSet(self,nodeNumberList,xCoordList,yCoordList,zCoordList):

              maxX = max(xCoordList)              
              minX = min(xCoordList)
              maxY = max(yCoordList)              
              minY = min(yCoordList)              
              maxZ = max(zCoordList)
              minZ = min(zCoordList)                  

              for i in range(len(nodeNumberList)):
                     
                     surface = False       

                     if xCoordList[i] == maxX:
                            surface = True
                     elif xCoordList[i] == minX:
                            surface = True
                     elif yCoordList[i] == maxY:
                            surface = True
                     elif yCoordList[i] == minY:
                            surface = True
                     elif zCoordList[i] == minZ:
                            surface = True
                     elif zCoordList[i] == maxZ:
                            surface = True

                     if surface == True:

                            self.nodeNumber.append(nodeNumberList[i])
                            self.xCoord.append(xCoordList[i])
                            self.yCoord.append(yCoordList[i])
                            self.zCoord.append(zCoordList[i])


# NodeSetClass generates where NodeSet Numbers are stored in
class NodeClass(object):

       def __init__(self):
              self.nodeSet=[]                      


# Creates a List of Sorted Nodes with respect to sortCoord1 and sortCoord2 for faces
# Input: Nodeset of surface nodes
def CreatePeriodicNodeSets(Nodes,sortCoord1,sortCoord2,NodeSet):

       startList=[]
       sortedList=[]

       #print NodeSet
       for number in NodeSet:
              startList.append([number,sortCoord1[Nodes.index(number)],sortCoord2[Nodes.index(number)]])
       #print startList
       import operator

       startList.sort(key=operator.itemgetter(1, 2))
       for item in range(len(startList)):
              sortedList.append(startList[item][0])

       return sortedList


# Creates a List of Sorted Nodes with respect to sortCoord for edges
def CreatePeriodicEdgeSets(Nodes,sortCoord,NodeSet):

       startList=[]
       sortedList=[]

       #print NodeSet
       for number in NodeSet:
              startList.append([number,sortCoord[Nodes.index(number)]])

       import operator

       startList.sort(key=operator.itemgetter(1))
       for item in range(len(startList)):
              sortedList.append(startList[item][0])

       return sortedList


def RemoveItemInList(NodeList,ReplaceNodeList):
 
       for item in ReplaceNodeList:
              try:            
                     NodeList.remove(item)
              except ValueError:
                     pass

       return NodeList


def F_CreateSurfaceSet(nodeNumber,xCoord,yCoord,zCoord,SurfaceName):

       if SurfaceName=='LEFT':
              coordCrit=min(xCoord)
              coord=xCoord
       if SurfaceName=='RIGHT':
              coordCrit=max(xCoord)
              coord=xCoord
       if SurfaceName=='BOTTOM':
              coordCrit=min(yCoord)
              coord=yCoord
       if SurfaceName=='TOP':
              coordCrit=max(yCoord)
              coord=yCoord
       if SurfaceName=='FRONT':
              coordCrit=max(zCoord)
              coord=zCoord
       if SurfaceName=='REAR':
              coordCrit=min(zCoord)
              coord=zCoord

       #Nodes in Set
       surfaceNodes=[]
       
       for i in range(len(nodeNumber)):
              if coord[i]==coordCrit:
                     surfaceNodes.append(nodeNumber[i])

       return surfaceNodes

# -------------------------------------------------------------------
ta=time.perf_counter()
inputFile=geom_name
outputDirectory='PeriodicData/'
if not os.path.exists(str(outputDirectory)) == True:
       os.system('mkdir -p '+str(outputDirectory))

Input=InputFile()
Input.F_AssignNodesToCoordinates(inputFile)
print( 'Total Node Number: ', len(Input.nodeNumber))


# Create node sets for boundary conditions
leftSet=F_CreateSurfaceSet(Input.nodeNumber,Input.xCoord,Input.yCoord,Input.zCoord,'LEFT')
rightSet=F_CreateSurfaceSet(Input.nodeNumber,Input.xCoord,Input.yCoord,Input.zCoord,'RIGHT')
bottomSet=F_CreateSurfaceSet(Input.nodeNumber,Input.xCoord,Input.yCoord,Input.zCoord,'BOTTOM')
topSet=F_CreateSurfaceSet(Input.nodeNumber,Input.xCoord,Input.yCoord,Input.zCoord,'TOP')
frontSet=F_CreateSurfaceSet(Input.nodeNumber,Input.xCoord,Input.yCoord,Input.zCoord,'FRONT')
rearSet=F_CreateSurfaceSet(Input.nodeNumber,Input.xCoord,Input.yCoord,Input.zCoord,'REAR')

Surface=Set()
Surface.F_SurfaceSet(Input.nodeNumber,Input.xCoord,Input.yCoord,Input.zCoord)
print ('Surface nodes: ', len(Surface.nodeNumber))

######### Reads the nodesets from the input file
# NODE SETS
print ('Nodes before deleting double nodes')
Left=NodeClass()
Left.nodeSet=leftSet
print ('L: ', len(Left.nodeSet))

Right=NodeClass()
Right.nodeSet=rightSet
print ('R: ', len(Right.nodeSet))

Bottom=NodeClass()
Bottom.nodeSet=bottomSet
print ('B: ', len(Bottom.nodeSet))

Top=NodeClass()
Top.nodeSet=topSet
print ('T: ', len(Top.nodeSet))

Front=NodeClass()
Front.nodeSet=frontSet
print ('F: ', len(Front.nodeSet))

Rear=NodeClass()
Rear.nodeSet=rearSet
print ('R: ', len(Rear.nodeSet))

t1=time.perf_counter()
##################
# EDGES
##################
# top front edge
E_T1 = np.intersect1d(Front.nodeSet,Top.nodeSet)
E_T1 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.xCoord,E_T1)

# top back edge
E_T3 = np.intersect1d(Rear.nodeSet,Top.nodeSet)
E_T3 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.xCoord,E_T3)

# top left edge
E_T4 = np.intersect1d(Left.nodeSet,Top.nodeSet)
E_T4 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.zCoord,E_T4)

# top right edge
E_T2 = np.intersect1d(Right.nodeSet,Top.nodeSet)
E_T2 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.zCoord,E_T2)

##
# bottm front edge
E_B1 = np.intersect1d(Front.nodeSet,Bottom.nodeSet)
E_B1 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.xCoord,E_B1)

# bottm back edge
E_B3 = np.intersect1d(Rear.nodeSet,Bottom.nodeSet)
E_B3 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.xCoord,E_B3)

# bottm left edge
E_B4 = np.intersect1d(Left.nodeSet,Bottom.nodeSet)
E_B4 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.zCoord,E_B4)

# bottm right edge
E_B2 = np.intersect1d(Right.nodeSet,Bottom.nodeSet)
E_B2 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.zCoord,E_B2)

##
# left front edge
E_M1 = np.intersect1d(Left.nodeSet,Front.nodeSet)
E_M1 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.yCoord,E_M1)

# right front edge
E_M2 = np.intersect1d(Front.nodeSet,Right.nodeSet)
E_M2 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.yCoord,E_M2)

# left rear edge
E_M4 = np.intersect1d(Rear.nodeSet,Left.nodeSet)
E_M4 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.yCoord,E_M4)

# right rear edge
E_M3 = np.intersect1d(Rear.nodeSet,Right.nodeSet)
E_M3 = CreatePeriodicEdgeSets(Surface.nodeNumber,Surface.yCoord,E_M3)

t2=time.perf_counter()
##########
# VERTICES
##########
V1 = np.intersect1d(E_B1,E_B4)
V1 = V1[0]
V2 = np.intersect1d(E_B1,E_B2)
V2 = V2[0]
H1 = np.intersect1d(E_B3,E_B4)
H1 = H1[0]
H2 = np.intersect1d(E_B2,E_B3)
H2 = H2[0]
V3 = np.intersect1d(E_T1,E_T2)
V3 = V3[0]
H3 = np.intersect1d(E_T2,E_T3)
H3 = H3[0]
V4 = np.intersect1d(E_T1,E_T4)
V4 = V4[0]
H4 = np.intersect1d(E_T3,E_T4)
H4 = H4[0]

CORNERNODES=[V1,V2,V3,V4,H1,H2,H3,H4]

##########
# CreateEdgeNodeset
##########
EdgeNodes=[]
EdgeNodes.extend(E_T1)
EdgeNodes.extend(E_T3)
EdgeNodes.extend(E_T2)
EdgeNodes.extend(E_T4)
EdgeNodes.extend(E_B1)
EdgeNodes.extend(E_B3)
EdgeNodes.extend(E_B2)
EdgeNodes.extend(E_B4)
EdgeNodes.extend(E_M1)
EdgeNodes.extend(E_M2)
EdgeNodes.extend(E_M4)
EdgeNodes.extend(E_M3)

# Remove Corner Nodes from Edge Nodesets

E_T1 = RemoveItemInList(E_T1,CORNERNODES)
E_T3 = RemoveItemInList(E_T3,CORNERNODES)
E_T4 = RemoveItemInList(E_T4,CORNERNODES)
E_T2 = RemoveItemInList(E_T2,CORNERNODES)
E_B1 = RemoveItemInList(E_B1,CORNERNODES)
E_B3 = RemoveItemInList(E_B3,CORNERNODES)
E_B4 = RemoveItemInList(E_B4,CORNERNODES)
E_B2 = RemoveItemInList(E_B2,CORNERNODES)
E_M1 = RemoveItemInList(E_M1,CORNERNODES)
E_M2 = RemoveItemInList(E_M2,CORNERNODES)
E_M4 = RemoveItemInList(E_M4,CORNERNODES)
E_M3 = RemoveItemInList(E_M3,CORNERNODES)

#print len(E_M1),len(E_M2),len(E_M3),len(E_M4)
#print len(E_T1),len(E_T3),len(E_B1),len(E_B3)

######### Sorts the nodesets with respect to their coordinates
print ('Pure surface nodes after deleting edge and corner nodes')
t3=time.perf_counter()
BottomSet = Bottom.nodeSet
BottomSet = RemoveItemInList(BottomSet,CORNERNODES)
BottomSet = RemoveItemInList(BottomSet,EdgeNodes)
print  ('B: ', len(BottomSet))

TopSet = Top.nodeSet
TopSet = RemoveItemInList(TopSet,CORNERNODES)
TopSet = RemoveItemInList(TopSet,EdgeNodes)
print( 'T: ', len(TopSet))

LeftSet = Left.nodeSet
LeftSet = RemoveItemInList(LeftSet,CORNERNODES)
LeftSet = RemoveItemInList(LeftSet,TopSet)
LeftSet = RemoveItemInList(LeftSet,BottomSet)
LeftSet = RemoveItemInList(LeftSet,EdgeNodes)
print ('L: ', len(LeftSet))

RightSet = Right.nodeSet
RightSet = RemoveItemInList(RightSet,CORNERNODES)
RightSet = RemoveItemInList(RightSet,TopSet)
RightSet = RemoveItemInList(RightSet,BottomSet)
RightSet = RemoveItemInList(RightSet,EdgeNodes)
print ('R: ', len(RightSet))

FrontSet = Front.nodeSet
FrontSet = RemoveItemInList(FrontSet,CORNERNODES)
FrontSet = RemoveItemInList(FrontSet,TopSet)
FrontSet = RemoveItemInList(FrontSet,BottomSet)
FrontSet = RemoveItemInList(FrontSet,LeftSet)
FrontSet = RemoveItemInList(FrontSet,RightSet)
FrontSet = RemoveItemInList(FrontSet,EdgeNodes)
print ('F: ', len(FrontSet))

RearSet = Rear.nodeSet
RearSet = RemoveItemInList(RearSet,CORNERNODES)
RearSet = RemoveItemInList(RearSet,TopSet)
RearSet = RemoveItemInList(RearSet,BottomSet)
RearSet = RemoveItemInList(RearSet,LeftSet)
RearSet = RemoveItemInList(RearSet,RightSet)
RearSet = RemoveItemInList(RearSet,EdgeNodes)
print ('R: ', len(RearSet))

t4=time.perf_counter()
# Order opposite surfaces in the same way so that corresponding nodes are directly at same position in nodeSet
BottomSet = CreatePeriodicNodeSets(Surface.nodeNumber,Surface.xCoord,Surface.zCoord,BottomSet)

TopSet = CreatePeriodicNodeSets(Surface.nodeNumber,Surface.xCoord,Surface.zCoord,TopSet)

LeftSet = CreatePeriodicNodeSets(Surface.nodeNumber,Surface.yCoord,Surface.zCoord,LeftSet)

RightSet = CreatePeriodicNodeSets(Surface.nodeNumber,Surface.yCoord,Surface.zCoord,RightSet)

FrontSet = CreatePeriodicNodeSets(Surface.nodeNumber,Surface.xCoord,Surface.yCoord,FrontSet)

RearSet = CreatePeriodicNodeSets(Surface.nodeNumber,Surface.xCoord,Surface.yCoord,RearSet)




OutPutFile=open(str(outputDirectory)+'LeftToRight.inp','w')

OutPutFile.write('**** X-DIR \n')
for i in range(len(LeftSet)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(RightSet[i])+',1, 1 \n')
       OutPutFile.write(str(LeftSet[i])+',1,-1 \n')
       OutPutFile.write(str(V2)+',1,-1 \n')
       OutPutFile.write(str(V1)+',1, 1 \n')

OutPutFile.write('**** \n')
OutPutFile.write('**** Y-DIR \n')
for i in range(len(LeftSet)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(RightSet[i])+',2, 1 \n')
       OutPutFile.write(str(LeftSet[i])+',2,-1 \n')
       OutPutFile.write(str(V2)+',2,-1 \n')
       OutPutFile.write(str(V1)+',2, 1 \n')

OutPutFile.write('**** \n')
OutPutFile.write('**** Z-DIR \n')
for i in range(len(LeftSet)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(RightSet[i])+',3, 1 \n')
       OutPutFile.write(str(LeftSet[i])+',3,-1 \n')
       OutPutFile.write(str(V2)+',3,-1 \n')
       OutPutFile.write(str(V1)+',3, 1 \n')



OutPutFile=open(str(outputDirectory)+'BottomToTop.inp','w')

OutPutFile.write('**** X-DIR \n')
for i in range(len(BottomSet)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(BottomSet[i])+',1,1 \n')
       OutPutFile.write(str(TopSet[i])+',1,-1 \n')
       OutPutFile.write(str(V1)+',1,-1 \n')
       OutPutFile.write(str(V4)+',1,1 \n')


OutPutFile.write('**** \n')
OutPutFile.write('**** Y-DIR \n')
for i in range(len(BottomSet)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(BottomSet[i])+',2, 1 \n')
       OutPutFile.write(str(TopSet[i])+',2,-1 \n')
       OutPutFile.write(str(V1)+',2,-1 \n')
       OutPutFile.write(str(V4)+',2, 1 \n')


OutPutFile.write('**** \n')
OutPutFile.write('**** Z-DIR \n')
for i in range(len(BottomSet)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(BottomSet[i])+',3, 1 \n')
       OutPutFile.write(str(TopSet[i])+',3,-1 \n')
       OutPutFile.write(str(V1)+',3,-1 \n')
       OutPutFile.write(str(V4)+',3, 1 \n')


OutPutFile=open(str(outputDirectory)+'FrontToRear.inp','w')

OutPutFile.write('**** X-DIR \n')
for i in range(len(RearSet)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(RearSet[i])+',1,-1 \n')
       OutPutFile.write(str(FrontSet[i])+',1,1 \n')
       OutPutFile.write(str(V1)+',1,-1 \n')
       OutPutFile.write(str(H1)+',1,1 \n')


OutPutFile.write('**** \n')
OutPutFile.write('**** Y-DIR \n')
for i in range(len(RearSet)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(RearSet[i])+',2,-1 \n')
       OutPutFile.write(str(FrontSet[i])+',2,1 \n')
       OutPutFile.write(str(V1)+',2,-1 \n')
       OutPutFile.write(str(H1)+',2,1 \n')


OutPutFile.write('**** \n')
OutPutFile.write('**** Z-DIR \n')
for i in range(len(RearSet)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(RearSet[i])+',3,-1 \n')
       OutPutFile.write(str(FrontSet[i])+',3,1 \n')
       OutPutFile.write(str(V1)+',3,-1 \n')
       OutPutFile.write(str(H1)+',3,1 \n')


OutPutFile=open(str(outputDirectory)+'Edges.inp','w')


# Edges in x-y Plane
#right top edge to left top edge
OutPutFile.write('**** X-DIR \n')
for i in range(len(E_T2)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T2[i])+',1, 1 \n')
       OutPutFile.write(str(E_T4[i])+',1,-1 \n')
       OutPutFile.write(str(V2)+',1,-1 \n')
       OutPutFile.write(str(V1)+',1, 1 \n')

OutPutFile.write('**** Y-DIR \n')
for i in range(len(E_T2)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T2[i])+',2, 1 \n')
       OutPutFile.write(str(E_T4[i])+',2,-1 \n')
       OutPutFile.write(str(V2)+',2,-1 \n')
       OutPutFile.write(str(V1)+',2, 1 \n')

OutPutFile.write('**** Z-DIR \n')
for i in range(len(E_T2)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T2[i])+',3, 1 \n')
       OutPutFile.write(str(E_T4[i])+',3,-1 \n')
       OutPutFile.write(str(V2)+',3,-1 \n')
       OutPutFile.write(str(V1)+',3, 1 \n')


#right bottom edge to left bottom edge
OutPutFile.write('**** X-DIR \n')
for i in range(len(E_B2)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_B2[i])+',1, 1 \n')
       OutPutFile.write(str(E_B4[i])+',1,-1 \n')
       OutPutFile.write(str(V2)+',1,-1 \n')
       OutPutFile.write(str(V1)+',1, 1 \n')

OutPutFile.write('**** Y-DIR \n')
for i in range(len(E_B2)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_B2[i])+',2, 1 \n')
       OutPutFile.write(str(E_B4[i])+',2,-1 \n')
       OutPutFile.write(str(V2)+',2,-1 \n')
       OutPutFile.write(str(V1)+',2, 1 \n')

OutPutFile.write('**** Z-DIR \n')
for i in range(len(E_B2)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_B2[i])+',3, 1 \n')
       OutPutFile.write(str(E_B4[i])+',3,-1 \n')
       OutPutFile.write(str(V2)+',3,-1 \n')
       OutPutFile.write(str(V1)+',3, 1 \n')


#left top edge to left bottom edge
OutPutFile.write('**** X-DIR \n')
for i in range(len(E_T4)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T4[i])+',1, 1 \n')
       OutPutFile.write(str(E_B4[i])+',1,-1 \n')
       OutPutFile.write(str(V4)+',1,-1 \n')
       OutPutFile.write(str(V1)+',1, 1 \n')

OutPutFile.write('**** Y-DIR \n')
for i in range(len(E_T4)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T4[i])+',2, 1 \n')
       OutPutFile.write(str(E_B4[i])+',2,-1 \n')
       OutPutFile.write(str(V4)+',2,-1 \n')
       OutPutFile.write(str(V1)+',2, 1 \n')

OutPutFile.write('**** Z-DIR \n')
for i in range(len(E_T4)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T4[i])+',3, 1 \n')
       OutPutFile.write(str(E_B4[i])+',3,-1 \n')
       OutPutFile.write(str(V4)+',3,-1 \n')
       OutPutFile.write(str(V1)+',3, 1 \n')


# Edges in y-z Plane
# top back edge to top front edge
OutPutFile.write('**** X-DIR \n')
for i in range(len(E_T3)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T3[i])+',1, 1 \n')
       OutPutFile.write(str(E_T1[i])+',1,-1 \n')
       OutPutFile.write(str(H1)+',1,-1 \n')
       OutPutFile.write(str(V1)+',1, 1 \n')

OutPutFile.write('**** Y-DIR \n')
for i in range(len(E_T3)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T3[i])+',2, 1 \n')
       OutPutFile.write(str(E_T1[i])+',2,-1 \n')
       OutPutFile.write(str(H1)+',2,-1 \n')
       OutPutFile.write(str(V1)+',2, 1 \n')

OutPutFile.write('**** Z-DIR \n')
for i in range(len(E_T3)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T3[i])+',3, 1 \n')
       OutPutFile.write(str(E_T1[i])+',3,-1 \n')
       OutPutFile.write(str(H1)+',3,-1 \n')
       OutPutFile.write(str(V1)+',3, 1 \n')


# Botom back edge to bottom front edge
OutPutFile.write('**** X-DIR \n')
for i in range(len(E_B3)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_B3[i])+',1, 1 \n')
       OutPutFile.write(str(E_B1[i])+',1,-1 \n')
       OutPutFile.write(str(H1)+',1,-1 \n')
       OutPutFile.write(str(V1)+',1, 1 \n')

OutPutFile.write('**** Y-DIR \n')
for i in range(len(E_B3)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_B3[i])+',2, 1 \n')
       OutPutFile.write(str(E_B1[i])+',2,-1 \n')
       OutPutFile.write(str(H1)+',2,-1 \n')
       OutPutFile.write(str(V1)+',2, 1 \n')

OutPutFile.write('**** Z-DIR \n')
for i in range(len(E_B3)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_B3[i])+',3, 1 \n')
       OutPutFile.write(str(E_B1[i])+',3,-1 \n')
       OutPutFile.write(str(H1)+',3,-1 \n')
       OutPutFile.write(str(V1)+',3, 1 \n')


#top front edge to bottom front edge
OutPutFile.write('**** X-DIR \n')
for i in range(len(E_T1)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T1[i])+',1, 1 \n')
       OutPutFile.write(str(E_B1[i])+',1,-1 \n')
       OutPutFile.write(str(V4)+',1,-1 \n')
       OutPutFile.write(str(V1)+',1, 1 \n')

OutPutFile.write('**** Y-DIR \n')
for i in range(len(E_T1)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T1[i])+',2, 1 \n')
       OutPutFile.write(str(E_B1[i])+',2,-1 \n')
       OutPutFile.write(str(V4)+',2,-1 \n')
       OutPutFile.write(str(V1)+',2, 1 \n')

OutPutFile.write('**** Z-DIR \n')
for i in range(len(E_T1)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_T1[i])+',3, 1 \n')
       OutPutFile.write(str(E_B1[i])+',3,-1 \n')
       OutPutFile.write(str(V4)+',3,-1 \n')
       OutPutFile.write(str(V1)+',3, 1 \n')


# Edges in x-z Plane
# rear right edge to rear left edge
OutPutFile.write('**** X-DIR \n')
for i in range(len(E_M3)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_M3[i])+',1, 1 \n')
       OutPutFile.write(str(E_M4[i])+',1,-1 \n')
       OutPutFile.write(str(V2)+',1,-1 \n')
       OutPutFile.write(str(V1)+',1, 1 \n')

OutPutFile.write('**** Y-DIR \n')
for i in range(len(E_M3)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_M3[i])+',2, 1 \n')
       OutPutFile.write(str(E_M4[i])+',2,-1 \n')
       OutPutFile.write(str(V2)+',2,-1 \n')
       OutPutFile.write(str(V1)+',2, 1 \n')

OutPutFile.write('**** Z-DIR \n')
for i in range(len(E_M3)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_M3[i])+',3, 1 \n')
       OutPutFile.write(str(E_M4[i])+',3,-1 \n')
       OutPutFile.write(str(V2)+',3,-1 \n')
       OutPutFile.write(str(V1)+',3, 1 \n')


# front right edge to front left edge
OutPutFile.write('**** X-DIR \n')
for i in range(len(E_M2)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_M2[i])+',1, 1 \n')
       OutPutFile.write(str(E_M1[i])+',1,-1 \n')
       OutPutFile.write(str(V2)+',1,-1 \n')
       OutPutFile.write(str(V1)+',1, 1 \n')

OutPutFile.write('**** Y-DIR \n')
for i in range(len(E_M2)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_M2[i])+',2, 1 \n')
       OutPutFile.write(str(E_M1[i])+',2,-1 \n')
       OutPutFile.write(str(V2)+',2,-1 \n')
       OutPutFile.write(str(V1)+',2, 1 \n')

OutPutFile.write('**** Z-DIR \n')
for i in range(len(E_M2)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_M2[i])+',3, 1 \n')
       OutPutFile.write(str(E_M1[i])+',3,-1 \n')
       OutPutFile.write(str(V2)+',3,-1 \n')
       OutPutFile.write(str(V1)+',3, 1 \n')


#top front edge to bottom front edge
OutPutFile.write('**** X-DIR \n')
for i in range(len(E_M4)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_M4[i])+',1, 1 \n')
       OutPutFile.write(str(E_M1[i])+',1,-1 \n')
       OutPutFile.write(str(H1)+',1,-1 \n')
       OutPutFile.write(str(V1)+',1, 1 \n')

OutPutFile.write('**** Y-DIR \n')
for i in range(len(E_M4)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_M4[i])+',2, 1 \n')
       OutPutFile.write(str(E_M1[i])+',2,-1 \n')
       OutPutFile.write(str(H1)+',2,-1 \n')
       OutPutFile.write(str(V1)+',2, 1 \n')

OutPutFile.write('**** Z-DIR \n')
for i in range(len(E_M4)):
       #print item
       OutPutFile.write('*Equation \n')
       OutPutFile.write('4 \n')
       OutPutFile.write(str(E_M4[i])+',3, 1 \n')
       OutPutFile.write(str(E_M1[i])+',3,-1 \n')
       OutPutFile.write(str(H1)+',3,-1 \n')
       OutPutFile.write(str(V1)+',3, 1 \n')


OutPutFile=open(str(outputDirectory)+'Corners.inp','w')

#V3 zu V4
OutPutFile.write('**** X-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(V3)+',1, 1 \n')
OutPutFile.write(str(V4)+',1,-1 \n')
OutPutFile.write(str(V2)+',1,-1 \n')
OutPutFile.write(str(V1)+',1, 1 \n')
OutPutFile.write('**** y-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(V3)+',2, 1 \n')
OutPutFile.write(str(V4)+',2,-1 \n')
OutPutFile.write(str(V2)+',2,-1 \n')
OutPutFile.write(str(V1)+',2, 1 \n')
OutPutFile.write('**** z-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(V3)+',3, 1 \n')
OutPutFile.write(str(V4)+',3,-1 \n')
OutPutFile.write(str(V2)+',3,-1 \n')
OutPutFile.write(str(V1)+',3, 1 \n')

#H4 zu V4
OutPutFile.write('**** X-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(H4)+',1, 1 \n')
OutPutFile.write(str(V4)+',1,-1 \n')
OutPutFile.write(str(H1)+',1,-1 \n')
OutPutFile.write(str(V1)+',1, 1 \n')
OutPutFile.write('**** y-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(H4)+',2, 1 \n')
OutPutFile.write(str(V4)+',2,-1 \n')
OutPutFile.write(str(H1)+',2,-1 \n')
OutPutFile.write(str(V1)+',2, 1 \n')
OutPutFile.write('**** z-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(H4)+',3, 1 \n')
OutPutFile.write(str(V4)+',3,-1 \n')
OutPutFile.write(str(H1)+',3,-1 \n')
OutPutFile.write(str(V1)+',3, 1 \n')

#H3 zu V3
OutPutFile.write('**** X-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(H3)+',1, 1 \n')
OutPutFile.write(str(V3)+',1,-1 \n')
OutPutFile.write(str(H1)+',1,-1 \n')
OutPutFile.write(str(V1)+',1, 1 \n')
OutPutFile.write('**** y-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(H3)+',2, 1 \n')
OutPutFile.write(str(V3)+',2,-1 \n')
OutPutFile.write(str(H1)+',2,-1 \n')
OutPutFile.write(str(V1)+',2, 1 \n')
OutPutFile.write('**** z-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(H3)+',3, 1 \n')
OutPutFile.write(str(V3)+',3,-1 \n')
OutPutFile.write(str(H1)+',3,-1 \n')
OutPutFile.write(str(V1)+',3, 1 \n')

#H2 zu V2
OutPutFile.write('**** X-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(H2)+',1, 1 \n')
OutPutFile.write(str(V2)+',1,-1 \n')
OutPutFile.write(str(H1)+',1,-1 \n')
OutPutFile.write(str(V1)+',1, 1 \n')
OutPutFile.write('**** y-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(H2)+',2, 1 \n')
OutPutFile.write(str(V2)+',2,-1 \n')
OutPutFile.write(str(H1)+',2,-1 \n')
OutPutFile.write(str(V1)+',2, 1 \n')
OutPutFile.write('**** z-DIR \n')
OutPutFile.write('*Equation \n')
OutPutFile.write('4 \n')
OutPutFile.write(str(H2)+',3, 1 \n')
OutPutFile.write(str(V2)+',3,-1 \n')
OutPutFile.write(str(H1)+',3,-1 \n')
OutPutFile.write(str(V1)+',3, 1 \n')

OutPutFile=open(str(outputDirectory)+'VerticeSets.inp','w')
OutPutFile.write('*Nset, nset=V1, instance=PART-1-1 \n')
OutPutFile.write(str(V1)+'\n')
OutPutFile.write('*Nset, nset=V2, instance=PART-1-1 \n')
OutPutFile.write(str(V2)+'\n')
OutPutFile.write('*Nset, nset=V3, instance=PART-1-1 \n')
OutPutFile.write(str(V3)+'\n')
OutPutFile.write('*Nset, nset=V4, instance=PART-1-1 \n')
OutPutFile.write(str(V4)+'\n')
OutPutFile.write('*Nset, nset=H1, instance=PART-1-1 \n')
OutPutFile.write(str(H1)+'\n')
OutPutFile.write('*Nset, nset=H2, instance=PART-1-1 \n')
OutPutFile.write(str(H2)+'\n')
OutPutFile.write('*Nset, nset=H3, instance=PART-1-1 \n')
OutPutFile.write(str(H3)+'\n')
OutPutFile.write('*Nset, nset=H4, instance=PART-1-1 \n')
OutPutFile.write(str(H4)+'\n')
OutPutFile.write('*Nset, nset=Vertices \n')
OutPutFile.write('V1,V2,V4,H1\n')
# -------------------------------------------------------------------
# Rewrite Input File and Include Equation Files

def WriteInputFile(fileName):
       daten=open(fileName,'r')
       jumpParameter=0
       content=[]
       for zeile in daten:
       # Cut of the Formatting Parameters
              zeile=zeile.strip('\r\n').split('\n')
              #print zeile
              if zeile[0] == '*End Part':
                     content.append(['*Include, input=PeriodicData/LeftToRight.inp'])
                     content.append(['*Include, input=PeriodicData/BottomToTop.inp'])
                     content.append(['*Include, input=PeriodicData/FrontToRear.inp'])
                     content.append(['*Include, input=PeriodicData/Edges.inp'])
                     content.append(['*Include, input=PeriodicData/Corners.inp'])
                     content.append(['*End Part'])

              elif zeile[0] == '*End Assembly':
                     content.append(['*Include, input=PeriodicData/VerticeSets.inp'])
                     content.append(['*End Assembly'])
              else:
                     content.append(zeile)
       

       OutPutFile=open(str(outputDirectory)+'geometry_Periodic.inp','w')
       for item in content:
              OutPutFile.write(str(item[0])+'\n')
       

WriteInputFile(inputFile)

te=time.perf_counter()
print ('******')
print ('Periodicity Done!')
print ('******')

print (te-ta)
print (t2-t1)
print (t4-t3)

