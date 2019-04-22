# -*- coding: mbcs -*-
# Do not delete thefollowing import lines
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF,

session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(

p1 = mdb.models['Model-1'].parts['Backboard']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['Backboard']
p.features['Solidextrude-1'].setValues(depth=0.5) 
p = mdb.models['Model-1'].parts['Backboard']
p.regenerate()
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].view.setValues(nearPlane=129.419,
 farPlane=209.645, width=7.92903, height=4.22947,

a1 = mdb.models['Model-1'].rootAssembly
a1.translate(instanceList=('Backboard-1',),vector=(0.0,0.0,0.1)) 
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(

a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['Backboard-1'], )
a.generateMesh(regions=partInstances)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(

mdb.Job(name='Job-38', model='Model-1', description='', type=ANALYSIS,
 atTime=None, waitMinutes=0, waitHours=0, queue=None,
 memoryUnits=PERCENTAGE,
 explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
 modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
 scratch='', multiprocessingMode=DEFAULT, numCpus=1,
mdb.jobs['NewJob'].writeInput(consistencyChecking=OFF)
