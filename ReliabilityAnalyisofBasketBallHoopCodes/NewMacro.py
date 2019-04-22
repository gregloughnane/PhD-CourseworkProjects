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
o1=session.openOdb(name='PertvsNew.odb') 
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
odb=session.odbs['PertvsNew.odb'] 
session.writeFieldReport(fileName='Stresses_vs.rpt',append=OFF, 
sortItem='Element Label', odb=odb, step=0, frame=1,
outputPosition=INTEGRATION_POINT, variable=(('S', INTEGRATION_POINT, ((INVARIANT, 'Mises'),)), ))
