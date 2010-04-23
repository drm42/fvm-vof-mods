#!/usr/bin/env python
import sys

import fvm
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers

import time
from numpy import *
import tecplotESBGK
import tecplot

fvm.set_atype('double')

from FluentCase import FluentCase

#fvmbaseExt.enableDebug("cdtor")

fileBase = None
numIterations = 10
fileBase = "/home/aerosun/a/schigull/FVM-trial/src/fvm/test/testKineticFlowModel"
def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)


# change as needed

outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"-tecplt.dat"

#import debug

   
reader = FluentCase(fileBase+".cas")

reader.read();
import sys

#fluent_meshes = reader.getMeshList()

import time
t0 = time.time()


meshes = reader.getMeshList()

geomFields =  fvm.models.GeomFields('geom')
metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()
#!!!!!!!!!!!!!!!!source code (metricsMeshCalcultaro)::init(), volumeField.synLocal, coordField.syncLocal() turned off

flowFields =  fvm.models.FlowFields('flow')

fmodel = fvm.models.FlowModelA(geomFields,flowFields,meshes)

## set bc for top to be a wall with x velocity
bcMap = fmodel.getBCMap()
if 3 in bcMap:
   bc3 = fmodel.getBCMap()[3]
   bc3.bcType = 'NoSlipWall'
   bc3.setVar('specifiedXVelocity',1)
				

## set viscosity and density, this is done per mesh since each mesh has its own VC object
vcMap = fmodel.getVCMap()
for vc in vcMap.values():
    vc.setVar('density',1.0)
    vc.setVar('viscosity',1.0)


momSolver = fvmbaseExt.AMG()
contSolver = fvmbaseExt.AMG()

# to test quadrature


import fvm.esbgk_atyped_double as esbgk

#import fvm.MacroParameters as macropr
#import fvm.DistFunctFields as f

foptions = fmodel.getOptions()
fmodel.init()

#kineticmodel=fvm.models.KineticModelA(meshes,flowFields,macroFields,quad)
#kineticmodel.init()

#cartesian
#import ddd
quad0=esbgk.QuadratureD(10,10,10,5.5,1.0) #cartesian
quad1=esbgk.QuadratureD(0,16,0,32,0,16) #spherical
quad2=esbgk.QuadratureD(16,16,1,8,1,4) #gauss-hermit quadrature and 3/8th rule

macroFields=esbgk.MacroFields('flow')
esbgk1=esbgk.KineticModelD(meshes,geomFields,macroFields,quad0)

#dsff=esbgk.DistFunctFields(
#esbgk.KineticModelD.InitializeMacroparameters()
esbgk1.InitializeMacroparameters()
esbgk1.advance(1)

#numIterations=1
#def advance(niter):
#    for i in range(0,niter):
#        esbgk1.advance(1)
#       
#advance(numIterations)        


cellSite = meshes[0].getCells()
velField = flowFields.velocity[ cellSite ].asNumPyArray() 
#print  flowFields.velocity[ cellSite ].asNumPyArray()

#import debug
#fmodel.advance(100)
tecplotESBGK.esbgkTecplotFile(meshes,macroFields)
tecplot.dumpTecplotFile(1,meshes,macroFields)
