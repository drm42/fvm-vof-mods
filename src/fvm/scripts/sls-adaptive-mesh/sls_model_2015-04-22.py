import sys
import numpy as np
import math
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
import pdb
import os


class Point:
  def __init__(self,x,y):
    self.x = x
    self.y = y

class Rectangle:
  def __init__(self,p1,p2,p3,p4):
    self.verts = [p1,p2,p3,p4]
    self.findMinMax()

  def findMinMax(self):
    self.max_x = self.verts[0].x
    self.max_y = self.verts[0].y
    self.min_x = self.verts[0].x
    self.min_y = self.verts[0].y
    for i in range(1,4):
      if self.verts[i].x > self.max_x:
        self.max_x = self.verts[i].x
      if self.verts[i].x < self.min_x:
        self.min_x = self.verts[i].x
      if self.verts[i].y > self.max_y:
        self.max_y = self.verts[i].y
      if self.verts[i].y < self.min_y:
        self.min_y = self.verts[i].y

  def rotate(self,theta):
    cx = (self.max_x + self.min_x)/2.0
    cy = (self.max_y + self.min_y)/2.0

    for i in range(4):
        tempX = self.verts[i].x - cx
        tempY = self.verts[i].y - cy

        rotatedX = tempX*math.cos(theta) - tempY*math.sin(theta)
        rotatedY = tempY*math.cos(theta) + tempX*math.sin(theta)
        
        self.verts[i] = Point(rotatedX + cx, rotatedY + cy)
    self.findMinMax()

#Laser diameter
laser_diam = 0.00008
#Powder bed reflectivity
R = 0.82
#Min and max timestep
minTimeStep = 1.e-5
maxTimeStep = 1.e-4
#Initial bed temperature
initTemp = 298.0
#Base filename for output vtk files
outfilebase = 'sls_model_new'
#Max number of iterations per linear solve
maxIterPerTimeStep = 100
#File with specific heat function data
specificHeatFile = '/home/drmoser/Documents/sls_model-2015-02-18/cp_model'
#Latenet heat of fusion
latentHeat = 321000.0
#Solid and liquidus temperatures
solidTemp = 855.0
liquidTemp = 925.0
#Solid thermal conductivity
thermalConductivity = 170.0
#Powder density
density = 1350.0
#Powder extinction coefficient
beta = 9914.0
#Initial Domain Size
nx = 1000
ny = 1000
nz = 2
xlen = 0.10
ylen = 0.10
zlen = 0.0002
#Layer thickness
zlayer = 0.0001
nlayer = 1
#Output vtk file frequency
outFreq = 5


#Set defaults for system variables
laser_pos = Point(0,0)
laser_pos_next = Point(0,0)
v_laser = [0, 0]
laser_power = 0.0
laser_power_next = 0.0
sqrt2 = math.sqrt(2)
time = 0
timeStep = minTimeStep
numTimeSteps = 1
curr_cycle = -1
v_laser = []
int_time = -1
###############################################################################
##Define how conductivity varies with temperature
def condFunction(T):
  if(T < solidTemp):
    cond = 0.07*(T-273.0)+190
  else:
    cond = 0.07*(solidTemp-273.0)+190
  return cond
###############################################################################

###############################################################################
##Define how conductivity varies with melt fraction
def condMeltFunction(k_s,frac):
  return k_s*0.05+0.95*k_s*frac

###############################################################################
##Determine min and max extents for a cell
def getNodeCoords(c,CellCoord,NodeCoord,nodes):
  x = CellCoord[c,0]
  y = CellCoord[c,1]
  z = CellCoord[c,2]
  lim = [x,x,y,y,z,z]
  for nindex in range(nodes.getCount(c)):
    n = nodes(c,nindex)
    if NodeCoord[n,0] > lim[1]:
      lim[1] = NodeCoord[n,0]
    if NodeCoord[n,0] < lim[0]:
      lim[0] = NodeCoord[n,0]
    if NodeCoord[n,1] > lim[3]:
      lim[3] = NodeCoord[n,1]
    if NodeCoord[n,1] < lim[2]:
      lim[2] = NodeCoord[n,1]
    if NodeCoord[n,2] > lim[5]:
      lim[5] = NodeCoord[n,2]
    if NodeCoord[n,2] < lim[4]:
      lim[4] = NodeCoord[n,2]
  return lim
###############################################################################

###############################################################################
##Define the retangle made by the laser path
def getLaserRect(p1,p2,v,r):
  x0 = p1.x
  y0 = p1.y
  x1 = p2.x
  y1 = p2.y
  theta = math.atan2(v[1],v[0])
  p1 = Point(x0-r*math.sin(theta),y0+r*math.cos(theta))
  p2 = Point(x0+r*math.sin(theta),y0-r*math.cos(theta))
  p3 = Point(x1-r*math.sin(theta),y1+r*math.cos(theta))
  p4 = Point(x1+r*math.sin(theta),y1-r*math.cos(theta))
  rect = Rectangle(p1,p3,p4,p2)
  return rect
###############################################################################

###############################################################################
#Determine if a cell and a rectangle overlap
def rectIntersection(lim,rect):
  result = True
  test = [0,0,0,0,0]
  if (rect.max_x < lim[0]) or (rect.min_x > lim[1]):
    result = False
  elif (rect.max_y < lim[2]) or (rect.min_y > lim[3]):
    result = False
  else:
    for i in range(4):
      j = (i + 1)%4
      test[0] = ((rect.verts[j].x-rect.verts[i].x)*(lim[2]-rect.verts[i].y)-
          (rect.verts[j].y-rect.verts[i].y)*(lim[0]-rect.verts[i].x))
      if test[0] == 0:
        continue
      test[1] = ((rect.verts[j].x-rect.verts[i].x)*(lim[2]-rect.verts[i].y)-
          (rect.verts[j].y-rect.verts[i].y)*(lim[1]-rect.verts[i].x))
      if (test[1] == 0) or ((test[0]<0)!=(test[1]<0)):
        continue
      test[2] = ((rect.verts[j].x-rect.verts[i].x)*(lim[3]-rect.verts[i].y)-
          (rect.verts[j].y-rect.verts[i].y)*(lim[0]-rect.verts[i].x))
      if (test[2] == 0) or ((test[1]<0)!=(test[2]<0)):
        continue
      test[3] = ((rect.verts[j].x-rect.verts[i].x)*(lim[3]-rect.verts[i].y)-
          (rect.verts[j].y-rect.verts[i].y)*(lim[1]-rect.verts[i].x))
      if (test[3] == 0) or ((test[2]<0)!=(test[3]<0)):
        continue
      for k in range(4):
        if k == i or k == j:
          continue
        else:
          test[4] = ((rect.verts[j].x-rect.verts[i].x)*(rect.verts[k].y-
                        rect.verts[i].y)-(rect.verts[j].y-rect.verts[i].y)*
                        (rect.verts[k].x-rect.verts[i].x))
          break
      if ((test[4]<0)!=(test[3]<0)):
        result = False
        break 
  return result
###############################################################################

###############################################################################
#Determine if a cell and a circle overlap
def circleIntersection(lim,center,r):
  rcenter = Point((lim[0] + lim[1])/2.0, (lim[2] + lim[3])/2.0)

  dist_x = abs(center.x-rcenter.x)
  dist_y = abs(center.y-rcenter.y)

  if dist_x > ((lim[1]-lim[0])/2.0 + r):
    return False
  if dist_y > ((lim[3]-lim[2])/2.0 + r):
    return False

  if (dist_x <= ((lim[1]-lim[0])/2.0)):
    return True
  if (dist_y <= ((lim[3]-lim[2])/2.0)):
    return True

  cornerDist = ((dist_x - ((lim[1]-lim[0])/2.0))**2+(dist_y - 
               ((lim[3]-lim[2])/2.0))**2)
  return (cornerDist <= r**2)
###############################################################################

###############################################################################
#Determine if the laser will effect a cell
def cellLaserOverlap(lim,p1,p2,r,rect):
  result = False
  if(rectIntersection(lim,rect)):
    result = True
  elif(circleIntersection(lim,p1,r)):
    result = True
  elif(circleIntersection(lim,p2,r)):
    result = True
  return result
###############################################################################
##Calculate the times at which a laser will overlap a cell during a timestep
def getTimes(lim, center, r, v, stepSize):
  cx = center.x
  cy = center.y
  vx = v[0]
  vy = v[1]
  xmin = lim[0] - r
  xmax = lim[1] + r
  ymin = lim[2] - r
  ymax = lim[3] + r
  result = [0.0, stepSize]

  if (vx > 0):
    tmin = (xmin - cx)/vx
    tmax = (xmax - cx)/vx
  elif (vx < 0):
    tmin = (xmax - cx)/vx
    tmax = (xmin - cx)/vx
  else:
    if((cx > xmin) and (cx < xmax)):
      tmin = float('-inf')
      tmax = float('inf')
    else:
      tmin = float('inf')
      tmax = float('-inf')

  if (vy > 0):
    tymin = (ymin - cy)/vy
    tymax = (ymax - cy)/vy
  elif (vy < 0):
    tymin = (ymax - cy)/vy
    tymax = (ymin - cy)/vy
  else:
    if((cy > ymin) and (cy < ymax)):
      tymin = float('-inf')
      tymax = float('inf')
    else:
      tymin = float('inf')
      tymax = float('-inf')

  if ((tmin > tymax) or (tymin > tmax)):
    return []
  
  tmin = max(tmin, tymin)
  tmax = min(tmax, tymax)

  if (tmin >= stepSize):
    return []
  if (tmax <= 0.0):
    return []
  if (tmin > 0.0):
    result[0] = tmin
  if (tmax < stepSize):
    result[1] = tmax

  return result

###############################################################################
# Recursive generation of the Legendre polynomial of order n
def Legendre(n,x):
  if (n==0):
    return x*0+1.0
  elif (n==1):
    return x
  else:
    return ((2.0*n-1.0)*x*Legendre(n-1,x)-(n-1)*Legendre(n-2,x))/n
###############################################################################
 
###############################################################################
# Derivative of the Legendre polynomials
def DLegendre(n,x):
  if (n==0):
    return x*0
  elif (n==1):
    return x*0+1.0
  else:
    return (n/(x**2-1.0))*(x*Legendre(n,x)-Legendre(n-1,x))
###############################################################################

###############################################################################
# Roots of the polynomial obtained using Newton-Raphson method
def LegendreRoots(polyorder,tolerance=1e-20):
  if polyorder<2:
    err=1 # bad polyorder no roots can be found
  else:
    roots=[]
    # The polynomials are alternately even and odd functions. 
    #So we evaluate only half the number of roots. 
    for i in range(1,polyorder+1):
      x=math.cos(math.pi*(i-0.25)/(polyorder+0.5))
      error=10*tolerance
      iters=0
      while (error>tolerance) and (iters<1000):
        dx=-Legendre(polyorder,x)/DLegendre(polyorder,x)
	x=x+dx
	iters=iters+1
	error=abs(dx)
      roots.append(x)
  return roots
###############################################################################

###############################################################################
# Weight coefficients
def GaussLegendreWeights(polyorder,xis):
  weights=[]
  for x in xis:
    W=2.0/( (1.0-x**2)*(DLegendre(polyorder,x)**2) )
    weights.append(W)
  return weights
###############################################################################

###############################################################################
# The integral value 
# func 		   : the integrand
# a, b 		   : lower and upper limits of the integral
# xs, Ws           : the roots and weights for Gauss Legrendre
# lim, laser, P, v : values needed by the laser source function
def GaussLegendreQuadrature(func, xs, Ws, a, b, lim, laser, P, v):
  ans = 0
  for i in range(len(xs)):
    ans=ans + Ws[i]*func((b-a)*0.5*xs[i]+ (b+a)*0.5,lim,laser,P,v) 
  ans = ans * (b-a)*0.5
  return ans
###############################################################################

###############################################################################
#NOTE: This method of calculating the source term assumes the surface of the 
#powder bed is at z=0
def laserSource(t,lim,laser,P,v):
  x1 = lim[0]
  x2 = lim[1]
  y1 = lim[2]
  y2 = lim[3]
  z1 = lim[4]-zlen
  z2 = lim[5]-zlen
  x_l = laser.x
  y_l = laser.y
  vol = (x2-x1)*(y2-y1)*(z2-z1)
  w = laser_diam/2/2.146
  I = 2*P/math.pi/w**2
  source = (1-R)*I*(w**2)*math.pi/8*((math.erf(sqrt2*(x2-x_l-v[0]*t)/w)-
                                    math.erf(sqrt2*(x1-x_l-v[0]*t)/w))*
                                    (math.erf(sqrt2*(y2-y_l-v[1]*t)/w)-
                                    math.erf(sqrt2*(y1-y_l-v[1]*t)/w))*
                                    (math.exp(beta*z2)-math.exp(beta*z1)))
  source = source/vol
  return source
###############################################################################

###############################################################################
def generateBoxMesh(xmax, ymax, zmax, imax, jmax, kmax):


    xcoord = np.linspace(0.0, xmax, num=(imax+1))
    xcoord = xcoord.astype('double')
    ycoord = np.linspace(0.0, ymax, num=(jmax+1))
    ycoord = ycoord.astype('double')
    zcoord = np.linspace(0.0, zmax, num=(kmax+1))
    zcoord = zcoord.astype('double')

    nCells = imax*jmax*kmax

    nodeIndex = np.zeros(shape=(imax+1,jmax+1,kmax+1), dtype='int')
    nNodes = 0

    nodeCoordsN = models.newVec3Array((imax+1)*(jmax+1)*(kmax+1))
    nodeCoordsA = nodeCoordsN.asNumPyArray()

    for k in range(kmax+1):
        for j in range(jmax+1):
            for i in range(imax+1):
                nodeIndex[i,j,k] = nNodes
                nodeCoordsA[nNodes,0] = xcoord[i]
                nodeCoordsA[nNodes,1] = ycoord[j]
                nodeCoordsA[nNodes,2] = zcoord[k]
                nNodes = nNodes + 1

    cellIndex = np.zeros(shape=(imax,jmax,kmax), dtype='int')
    nCells = 0

    for k in range(kmax):
        for j in range(jmax):
            for i in range(imax):
                cellIndex[i,j,k] = nCells
                nCells = nCells + 1
                
    nFacesInterior = imax*jmax*(kmax-1)+imax*(jmax-1)*kmax+(imax-1)*jmax*kmax

    nFaceZones = 7

    faceGroupCountN = fvmbaseExt.newIntArray(nFaceZones)
    faceGroupCount = faceGroupCountN.asNumPyArray()

    #interior faces
    faceGroupCount[0] = nFacesInterior 
    #xmin faces
    faceGroupCount[1] = jmax*kmax
    #xmax faces
    faceGroupCount[2] = jmax*kmax
    #ymin faces
    faceGroupCount[3] = imax*kmax
    #ymax faces
    faceGroupCount[4] = imax*kmax
    #zmin faces
    faceGroupCount[5] = imax*jmax
    #zmax faces
    faceGroupCount[6] = imax*jmax

    nFaces = int(faceGroupCount.sum())

    ## allocate arrays for face nodes and cells
    faceNodeCountN = fvmbaseExt.newIntArray(nFaces)
    faceNodesN = fvmbaseExt.newIntArray(nFaces*4)
    faceCellsN = fvmbaseExt.newIntArray(nFaces*2)

    faceNodesA = faceNodesN.asNumPyArray()
    faceCellsA = faceCellsN.asNumPyArray()

    faceNodeCountA = faceNodeCountN.asNumPyArray()
    faceNodeCountA[:] = 4
    
    ## reshape for convenience

    faceNodes = faceNodesA.reshape((nFaces,4))
    faceCells = faceCellsA.reshape((nFaces,2))
    
    nf = 0
    # interior x faces
    for i in range(1,imax):
        for j in range(jmax):
            for k in range(kmax):
                faceNodes[nf,0] = nodeIndex[i,j,k]
                faceNodes[nf,1] = nodeIndex[i,j+1,k]
                faceNodes[nf,2] = nodeIndex[i,j+1,k+1]
                faceNodes[nf,3] = nodeIndex[i,j,k+1]

                faceCells[nf,0] = cellIndex[i,j,k]
                faceCells[nf,1] = cellIndex[i-1,j,k]

                nf = nf + 1

    # interior y faces
    for j in range(1,jmax):
        for i in range(imax):
            for k in range(kmax):
                faceNodes[nf,0] = nodeIndex[i,j,k]
                faceNodes[nf,1] = nodeIndex[i+1,j,k]
                faceNodes[nf,2] = nodeIndex[i+1,j,k+1]
                faceNodes[nf,3] = nodeIndex[i,j,k+1]

                faceCells[nf,0] = cellIndex[i,j-1,k]
                faceCells[nf,1] = cellIndex[i,j,k]

                nf = nf + 1

    # interior z faces
    for k in range(1,kmax):
        for i in range(imax):
            for j in range(jmax):
                faceNodes[nf,0] = nodeIndex[i,j,k]
                faceNodes[nf,1] = nodeIndex[i+1,j,k]
                faceNodes[nf,2] = nodeIndex[i+1,j+1,k]
                faceNodes[nf,3] = nodeIndex[i,j+1,k]

                faceCells[nf,0] = cellIndex[i,j,k]
                faceCells[nf,1] = cellIndex[i,j,k-1]

                nf = nf + 1
        
    
    
    nb = 0
    #x-boundaries
    i = 0
    for j in range(jmax):
        for k in range(kmax):
            faceNodes[nf,0] = nodeIndex[i,j,k]
            faceNodes[nf,1] = nodeIndex[i,j+1,k]
            faceNodes[nf,2] = nodeIndex[i,j+1,k+1]
            faceNodes[nf,3] = nodeIndex[i,j,k+1]

            faceCells[nf,0] = cellIndex[i,j,k]
            faceCells[nf,1] = nCells + nb

            nf = nf + 1
            nb = nb + 1

    i = imax
    for j in range(jmax):
        for k in range(kmax):
            faceNodes[nf,0] = nodeIndex[i,j,k]
            faceNodes[nf,1] = nodeIndex[i,j,k+1]
            faceNodes[nf,2] = nodeIndex[i,j+1,k+1]
            faceNodes[nf,3] = nodeIndex[i,j+1,k]

            faceCells[nf,0] = cellIndex[i-1,j,k]
            faceCells[nf,1] = nCells + nb

            nf = nf + 1
            nb = nb + 1

    #y-boundaries
    j = 0
    for i in range(imax):
        for k in range(kmax):
            faceNodes[nf,0] = nodeIndex[i,j,k]
            faceNodes[nf,1] = nodeIndex[i,j,k+1]
            faceNodes[nf,2] = nodeIndex[i+1,j,k+1]
            faceNodes[nf,3] = nodeIndex[i+1,j,k]

            faceCells[nf,0] = cellIndex[i,j,k]
            faceCells[nf,1] = nCells + nb

            nf = nf + 1
            nb = nb + 1

    j = jmax
    for i in range(imax):
        for k in range(kmax):
            faceNodes[nf,0] = nodeIndex[i,j,k]
            faceNodes[nf,1] = nodeIndex[i+1,j,k]
            faceNodes[nf,2] = nodeIndex[i+1,j,k+1]
            faceNodes[nf,3] = nodeIndex[i,j,k+1]

            faceCells[nf,0] = cellIndex[i,j-1,k]
            faceCells[nf,1] = nCells + nb

            nf = nf + 1
            nb = nb + 1

    #z-boundaries
    k = 0
    for i in range(imax):
        for j in range(jmax):
            faceNodes[nf,0] = nodeIndex[i,j,k]
            faceNodes[nf,1] = nodeIndex[i+1,j,k]
            faceNodes[nf,2] = nodeIndex[i+1,j+1,k]
            faceNodes[nf,3] = nodeIndex[i,j+1,k]

            faceCells[nf,0] = cellIndex[i,j,k]
            faceCells[nf,1] = nCells + nb

            nf = nf + 1
            nb = nb + 1

    k = kmax
    for i in range(imax):
        for j in range(jmax):
            faceNodes[nf,0] = nodeIndex[i,j,k]
            faceNodes[nf,1] = nodeIndex[i,j+1,k]
            faceNodes[nf,2] = nodeIndex[i+1,j+1,k]
            faceNodes[nf,3] = nodeIndex[i+1,j,k]

            faceCells[nf,0] = cellIndex[i,j,k-1]
            faceCells[nf,1] = nCells + nb

            nf = nf + 1
            nb = nb + 1
    

    mesh = fvmbaseExt.Mesh(3, nCells, nodeCoordsN, faceCellsN,
                           faceNodesN, faceNodeCountN, faceGroupCountN)

    return mesh
##########################################################################


####################################################################
def write_restart_file(values):
  old_ncells, old_temp, old_tempN1, old_enth, old_enthN1, old_enthInv, old_dHdT, old_meltFrac = values
  if(os.path.exists('restart_file')):
    os.rename('restart_file', 'restart_file.bak')
  f = open('restart_file', 'w')
  f.write(str(xlen) + ' ' + str(ylen) + ' ' + str(zlen) + '\n')
  f.write(str(nx) + ' ' + str(ny) + ' ' + str(nz) + '\n')
  f.write(str(time) + ' ' + str(numLayers) + ' ' + str(numTimeSteps) + ' ' + str(laser_pos.x) + ' ' + str(laser_pos.y) + '\n')
  f.write(str(old_ncells) + '\n')
  f.write(str(int_time) + ' ' + str(curr_cycle) + '\n')
  for c in range(old_ncells):
    f.write(str(old_temp[c]) + ' ' + str(old_tempN1[c]) + ' ' +
            str(old_enth[c]) + ' ' + str(old_enthN1[c]) + ' ' +
            str(old_enthInv[c]) + ' ' + str(old_dHdT[c]) + ' ' +
            str(old_meltFrac[c]) + '\n')
  f.close()
####################################################################
def read_restart_file():
  global xlen
  global ylen
  global zlen
  global nx
  global ny
  global nz
  global time
  global laser_pos
  global numLayers
  global numTimeSteps
  global int_time
  global curr_cycle

  f = open('restart_file', 'r')
  lengths = f.readline()
  lengths = lengths.split()
  xlen = float(lengths[0])
  ylen = float(lengths[1])
  zlen = float(lengths[2])
  n = f.readline()
  n = n.split()
  nx = int(n[0])
  ny = int(n[1])
  nz = int(n[2])
  t_data = f.readline()
  t_data = t_data.split()
  time = float(t_data[0])
  numLayers = int(t_data[1])
  numTimeSteps = int(t_data[2]) + 1
  laser_pos.x = float(t_data[3])
  laser_pos.y = float(t_data[4])

  ncells = int(f.readline())
  angle = f.readline()
  angle = angle.split()
  int_time = float(angle[0])
  curr_cycle = int(angle[1])
  
  temp = np.zeros(ncells)
  tempN1 = np.zeros(ncells)
  enth = np.zeros(ncells)
  enthN1 = np.zeros(ncells)
  enthInv = np.zeros(ncells)
  dHdT = np.zeros(ncells)
  meltFrac = np.zeros(ncells)
  c = 0

  for line in f:
    words = line.split()
    temp[c] = float(words[0])
    tempN1[c] = float(words[1])
    enth[c] = float(words[2])
    enthN1[c] = float(words[3])
    enthInv[c] = float(words[4])
    dHdT[c] = float(words[5])
    meltFrac[c] = float(words[6])
    c = c + 1

  f.close()
  return (ncells, temp, tempN1, enth, enthN1, 
          enthInv, dHdT, meltFrac)
####################################################################

##########################################################################
def RunSim(xmax,ymax,zmax,imax,jmax,kmax,old_values,laser_file):
#pdb.set_trace()
    mesh = generateBoxMesh(xmax, ymax, zmax, imax, jmax, kmax)
    meshes = (mesh,)
    geomFields = models.GeomFields('geom')
    metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)
    metricsCalculator.init()
    thermalFields = models.ThermalFields('temperature') 
    tmodel = models.ThermalModelA(geomFields,thermalFields,meshes) 

    bcmap = tmodel.getBCMap()
    vcmap = tmodel.getVCMap()
    
    global time
    global numTimeSteps
    global curr_cycle
    global v_laser
    global laser_pos
    global laser_pos_next
    global int_time

    vcZone = vcmap[mesh.getID()]
    # Set boundary conditions and material properties
    vcZone['thermalConductivity'] = thermalConductivity
    vcZone['density'] = density
#    vcZone['specificHeat'] = 1546

    xmin_id = 1
    xmax_id = 2
    ymin_id = 3
    ymax_id = 4
    zmin_id = 5
    zmax_id = 6

    bc_xmin = bcmap[xmin_id]
    bc_xmax = bcmap[xmax_id]
    bc_ymin = bcmap[ymin_id]
    bc_ymax = bcmap[ymax_id]
    bc_zmin = bcmap[zmin_id]
    bc_zmax = bcmap[zmax_id]

    bc_zmax.bcType = 'Mixed'
    bc_zmax.setVar('farFieldTemperature',initTemp)
    bc_zmax.setVar('surfaceEmissivity',-1.0)
    bc_zmax.setVar('convectiveCoefficient',-30)

    bc_zmin.bcType = 'Convective'
    bc_zmin.setVar('farFieldTemperature',initTemp)
    bc_zmin.setVar('convectiveCoefficient',-1332)
#    bc_zmin.setVar('specifiedHeatFlux',0)

    bc_ymin.bcType = 'SpecifiedTemperature'
    bc_ymin.setVar('specifiedTemperature',298.0)

    bc_ymax.bcType = 'SpecifiedTemperature'
    bc_ymax.setVar('specifiedTemperature',298.0)

    bc_xmin.bcType = 'SpecifiedTemperature'
    bc_xmin.setVar('specifiedTemperature',298.0)

    bc_xmax.bcType = 'SpecifiedTemperature'
    bc_xmax.setVar('specifiedTemperature',298.0)

    # Set solver options
    toptions = tmodel.getOptions()
    toptions.transient = True
    toptions.enthalpyModel = True
    toptions.polynomialCp = True
    toptions['latentHeat'] = latentHeat
    toptions['solidTemp'] = solidTemp
    toptions['liquidTemp'] = liquidTemp
    toptions.cpFile = specificHeatFile
    toptions.setVar('initialTemperature', initTemp)
    pc = fvmbaseExt.AMG()
    pc.verbosity = 0
    defSolver = fvmbaseExt.BCGStab()
    defSolver.preconditioner = pc
    defSolver.relativeTolerance = 1.e-10
    defSolver.absoluteTolerance = 1.e-10
    defSolver.nMaxIteractions = 1000
    defSolver.verbosity = 0
    toptions.linearSolver = defSolver
    tmodel.init()

    cells = mesh.getCells()
    ncells = cells.getSelfCount()

    old_ncells, old_temp, old_tempN1, old_enth, old_enthN1, old_enthInv, old_dHdT, old_meltFrac = old_values

    temp = thermalFields.temperature
    tempN1 = thermalFields.temperatureN1
    enth = thermalFields.enthalpy
    enthN1 = thermalFields.enthalpyN1
    enthInv = thermalFields.enthalpyInverse
    dHdT = thermalFields.dHdT
    meltFrac = thermalFields.meltFrac
    temp = temp[cells].asNumPyArray()
    tempN1 = tempN1[cells].asNumPyArray()
    enth = enth[cells].asNumPyArray()
    enthN1 = enthN1[cells].asNumPyArray()
    enthInv = enthInv[cells].asNumPyArray()
    dHdT = dHdT[cells].asNumPyArray()
    meltFrac = meltFrac[cells].asNumPyArray()
    ##recover data from previous mesh if it exists
    if(old_ncells != 0):
      for c in range(old_ncells):
        temp[c] = old_temp[c]
        tempN1[c] = old_tempN1[c]
        enth[c] = old_enth[c]
        enthN1[c] = old_enthN1[c]
        enthInv[c] = old_enthInv[c]
        dHdT[c] = old_dHdT[c]
        meltFrac[c] = old_meltFrac[c]

    
    nodes=mesh.getCellNodes()
    NodeCoord=mesh.getNodeCoordinates().asNumPyArray()
    CellCoord=geomFields.coordinate[cells].asNumPyArray()
    laser_width = 3*laser_diam/4/2.146
    laser_path = open(laser_file,'r')
    if (curr_cycle == -1): 
      curr_cycle = 1
      line = laser_path.readline()
      words = line.split()
      if(time != float(words[0])):
         raise NameError('Time does not match start time in laser_file')
      laser_pos.x = float(words[1])
      laser_pos.y = float(words[2])
      laser_power = float(words[3])
      line = laser_path.readline()
      words = line.split()
      int_time = float(words[0])-time
      v_laser = [(float(words[1])-laser_pos.x)/int_time,(float(words[2])-laser_pos.y)/int_time]
      laser_power_next = float(words[3])
    else:
      line = laser_path.readline()
      words = line.split()
      laser_power = float(words[3])
      if (float(words[0]) > time):
        raise NameError('Current time not within laser_file for restart')
      else:
        while(True):
          line = laser_path.readline()
          words = line.split()
          if (float(words[0]) > time):
            int_time = float(words[0])-time
            v_laser = [(float(words[1])-laser_pos.x)/int_time,(float(words[2])-laser_pos.y)/int_time]
            laser_power_next = float(words[3])
            break
          else:
            laser_power = float(words[3])

    f = open('laser_position','a')
    f.write(str(time) + ' ' + str(numTimeSteps) + ' ' + str(laser_pos.x) + ' ' + str(laser_pos.y) + '\n')

    MeltFracField = thermalFields.meltFrac
    meltFracArray = MeltFracField[cells].asNumPyArray()
    for c in range(ncells):
      lim = getNodeCoords(c,CellCoord,NodeCoord,nodes)
      if(lim[0] >= 0.01) and (lim[1] <= 0.09) and (lim[2] >= 0.01) and (lim[3] <= 0.09) and (lim[5] <= 0.0001):
        meltFracArray[c] = 1.0
      else:
        meltFracArray[c] = 0.0

    while(True):
      #find the timestep.
      if (int_time <= maxTimeStep):
        timeStep = int_time
      else:
        timeStep = maxTimeStep
      laser_pos_next.x = laser_pos.x+v_laser[0]*timeStep
      laser_pos_next.y = laser_pos.y+v_laser[1]*timeStep  
      int_time = int_time - timeStep

      toptions.setVar('timeStep',timeStep)
      print "Current Time: " + str(time)
      print "Time Step: " + str(timeStep)

      laser_rect = getLaserRect(laser_pos,laser_pos_next,v_laser,laser_width)

      SourceField = thermalFields.source
      sourceArray = SourceField[cells].asNumPyArray()
      CondField = thermalFields.conductivity
      condArray = CondField[cells].asNumPyArray()
      TempField = thermalFields.temperature
      tempArray = TempField[cells].asNumPyArray()
      MeltFracField = thermalFields.meltFrac
      meltFracArray = MeltFracField[cells].asNumPyArray()

      for c in range(ncells):
        #set source term
        lim = getNodeCoords(c,CellCoord,NodeCoord,nodes)
        if ((lim[5]-zlen) < -1/beta*math.log(0.01*beta/(1-R))):
          source = 0
        elif (cellLaserOverlap(lim,laser_pos,laser_pos_next,laser_width,
                               laser_rect)):
          ts = getTimes(lim,laser_pos,laser_width,v_laser,timeStep)
          if (len(ts) == 2):
            source = GaussLegendreQuadrature(laserSource,xs,Ws,ts[0],ts[1],lim,laser_pos,laser_power,v_laser)/timeStep
          elif (len(ts) == 0):
            source = 0.0
          else:
            raise NameError('Bad number of solutions to laser overlap cell ' + str(c))
        else:
          source = 0
        sourceArray[c] = source 

        ##Set conductivity based on previous temperature
        solidCond = condFunction(tempArray[c])
        condArray[c] = condMeltFunction(solidCond,meltFracArray[c])
      
      tmodel.advance(maxIterPerTimeStep)

      ##Picard iterate out any temperature dependent conductivity
      picardIters = 1
      urf = 0.9
      while(picardIters < 200):
        for c in range(ncells):
          solidCond = condFunction(tempArray[c])
          if((meltFracArray[c] >= 1.0) or (tempArray[c] >= liquidTemp)):
            condArray[c] = (1.0-urf)*condArray[c] + urf*solidCond
          else:
            condArray[c] = (1.0-urf)*condArray[c] + urf*condMeltFunction(solidCond,
                                                                   max(meltFracArray[c],(tempArray[c]-solidTemp)/
                                                                       (liquidTemp-solidTemp)))
        picardIters = picardIters + 1
        prev_temp = np.copy(tempArray)
        tmodel.advance(maxIterPerTimeStep)
        err_norm = np.linalg.norm(prev_temp-tempArray)/np.linalg.norm(tempArray)
        print err_norm
        if (err_norm < 1e-4*urf):
          break
        if (picardIters == 20):
          print 'Picard Iterations did not converge after 20 iterations. Decreasing urf to 0.8'
          urf = 0.8
        if (picardIters == 60):
          print 'Picard Iterations did not converge after 60 iterations. Decreasing urf to 0.7'
          urf = 0.7
        if (picardIters == 120):
          print 'Picard Iterations did not converge after 120 iterations. Decreasing urf to 0.6'
          urf = 0.6
      else:
        raise NameError('Picard iterations did not converge after 200 iterations')
      ##Update time
      tmodel.updateTime()
        

      time = time + timeStep
      # Update laser position
      laser_pos.x = laser_pos_next.x
      laser_pos.y = laser_pos_next.y

      f.write(str(time) + ' ' + str(numTimeSteps) + ' ' + str(laser_pos.x) + ' ' + str(laser_pos.y) + '\n')
      # Output data for current timestep
      if((numTimeSteps % outFreq) == 0):
        filename = outfilebase + '_' + str(numTimeSteps) + '.vtk'
        writer = exporters.VTKWriterA(geomFields,meshes,filename,
                                      "TestTemperature",False,0)
        writer.init()
        writer.writeScalarField(thermalFields.temperature, "Temperature")
        writer.writeScalarField(thermalFields.meltFrac, "MeltFraction")
        #writer.writeScalarField(thermalFields.conductivity, "Conductivity")
        writer.finish()

      write_restart_file((ncells, temp, tempN1, enth, enthN1, enthInv, dHdT, meltFrac))
      numTimeSteps = numTimeSteps + 1

      if(int_time == 0.0):
        laser_power = laser_power_next
        line = laser_path.readline()
        if (line == ''):
          break
        else:
          words = line.split()
          if (len(words) == 0):
            break
          else:
            int_time = float(words[0])-time
            v_laser = [(float(words[1])-laser_pos.x)/int_time,(float(words[2])-laser_pos.y)/int_time]
            laser_power_next = float(words[3])
          
        
    f.close()                
    return_temp = np.copy(temp)
    return_tempN1 = np.copy(tempN1)
    return_enth = np.copy(enth)
    return_enthN1 = np.copy(enthN1)
    return_enthInv = np.copy(enthInv)
    return_dHdT = np.copy(dHdT)
    return_meltFrac = np.copy(meltFrac)
    
    return (ncells, return_temp, return_tempN1, return_enth, 
            return_enthN1, return_enthInv, return_dHdT, return_meltFrac)
####################################################################

raw_input("Press Enter to continue.")
#pdb.set_trace()
values_old = (0,[],[],[],[],[],[],[])
#Calculate roots and weights for 5th order Gauss Quardature integration
xs = LegendreRoots(5)
Ws = GaussLegendreWeights(5,xs)
numLayers = 0

if(os.path.exists('restart_file')):
  values_old = read_restart_file()

while numLayers < 1:
  values_old = RunSim(xlen,ylen,zlen,nx,ny,nz,values_old,'laser_path-' + str(numLayers))
  zlen = zlen + 0.0001
  nz = nz + 1
  numLayers = numLayers + 1
  curr_cycle = -1
