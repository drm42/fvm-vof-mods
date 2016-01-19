import numpy as np
import fvm.fvmbaseExt as fvmbaseExt
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
import pdb

class Grid:
  def __init__(self,xmax,ymax,zmax,imax,jmax,kmax):
    #set up extents in grid space and coordinate space
    self.xmax = xmax
    self.ymax = ymax
    self.zmax = zmax

    self.imax = imax
    self.jmax = jmax
    self.kmax = kmax

    #set up node coordinates
    self.xcoord = np.linspace(0.0, self.xmax, num=(self.imax+1))
    self.ycoord = np.linspace(0.0, self.ymax, num=(self.jmax+1))
    self.zcoord = np.linspace(0.0, self.zmax, num=(self.kmax+1))

    #set up node and cell arrays
    self.nNodes = (self.imax+1)*(self.jmax+1)*(self.kmax+1)
    self.nodes = [[[None for k in range(self.kmax+1)] 
                   for j in range(self.jmax+1)] 
                  for i in range(self.imax+1)]

    self.nCells = self.imax*self.jmax*self.kmax
    self.cells = [[[None for k in range(self.kmax)] 
                   for j in range(self.jmax)] 
                  for i in range(self.imax)]

    #set up face arrays
    self.xfaces = [[[None for k in range(self.kmax)] 
                   for j in range(self.jmax)] 
                  for i in range(self.imax+1)]
    self.yfaces = [[[None for k in range(self.kmax)] 
                   for j in range(self.jmax+1)] 
                  for i in range(self.imax)]
    self.zfaces = [[[None for k in range(self.kmax+1)] 
                   for j in range(self.jmax)] 
                  for i in range(self.imax)]


class Index: 
  def __init__(self,i,j,k):
    self.i = i
    self.j = j
    self.k = k

class Cell:
  def __init__(self,level,index,cID=-1,isGhost=False):
    self.level = level
    self.index = index
    self.faces = []
    self.cID = cID
    self.isGhost = isGhost
    self.parent = None
    self.children = []

class Face:
  def __init__(self,level,kind,index,nodes,fID=-1,isBoundary=False):
    self.level = level
    self.index = index
    self.kind = kind
    self.cell0 = None
    self.cell1 = None
    self.nodes = nodes
    self.fID = fID
    self.faceList = None
    self.isBoundary = isBoundary

class Node:
  def __init__(self,indicies,x,y,z,nID=-1):
    self.indicies = indicies
    self.nID = nID
    self.x = x
    self.y = y
    self.z = z

class Mesh:
  def __init__(self,grid):
    #add input grid as level 0 (coarsest) grid
    self.grids = [grid]
    #start mesh node count
    self.nNodes = 0
    self.nodes = []
    #setup nodes in the grid and node array
    for k in range(grid.kmax+1):
        for j in range(grid.jmax+1):
            for i in range(grid.imax+1):
                node = Node([Index(i,j,k)], grid.xcoord[i], grid.ycoord[j], 
                            grid.zcoord[k], self.nNodes)
                self.grids[0].nodes[i][j][k] = node
                self.nodes.append(node)
                self.nNodes = self.nNodes + 1
    #set cells ids in the grid cell id array
    self.nCells = 0
    self.cells = []
    for k in range(grid.kmax):
        for j in range(grid.jmax):
            for i in range(grid.imax):
                cell = Cell(level=0,index=Index(i,j,k),
                            cID=self.nCells)
                self.grids[0].cells[i][j][k] = cell
                self.cells.append(cell)
                self.nCells = self.nCells + 1

      
    self.interiorFaces = []
    self.nFaces = 0
    # interior x faces
    for i in range(1,grid.imax):
        for j in range(grid.jmax):
            for k in range(grid.kmax):
                self.makeXface(0,i,j,k)

    # interior y faces
    for j in range(1,grid.jmax):
        for i in range(grid.imax):
            for k in range(grid.kmax):
                self.makeYface(0,i,j,k)

    # interior z faces
    for k in range(1,grid.kmax):
        for i in range(grid.imax):
            for j in range(grid.jmax):
                self.makeZface(0,i,j,k)

    self.BoundaryCells = []
    #x-boundaries
    self.xminFaces = []
    i = 0
    for j in range(grid.jmax):
        for k in range(grid.kmax):
            self.makeXface(0,i,j,k)
    self.xmaxFaces = []
    i = grid.imax
    for j in range(grid.jmax):
        for k in range(grid.kmax):
            self.makeXface(0,i,j,k)

    #y-boundaries
    self.yminFaces = []
    j = 0
    for i in range(grid.imax):
        for k in range(grid.kmax):
            self.makeYface(0,i,j,k)
    self.ymaxFaces = []
    j = grid.jmax
    for i in range(grid.imax):
        for k in range(grid.kmax):
            self.makeYface(0,i,j,k)

    #z-boundaries
    self.zminFaces = []
    k = 0
    for i in range(grid.imax):
        for j in range(grid.jmax):
            self.makeZface(0,i,j,k)
    self.zmaxFaces = []
    k = grid.kmax
    for i in range(grid.imax):
        for j in range(grid.jmax):
            self.makeZface(0,i,j,k)

#######################################################################################
  def numberBoundaryCells(self):
    nb = 0
    for c in self.BoundaryCells:
      c.cID = self.nCells + nb
      nb = nb + 1
#######################################################################################
  def countFaces(self):
    self.nFaces = (len(self.interiorFaces) + len(self.xminFaces) + len(self.xmaxFaces) +
                   len(self.yminFaces) + len(self.ymaxFaces) + len(self.zminFaces) +
                   len(self.zmaxFaces))
#######################################################################################
  def deactivateFace(self,face):
    if(not face.faceList):
      raise NameError('Active Face is not in a Face List')
    face.faceList.remove(face)
    face.faceList = None
    if(face.isBoundary):
      self.BoundaryCells.remove(face.cell1)
    face.cell0.faces.remove(face)
    face.cell1.faces.remove(face)
    f.cell0 = None
    f.cell1 = None
#######################################################################################
  def makeXface(self,level,i,j,k):
    #xmin face
    if(i == 0):
      if(self.grids[level].xfaces[i][j][k] != None):
        #face already exists, get it
        face = self.grids[level].xfaces[i][j][k]
        if(face.faceList):
          raise NameError('Face to be created is already in a Face List')
      else:
        #face does not exist, make it
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i][j+1][k]
        node2 = self.grids[level].nodes[i][j+1][k+1]
        node3 = self.grids[level].nodes[i][j][k+1]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind='x',index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3],isBoundary=True)
        self.grids[level].xfaces[i][j][k] = face
      if(self.grids[level].cells[i][j][k] == None):
        raise NameError('Cell needed to create face does not exist on desired level')
      #activate face
      self.xminFaces.append(face)
      face.faceList = self.xminFaces
      face.cell0 = self.grids[level].cells[i][j][k]
      boundaryCell = Cell(level,Index(i-1,j,k),isGhost=True)
      self.BoundaryCells.append(boundaryCell)
      face.cell1 = boundaryCell
      face.cell0.faces.append(face)
      face.cell1.faces.append(face)
      return face

    #xmax face
    elif(i == self.grids[level].imax):
      if(self.grids[level].xfaces[i][j][k] != None):
        #face already exists, get it
        face = self.grids[level].xfaces[i][j][k]
        if(face.faceList):
          raise NameError('Face to be created is already in a Face List')
      else:
        #face does not exist, make it
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i][j][k+1]
        node2 = self.grids[level].nodes[i][j+1][k+1]
        node3 = self.grids[level].nodes[i][j+1][k]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind='x',index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3],isBoundary=True)
        self.grids[level].xfaces[i][j][k] = face
      if(self.grids[level].cells[i-1][j][k] == None):
        raise NameError('Cell needed to create face does not exist on desired level')
      #activate face
      self.xmaxFaces.append(face)
      face.faceList = self.xmaxFaces
      face.cell0 = self.grids[level].cells[i-1][j][k]
      boundaryCell = Cell(level,Index(i,j,k),isGhost=True)
      self.BoundaryCells.append(boundaryCell)
      face.cell1 = boundaryCell
      face.cell0.faces.append(face)
      face.cell1.faces.append(face)
      return face

    #interior face
    else:
      if(self.grids[level].xfaces[i][j][k] != None):
        #face already exists, get it
        face = self.grids[level].xfaces[i][j][k]
        if(face.faceList):
          raise NameError('Face to be created is already in a Face List')
      else:
        #face does not exist, make it
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i][j+1][k]
        node2 = self.grids[level].nodes[i][j+1][k+1]
        node3 = self.grids[level].nodes[i][j][k+1]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind='x',index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        self.grids[level].xfaces[i][j][k] = face
      if(self.grids[level].cells[i][j][k] == None or
         self.grids[level].cells[i-1][j][k] == None):
        raise NameError('Cell needed to create face does not exist on desired level')
      self.interiorFaces.append(face)
      face.faceList = self.interiorFaces
      face.cell0 = self.grids[level].cells[i][j][k]
      face.cell1 = self.grids[level].cells[i-1][j][k]
      face.cell0.faces.append(face)
      face.cell1.faces.append(face)
      return face
#######################################################################################
  def makeYface(self,level,i,j,k):
    #ymin face
    if(j == 0):
      if(self.grids[level].yfaces[i][j][k] != None):
        #face already exists, get it
        face = self.grids[level].yfaces[i][j][k]
        if(face.faceList):
          raise NameError('Face to be created is already in a Face List')
      else:
        #face does not exist, make it
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i][j][k+1]
        node2 = self.grids[level].nodes[i+1][j][k+1]
        node3 = self.grids[level].nodes[i+1][j][k]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind='y',index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3],isBoundary=True)
        self.grids[level].yfaces[i][j][k] = face
      if(self.grids[level].cells[i][j][k] == None):
        raise NameError('Cell needed to create face does not exist on desired level')
      self.yminFaces.append(face)
      face.faceList = self.yminFaces
      face.cell0 = self.grids[level].cells[i][j][k]
      boundaryCell = Cell(level,Index(i,j-1,k),isGhost=True)
      self.BoundaryCells.append(boundaryCell)
      face.cell1 = boundaryCell
      face.cell0.faces.append(face)
      face.cell1.faces.append(face)
      return face

    #ymax face
    elif(j == self.grids[level].jmax):
      if(self.grids[level].yfaces[i][j][k] != None):
        #face already exists, get it
        face = self.grids[level].yfaces[i][j][k]
        if(face.faceList):
          raise NameError('Face to be created is already in a Face List')
      else:
        #face does not exist, make it
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i+1][j][k]
        node2 = self.grids[level].nodes[i+1][j][k+1]
        node3 = self.grids[level].nodes[i][j][k+1]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind='y',index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3],isBoundary=True)
        self.grids[level].yfaces[i][j][k] = face
      if(self.grids[level].cells[i][j-1][k] == None):
        raise NameError('Cell needed to create face does not exist on desired level')
      self.ymaxFaces.append(face)
      face.faceList = self.ymaxFaces
      face.cell0 = self.grids[level].cells[i][j-1][k]
      boundaryCell = Cell(level,Index(i,j,k),isGhost=True)
      self.BoundaryCells.append(boundaryCell)
      face.cell1 = boundaryCell
      face.cell0.faces.append(face)
      face.cell1.faces.append(face)
      return face
    
    #interior face
    else:
      if(self.grids[level].yfaces[i][j][k] != None):
        #face already exists, get it
        face = self.grids[level].yfaces[i][j][k]
        if(face.faceList):
          raise NameError('Face to be created is already in a Face List')
      else:
        #face does not exist, make it
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i+1][j][k]
        node2 = self.grids[level].nodes[i+1][j][k+1]
        node3 = self.grids[level].nodes[i][j][k+1]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind='y',index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        self.grids[level].yfaces[i][j][k] = face
      if(self.grids[level].cells[i][j-1][k] == None or
         self.grids[level].cells[i][j][k] == None):
        raise NameError('Cell needed to create face does not exist on desired level')
      self.interiorFaces.append(face)
      face.faceList = self.interiorFaces
      face.cell0 = self.grids[level].cells[i][j-1][k]
      face.cell1 = self.grids[level].cells[i][j][k]
      face.cell0.faces.append(face)
      face.cell1.faces.append(face)
      return face
#######################################################################################
  def makeZface(self,level,i,j,k):
    #zmin face
    if(k == 0):
      if(self.grids[level].zfaces[i][j][k] != None):
        #face already exists, get it
        face = self.grids[level].zfaces[i][j][k]
        if(face.faceList):
          raise NameError('Face to be created is already in a Face List')
      else:
        #face does not exist, make it
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i+1][j][k]
        node2 = self.grids[level].nodes[i+1][j+1][k]
        node3 = self.grids[level].nodes[i][j+1][k]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind='z',index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3],isBoundary=True)
        self.grids[level].zfaces[i][j][k] = face
      if(self.grids[level].cells[i][j][k] == None):
        raise NameError('Cell needed to create face does not exist on desired level')
      self.zminFaces.append(face)
      face.faceList = self.zminFaces
      face.cell0 = self.grids[level].cells[i][j][k]
      boundaryCell = Cell(level,Index(i,j,k-1),isGhost=True)
      self.BoundaryCells.append(boundaryCell)
      face.cell1 = boundaryCell
      face.cell0.faces.append(face)
      face.cell1.faces.append(face)
      return face

    #zmax face
    elif(k == self.grids[level].kmax):
      if(self.grids[level].zfaces[i][j][k] != None):
        #face already exists, get it
        face = self.grids[level].zfaces[i][j][k]
        if(face.faceList):
          raise NameError('Face to be created is already in a Face List')
      else:
        #face does not exist, make it
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i][j+1][k]
        node2 = self.grids[level].nodes[i+1][j+1][k]
        node3 = self.grids[level].nodes[i+1][j][k]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind='z',index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3],isBoundary=True)
        self.grids[level].zfaces[i][j][k] = face
      if(self.grids[level].cells[i][j][k-1] == None):
        raise NameError('Cell needed to create face does not exist on desired level')
      self.zmaxFaces.append(face)
      face.faceList = self.zmaxFaces
      face.cell0 = self.grids[level].cells[i][j][k-1]
      boundaryCell = Cell(level,Index(i,j,k),isGhost=True)
      self.BoundaryCells.append(boundaryCell)
      face.cell1 = boundaryCell
      face.cell0.faces.append(face)
      face.cell1.faces.append(face)
      return face
    
    #interior face
    else:
      if(self.grids[level].zfaces[i][j][k] != None):
        #face already exists, get it
        face = self.grids[level].zfaces[i][j][k]
        if(face.faceList):
          raise NameError('Face to be created is already in a Face List')
      else:
        #face does not exist, make it
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i+1][j][k]
        node2 = self.grids[level].nodes[i+1][j+1][k]
        node3 = self.grids[level].nodes[i][j+1][k]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind='z',index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        self.grids[level].zfaces[i][j][k] = face
      if(self.grids[level].cells[i][j][k] == None or
         self.grids[level].cells[i][j][k-1] == None):
        raise NameError('Cell needed to create face does not exist on desired level')
      self.interiorFaces.append(face)
      face.faceList = self.interiorFaces
      face.cell0 = self.grids[level].cells[i][j][k]
      face.cell1 = self.grids[level].cells[i][j][k-1]
      face.cell0.faces.append(face)
      face.cell1.faces.append(face)
      return face
#######################################################################################
  def makeMesh(self):
    nodeCoordsN = models.newVec3Array(self.nNodes)
    nodeCoordsA = nodeCoordsN.asNumPyArray()

    for n in range(self.nNodes):
      nodeCoordsA[n,0] = self.nodes[n].x
      nodeCoordsA[n,1] = self.nodes[n].y
      nodeCoordsA[n,2] = self.nodes[n].z

    nFaceZones = 7
    faceGroupCountN = fvmbaseExt.newIntArray(nFaceZones)
    faceGroupCount = faceGroupCountN.asNumPyArray()

    #interior faces
    faceGroupCount[0] = len(self.interiorFaces)
    #xmin faces
    faceGroupCount[1] = len(self.xminFaces)
    #xmax faces
    faceGroupCount[2] = len(self.xmaxFaces)
    #ymin faces
    faceGroupCount[3] = len(self.yminFaces)
    #ymax faces
    faceGroupCount[4] = len(self.ymaxFaces)
    #zmin faces
    faceGroupCount[5] = len(self.zminFaces)
    #zmax faces
    faceGroupCount[6] = len(self.zmaxFaces)
    
    ## allocate arrays for face nodes and cells
    faceNodeCountN = fvmbaseExt.newIntArray(self.nFaces)
    faceNodesN = fvmbaseExt.newIntArray(self.nFaces*4)
    faceCellsN = fvmbaseExt.newIntArray(self.nFaces*2)
    faceNodesA = faceNodesN.asNumPyArray()
    faceCellsA = faceCellsN.asNumPyArray()

    faceNodeCountA = faceNodeCountN.asNumPyArray()
    faceNodeCountA[:] = 4
    
    ## reshape for convenience
    faceNodes = faceNodesA.reshape((self.nFaces,4))
    faceCells = faceCellsA.reshape((self.nFaces,2))

    nf = 0
    for f in self.interiorFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      f.fID = nf
      nf = nf + 1

    #number boundary cells for use with boundary faces
    self.numberBoundaryCells()

    for f in self.xminFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      f.fID = nf
      nf = nf + 1

    for f in self.xmaxFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      f.fID = nf
      nf = nf + 1

    for f in self.yminFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      f.fID = nf
      nf = nf + 1

    for f in self.ymaxFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      f.fID = nf
      nf = nf + 1

    for f in self.zminFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      f.fID = nf
      nf = nf + 1

    for f in self.zmaxFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      f.fID = nf
      nf = nf + 1
    return fvmbaseExt.Mesh(3, self.nCells, nodeCoordsN, faceCellsN, 
                           faceNodesN, faceNodeCountN, faceGroupCountN)
    
#######################################################################################
  def createNewGridLevel(self):
    currLevel = len(self.grids) - 1
    currGrid = self.grids[currLevel]
    newGrid = Grid(currGrid.xmax,currGrid.ymax,currGrid.zmax,2*currGrid.imax,
                   2*currGrid.jmax,2*currGrid.kmax)
    self.grids.append(newGrid)
    #copy current level nodes down to new grid
    for k in range(currGrid.kmax+1):
      for j in range(currGrid.jmax+1):
        for i in range(currGrid.imax+1):
          node = currGrid.nodes[i][j][k]
          if(node != None):
            if(newGrid.nodes[2*i][2*j][2*k] != None):
              raise NameError('Node already exists on new grid')
            newGrid.nodes[2*i][2*j][2*k] = node
            node.indicies.append(Index(2*i,2*j,2*k))
    #copy current level cells down to new grid
    for k in range(currGrid.kmax):
      for j in range(currGrid.jmax):
        for i in range(currGrid.imax):
          cell = currGrid.cells[i][j][k]
          newGrid.cells[2*i][2*j][2*k] = cell
          newGrid.cells[2*i][2*j][2*k+1] = cell
          newGrid.cells[2*i][2*j+1][2*k] = cell
          newGrid.cells[2*i][2*j+1][2*k+1] = cell
          newGrid.cells[2*i+1][2*j][2*k] = cell
          newGrid.cells[2*i+1][2*j][2*k+1] = cell
          newGrid.cells[2*i+1][2*j+1][2*k] = cell
          newGrid.cells[2*i+1][2*j+1][2*k+1] = cell

#######################################################################################
  def refineCell(self,currCell):
    if(currCell.children):
      raise NameError('Cell is already refined')
    #get current cell level
    currLevel = currCell.level
    newLevel = currLevel + 1

    #check and see if new grid level exists
    if(len(self.grids) < newLevel + 1):
      #create it if not
      self.createNewGridLevel()

    #starting index on new level
    ni = currCell.index.i*2
    nj = currCell.index.j*2
    nk = currCell.index.k*2
    #create new nodes
    for i in range(3):
      for j in range(3):
        for k in range(3):
          #check if node already exists, if not, create it
          n = self.grids[newLevel].nodes[ni+i][nj+j][nk+k]
          #node does not exist
          if(n == None):
            node = Node([Index(ni+i,nj+j,nk+k)], self.grids[newLevel].xcoord[ni+i], 
                        self.grids[newLevel].ycoord[nj+j], 
                        self.grids[newLevel].zcoord[nk+k], self.nNodes)
            self.grids[newLevel].nodes[ni+i][nj+j][nk+k] = node
            self.nodes.append(node)
            self.nNodes = self.nNodes + 1
          #node exists, but is not active
          elif(n.nID == -1):
            self.nodes.append(node)
            self.nID = self.nNodes
            self.nNodes = self.nNodes + 1
            
    
    
    #create new cells
    for i in range(2):
      for j in range(2):
        for k in range(2):
          newCell = Cell(level=newLevel,index=Index(ni+i,nj+j,nk+k))
          self.cells.append(newCell)
          self.nCells = self.nCells + 1
          currCell.children.append(newCell)
          newCell.parent = currCell
          self.grids[newLevel].cells[ni+i][nj+j][nk+k] = newCell
    #remove old cell
    currCell.cID = -1
    self.cells.remove(currCell)
    self.nCells = self.nCells - 1

    #deep copy current cell faces
    cellFaces = []
    for f in currCell.faces:
      cellFaces.append(f)

    #deactivate cell faces
    for f in cellFaces:
      self.deactivateFace(f)

    #create/update faces
    for i in range(3):
      for j in range(2):
        for k in range(2):
          self.makeXface(newLevel,ni+i,nj+j,nk+k)
    
    for j in range(3):
      for i in range(2):
        for k in range(2):
          self.makeYface(newLevel,ni+i,nj+j,nk+k)

    for k in range(3):
      for i in range(2):
        for j in range(2):
          self.makeZface(newLevel,ni+i,nj+j,nk+k)
    self.countFaces()

###############################################################################
  def coarsenCell(self,currCell):
    if(not currCell.children):
      raise NameError('Cell has no children')
    for c in currCell.children:
      if(currCell.children):
        raise NameError('Cell has grandchildren')

    level = currCell.level
    fineLevel = level + 1

    #starting index of child cells
    ni = currCell.index.i*2
    nj = currCell.index.j*2
    nk = currCell.index.k*2

    #copy cell down to fine level grid
    for i in range(2):
      for j in range(2):
        for k in range(2):
          self.grids[fineLevel].cells[ni+i][nj+j][nk+k] = currCell

    #remove center node
    n = self.grids[fineLevel].nodes[ni+1][nj+1][nk+1]
    self.nodes.remove(n)
    self.nNodes = self.nNodes - 1
    self.nID = -1

    #clear interior faces
    i = 1
    for j in range(2):
      for k in range(2):
        f = self.grids[level].xfaces[i][j][k]
        self.deactivateFace(f)
    j = 1
    for i in range(2):
      for k in range(2):
        f = self.grids[level].yfaces[i][j][k]
        self.deactivateFace(f)
    k = 1
    for i in range(2):
      for j in range(2):
        f = self.grids[level].zfaces[i][j][k]
        self.deactivateFace(f)

    #clear child cells
    for c in currCell.children:
      c.faces = []
      c.cID = -1
      c.parent = None
      self.cells.remove(c)
      self.nCells = self.nCells - 1
    #activate parent cell
    self.cells.append(currCell)
    self.nCells = self.nCells + 1
    currCell.children = []

    #coarsen cell faces
    for i in range(2):
      self.coarsenFace(level, currCell.index.i + i, 
                       currCell.index.j, currCell.index.k, 'x')
      self.coarsenFace(level, currCell.index.i, 
                       currCell.index.j + i, currCell.index.k, 'y')
      self.coarsenFace(level, currCell.index.i, 
                       currCell.index.j, currCell.index.k + i, 'z')
###############################################################################
  def coarsenFace(level, i, j, k, kind):
    fineLevel = Level + 1
    ni = i*2
    nj = j*2
    nk = k*2
    if(kind == 'x'):
      if(i == 0 or i == self.grids[level].imax):
        #boundary face, remove extra nodes and activate the coarse face
        n = self.grids[fineLevel].nodes[ni][nj+1][nk+1]
        self.nodes.remove(n)
        self.nNodes = self.nNodes - 1
        self.nID = -1
        self.makeXface(level,i,j,k)
      else:
        #interior face, check if either neighbor cell is refined
        cell1 = self.grids[level].cells[i][j][k]
        cell0 = self.grids[level].cells[i-1][j][k]
        if(cell0.children or cell1.children):
          for j in range(2):
            for k in range(2):
              self.makeXface(fineLevel,ni,nj+j,nk+k)

###############################################################################
xmax = 1.0
ymax = 1.0
zmax = 1.0

imax = 2
jmax = 2
kmax = 2

coarseGrid = Grid(xmax,ymax,zmax,imax,jmax,kmax)
coarseMesh = Mesh(coarseGrid)
coarseMesh.refineCell(coarseMesh.grids[0].cells[0][0][0])
coarseMesh.refineCell(coarseMesh.grids[0].cells[1][1][1])
#coarseMesh.refineCell(coarseMesh.grids[0].cells[0][0][1])
#coarseMesh.refineCell(coarseMesh.grids[0].cells[0][1][0])
coarseMesh.refineCell(coarseMesh.grids[0].cells[0][1][1])
#coarseMesh.refineCell(coarseMesh.grids[0].cells[1][0][0])
#coarseMesh.refineCell(coarseMesh.grids[0].cells[1][0][1])
#coarseMesh.refineCell(coarseMesh.grids[0].cells[1][1][0])
coarseMesh.refineCell(coarseMesh.grids[1].cells[0][0][0])
coarseMesh.refineCell(coarseMesh.grids[1].cells[0][3][3])
coarseMesh.refineCell(coarseMesh.grids[2].cells[0][0][0])
fvmMesh = coarseMesh.makeMesh()

meshes = (fvmMesh,)

geomFields = models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()
thermalFields = models.ThermalFields('temperature') 
tmodel = models.ThermalModelA(geomFields,thermalFields,meshes) 

bcmap = tmodel.getBCMap()
vcmap = tmodel.getVCMap()

vcZone = vcmap[fvmMesh.getID()]
vcZone['thermalConductivity'] = 1.0
vcZone['density'] = 1.0
vcZone['specificHeat'] = 1.0

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

bc_zmax.bcType = 'SpecifiedHeatFlux'
bc_zmax.setVar('specifiedHeatFlux',0)

bc_zmin.bcType = 'SpecifiedHeatFlux'
bc_zmin.setVar('specifiedHeatFlux',0)

bc_ymin.bcType = 'SpecifiedHeatFlux'
bc_ymin.setVar('specifiedHeatFlux',0)

bc_ymax.bcType = 'SpecifiedHeatFlux'
bc_ymax.setVar('specifiedHeatFlux',0)

bc_xmin.bcType = 'SpecifiedTemperature'
bc_xmin.setVar('specifiedTemperature',300.0)

bc_xmax.bcType = 'SpecifiedTemperature'
bc_xmax.setVar('specifiedTemperature',400.0)

# Set solver options
toptions = tmodel.getOptions()
toptions.transient = False
toptions.enthalpyModel = False
toptions.polynomialCp = False
#toptions['latentHeat'] = latentHeat
#toptions['solidTemp'] = solidTemp
#toptions['liquidTemp'] = liquidTemp
#toptions.cpFile = specificHeatFile
toptions.setVar('initialTemperature', 298.0)
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
tmodel.advance(100)

writer = exporters.VTKWriterA(geomFields,meshes,"adaptive_test.vtk",
                              "TestTemperature",False,0)
writer.init()
writer.writeScalarField(thermalFields.temperature, "Temperature")
#writer.writeScalarField(thermalFields.meltFrac, "MeltFraction")
#writer.writeScalarField(thermalFields.conductivity, "Conductivity")
writer.finish()
