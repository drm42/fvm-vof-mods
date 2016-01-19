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
    self.nfaces = 0
    self.faces = []
    self.cID = cID
    self.isGhost = isGhost
    self.parent = None
    self.children = []

class Face:
  def __init__(self,level,kind,index,nodes,fID=-1):
    self.level = level
    self.index = index
    self.kind = kind
    self.cell0 = None
    self.cell1 = None
    self.nodes = nodes
    self.fID = fID
    self.faceGroup = -1
    self.parent = None
    self.children = []

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
                face = Face(level=0,kind=0,index=Index(i,j,k),
                            nodes=[self.grids[0].nodes[i][j][k],
                                   self.grids[0].nodes[i][j+1][k],
                                   self.grids[0].nodes[i][j+1][k+1],
                                   self.grids[0].nodes[i][j][k+1]],
                            fID=self.nFaces)
                self.interiorFaces.append(face)

                face.cell0 = self.grids[0].cells[i][j][k]
                face.cell1 = self.grids[0].cells[i-1][j][k]
                face.cell0.faces.append(face)
                face.cell1.faces.append(face)

                self.nFaces = self.nFaces + 1

    # interior y faces
    for j in range(1,grid.jmax):
        for i in range(grid.imax):
            for k in range(grid.kmax):
                face = Face(level=0,kind=1,index=Index(i,j,k),
                            nodes=[self.grids[0].nodes[i][j][k],
                                   self.grids[0].nodes[i+1][j][k],
                                   self.grids[0].nodes[i+1][j][k+1],
                                   self.grids[0].nodes[i][j][k+1]], 
                            fID=self.nFaces)
                self.interiorFaces.append(face)

                face.cell0 = self.grids[0].cells[i][j-1][k]
                face.cell1 = self.grids[0].cells[i][j][k]
                face.cell0.faces.append(face)
                face.cell1.faces.append(face)

                self.nFaces = self.nFaces + 1

    # interior z faces
    for k in range(1,kmax):
        for i in range(imax):
            for j in range(jmax):
                face = Face(level=0,kind=2,index=Index(i,j,k),
                            nodes=[self.grids[0].nodes[i][j][k],
                                   self.grids[0].nodes[i+1][j][k],
                                   self.grids[0].nodes[i+1][j+1][k],
                                   self.grids[0].nodes[i][j+1][k]], 
                            fID=self.nFaces)
                self.interiorFaces.append(face)

                face.cell0 = self.grids[0].cells[i][j][k]
                face.cell1 = self.grids[0].cells[i][j][k-1]
                face.cell0.faces.append(face)
                face.cell1.faces.append(face)

                self.nFaces = self.nFaces + 1
    nb = 0
    self.BoundaryCells = []
    #x-boundaries
    self.xminFaces = []
    i = 0
    for j in range(jmax):
        for k in range(kmax):
            face = Face(level=0,kind=0,index=Index(i,j,k),
                        nodes=[self.grids[0].nodes[i][j][k],
                               self.grids[0].nodes[i][j+1][k],
                               self.grids[0].nodes[i][j+1][k+1],
                               self.grids[0].nodes[i][j][k+1]],
                        fID=self.nFaces)
            self.xminFaces.append(face) 

            face.cell0 = self.grids[0].cells[i][j][k]
            self.BoundaryCells.append(Cell(0,Index(i-1,j,k),self.nCells + nb,True))
            face.cell1 = self.BoundaryCells[nb]
            face.cell0.faces.append(face)
            face.cell1.faces.append(face)

            self.nFaces = self.nFaces + 1
            nb = nb + 1
    self.xmaxFaces = []
    i = imax
    for j in range(jmax):
        for k in range(kmax):
            face = Face(level=0,kind=0,index=Index(i,j,k),
                        nodes=[self.grids[0].nodes[i][j][k],
                               self.grids[0].nodes[i][j][k+1],
                               self.grids[0].nodes[i][j+1][k+1],
                               self.grids[0].nodes[i][j+1][k]],
                        fID=self.nFaces)
            self.xmaxFaces.append(face)

            face.cell0 = self.grids[0].cells[i-1][j][k]
            self.BoundaryCells.append(Cell(0,Index(i,j,k),self.nCells + nb,True))
            face.cell1 = self.BoundaryCells[nb]
            face.cell0.faces.append(face)
            face.cell1.faces.append(face)

            self.nFaces = self.nFaces + 1
            nb = nb + 1

    #y-boundaries
    self.yminFaces = []
    j = 0
    for i in range(imax):
        for k in range(kmax):
            face = Face(level=0,kind=1,index=Index(i,j,k),
                        nodes=[self.grids[0].nodes[i][j][k],
                               self.grids[0].nodes[i][j][k+1],
                               self.grids[0].nodes[i+1][j][k+1],
                               self.grids[0].nodes[i+1][j][k]],
                        fID=self.nFaces)
            self.yminFaces.append(face)

            face.cell0 = self.grids[0].cells[i][j][k]
            self.BoundaryCells.append(Cell(0,Index(i,j-1,k),self.nCells + nb,True))
            face.cell1 = self.BoundaryCells[nb]
            face.cell0.faces.append(face)
            face.cell1.faces.append(face)

            self.nFaces = self.nFaces + 1
            nb = nb + 1
            
    self.ymaxFaces = []
    j = jmax
    for i in range(imax):
        for k in range(kmax):
            face = Face(level=0,kind=1,index=Index(i,j,k),
                        nodes=[self.grids[0].nodes[i][j][k],
                               self.grids[0].nodes[i+1][j][k],
                               self.grids[0].nodes[i+1][j][k+1],
                               self.grids[0].nodes[i][j][k+1]],
                        fID=self.nFaces)
            self.ymaxFaces.append(face)

            face.cell0 = self.grids[0].cells[i][j-1][k]
            self.BoundaryCells.append(Cell(0,Index(i,j,k),self.nCells + nb,True))
            face.cell1 = self.BoundaryCells[nb]
            face.cell0.faces.append(face)
            face.cell1.faces.append(face)

            self.nFaces = self.nFaces + 1
            nb = nb + 1

    #z-boundaries
    self.zminFaces = []
    k = 0
    for i in range(imax):
        for j in range(jmax):
            face = Face(level=0,kind=2,index=Index(i,j,k),
                        nodes=[self.grids[0].nodes[i][j][k],
                               self.grids[0].nodes[i+1][j][k],
                               self.grids[0].nodes[i+1][j+1][k],
                               self.grids[0].nodes[i][j+1][k]],
                        fID=self.nFaces)
            self.zminFaces.append(face) 
            
            face.cell0 = self.grids[0].cells[i][j][k]
            self.BoundaryCells.append(Cell(0,Index(i,j,k-1),self.nCells + nb,True))
            face.cell1 = self.BoundaryCells[nb]
            face.cell0.faces.append(face)
            face.cell1.faces.append(face)

            self.nFaces = self.nFaces + 1
            nb = nb + 1

    self.zmaxFaces = []
    k = kmax
    for i in range(imax):
        for j in range(jmax):
            face = Face(level=0,kind=2,index=Index(i,j,k),
                        nodes=[self.grids[0].nodes[i][j][k],
                               self.grids[0].nodes[i][j+1][k],
                               self.grids[0].nodes[i+1][j+1][k],
                               self.grids[0].nodes[i+1][j][k]],
                        fID=self.nFaces)
            self.zmaxFaces.append(face)

            face.cell0 = self.grids[0].cells[i][j][k-1]
            self.BoundaryCells.append(Cell(0,Index(i,j,k),self.nCells + nb,True))
            face.cell1 = self.BoundaryCells[nb]
            face.cell0.faces.append(face)
            face.cell1.faces.append(face)

            self.nFaces = self.nFaces + 1
            nb = nb + 1
#######################################################################################
  def makeXface(self,level,i,j,k):
    if(self.grids[level].xfaces[i][j][k] != None):
      #face already exists, get it
      face = self.grids[level].xfaces[i][j][k]
      if(i == 0 or i == self.grids[level].imax):
        #face already exists and is a boundary face, no updates possible
        return face
      #otherwise see if the cell relationships need to be updated
      if(self.grids[level].cells[i][j][k] != None and 
         face.cell0 is not self.grids[level].cells[i][j][k]):
        #disconnect face from current cell
        face.cell0.faces.remove(face)
        #attach face to new cell
        face.cell0 = self.grids[level].cells[i][j][k]
        face.cell0.faces.append(face)
      if(self.grids[level].cells[i-1][j][k] != None and
           face.cell1 is not self.grids[level].cells[i-1][j][k]):
        #disconnect face from current cell
        face.cell1.faces.remove(face)
        #attach face to new cell
        face.cell1 = self.grids[level].cells[i-1][j][k]
        face.cell1.faces.append(face)
      return face        
                
    else:
      if(i == 0):
        #xmin face
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i][j+1][k],
        node2 = self.grids[level].nodes[i][j+1][k+1],
        node3 = self.grids[level].nodes[i][j][k+1]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind=0,index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        if(self.grids[level].cells[i][j][k] == None):
          raise NameError('Cell needed to create face does not exist on desired level')
        self.xminFaces.append(face)
        face.cell0 = self.grids[level].cells[i][j][k]
        boundaryCell = Cell(level,Index(i-1,j,k),isGhost=True)
        self.BoundaryCells.append(boundaryCell)
        face.cell1 = boundaryCell
        self.grids[level].xfaces[i][j][k] = face
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)
        return face
      elif(i == self.grids[level].imax):
        #xmax face
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i][j][k+1]
        node2 = self.grids[level].nodes[i][j+1][k+1]
        node3 = self.grids[level].nodes[i][j+1][k]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind=0,index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        if(self.grids[level].cells[i-1][j][k] == None):
          raise NameError('Cell needed to create face does not exist on desired level')
        self.xmaxFaces.append(face)
        face.cell0 = self.grids[level].cells[i-1][j][k]
        boundaryCell = Cell(level,Index(i,j,k),isGhost=True)
        self.BoundaryCells.append(boundaryCell)
        face.cell1 = boundaryCell
        self.grids[level].xfaces[i][j][k] = face
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)
        return face
      else:
        #interior face
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i][j+1][k]
        node2 = self.grids[level].nodes[i][j+1][k+1]
        node3 = self.grids[level].nodes[i][j][k+1]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind=0,index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        if(self.grids[level].cells[i][j][k] == None or
           self.grids[level].cells[i-1][j][k] == None):
          raise NameError('Cell needed to create face does not exist on desired level')
        self.interiorFaces.append(face)
        face.cell0 = self.grids[level].cells[i][j][k]
        face.cell1 = self.grids[level].cells[i-1][j][k]
        self.grids[level].xfaces[i][j][k] = face
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)
        return face

#######################################################################################
  def makeYface(self,level,i,j,k):
    if(self.grids[level].yfaces[i][j][k] != None):
      #face already exists, get it
      face = self.grids[level].yfaces[i][j][k]
      if(j == 0 or j == self.grids[level].jmax):
        #face already exists and is a boundary face, no updates possible
        return face
      #otherwise see if the cell relationships need to be updated
      if(self.grids[level].cells[i][j-1][k] != None and 
         face.cell0 is not self.grids[level].cells[i][j-1][k]):
        #disconnect face from current cell
        face.cell0.faces.remove(face)
        #attach face to new cell
        face.cell0 = self.grids[level].cells[i][j-1][k]
        face.cell0.faces.append(face)
      if(self.grids[level].cells[i][j][k] != None and
           face.cell1 is not self.grids[level].cells[i][j][k]):
        #disconnect face from current cell
        face.cell1.faces.remove(face)
        #attach face to new cell
        face.cell1 = self.grids[level].cells[i][j][k]
        face.cell1.faces.append(face)
      return face        
                
    else:
      if(j == 0):
        #ymin face
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i][j][k+1],
        node2 = self.grids[level].nodes[i+1][j][k+1],
        node3 = self.grids[level].nodes[i+1][j][k]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind=1,index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        if(self.grids[level].cells[i][j][k] == None):
          raise NameError('Cell needed to create face does not exist on desired level')
        self.yminFaces.append(face)
        face.cell0 = self.grids[level].cells[i][j][k]
        boundaryCell = Cell(level,Index(i,j-1,k),isGhost=True)
        self.BoundaryCells.append(boundaryCell)
        face.cell1 = boundaryCell
        self.grids[level].yfaces[i][j][k] = face
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)
        return face
      elif(j == self.grids[level].jmax):
        #ymax face
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i+1][j][k]
        node2 = self.grids[level].nodes[i+1][j][k+1]
        node3 = self.grids[level].nodes[i][j][k+1]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind=1,index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        if(self.grids[level].cells[i][j-1][k] == None):
          raise NameError('Cell needed to create face does not exist on desired level')
        self.ymaxFaces.append(face)
        face.cell0 = self.grids[level].cells[i][j-1][k]
        boundaryCell = Cell(level,Index(i,j,k),isGhost=True)
        self.BoundaryCells.append(boundaryCell)
        face.cell1 = boundaryCell
        self.grids[level].yfaces[i][j][k] = face
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)
        return face
      else:
        #interior face
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i+1][j][k]
        node2 = self.grids[level].nodes[i+1][j][k+1]
        node3 = self.grids[level].nodes[i][j][k+1]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind=1,index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        if(self.grids[level].cells[i][j-1][k] == None or
           self.grids[level].cells[i][j][k] == None):
          raise NameError('Cell needed to create face does not exist on desired level')
        self.interiorFaces.append(face)
        face.cell0 = self.grids[level].cells[i][j-1][k]
        face.cell1 = self.grids[level].cells[i][j][k]
        self.grids[level].yfaces[i][j][k] = face
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)
        return face
    
#######################################################################################
  def makeZface(self,level,i,j,k):
    if(self.grids[level].zfaces[i][j][k] != None):
      #face already exists, get it
      face = self.grids[level].zfaces[i][j][k]
      if(k == 0 or k == self.grids[level].kmax):
        #face already exists and is a boundary face, no updates possible
        return face
      #otherwise see if the cell relationships need to be updated
      if(self.grids[level].cells[i][j][k] != None and 
         face.cell0 is not self.grids[level].cells[i][j][k]):
        #disconnect face from current cell
        face.cell0.faces.remove(face)
        #attach face to new cell
        face.cell0 = self.grids[level].cells[i][j][k]
        face.cell0.faces.append(face)
      if(self.grids[level].cells[i][j][k-1] != None and
           face.cell1 is not self.grids[level].cells[i][j][k-1]):
        #disconnect face from current cell
        face.cell1.faces.remove(face)
        #attach face to new cell
        face.cell1 = self.grids[level].cells[i][j][k-1]
        face.cell1.faces.append(face)
      return face        
                
    else:
      if(k == 0):
        #zmin face
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i+1][j][k],
        node2 = self.grids[level].nodes[i+1][j+1][k],
        node3 = self.grids[level].nodes[i][j+1][k]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind=2,index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        if(self.grids[level].cells[i][j][k] == None):
          raise NameError('Cell needed to create face does not exist on desired level')
        self.zminFaces.append(face)
        face.cell0 = self.grids[level].cells[i][j][k]
        boundaryCell = Cell(level,Index(i,j,k-1),isGhost=True)
        self.BoundaryCells.append(boundaryCell)
        face.cell1 = boundaryCell
        self.grids[level].zfaces[i][j][k] = face
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)
        return face
      elif(k == self.grids[level].kmax):
        #zmax face
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i][j+1][k]
        node2 = self.grids[level].nodes[i+1][j+1][k]
        node3 = self.grids[level].nodes[i+1][j][k]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind=2,index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        if(self.grids[level].cells[i][j][k-1] == None):
          raise NameError('Cell needed to create face does not exist on desired level')
        self.zmaxFaces.append(face)
        face.cell0 = self.grids[level].cells[i][j][k-1]
        boundaryCell = Cell(level,Index(i,j,k),isGhost=True)
        self.BoundaryCells.append(boundaryCell)
        face.cell1 = boundaryCell
        self.grids[level].zfaces[i][j][k] = face
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)
        return face
      else:
        #interior face
        node0 = self.grids[level].nodes[i][j][k]
        node1 = self.grids[level].nodes[i+1][j][k]
        node2 = self.grids[level].nodes[i+1][j+1][k]
        node3 = self.grids[level].nodes[i][j+1][k]
        if(node0 == None or node1 == None or 
           node2 == None or node3 == None):
          raise NameError('Node needed to create face does not exist on desired level')
        face = Face(level=level,kind=2,index=Index(i,j,k),
                    nodes=[node0,node1,node2,node3])
        if(self.grids[level].cells[i][j][k] == None or
           self.grids[level].cells[i][j][k-1] == None):
          raise NameError('Cell needed to create face does not exist on desired level')
        self.interiorFaces.append(face)
        face.cell0 = self.grids[level].cells[i][j][k]
        face.cell1 = self.grids[level].cells[i][j][k-1]
        self.grids[level].zfaces[i][j][k] = face
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

      nf = nf + 1

    for f in self.xminFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      nf = nf + 1

    for f in self.xmaxFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      nf = nf + 1

    for f in self.yminFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      nf = nf + 1

    for f in self.ymaxFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      nf = nf + 1

    for f in self.zminFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      nf = nf + 1

    for f in self.zmaxFaces:
      faceNodes[nf,0] = f.nodes[0].nID
      faceNodes[nf,1] = f.nodes[1].nID
      faceNodes[nf,2] = f.nodes[2].nID
      faceNodes[nf,3] = f.nodes[3].nID

      faceCells[nf,0] = f.cell0.cID
      faceCells[nf,1] = f.cell1.cID

      nf = nf + 1
    
    pdb.set_trace()
    return fvmbaseExt.Mesh(3, self.nCells, nodeCoordsN, faceCellsN, 
                           faceNodesN, faceNodeCountN, faceGroupCountN)
    
    


#######################################################################################
  def refineCell(self,cID):
    #get current cell
    currCell = self.cells[cID]
    #get current cell level
    currLevel = currCell.level
    newLevel = currLevel + 1

    #check and see if new grid level exists
    if(len(self.grids) < newLevel + 1):
      #create it if not
      newGrid = Grid(self.grids[currLevel].xmax, self.grids[currLevel].ymax,
                     self.grids[currLevel].zmax, 2*self.grids[currLevel].imax,
                     2*self.grids[currLevel].jmax, 2*self.grids[currLevel].kmax)
      self.grids.append(newGrid)
      #copy current level nodes down to new grid
      for k in range(self.grids[currLevel].kmax+1):
        for j in range(self.grids[currLevel].jmax+1):
          for i in range(self.grids[currLevel].imax+1):
            node = self.grids[currLevel].nodes[i][j][k]
            if(self.grids[newLevel].nodes[2*i][2*j][2*k] != None):
              raise NameError('Node already exists on new grid')
            self.grids[newLevel].nodes[2*i][2*j][2*k] = node
            node.indicies.append(Index(2*i,2*j,2*k))
      #copy current level cells down to new grid
      for k in range(self.grids[currLevel].kmax):
        for j in range(self.grids[currLevel].jmax):
          for i in range(self.grids[currLevel].imax):
            cell = self.grids[currLevel].cells[i][j][k]
            self.grids[newLevel].cells[2*i][2*j][2*k] = cell
            self.grids[newLevel].cells[2*i][2*j][2*k+1] = cell
            self.grids[newLevel].cells[2*i][2*j+1][2*k] = cell
            self.grids[newLevel].cells[2*i][2*j+1][2*k+1] = cell
            self.grids[newLevel].cells[2*i+1][2*j][2*k] = cell
            self.grids[newLevel].cells[2*i+1][2*j][2*k+1] = cell
            self.grids[newLevel].cells[2*i+1][2*j+1][2*k] = cell
            self.grids[newLevel].cells[2*i+1][2*j+1][2*k+1] = cell

    #starting index on new level
    ni = currCell.index.i
    nj = currCell.index.j
    nk = currCell.index.k
    #create new nodes
    for i in range(3):
      for j in range(3):
        for k in range(3):
          #check if node already exists, if not, create it
          if(self.grids[newLevel].nodes[ni+i][nj+j][nk+k] == None):
            node = Node([Index(ni+i,nj+j,nk+k)], self.grids[newLevel].xcoord[ni+i], 
                        self.grids[newLevel].ycoord[nj+j], 
                        self.grids[newLevel].zcoord[nk+k], self.nNodes)
            self.grids[newLevel].node[nI+i][nJ+j][nK+k] = node
            self.nodes.append(node)
            self.nNodes = self.nNodes + 1
    
    
    #create new cells
    for i in range(2):
      for j in range(2):
        for k in range(2):
          newCell = Cell(level=newLevel,index=Index(ni+i,nj+j,nk+k))
          if(i == 0 and j ==0 and k == 0):
            #reuse existing slot and id for first cell
            newCell.cID = currCell.cID
            currCell.cID = -1
            self.cells[cID] = newCell
          else:
            #create new id and append
            newCell.cID = self.nCells
            self.cells.append(newCell)
            self.nCells = self.nCells + 1
          currCell.children.append(newCell)
          newCell.parent = currCell
          self.grids[newLevel].cells[ni+i][nj+j][nk+k] = cell


    #create new interior faces
    # interior x faces
    i = 1
    for j in range(2):
      for k in range(2):
        face = Face(level=newLevel,kind=0,index=Index(ni+i,nj+j,nk+k),
                    nodes=[self.grids[newLevel].nodes[ni+i][nj+j][nk+k],
                           self.grids[newLevel].nodes[ni+i][nj+j+1][nk+k],
                           self.grids[newLevel].nodes[ni+i][nj+j+1][nk+k+1],
                           self.grids[newLevel].nodes[ni+i][nj+j][nk+k+1]])
        self.interiorFaces.append(face)

        face.cell0 = self.grids[newLevel].cells[ni+i][nj+j][nk+k]
        face.cell1 = self.grids[newLevel].cells[ni][nj+j][nk+k]
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)

    # interior y faces
    j = 1
    for i in range(2):
      for k in range(2):
        face = Face(level=newLevel,kind=1,index=Index(ni+i,nj+j,nk+k),
                    nodes=[self.grids[newLevel].nodes[ni+i][nj+j][nk+k],
                           self.grids[newLevel].nodes[ni+i+1][nj+j][nk+k],
                           self.grids[newLevel].nodes[ni+i+1][nj+j][nk+k+1],
                           self.grids[newLevel].nodes[ni+i][nj+j][nk+k+1]])
        self.interiorFaces.append(face)

        face.cell0 = self.grids[newLevel].cells[ni+i][nj][nk+k]
        face.cell1 = self.grids[newLevel].cells[ni+i][nj+j][nk+k]
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)

    # interior z faces
    k = 1
    for i in range(2):
      for j in range(2):
        face = Face(level=newLevel,kind=1,index=Index(ni+i,nj+j,nk+k),
                    nodes=[self.grids[newLevel].nodes[ni+i][nj+j][nk+k],
                           self.grids[newLevel].nodes[ni+i+1][nj+j][nk+k],
                           self.grids[newLevel].nodes[ni+i+1][nj+j+1][nk+k],
                           self.grids[newLevel].nodes[ni+i][nj+j+1][nk+k]])
        self.interiorFaces.append(face)

        face.cell0 = self.grids[newLevel].cells[ni+i][nj][nk+k]
        face.cell1 = self.grids[newLevel].cells[ni+i][nj+j][nk]
        face.cell0.faces.append(face)
        face.cell1.faces.append(face)

    #boundary faces
    cellFaces = currCell.faces
    for f in cellFaces:
      #if face is already at new level, just need to update cell references
      if(f.level == newLevel):
        #x-face
        if(f.kind == 0):
          if(f.index.i == ni):
            newFaceCell = self.grids[newLevel].cells[ni][f.index.j][f.index.k]
          else:
            newFaceCell = self.grids[newLevel].cells[ni+1][f.index.j][f.index.k]
        #y-face
        elif(f.kind == 1):
          if(f.index.j == nj):
            newFaceCell = self.grids[newLevel].cells[f.index.i][nj][f.index.k]
          else:
            newFaceCell = self.grids[newLevel].cells[f.index.i][nj+1][f.index.k]
        #z-face
        elif(f.kind == 2):
          if(f.index.k == nk):
            newFaceCell = self.grids[newLevel].cells[f.index.i][f.index.j][nk]
          else:
            newFaceCell = self.grids[newLevel].cells[f.index.i][f.index.j][nk+1]
        if(f.cell0 is currCell):
          f.cell0 = newFaceCell
        else:
          f.cell1 = newFaceCell
        newFaceCell.faces.append(f)
      #if face is at coarser level, need to create new faces
      elif(f.level < newLevel):
        #x-face
        if(f.kind == 0):
          #xmin boundary
          if(f.index.i == 0):
            self.xminFaces.remove(f)
            self.BoundaryCells.remove(f.cell1)
            f.cell1.faces = []
            for j in range(2):
              for k in range(2):
                face = Face(level=newLevel,kind=0,index=Index(0,j,k),
                            nodes=[self.grids[newLevel].nodes[0][nj+j][nk+k],
                                   self.grids[newLevel].nodes[0][nj+j+1][nk+k],
                                   self.grids[newLevel].nodes[0][nj+j+1][nk+k+1],
                                   self.grids[newLevel].nodes[0][nj+j][nk+k+1]])
                boundaryCell = Cell(newLevel,Index(-1,nj+j,nk+k),isGhost=True)
                self.xminFaces.append(face)
                face.cell0 = self.grids[newLevel].cells[0][nj+j][nk+k]
                self.BoundaryCells.append(boundaryCell)
                face.cell1 = boundaryCell
                face.cell0.faces.append(face)
                face.cell1.faces.append(face)
                face.parent = f
                f.children.append(face)
          #xmax boundary
          elif(f.index.i == self.grids[f.level].imax):
            self.xmaxFaces.remove(f)
            self.BoundaryCells.remove(f.cell1)
            f.cell1.faces = []
            i = self.grids[newLevel].imax
            for j in range(2):
              for k in range(2):
                face = Face(level=newLevel,kind=0,index=Index(i,j,k),
                            nodes=[self.grids[0].nodes[i][nj+j][nk+k],
                                   self.grids[0].nodes[i][nj+j][nk+k+1],
                                   self.grids[0].nodes[i][nj+j+1][nk+k+1],
                                   self.grids[0].nodes[i][nj+j+1][nk+k]])
                boundaryCell = Cell(newLevel,Index(i,nj+j,nk+k),isGhost=True)
                self.xmaxFaces.append(face)
                face.cell0 = self.grids[newLevel].cells[i-1][nj+j][nk+k]
                self.BoundaryCells.append(boundaryCell)
                face.cell1 = boundaryCell
                face.cell0.faces.append(face)
                face.cell1.faces.append(face)
                face.parent = f
                f.children.append(face)
            
        
        
    currCell.faces = []

    #xmin face
    nb = 0
    if(nI == 0):
      #face is boundary face
      for j in range(2):
        for k in range(2):
          newFaceNodes[1].append([self.grids[newLevel].nodeIndex[0,nJ+j,nK+k],
                                  self.grids[newLevel].nodeIndex[0,nJ+j+1,nK+k],
                                  self.grids[newLevel].nodeIndex[0,nJ+j+1,nK+k+1],
                                  self.grids[newLevel].nodeIndex[0,nJ+j,nK+k+1]])
          newFaceCells[1].append([self.grids[newLevel].cellIndex[0,nJ+j,nK+k],
                                  self.nCells + nb])
          nb = nb + 1
    elif(self.grids[newLevel].cellIndex[nI-1,nJ,nK] == -2):
      #face borders cell that is at coarser level
      coarseI = cI
      coarseJ = cJ
      coarseK = cK
      coarseLevel = currLevel
      while(self.grids[coarseLevel].cellIndex[coarseI-1,coarseJ,coarseK] == -2):
        coarseI = coarseI/2
        coarseJ = coarseJ/2
        coarseK = coarseK/2
        coarseLevel = coarseLevel - 1
    elif(self.grids[newLevel].cellIndex[nI-1,nJ,nK] == -1):
      #face borders cells that are at finer level
      x = 2
    else:
      #face borders cells that are at current level 
      x = 1
  
  def getCellsW(self,i,j,k,level):
    l = level
    li = i
    lj = j
    lk = k
    #cell is at coarser level
    if(self.grids[level].cellIndex[i-1,j,k] == -2):
      while(self.grids[l].cellIndex[li-1,lj,lk] == -2 and l > 0):
        l = l - 1
        li = li/2
        lj = lj/2
        lk = lk/2
      return None


###############################################################################

xmax = 1.0
ymax = 1.0
zmax = 1.0

imax = 2
jmax = 2
kmax = 2

coarseGrid = Grid(xmax,ymax,zmax,imax,jmax,kmax)
coarseMesh = Mesh(coarseGrid)
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
