
# Jacqueline Lewis
# moleculefit.py


# This file computes the general shape of the iridium complex formed
# from coordinating pyridines and bipyrides to an iridium core. The
# computed complex is printed to an xyz file for input into avogadro.


import math
from sympy import nsolve, Symbol
import mpmath
from Tkinter import *
import Tkinter, Tkconstants, tkFileDialog
import os
import subprocess


# File I/O functions, from 15-112
def readFile(path):
    with open(path, "rt") as f:
        return f.read()

def writeFile(path, contents):
    with open(path, "wt") as f:
        f.write(contents)

# Simple math equations
def eq(x,y):
    return abs(x-y) < 0.01 # high tolerance due to variance in atoms

def distance(p1,p2):
    x1,y1,z1 = p1
    x2,y2,z2 = p2
    return math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def vector(p1,p2): # returns direction of vector from p1 to p2
    x1,y1,z1 = p1
    x2,y2,z2 = p2
    return x2-x1,y2-y1,z2-z1

def crossProduct(p1,p2):
    x = p1[1]*p2[2]-p1[2]*p2[1]
    y = -(p1[0]*p2[2]-p1[2]*p2[0])
    z = p1[0]*p2[1]-p1[1]*p2[0]
    return(x,y,z)

def push(v,coeff,init): # travels down vector
    x0,y0,z0 = init
    x,y,z = v
    return coeff*x+x0,coeff*y+y0,coeff*z+z0


# This class defines the necessary information for an atom and its
# movement through 3D space.

class Atom(object): 

    def __init__(self, x, y, z, name):
        self.x = x
        self.y = y
        self.z = z
        self.p = (self.x,self.y,self.z)
        self.name = name

    def __eq__(self,other):
        return (isinstance(other, Atom) and self.x == other.x 
            and self.y == other.y and self.z == other.z 
                and self.name == other.name)

    def __repr__(self):
        return "Atom(%s,%f,%f,%f)" % (self.name,self.x,self.y,self.z)

    def atomDistance(self,other):
        return distance(self.p,other.p)

    # This function changes the position or name of the atom as defined.
    def update(self,p=None,name=None):
        if p != None: 
            self.x,self.y,self.z = p
            self.p = p
        if name != None: self.name = name

    # This function takes an atom not in the original plane and orients
    # it to be at the same position to the new plane as it was to the 
    # original plane
    def rotate(self,origPlane,plane,origPoints,tPoints):

        orig1,orig2 = origPoints
        p1,p2 = tPoints
        dist = origPlane.distToPlane(self.p)
        above = origPlane.above(self.p)
        
        # moves atom onto original plane
        if above: along = -1
        else: along = 1
        self.update(origPlane.travel(self.p,dist,along))
        
        # translates position to new plane
        self.update(plane.findPlace(p1,p2,
            distance(self.p,orig1),distance(self.p,orig2)))
        
        # shifts atom away from new plane to proper distance
        if above: along = 1
        else: along = -1
        self.update(plane.travel(self.p,dist,along))


# This class defines a grouping of atoms that make up a chemical 
# complex. Functions to change the position of the complex in 
# cartesian coordinates are provided. 

class Complex(object):

    def __init__(self,atoms):
        self.atoms = atoms
        self.size = len(atoms)

    def __repr__(self):
        return "Molecule("+str(self.atoms)+")"

    # This function removes the specified atom from the complex, should
    # it exist.
    def deleteAtom(self,atom):
        if atom not in self.atoms: pass
        for i in range(self.size):
            if self.atoms[i] == atom:
                index = i
                break
        # ensures there are no array access errors
        if index == self.size-1: 
            self.atoms = self.atoms[:-1]
            self.size -= 1
        else: 
            self.atoms = self.atoms[:index]+self.atoms[index+1:]
            self.size -= 1

    # This function finds the nitrogens in the molecule.
    def getNs(self): return getNitrogens(self.atoms)

    # This function finds all atoms of the given type in the complex.
    def getAtom(self,name):
        res = []
        for a in self.atoms:
            if a.name == name: res.append(a)
        return res

    # This function finds all bound neighbors of the atom within tolerance.
    def getNeighbors(self,atom,tol=1.6):
        if atom not in self.atoms: return False
        neighbors = []
        for a in self.atoms:
            if a == atom: continue
            # checks if atom is in bonding range (this number may vary)
            if atom.atomDistance(a) < tol: 
                neighbors.append(a)
        return neighbors

    # This function flips all atoms across a plane, creating a reflection.
    def flip(self,plane):
        for atom in self.atoms:
            dist = plane.distToPlane(atom.p)
            atom.update(plane.travel(atom.p,2*dist,1))


# Ligands are complexes containing coordinating locations marked by Ag/Au.

class Ligand(Complex):

    def __init__(self,atoms):
        self.atoms = atoms
        self.size = len(atoms)
        self.nitrogens = getNitrogens(self.atoms)
        self.pivots = self.getPivots()

    # This function finds the Coordinating position marked by the Au and Ag.
    def findCoord(self):
        elems = self.getAtom("Au")+self.getAtom("Ag")
        return elems[0],elems[1]

    # This function finds the points that can be used to make the pyridinal 
    # plane using the Au/Ag and a planar neighbor.
    def getPivots(self):
        one,two = self.findCoord()
        others = self.getNeighbors(one)
        # ensures that the last pivot chosen is not in line with Au/Ag
        if two.atomDistance(others[0]) > two.atomDistance(others[1]):
            three = others[1]
        else: three = others[0]
        return (one,two,three)

    # This function removes any H bound to the coordinating atoms.
    def deleteExtraH(self):
        for atom in self.findCoord():
            for neigh in self.getNeighbors(atom):
                if neigh.name == "H":
                    self.deleteAtom(neigh)

    # This function moves the complex to a specified location and plane.
    def moveComplex(self,p1,p2,c):

        # calculates starting and ending position
        plane = Plane((p1,p2,c),"points")
        orig1,orig2 = self.pivots[0].p,self.pivots[1].p
        origPlane = Plane((self.pivots[0].p,self.pivots[1].p,self.pivots[2].p),"points")
        p1,p2 = updatePoints((p1,p2),(self.pivots[0].p,self.pivots[1].p))
        upper = []
        
        # moves each atom to its equivalent position to the new plane
        for atom in self.atoms:
            if atom == self.pivots[0]: atom.update(p1)
            elif atom == self.pivots[1]: atom.update(p2)
            elif not origPlane.onPlane(atom.p):
                atom.rotate(origPlane,plane,(orig1,orig2),(p1,p2))
            else:
                orthPlane = origPlane.orthPlane(orig1,orig2)
                # atoms above the orthogonal plane sometimes are placed 
                # on the wrong side of the new orthogonal plane and must
                # be fixed
                if orthPlane.above(atom.p): 
                    dist = orthPlane.distToPlane(atom.p)
                    upper.append((atom,dist))
                atom.update(plane.findPlace(p1,p2,
                    distance(atom.p,orig1),distance(atom.p,orig2)))
        
        # determines if the complex sits in the correct position relative 
        # to the center
        newOrthPlane = plane.orthPlane(self.pivots[0].p,self.pivots[1].p)
        if distance(self.pivots[2].p,c) < distance(self.pivots[0].p,c):
            self.flip(newOrthPlane)
        
        # fixes the above plane atoms to be in the correct position
        for (atom,dist) in upper:
            atom.update(newOrthPlane.travel(atom.p,2*dist,1))
            if newOrthPlane.distToPlane(atom.p) > dist:
                atom.update(newOrthPlane.travel(atom.p,4*dist,-1))


# This class defines a plane and its uses in translating points. The
# plane can be defined in multiple ways depending upon the information
# available.

class Plane(object):

    def __init__(self,args,method):
        # uses 3 points to define the plane
        if method == "points": 
            p1,p2,p3 = args
            self.p = self.getPlane(p1,p2,p3)
            self.points = [p1,p2,p3]
            self.a,self.b,self.c,self.d = self.p
        # uses a point and the normal vector to define the plane
        if method == "vector":
            p,v = args
            self.a,self.b,self.c = v
            self.d = self.getD(p)
            self.points = [p]
            self.p = (self.a,self.b,self.c,self.d)
        self.normal = (self.a,self.b,self.c)
        self.magnitude = math.sqrt(self.a**2+self.b**2+self.c**2)

    def __repr__(self):
        return "(%f,%f,%f,%f)" % (self.a,self.b,self.c,self.d)

    # This function determines the distance form a point to a plane
    def distToPlane(self,p1): 
        x1,y1,z1 = p1
        x0,y0,z0 = self.points[0]
        num = abs(self.a*(x1-x0)+self.b*(y1-y0)+self.c*(z1-z0))
        den = self.magnitude
        return(num/den)

    # This function finds the D value in the plane equation
    # Ax + By + Cz + D = 0.
    def getD(self,p):
        x,y,z = p
        return -(self.a*x+self.b*y+self.c*z)

    # This function moves a point along a normal vector.
    def travel(self,p,dist,along):
        coeff = along*dist/self.magnitude
        return push(self.normal,coeff,p)

    # This function checks if a point is on a plane.
    def onPlane(self,p1):
        return eq(self.distToPlane(p1),0)

    # This function finds the location of a point on a new plane
    # based on its location on a different plane.
    def findPlace(self,p1,p2,dist1,dist2):
        a,b,c,d = self.p
        x1,y1,z1 = p1
        x2,y2,z2 = p2
        x,y,z = Symbol('x'),Symbol('y'),Symbol('z')
        # declares the functions defining the location of the point
        f1 = a*x+b*y+c*z+d # on plane 
        f2 = (x-x1)**2+(y-y1)**2+(z-z1)**2-dist1**2 # dist from pivot
        f3 = (x-x2)**2+(y-y2)**2+(z-z2)**2-dist2**2 # dist from pivot
        # solves the system of nonlinear equations
        x,y,z = nsolve((f1,f2,f3),(x,y,z),(1,1,1),tol=1*10**-7)
        return(round(x,5),round(y,5),round(z,5))

    # This function determines the equation of a plane from 3 points.
    def getPlane(self,p1,p2,p3):
        x1,y1,z1 = p1
        x2,y2,z2 = p2
        x3,y3,z3 = p3
        # finds the vectors between points
        AB = vector(p1,p2)
        AC = vector(p1,p3)
        a,b,c = crossProduct(AB,AC)
        # finds the d value for the plane eqn Ax + By + Cz + D = 0
        d = -(a*x1+b*y1+c*z1)
        return(a,b,c,d)

    # This function determines if the point is above the plane.
    def above(self,p1):
        x1,y1,z1 = p1
        return x1*self.a+y1*self.b+z1*self.c+self.d > 0

    # This function finds the angle between the point and the plane
    # from a specified center of the plane.
    def getAngle(self,p,c):
        x1,y1,z1 = p
        x2,y2,z2 = c
        v1,v2,v3 = (x2-x1,y2-y1,z2-z1)
        p1,p2,p3 = self.a,self.b,self.c
        magn1 = self.magnitude
        magn2 = math.sqrt(v1**2+v2**2+v3**2)
        angle = math.degrees(math.asin(abs(p1*v1+p2*v2+p3*v3)/(magn1*magn2)))
        return angle

    # This function creates a plane orthogonal to the given plane, 
    # from two points on both.
    def orthPlane(self,p1,p2):
        x1,y1,z1 = p1
        x2,y2,z2 = p2
        return Plane((p1,crossProduct(self.normal,
            (x1-x2,y1-y2,z1-z2))),"vector")


# This function updates the pivot points to be the correct distance
# from one another, as observed in the molecule.
def updatePoints(target,given):
    t1,t2 = target
    g1,g2 = given
    tDist = distance(t1,t2)
    gDist = distance(g1,g2)
    dif = (gDist-tDist)/2
    v = vector(t1,t2)
    v1,v2,v3 = v
    # moves the two points apart or together as needed
    if dif > 0:
        t1 = push(v,-dif/tDist,t1)
        t2 = push(v,dif/tDist,t2)
    else:
        t1 = push(v,dif/tDist,t1)
        t2 = push(v,-dif/tDist,t2)
    return (t1,t2)


# This function finds all the nitrogens in a group of atoms.
def getNitrogens(atoms):
    nitros = []
    for atom in atoms:
        if atom.name == "N":
            nitros.append(atom)
    return nitros


# This function creates a complex from a file containing the xyz 
# coordinates of the compositional atoms.
def makeComplex(path):
    raw = readFile(path).replace("\n"," ").replace("\r"," ").split(" ")
    raw = [x for x in raw if x != ""]
    numAtoms = int(raw[0])
    index = 1
    # if extraneous information is included, it is removed
    for i in range(len(raw)):
        if "Energy" in raw[i]:
            index = i+2
            break
    table = raw[index:]
    atoms = []
    # translates strings to atoms
    for i in range(numAtoms):
        name = table[i*4]
        x = float(table[i*4+1])
        y = float(table[i*4+2])
        z = float(table[i*4+3])
        atoms.append(Atom(x,y,z,name))
    return Complex(atoms)


# This function creates a complex that has the coordinating Ag/Au.
def makeLigand(mol):
    return Ligand(mol.atoms)


# This function creates a complex from the specified file and moves
# it to the correct position and orientation in space.
def compute(file,p1,p2,c,name):
    mol = makeLigand(makeComplex(file))
    mol.moveComplex(p1,p2,c)
    mol.deleteExtraH()
    mol.pivots[0].update(name="N")
    if name == "NN": mol.pivots[1].update(name="N")
    else: mol.pivots[1].update(name="C")
    return mol


# This function details the coordinates of the iridium center being used 
# for these computations. If a different file is used, the coordinates must
# be changed.
def iridiumCoordinates():

    coord = [0]*6
    # bpy binding site
    coord[0],coord[1] = (0.70272,1.29041,0.34899),(0.70272,-1.29041,-0.34899)
    # ppy binding sites (N first)
    coord[2],coord[3] = (-1.12624,-0.77394,1.91881),(-2.44763,-1.40246,-0.28744)
    coord[4],coord[5] = (-1.12624,0.77394,-1.91881),(-2.44763,1.40246,0.28744)
    return coord


# This function combines multiple complexes into the iridium complex
# formed through conjugation.
def joinComplex(mols,IrCoord):
    atoms = [Atom(IrCoord[0],IrCoord[1],IrCoord[2],"Ir")]
    # puts all atoms into one complex file
    for i in range(len(mols)): atoms += mols[i].atoms
    return atoms


####################################
# UI #
####################################

def init(data): 

    data.CNfiles = [""]*16
    data.NNfiles = [""]*16
    data.complete = False
    data.output = "output"
    data.header = getDefault()

    # editing tools
    data.folder = "output"
    data.hEdit = getDefault()
    data.editing = [False,False]
    data.pipe = [False,False]
    data.time = 0
    data.index = (0,0)

# This function opens a file dialogue and adds the chosen file to the 
# corresponding location in the interface.
def fileExplorer(data,lst,index):
    location = os.getcwd()
    name = tkFileDialog.askopenfilename(initialdir = location,
        title = "Select file",filetypes = (("XYZ","*.xyz"),("all files","*.*")))
    if name == (): name = ""
    lst[index-9] = name

# This function returns a list with no empty string elements.
def clearFluff(lst):
    return [x for x in lst if x != ""]

# This function determines the size specifications for the UI.
def getSizeSpecs(data):
    center = data.width/2
    bheight = data.height/26
    bwidth = data.width/3
    left = bwidth/3
    right = center+bwidth/6
    return [center,bheight,bwidth,left,right]

# This function reacts to mouse clicks and reacts if a button is pressed.
def mousePressed(event, data): 

    # size specs
    center,bheight,bwidth,left,right = getSizeSpecs(data)
    index = int((event.y-bheight)//bheight)

    # checks if the left column is clicked on
    if event.x > left and event.x < left+bwidth:
        if index > 8 and index < 22: 
            data.complete = False
            fileExplorer(data,data.CNfiles,index)
    # checks if the right column is clicked on
    if event.x > right and event.x < right+bwidth:
        if index > 8 and index < 22: 
            data.complete = False
            fileExplorer(data,data.NNfiles,index)
    # checks if the bottom button is clicked on
    if event.x > left+bwidth/2 and event.x < right+bwidth/2:
        if index == 23: 
            mix(data,"output")
            data.complete = True

    if event.x > left and event.x < right+bwidth:
        # checks if the top button is clicked on
        if index == 0:
            if not data.editing[1]: data.editing[0] = True
        # checks if the header area is clicked on
        if index > 1 and index < 7:
            if not data.editing[0]: 
                # clicks on the specific line
                data.editing[1] = True
                idx = int((event.y-bheight*3)//(bheight*5/10))
                lines = data.header.split("\n")
                if len(lines) <= idx: idx = len(lines)-1
                data.index = (idx,len(lines[idx]))

# This function checks if a folder exists.
def isValidFolder(folder):
    return os.path.isdir(folder)

# This function backspaces out the character before the cursor in the header.
def removeChar(data):
    # determines which line of the multiline header is being edited
    lines = data.hEdit.split("\n")
    idx1,idx2 = data.index
    line = lines[idx1]
    # ensures no array access errors, removes elem from end of string
    if idx2 >= len(line) and len(line) > 0: 
        line = line[:-1]
        data.index = (idx1,idx2-1)
        lines[idx1] = line
    # removes element from middle of string
    elif idx2 > 0: 
        line = line[:idx2-1]+line[idx2:]
        data.index = (idx1,idx2-1)
        lines[idx1] = line
    # deletes newline separation between line before and current
    elif idx1 != 0: 
        data.index = (idx1-1,len(data.hEdit.split("\n")[idx1-1]))
        newLine = lines[idx1-1]+lines[idx1]
        # joins the two lines into one
        if idx1+1 < len(lines):
            lines = lines[:idx1-1]+[newLine]+lines[idx1+1:]
        else: lines = lines[:idx1-1]+[newLine]
    return "\n".join(lines)

# This function adds a character to the header.
def addChar(data,event):
    lines = data.hEdit.split("\n")
    idx1,idx2 = data.index
    line = lines[idx1]
    # ensures no array access errors
    if idx2 >= len(line): line = line + event.char
    elif idx2 >= 0: line = line[:idx2]+event.char+line[idx2:]
    # ensures all newline removing gets all carriage returns
    lines[idx1] = line.replace("\r","\n")
    # adding a newline changes the line number
    if event.keysym == "Return": data.index = (idx1+1,0)
    else: data.index = (data.index[0],data.index[1]+1)
    return "\n".join(lines)

# This function ends the editing operation for header and folder.
def finishEditing(data,textType):
    # a folder must be made and the overall variable changed
    if textType == "folder":
        makeFolder(data.folder)
        data.output = data.folder
        index = 0
    # the index is reset and the overall variable is set as the edited
    elif textType == "header":
        data.header = data.hEdit
        data.index = (0,0)
        index = 1
    # the editing operation is over
    data.editing[index] = False
    data.pipe[index] = False

# This function moves the cursor horizontally in header editing mode.
def horiz(data,direction):
    lines = data.hEdit.split("\n")
    if direction == "Left":
        # moves cursor left with wraparound, if possible
        if data.index[1] != 0: data.index = (data.index[0],data.index[1]-1)
        elif data.index[0] != 0: 
            index = data.index[0]-1
            data.index = (index,len(lines[index]))
    elif direction == "Right":
        # moves cursor right with wraparound, if possible
        if data.index[1] != len(lines[data.index[0]]):
            data.index = (data.index[0],data.index[1]+1)
        elif data.index[0] < len(lines)-1: 
            index = data.index[0]+1
            data.index = (index,0)

# This function moves the cursor vertically in header editing mode.
def vert(data,direction):
    lines = data.hEdit.split("\n")
    if direction == "Up":
        # moves cursor up, if possible
        if data.index[0] != 0:
            idx = data.index[0]-1
            lineLength = len(lines[idx])
            if data.index[1] > lineLength: data.index = (idx,lineLength)
            else: data.index = (idx,data.index[1])
    elif direction == "Down":
        # moves cursor down, if possible
        if data.index[0] < len(lines)-1:
            idx = data.index[0]+1
            lineLength = len(lines[idx])
            if data.index[1] > lineLength: data.index = (idx,lineLength)
            else: data.index = (idx,data.index[1])

# This function reacts to user key strokes and responds accordingly.
def keyPressed(event, data): 
    # edits the folder name
    if data.editing[0]:
        if event.keysym == "Return": finishEditing(data,"folder")
        elif event.keysym == "BackSpace": data.folder = data.folder[:-1]
        else: data.folder += event.char
        # prevents excessive characters from ruining the formatting
        if len(data.folder) > 24: data.folder = data.folder[:-1]
    # edits the header
    elif data.editing[1]:
        if event.keysym == "Escape": finishEditing(data,"header")
        elif event.keysym == "BackSpace": data.hEdit = removeChar(data)
        elif event.keysym == "Left" or event.keysym == "Right":
            horiz(data,event.keysym)
        elif event.keysym == "Up" or event.keysym == "Down":
            vert(data,event.keysym)
        else: data.hEdit = addChar(data,event)

# This function maintains the timing of the system.
def timerFired(data):
    # time is used to blink cursor in editing mode
    if data.time % 5 == 0 and data.editing[0]:
        data.pipe[0] = not data.pipe[0]
    if data.time % 5 == 0 and data.editing[1]:
        data.pipe[1]= not data.pipe[1]
    data.time += 1

# This function returns the name of the file (not the whole path) without 
# its "." ending.
def trimFileName(fileName):
    return fileName.split("/")[-1].split(".")[0]

# This function adds the specified heading to the file and formats to 
# produce a gaussian-readable script.
def addHeader(file,contents,heading):
    # if the title is left alone in the UI, the ligand names go there
    pieces = heading.split("Title")
    title = trimFileName(file)
    res = ""
    for piece in pieces[:-1]:
        res += piece + title
    return res+pieces[-1]+"\n"+contents+"\r\r\r"

# This function produces the default heading of a gaussian-readable file.
def getDefault():
    line = []
    line.append("%NProcShared=3")
    line.append("%Chk=Title.chk")
    line.append(chr(37)+"mem=6GB")
    line.append("#n B3LYP/LANL2DZ Opt")
    line.append("")
    line.append(" Title")
    line.append("")
    line.append("1 1")
    return "\n".join(line)

# This function makes the folder as specified, if it does not exist.
def makeFolder(foldName):
    if not isValidFolder(foldName):
        subprocess.check_call(["mkdir",foldName])

# This function creates the output .xyz and .com files for a complex.
def makeOutput(data,atoms,file="new"):
    xyzStart = str(len(atoms)) + "\n\n"
    res = ""
    # collects all atom information
    for atom in atoms:
        res += "%-10s%-14f%-14f%f\n" % (atom.name,atom.x,atom.y,atom.z)
    comContents = addHeader(file,res,data.header)
    # places both files into their specific folder
    makeFolder(file)
    xyzFileName = file + "/" + trimFileName(file) + ".xyz"
    comFileName = file + "/" + trimFileName(file) + ".com"
    writeFile(xyzFileName,xyzStart+res)
    writeFile(comFileName,comContents)

# This function creates all possible combinations of NN and CN complexes 
# entered and produces the output files for each pairing.
def mix(data,folder="output"):
    # gets the specific coordinates of the iridium coordinating atoms
    coord,IrCoord = iridiumCoordinates(),(-0.94992,0,0)
    CNs,NNs = [],[]
    # runs through all inputs and computes new positions
    for file in clearFluff(data.NNfiles):
        NNs.append((file,compute(file,coord[0],coord[1],IrCoord,"NN")))
    for file in clearFluff(data.CNfiles):
        CNs.append((file,compute(file,coord[2],coord[3],IrCoord,"CN"),
            compute(file,coord[4],coord[5],IrCoord,"CN")))
    # creates output files for all pairings
    for (file1,CN1,CN2) in CNs:
        for (file2,NN) in NNs:
            path = folder+"/"+"Ir_CN"+trimFileName(file2)+"_2_NN"+trimFileName(file1)
            makeOutput(data,joinComplex((NN,CN1,CN2),IrCoord),path)

# This function draws the grid elements for the CN and NN inputs, as well as
# the compute button.
def drawTable(canvas, data):
    center,bheight,bwidth,left,right = getSizeSpecs(data)

    # draws the grid and text
    for i in range(14):
        canvas.create_rectangle(left, (i+9)*bheight, left+bwidth, (i+10)*bheight,fill="white")
        canvas.create_rectangle(right, (i+9)*bheight, right+bwidth, (i+10)*bheight,fill="white")
        if i == 0: continue
        canvas.create_text(left+bwidth/2,(i+9.5)*bheight,text=trimFileName(data.CNfiles[i-1]),font="Arial 15 bold")
        canvas.create_text(right+bwidth/2,(i+9.5)*bheight,text=trimFileName(data.NNfiles[i-1]),font="Arial 15 bold")
    canvas.create_text(left+bwidth/2,9.5*bheight,text="CN", font="Arial 20 bold")
    canvas.create_text(right+bwidth/2,9.5*bheight,text="NN",font="Arial 20 bold")
    # compute button
    canvas.create_rectangle(left+bwidth/2,24*bheight,right+bwidth/2,25*bheight,fill="white")
    canvas.create_text(center,24.5*bheight,text="Compute",font="Arial 20 bold")
    # computation is complete notification
    if data.complete: canvas.create_text(center,25.5*bheight,text="Complete",font="Arial 20 bold")
    
# This function draws the folder and header editing locations and progress.
def drawInput(canvas, data):
    center,bheight,bwidth,left,right = getSizeSpecs(data)

    text = data.hEdit
    # adds the blinking cursor to the correct space in the header
    if data.editing[1]: 
        canvas.create_text(center,2.5*bheight,text="Press escape to end edit",font="Arial 20 bold")
        lines = text.split("\n")
        line = lines[data.index[0]]
        if data.pipe[1]: pipe = "|"
        else: pipe = " "
        if data.index[1] < len(line):
            line = line[:data.index[1]]+pipe+line[data.index[1]:]
        else: line = line+pipe
        lines[data.index[0]] = line
        text = "\n".join(lines)

    canvas.create_rectangle(left,bheight*3,right+bwidth,bheight*8,fill="white")
    canvas.create_text(left+5,bheight*3,anchor="nw",text=text,font="Arial 10 bold")
    canvas.create_rectangle(left,1*bheight,right+bwidth,2*bheight,fill="white")
    # adds the blinking cursor to the end of the folder name
    if data.pipe[0]: pipe = "|"
    else: pipe = " "
    canvas.create_text(left+5,1.5*bheight,text="Output Folder Name: "+data.folder+pipe,font="Arial 15 bold",anchor="w")

# This function draws the UI every frame.
def redrawAll(canvas, data): 
    canvas.create_rectangle(0,0,data.width+5,data.height+5,fill="purple2")
    drawTable(canvas,data)
    drawInput(canvas,data)


####################################
# runUI function # from 15-112 #
####################################

def runUI(width=300, height=300):
    def redrawAllWrapper(canvas, data):
        canvas.delete(ALL)
        redrawAll(canvas, data)
        canvas.update()    

    def mousePressedWrapper(event, canvas, data):
        mousePressed(event, data)
        redrawAllWrapper(canvas, data)

    def keyPressedWrapper(event, canvas, data):
        keyPressed(event, data)
        redrawAllWrapper(canvas, data)

    def timerFiredWrapper(canvas, data):
        timerFired(data)
        redrawAllWrapper(canvas, data)
        # pause, then call timerFired again
        canvas.after(data.timerDelay, timerFiredWrapper, canvas, data)
    # Set up data and call init
    class Struct(object): pass
    data = Struct()
    data.width = width
    data.height = height
    data.timerDelay = 100 # milliseconds
    init(data)
    # create the root and the canvas
    root = Tk()
    canvas = Canvas(root, width=data.width, height=data.height)
    canvas.pack()
    # set up events
    root.bind("<Button-1>", lambda event:
                            mousePressedWrapper(event, canvas, data))
    root.bind("<Key>", lambda event:
                            keyPressedWrapper(event, canvas, data))
    timerFiredWrapper(canvas, data)
    # and launch the app
    root.mainloop()  # blocks until window is closed

runUI(600, 800)