import sys
import numpy
import math
import random


class Errors:
    def __init__(self, etext):
        print "Error msg:", etext
        sys.exit(1)


def warrning(wtext):
    sys.stderr.write('Warrning: ' + wtext + '\n')


class mathTools:
    @staticmethod
    def getRMatrix(rVec, ra):
        rVec = numpy.array(rVec)
        rVec = rVec/numpy.linalg.norm(rVec)
        rmatrix = [[math.cos(ra) + rVec[0]*rVec[0]*(1-math.cos(ra)),
                    rVec[0]*rVec[1]*(1-math.cos(ra)) - rVec[2]*math.sin(ra),
                    rVec[0]*rVec[2]*(1-math.cos(ra)) + rVec[1]*math.sin(ra)],
                    [rVec[1]*rVec[0]*(1-math.cos(ra)) + rVec[2]*math.sin(ra),
                    math.cos(ra) + rVec[1]*rVec[1]*(1-math.cos(ra)),
                    rVec[1]*rVec[2]*(1-math.cos(ra)) - rVec[0]*math.sin(ra)],
                    [rVec[2]*rVec[0]*(1-math.cos(ra)) - rVec[1]*math.sin(ra),
                    rVec[2]*rVec[1]*(1-math.cos(ra)) + rVec[0]*math.sin(ra),
                    math.cos(ra) + rVec[2]*rVec[2]*(1-math.cos(ra))]]

        return numpy.array(rmatrix)

    @staticmethod
    def getRandomUS(numpyArr=True):
        angle = random.random()*math.pi*2
        u = random.uniform(-1, 1)
        u2 = u*u

        x = math.sqrt(1-u2) * math.cos(angle)
        y = math.sqrt(1-u2) * math.sin(angle)
        z = u

        if numpyArr:
            return numpy.array([x, y, z], dtype=float)
        else:
            return [x, y, z]


class basicShapes:
    @staticmethod
    def test():
        pass

    @staticmethod
    def sphere(r, x=0., y=0., z=0.):
        r2 = r*r
        return lambda xc, yc, zc: (xc - x)**2 + (yc - y)**2 + (zc - z)**2 < r2


class Lattice:
    def __init__(self, boxSize=None, system=None, origin=None, lvec=None):
        # Setting up default parameters
        if boxSize:
            self.setBoxSize(boxSize)
        else:
            self.setBoxSize(numpy.array([[-10., 10.], [-10., 10.], [-10., 10.]]))

        if system:
            self.setSystem(system)
        else:
            self.setSystem("pcc")

        if origin:
            self.setOrigin(origin)
        else:
            self.setOrigin(numpy.array([0., 0., 0.]))

        if lvec:
            self.setLatticeVec(lvec)
        else:
            self.setLatticeVec(numpy.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]))

    def getNodes(self):
        baseNodes = self._createNodes()
        latVec = self._latticeVec

        if (self._system in ["pcc", "simple", "P", "p"]):
            return baseNodes

        elif (self._system in ["bcc", "B", "b"]):
            self._origin += numpy.sum(self.latVec * 0.5, axis=0)
            exNodes = self._createNodes()
            return numpy.concatenate((baseNodes, exNodes), axis=0)

        elif (self._system in ["fcc", "F", "f"]):
            baseOrigin = numpy.copy(self._origin)
            self._origin += (latVec[0] + latVec[1]) * 0.5
            exNodes0 = self._createNodes()

            self._origin = numpy.copy(baseOrigin)
            self._origin += (latVec[1] + latVec[2]) * 0.5
            exNodes1 = self._createNodes()

            self._origin = numpy.copy(baseOrigin)
            self._origin += (latVec[2] + latVec[0]) * 0.5
            exNodes2 = self._createNodes()

            return numpy.concatenate((baseNodes, exNodes0, exNodes1, exNodes2), axis=0)

    def setMiller(self, mil, orient=[0., 0., 1.], angle=0.0):
        pointVec = numpy.array(orient)
        lVec = self._latticeVec
        sumLatVec = numpy.sum(lVec, axis=0)
        mVec = numpy.array([mil[0], mil[1], mil[2]]) * sumLatVec

        rVec = numpy.cross(mVec, pointVec)
        rVec = rVec / numpy.linalg.norm(rVec)

        ra = math.acos(numpy.dot(mVec, pointVec) /
                       (numpy.linalg.norm(mVec) * numpy.linalg.norm(pointVec)))

        rmatrix = self._getRMatrix(rVec, ra)
        rmatrixFin = self._getRMatrix(pointVec, angle)

        tmp0 = numpy.dot(rmatrixFin, numpy.dot(rmatrix, self._latticeVec[0]))
        tmp1 = numpy.dot(rmatrixFin, numpy.dot(rmatrix, self._latticeVec[1]))
        tmp2 = numpy.dot(rmatrixFin, numpy.dot(rmatrix, self._latticeVec[2]))

        self._latticeVec = numpy.array([tmp0, tmp1, tmp2])

    def setBoxSize(self, boxSize):
        xp = numpy.array([float(xp_) for xp_ in boxSize[0]])
        yp = numpy.array([float(yp_) for yp_ in boxSize[1]])
        zp = numpy.array([float(zp_) for zp_ in boxSize[2]])

        if (xp[0] >= xp[1]):
            Errors("Ill-defined X demension for box")
        if (yp[0] >= yp[1]):
            Errors("Ill-defined Y demension for box")
        if (zp[0] >= zp[1]):
            Errors("Ill-defined Z demension for box")

        self._boxSize = numpy.array([xp, yp, zp])

        bs = self._boxSize
        self._boxCorners = numpy.array([
            [bs[0][0], bs[1][0], bs[2][0]],  # 0
            [bs[0][0], bs[1][0], bs[2][1]],  # 1
            [bs[0][0], bs[1][1], bs[2][0]],  # 2
            [bs[0][0], bs[1][1], bs[2][1]],  # 3
            [bs[0][1], bs[1][0], bs[2][0]],  # 4
            [bs[0][1], bs[1][0], bs[2][1]],  # 5
            [bs[0][1], bs[1][1], bs[2][0]],  # 6
            [bs[0][1], bs[1][1], bs[2][1]],  # 7
        ])

    def setSystem(self, lType="pcc"):
        # simple
        if (lType in ["pcc", "simple", "P", "p"]):
            self._system = "pcc"

        # Body Centered Cubic
        elif (lType in ["bcc", "B", "b"]):
            self._system = "bcc"

        # Face Centered Cubic
        elif (lType in ["fcc", "F", "f"]):
            self._system = "fcc"

        else:
            Errors("unknown lattice type")

    def setOrigin(self, origin):
        origin = [float(origin[0]), float(origin[1]), float(origin[2])]
        self._origin = numpy.array(origin)

    def setLatticeVec(self, lVectors):
        al = numpy.array([float(ap_) for ap_ in lVectors[0]])
        bl = numpy.array([float(bp_) for bp_ in lVectors[1]])
        cl = numpy.array([float(cp_) for cp_ in lVectors[2]])

        if(numpy.dot(al, al) == 0):
            Errors("A vector has zero lenght")
        if(numpy.dot(bl, bl) == 0):
            Errors("B vector has zero lenght")
        if(numpy.dot(cl, cl) == 0):
            Errors("C vector has zero lenght")

        crAB = numpy.cross(al, bl)
        if(numpy.dot(crAB, crAB) == 0):
            Errors("A x B = 0, vectors may by parallel")

        if(numpy.dot(crAB, cl) == 0):
            Errors("C lies on AB plane")

        if(numpy.dot(crAB, cl) < 0):
            Errors("A B C vectros don't form right-handed system")

        self._latticeVec = numpy.array([al, bl, cl])

    def _createNodes(self):
        a = self._latticeVec[0]
        b = self._latticeVec[1]
        c = self._latticeVec[2]

        bs = self._boxSize
        bsVec = numpy.array([[bs[0][0],bs[1][0],bs[2][0]],
                            [bs[0][1],bs[1][1],bs[2][1]]])

        # find jumps and starting point
        pointOnPlane = numpy.empty([3, 3])
        pNormVec = numpy.empty([3, 3])
        jumps = numpy.empty([3], dtype=int)
        jumps[0], pNormVec[0], pointOnPlane[0] = self._findJumpsVec(b, c, a)
        jumps[1], pNormVec[1], pointOnPlane[1] = self._findJumpsVec(c, a, b)
        jumps[2], pNormVec[2], pointOnPlane[2] = self._findJumpsVec(a, b, c)

        dPlaneFactor = numpy.empty([3])
        for i in xrange(3):
            dPlaneFactor[i] = (pNormVec[i][0]*pointOnPlane[i][0] + pNormVec[i][1]*pointOnPlane[i][1] + pNormVec[i][2]*pointOnPlane[i][2])




        #stPartLen = numpy.array([numpy.dot(partVec_, partVec_) for partVec_ in stPartVec])

        startVec = numpy.linalg.solve(pNormVec, dPlaneFactor)

        sizeMax = (int(jumps[0])+1)*(int(jumps[1])+1)*(int(jumps[2])+1)
        nodes = numpy.empty([sizeMax, 3])
        nodeId = 0

        cLvl = numpy.copy(startVec)
        for ci in xrange(int(jumps[2])+1):
            bLvl = numpy.copy(cLvl)
            for bi in xrange(int(jumps[1])+1):
                aLvl = numpy.copy(bLvl)
                for ai in xrange(int(jumps[0])+1):
                    if (aLvl > bsVec[0]).all() and (aLvl < bsVec[1]).all():
                        nodes[nodeId][0] = aLvl[0]
                        nodes[nodeId][1] = aLvl[1]
                        nodes[nodeId][2] = aLvl[2]
                        nodeId += 1
                    aLvl += a
                bLvl += b
            cLvl += c

        return nodes[:nodeId]

    def _findJumpsVec(self, v1, v2, v3):
        # 1: define plane vector
        pNormVec = numpy.cross(v1, v2)
        pNormVec = numpy.multiply(pNormVec, 1.0/numpy.linalg.norm(pNormVec))

        # 2: find max and min vectros
        corners = self._boxCorners

        minCr = numpy.dot(corners[0], pNormVec)
        maxCr = numpy.dot(corners[0], pNormVec)
        for crVec in corners:
            crLen = numpy.dot(crVec, pNormVec)
            if(crLen < minCr):
                minCr = crLen
            elif(crLen > maxCr):
                maxCr = crLen

        # calculate jumps and start point
        jumpSize = numpy.dot(pNormVec, v3)
        originShiftFull = numpy.dot(self._origin, pNormVec)
        originShift = jumpSize*math.modf(originShiftFull/jumpSize)[0]

        startLen = math.ceil((minCr + originShift)/jumpSize)*jumpSize + originShift
        jumps = math.floor((maxCr - startLen)/jumpSize)
        startVec = numpy.multiply(pNormVec, startLen)

        # return jumps, normal to plane (pNormVec), and point on plane (startVec)
        return (jumps, pNormVec, startVec)

    def _getRMatrix(self, rVec, ra):
        rVec = numpy.array(rVec)
        rVec = rVec/numpy.linalg.norm(rVec)
        rmatrix = [[math.cos(ra) + rVec[0]*rVec[0]*(1-math.cos(ra)),
                    rVec[0]*rVec[1]*(1-math.cos(ra)) - rVec[2]*math.sin(ra),
                    rVec[0]*rVec[2]*(1-math.cos(ra)) + rVec[1]*math.sin(ra)],
                   [rVec[1]*rVec[0]*(1-math.cos(ra)) + rVec[2]*math.sin(ra),
                    math.cos(ra) + rVec[1]*rVec[1]*(1-math.cos(ra)),
                    rVec[1]*rVec[2]*(1-math.cos(ra)) - rVec[0]*math.sin(ra)],
                   [rVec[2]*rVec[0]*(1-math.cos(ra)) - rVec[1]*math.sin(ra),
                    rVec[2]*rVec[1]*(1-math.cos(ra)) + rVec[0]*math.sin(ra),
                    math.cos(ra) + rVec[2]*rVec[2]*(1-math.cos(ra))]]

        return numpy.array(rmatrix)


class Box:
    def __init__(self, nodes, enableAll=False):
        self.nodes = numpy.copy(nodes)
        nodesLen = len(nodes)
        self.nodesLen = nodesLen

        if enableAll:
            self.activNodes = numpy.ones(nodesLen, dtype=bool)
        else:
            self.activNodes = numpy.zeros(nodesLen, dtype=bool)

    def add(self, fun):
        nodesLen = self.nodesLen
        nodes = self.nodes
        activNodes = self.activNodes
        for i in xrange(nodesLen):
            node = nodes[i]
            if fun(node[0], node[1], node[2]):
                activNodes[i] = True

    def subtract(self, fun):
        nodesLen = self.nodesLen
        nodes = self.nodes
        activNodes = self.activNodes
        for i in xrange(nodesLen):
            node = nodes[i]
            if fun(node[0], node[1], node[2]):
                activNodes[i] = False

    def getActiveNodes(self):
        return self.getNodes()

    def getNodes(self):
        nodesLen = self.nodesLen
        nodes = self.nodes
        activNodes = self.activNodes

        # True is equal to 1, so sum is equal to all trues in table
        outNodesLen = numpy.sum(activNodes)
        outNodes = numpy.empty([outNodesLen, 3])

        crOutNode = 0
        for i in xrange(nodesLen):
            if activNodes[i]:
                outNodes[crOutNode][0] = nodes[i][0]
                outNodes[crOutNode][1] = nodes[i][1]
                outNodes[crOutNode][2] = nodes[i][2]
                crOutNode += 1

        return outNodes


class Molecule:
    def __init__(self, inName, format="xyz"):
        self.comment = "Empty"
        if format == "xyz":
            self._readXYZ(inName)
        elif format == "atom":
            self._initSingleAtom(inName)
        else:
            self._formatError("Unknow format!")

    def move(self, mv):
        atomCoords = self.atomCoords
        mv = numpy.copy([float(el) for el in mv])

        atomsLen = len(atomCoords)
        for i in xrange(atomsLen):
            atomCoords[i] += mv

    # to do center!!!
    def center(self, atomId=0, center=None):
        atomCoords = self.atomCoords
        if center is None:
            refAtom = -1. * numpy.copy(atomCoords[atomId])
        else:
            refAtom = numpy.array(center, dtype=float)

        atomsLen = len(atomCoords)
        for i in xrange(atomsLen):
            atomCoords[i] += refAtom

    def rotate(self, rVec, ra, atomId=0, center=None):
        atomCoords = self.atomCoords
        rmatrix = mathTools.getRMatrix(rVec, ra)
        if center is None:
            rPoint = numpy.copy(atomCoords[atomId])
        else:
            rPoint = numpy.array(center, dtype=float)

        atomsLen = len(atomCoords)
        for i in xrange(atomsLen):
            atomPos = numpy.copy(atomCoords[i])
            atomPos -= rPoint
            atomPos = numpy.dot(rmatrix, atomPos)
            atomPos += rPoint
            atomCoords[i] = atomPos

    def getAtoms(self):
        return (self.atomCoords, self.atomTypes)

    def _formatError(self, msg):
        pass

    def _initSingleAtom(self, atom):
        self.atomCoords = numpy.zeros([1,3], dtype=float)
        self.atomTypes = [str(atom)]

    def _readXYZ(self, inName):
        inFile = open(inName, 'r')

        line = inFile.readline()
        atomsInFile = int(line)
        self.comment = inFile.readline()

        atomCoords = numpy.empty([atomsInFile, 3], dtype=float)
        atomTypes = [None]*atomsInFile

        for i in xrange(atomsInFile):
            line = inFile.readline()
            sp = line.split()
            atomTypes[i] = sp[0]
            atomCoords[i][0] = float(sp[1])
            atomCoords[i][1] = float(sp[2])
            atomCoords[i][2] = float(sp[3])

        self.atomCoords = atomCoords
        self.atomTypes = atomTypes

    def printIt(self):
        print len(self.atomCoords)
        print "kkk"
        for i in xrange(len(self.atomCoords)):
            print self.atomTypes[i], self.atomCoords[i][0], self.atomCoords[i][1], self.atomCoords[i][2]


class Crystal:
    def __init__(self, nodes, molecule):
        self.nodes = nodes
        self.molecule = None
        self.setBaseMol(molecule)

    def setBaseMol(self, molecule):
        self.molecule = molecule

    def getAtoms(self):
        self._genAtoms()
        return (self.atomCoords, self.atomTypes)

    def _genAtoms(self):
        nodes = self.nodes
        nodesLen = len(nodes)

        molecule = self.molecule
        molCoords = molecule.atomCoords
        molType = molecule.atomTypes
        molLen = len(molecule.atomCoords)

        atomTotal = nodesLen * molLen

        atomCoords = numpy.empty([atomTotal, 3], dtype=float)
        atomTypes = [None] * atomTotal

        localId = 0
        for i in xrange(nodesLen):
            for j in xrange(molLen):
                atomTypes[localId] = molType[j]
                atomCoords[localId][0] = nodes[i][0] + molCoords[j][0]
                atomCoords[localId][1] = nodes[i][1] + molCoords[j][1]
                atomCoords[localId][2] = nodes[i][2] + molCoords[j][2]
                localId += 1

        self.atomCoords = atomCoords
        self.atomTypes = atomTypes

    def writeXYZ(self, outName):
        nodes = self.nodes
        nodesLen = len(nodes)

        molecule = self.molecule
        molCo = molecule.atomCoords
        molType = molecule.atomTypes
        molLen = len(molecule.atomCoords)

        outFile = open(outName, 'w')
        outFile.write(str(nodesLen*molLen) + '\n')
        outFile.write(molecule.comment.rstrip('\n') + '\n')

        for i in xrange(nodesLen):
            for j in xrange(molLen):
                xl = nodes[i][0] + molCo[j][0]
                yl = nodes[i][1] + molCo[j][1]
                zl = nodes[i][2] + molCo[j][2]

                line = molType[j] + " " + str(xl) + " " + str(yl) + " " + str(zl) + "\n"
                outFile.write(line)

            #line = atomTypes[i] + str(atomCoords[i][0]) + str(atomCoords[i][1]) + str(atomCoords[i][2])
            #outFile.write(line)

class System:
    def __init__(self):
        self.atomCoords = numpy.empty([0, 3], dtype=float)
        self.atomTypes = []
        self.comment = "no comment"

    def addCrystal(self, crystal):
        atomCoords = self.atomCoords
        atomTypes = self.atomTypes

        crystalAtomCoords, crystalAtomTypes = crystal.getAtoms()

        atomCoords = numpy.append(atomCoords, crystalAtomCoords, 0)
        atomTypes.extend(crystalAtomTypes)

        self.atomCoords = atomCoords
        self.atomTypes = atomTypes

    def addMolecule(self, molecule):
        atomCoords = self.atomCoords
        atomTypes = self.atomTypes

        atomCoords = numpy.append(atomCoords, molecule.atomCoords, 0)
        atomTypes.extend(molecule.atomTypes)

        self.atomCoords = atomCoords
        self.atomTypes = atomTypes

    def addAtoms(self, atomCoords, atomTypes):
        self.atomCoords = numpy.append(self.atomCoords, atomCoords, 0)
        self.atomTypes.extend(atomTypes)

    def writeXYZ(self, outName):
        atomTypes = self.atomTypes
        atomCoords = self.atomCoords

        size = len(atomTypes)

        outFile = open(outName, 'w')
        outFile.write(str(size) + "\n")
        outFile.write(self.comment + "\n")

        for i in xrange(size):
            coord = atomCoords[i]
            line = atomTypes[i] + " " + str(coord[0]) + " " + str(coord[1]) + " " + str(coord[2]) + "\n"

            outFile.write(line)

    def _findBoxSize(self):
        return  [[-10, 10], [-10, 10], [-10, 10]]

    def writeLammps(self, outName, tDic):
        atomTypes = self.atomTypes
        atomCoords = self.atomCoords
        size = len(atomTypes)

        # Write head
        outFile = open(outName, 'w')
        outFile.write("# " + str(self.comment) + "\n\n")
        outFile.write(" " + str(size) + "  atoms\n")
        outFile.write(" " + str(len(tDic)) + "  atom types\n\n")

        boxSize = self._findBoxSize()
        outFile.write(" " + str(boxSize[0][0]) + " " + str(boxSize[0][1]) + "  xlo xhi\n")
        outFile.write(" " + str(boxSize[1][0]) + " " + str(boxSize[1][1]) + "  ylo yhi\n")
        outFile.write(" " + str(boxSize[2][0]) + " " + str(boxSize[2][1]) + "  zlo zhi\n\n")

        outFile.write("Masses\n\n")
        for key in xrange(len(tDic)+1):
            outFile.write(" " + str(key) + " \n")

        outFile.write("\nAtoms\n")
        for i in xrange(size):
            coord = atomCoords[i]
            lID = i+1
            line = "\n" + str(lID) + " " + str(tDic[atomTypes[i]]) + " 0 " + str(coord[0]) + " " + str(coord[1]) + " " + str(coord[2])
            outFile.write(line)
