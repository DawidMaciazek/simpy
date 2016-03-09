import numpy
import sys
import logging as log


log.basicConfig(format="%(levelname)s: %(message)s", level=10)

def getf(frame, name):
    fields = frame[0]['format']
    try:
        findex = fields.index(name)
        return frame[findex + 1]
    except ValueError:
        log.warning("Could not find field: " + str(name))
        return None


class ReadTraj:

    def __init__(self, inName, format="lammpstrj"):
        self.knownFormat = ["xyz", "lammpstrj", "slammpstrj"]
        self.inCounter = 0

        if format in self.knownFormat:
            self.format = format
        else:
            log.error("Unknown format: " + str(format))

        self.inName = inName

        # Open file for reading
        log.info("Opening file: " + str(inName))
        self.inFile = open(inName, 'r')

    def skipFrame(self):
        self.readFrame(skip=True)

    def readFrame(self, skip=False):
        self.inCounter += 1
        format = self.format
        if format == "xyz":
            log.info("Reading xyz frame: " + str(self.inCounter))
            return self._readFrameXYZ(skip)
        elif format == "lammpstrj":
            log.info("Reading lammpstrj frame: " + str(self.inCounter))
            return self._readLammpstrj(skip)
        elif format == "slammpstrj":
            log.info("Reading simple lammpstrj frame: " + str(self.inCounter))
            return self._readSimpleLammpstrj()

    # Read frames formats

    # Standart .xyz format:
    # (Element) (x) (y) (z)
    def _readFrameXYZ(self, skip=False):
        inFile = self.inFile
        addInfo = {}
        line = inFile.readline()
        if not line:
            return False

        frameAtoms = int(line)
        addInfo['atoms'] = int(frameAtoms)
        # !!! check if previously defined atoms is equal
        # frameComment = inFile.readline()
        line = inFile.readline()
        addInfo['comment'] = line
        addInfo['fields'] = ['elements', 'coords']

        # Allocate memory
        coords = numpy.empty([frameAtoms, 3], dtype=float)
        elements = []

        if skip:
            for i in xrange(frameAtoms):
                inFile.readline().split()
            return True

        for i in xrange(frameAtoms):
            sl = inFile.readline().split()
            elements.append(sl[0])
            coords[i][0] = sl[1]
            coords[i][1] = sl[2]
            coords[i][2] = sl[3]

        return [addInfo, elements, coords]

    # Simple lammps format (fast but limited use):
    # TO DO - retihink if it is even worth to keep this function...
    def _readSimpleLammpstrj(self):
        maxHead = 50
        inFile = self.inFile
        addInfo = {}

        line = inFile.readline()
        if not line:
            return False

        keyFieldName = line.split()[1]
        headLoops = 0
        while keyFieldName != "ATOMS":
            if keyFieldName == "TIMESTEP":
                timestep = int(inFile.readline())
                addInfo['timestep'] = int(timestep)

            elif keyFieldName == "NUMBER":
                frameAtoms = int(inFile.readline())
                addInfo['atoms'] = int(frameAtoms)

            elif keyFieldName == "BOX":
                boxSize = []
                sl = inFile.readline().split()
                boxSize.append([float(sl[0]), float(sl[1])])
                sl = inFile.readline().split()
                boxSize.append([float(sl[0]), float(sl[1])])
                sl = inFile.readline().split()
                boxSize.append([float(sl[0]), float(sl[1])])

                addInfo['box'] = boxSize

            elif keyFieldName == "TIME":
                realTime = float(inFile.readline().split()[0])
                addInfo['time'] = float(realTime)

            else:
                log.warn("Unknow keyFieldName: " + keyFieldName)

            line = inFile.readline()
            keyFieldName = line.split()[1]
            headLoops += 1
            if headLoops == maxHead:
                log.error("Max head iter (cant find ATOMS field)")


        # Header line
        # --- tmp pass that line ---

        # Alloc memory
        coords = numpy.zeros([frameAtoms, 3], dtype=float)
        elements = [None]*frameAtoms

        for i in xrange(frameAtoms):
            sp = inFile.readline().split()
            atomNum = int(sp[0]) - 1
            elements[atomNum] = sp[1]

            coords[atomNum][0] = sp[2]
            coords[atomNum][1] = sp[3]
            coords[atomNum][2] = sp[4]

        addInfo['format'] = ['elements', 'coords']
        return [addInfo, elements, coords]

    # Standard lammps dump format
    def _readLammpstrj(self, skip=False):
        maxHead = 50
        inFile = self.inFile
        addInfo = {}

        line = inFile.readline()
        if not line:
            return False

        keyFieldName = line.split()[1]
        headLoops = 0
        while keyFieldName != "ATOMS":
            if keyFieldName == "TIMESTEP":
                timestep = int(inFile.readline())
                addInfo['timestep'] = int(timestep)

            elif keyFieldName == "NUMBER":
                frameAtoms = int(inFile.readline())
                addInfo['atoms'] = int(frameAtoms)

            elif keyFieldName == "BOX":
                boxSize = []
                sl = inFile.readline().split()
                boxSize.append([float(sl[0]), float(sl[1])])
                sl = inFile.readline().split()
                boxSize.append([float(sl[0]), float(sl[1])])
                sl = inFile.readline().split()
                boxSize.append([float(sl[0]), float(sl[1])])

                addInfo['box'] = boxSize

            elif keyFieldName == "TIME":
                realTime = float(inFile.readline().split()[0])
                addInfo['time'] = float(realTime)

            else:
                log.warn("Unknow keyFieldName: " + keyFieldName)

            line = inFile.readline()
            keyFieldName = line.split()[1]
            headLoops += 1
            if headLoops == maxHead:
                log.error("Max head iter (cant find ATOMS field)")

        sp = line.split()
        keyFields = sp[2:]

        # Key field organize
        if 'x' in keyFields and 'y' in keyFields and 'z' in keyFields:
            xyzFlag = True
        else:
            xyzFlag = False

        if 'vx' in keyFields and 'vy' in keyFields and 'vz' in keyFields:
            vxyzFlag = True
        else:
            vxyzFlag = False

        addInfo['format'] = []
        # Convert index of input to index of output
        maxSize = 100
        keyProxy = [None]*maxSize
        # Pointer to list containging data
        fields = []
        keyIndex = 0
        proxyIndex = 0
        idField = -1
        coordIndex = maxSize
        velocityIndex = maxSize
        for key in keyFields:
            if xyzFlag and key in ['x', 'y', 'z']:
                if key == 'x':
                    fields.append(numpy.zeros([frameAtoms, 3], dtype=float))
                    addInfo['format'].append('coords')

                    coordIndex = keyIndex
                    keyProxy[proxyIndex] = [coordIndex, 0]
                    keyIndex += 1
                elif key == 'y':
                    keyProxy[proxyIndex] = [coordIndex, 1]
                elif key == 'z':
                    keyProxy[proxyIndex] = [coordIndex, 2]
                else:
                    log.error("Something strange happend while assigning fields")

            elif vxyzFlag and key in ['vx', 'vy', 'vz']:
                if key == 'vx':
                    fields.append(numpy.zeros([frameAtoms, 3], dtype=float))
                    addInfo['format'].append('velocity')

                    velocityIndex = keyIndex
                    keyProxy[proxyIndex] = [velocityIndex, 0]
                    keyIndex += 1
                elif key == 'vy':
                    keyProxy[proxyIndex] = [velocityIndex, 1]
                elif key == 'vz':
                    keyProxy[proxyIndex] = [velocityIndex, 2]
                else:
                    log.error("Something strange happend while assigning fields")

            elif key == 'type':
                fields.append(numpy.zeros(frameAtoms, dtype=int))
                addInfo['format'].append(key)
                keyProxy[proxyIndex] = keyIndex
                keyIndex += 1

            elif key == 'id':
                idField = keyIndex
            elif key == 'element':
                addInfo['format'].append(key)
                fields.append([None]*frameAtoms)
                keyProxy[proxyIndex] = keyIndex
                keyIndex += 1
            else:
                fields.append(numpy.zeros(frameAtoms, dtype=float))
                addInfo['format'].append(key)
                keyProxy[proxyIndex] = keyIndex
                keyIndex += 1

            proxyIndex += 1

        keyProxy = keyProxy[:proxyIndex]
        # Read body
        fNum = proxyIndex
        aId = -1

        # skipping part
        if skip:
            for i in xrange(frameAtoms):
                inFile.readline().split()
            return None

        sort = False;
        if sort:
            for i in xrange(frameAtoms):
                sp = inFile.readline().split()
                if idField != -1:
                    aId = int(sp[idField]) - 1
                else:
                    aId += 1

                for j in xrange(fNum):
                    # Proxy index for field
                    jprox = keyProxy[j]

                    if isinstance(jprox, list):
                        fields[jprox[0]][aId][jprox[1]] = sp[j]
                    elif jprox is not None:
                        fields[jprox][aId] = sp[j]

        else:
            aId = -1
            for i in xrange(frameAtoms):
                aId += 1
                sp = inFile.readline().split()

                for j in xrange(fNum):
                    # Proxy index for field
                    jprox = keyProxy[j]

                    if isinstance(jprox, list):
                        fields[jprox[0]][aId][jprox[1]] = sp[j]
                    elif jprox is not None:
                        fields[jprox][aId] = sp[j]


        output = [addInfo]
        output.extend(fields)

        return output

    knownFormat = ["xyz", "lammpsdat", "lammpsmol"]
    def writeFrame(frame, format='xyz'):
        pass

