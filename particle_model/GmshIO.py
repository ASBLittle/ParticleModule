""" Model containing class to read Gmsh files."""


class GmshMesh(object):
    """This is a class for storing nodes and elements.

    Members:
    nodes -- A dict of the form { nodeID: [ xcoord, ycoord, zcoord] }
    elements -- A dict of the form { elemID: (type, [tags], [nodeIDs]) }

    Methods:
    read(file) -- Parse a Gmsh version 1.0 or 2.0 mesh file
    write(file) -- Output a Gmsh version 2.0 mesh file
    """

    def __init__(self):
        self.nodes = {}
        self.elements = {}

    def read(self, filename):
        """Read a Gmsh .msh file.

        Reads Gmsh format 1.0 and 2.0 mesh files, storing the nodes and
        elements in the appropriate dicts.
        """
        import struct

        mshfile = open(filename, 'r')

        mode_dict = {'$NOD': 1, '$Nodes': 1,
                     '$ELM': 2,
                     '$Elements': 3,
                     '$MeshFormat': 4}

        readmode = 0
        line = 'a'
        while line:
            line = mshfile.readline()
            line = line.strip()
            if line.startswith('$'):
                readmode = mode_dict.get(line, 0)
            elif readmode:
                columns = line.split()
                if readmode == 4:
                    if len(columns) == 3:
                        vno, ftype, dsize = (float(columns[0]),
                                             int(columns[1]),
                                             int(columns[2]))
                        del vno
                    else:
                        endian = struct.unpack('i', columns[0])
                        del endian
                if readmode == 1:
                    # Version 1.0 or 2.0 Nodes
                    try:
                        if ftype == 0 and len(columns) == 4:
                            self.nodes[int(columns[0])] = [float(c) for c in columns[1:]]
                        elif ftype == 1:
                            nnods = int(columns[0])
                            for __ in range(nnods):
                                data = mshfile.read(4 + 3 * dsize)
                                values = struct.unpack('=i3d', data)
                                self.nodes[values[0]] = list(values[1:])
                            mshfile.read(1)
                    except ValueError:
                        readmode = 0
                elif ftype == 0 and readmode > 1 and len(columns) > 5:
                    # Version 1.0 or 2.0 Elements
                    try:
                        columns = [int(c) for c in columns]
                    except ValueError:
                        readmode = 0
                    else:
                        (ele_id, ele_type) = columns[0:2]
                        if readmode == 2:
                            # Version 1.0 Elements
                            tags = columns[2:4]
                            nodes = columns[5:]
                        else:
                            # Version 2.0 Elements
                            ntags = columns[2]
                            tags = columns[3:3 + ntags]
                            nodes = columns[3 + ntags:]
                        self.elements[ele_id] = (ele_type, tags, nodes)
                elif readmode == 3 and ftype == 1:
                    tdict = {1: 2, 2: 3, 3: 4, 4: 4, 5: 5, 6: 6, 7: 5, 8: 3, 9: 6, 10: 9, 11: 10}
                    try:
                        neles = int(columns[0])
                        k = 0
                        while k < neles:
                            ele_type, ntype, ntags = struct.unpack('=3i', mshfile.read(3 * 4))
                            k += ntype
                            for __ in range(ntype):
                                mysize = 1 + ntags + tdict[ele_type]
                                data = struct.unpack('=%di' % mysize,
                                                     mshfile.read(4 * mysize))
                                self.elements[data[0]] = (ele_type,
                                                          data[1:1 + ntags],
                                                          data[1 + ntags:])
                    except:
                        raise
                    mshfile.read(1)

        mshfile.close()

    def write(self, filename):
        """Dump the mesh out to a Gmsh 2.0 msh file."""

        with open(filename, 'w') as mshfile:

            mshfile.write('$MeshFormat\n2.0 0 8\n$EndMeshFormat\n')
            mshfile.write('$Nodes\n%d\n' % len(self.nodes))
            for node_id, coord in self.nodes.items():
                mshfile.write('%s %s %s\n' % (mshfile, node_id, ' '.join([str(c) for c in coord])))
            mshfile.write('$EndNodes\n')
            mshfile.write('$Elements\n%d\n' % len(self.elements))
            for ele_id, elem in self.elements.items():
                (ele_type, tags, nodes) = elem
                mshfile.write('%s %s %s ' % (ele_id, ele_type, len(tags)))
                mshfile.write('%s ' % ' '.join([str(c) for c in tags]))
                mshfile.write('%s\n' % ' '.join([str(c) for c in nodes]))
            mshfile.write('$EndElements\n')
            mshfile.close()
