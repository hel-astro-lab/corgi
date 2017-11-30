import unittest

import numpy as np
import sys
sys.path.append('pycorgi')

import corgi 


class Params:
    mins = None
    maxs = None
    lens = None


class Initialization(unittest.TestCase):
    Nx = 10
    Ny = 20

    xmin = 0.0
    xmax = 1.0
    ymin = 2.0
    ymax = 3.0

    def setUp(self):
        self.node = corgi.Node(self.Nx, self.Ny)

        self.node.setGridLims(self.xmin, self.xmax,
                              self.ymin, self.ymax
                              )

    def test_size(self):
        nx = self.node.getNx()
        ny = self.node.getNy()

        self.assertEqual(self.node.getNx(), self.Nx)
        self.assertEqual(self.node.getNy(), self.Ny)

    def test_physicalSize(self):
        self.assertEqual( self.node.getXmin(), self.xmin )
        self.assertEqual( self.node.getXmax(), self.xmax )

        self.assertEqual( self.node.getYmin(), self.ymin )
        self.assertEqual( self.node.getYmax(), self.ymax )


def cellID(i,j,Nx,Ny):
    return j*Nx + i


class Parallel(unittest.TestCase):
    
    Nx = 10
    Ny = 15

    xmin = 0.0
    xmax = 1.0
    ymin = 2.0
    ymax = 3.0


    def setUp(self):
        self.node = corgi.Node(self.Nx, self.Ny)
        self.node.setGridLims(self.xmin, self.xmax, self.ymin, self.ymax)

    def mpiInitialization(self):

        self.node.initMpi()

        self.refGrid = np.zeros((self.Nx, self.Ny), np.int)
        self.refGrid[0:5,   0:10] = 0
        self.refGrid[0:5,  10:15] = 1
        self.refGrid[5:10,  0:10] = 2
        self.refGrid[5:10, 10:15] = 3

        if self.node.master:
            for j in range(self.node.getNy()):
                for i in range(self.node.getNx()):
                    val = self.refGrid[i,j]
                    self.node.setMpiGrid(i, j, val )
        self.node.bcastMpiGrid()

        for j in range(self.node.getNy()):
            for i in range(self.node.getNx()):
                val = self.node.getMpiGrid(i,j)
                self.assertEqual(val, self.refGrid[i,j])

        self.node.finalizeMpi()


    def test_cid(self):
        for j in range(self.node.getNy()):
            for i in range(self.node.getNx()):
                cid = self.node.cellId(i, j)
                cidr = cellID( i, j, self.node.getNx(), self.node.getNy() )
                self.assertEqual(cid, cidr)


    def test_loading(self):

        #load cells
        k = 0
        for j in range(self.node.getNy()):
            for i in range(self.node.getNx()):
                c = corgi.Cell(i, j, 0, self.node.getNx(), self.node.getNy() )
                self.node.addCell(c) 
                k += 1
        self.assertEqual( k, self.Nx*self.Ny )

        cids = self.node.getCellIds() 
        self.assertEqual( len(cids), self.Nx*self.Ny )

        #now try and get then back
        for cid in cids:
            c = self.node.getCellPtr(cid)

            self.assertEqual(c.cid,   cid)
            self.assertEqual(c.owner, self.node.rank)
            self.assertEqual(c.local, True)

            #ci = c.i
            #cj = c.j
            #self.assertEqual(ci, ri)
            #self.assertEqual(cj, rj)





if __name__ == '__main__':
    unittest.main()


