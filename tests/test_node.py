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
        self.node.initMpi()

    def test_loading(self):

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
                val = self.node.mpiGrid(i,j)
                self.assertEqual(val, self.refGrid[i,j])
    

    def tearDown(self):
        self.node.finalizeMpi()




if __name__ == '__main__':
    unittest.main()


