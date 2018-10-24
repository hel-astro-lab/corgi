from mpi4py import MPI

import unittest

import numpy as np
import pycorgi


class Params:
    mins = None
    maxs = None
    lens = None


class Initialization(unittest.TestCase):
    Nx = 10
    Ny = 20
    Nz = 30

    xmin = 0.0
    xmax = 1.0
    ymin = 2.0
    ymax = 3.0
    zmin = 4.0
    zmax = 5.0

    def test_oneD(self):
        node = pycorgi.oneD.Node(self.Nx)
        self.assertEqual(node.getNx(), self.Nx)

        node.setGridLims(self.xmin, self.xmax)
        self.assertEqual( node.getXmin(), self.xmin )
        self.assertEqual( node.getXmax(), self.xmax )


    def test_oneD2(self):
        node = pycorgi.oneD.Node(self.Nx, self.Ny)
        self.assertEqual(node.getNx(), self.Nx)
        self.assertEqual(node.getNy(), 1)

        node.setGridLims(self.xmin, self.xmax, self.ymin, self.ymax)
        self.assertEqual( node.getXmin(), self.xmin )
        self.assertEqual( node.getXmax(), self.xmax )

        self.assertEqual( node.getYmin(), 0.0 )
        self.assertEqual( node.getYmax(), 1.0 )

    def test_oneD3(self):
        node = pycorgi.oneD.Node(self.Nx, self.Ny, self.Nz)
        self.assertEqual(node.getNx(), self.Nx)
        self.assertEqual(node.getNy(), 1)
        self.assertEqual(node.getNz(), 1)

        #TODO setGridLims
    
    def test_twoD2(self):
        node = pycorgi.twoD.Node(self.Nx, self.Ny)
        self.assertEqual(node.getNx(), self.Nx)
        self.assertEqual(node.getNy(), self.Ny)

        node.setGridLims(self.xmin, self.xmax, self.ymin, self.ymax)
        self.assertEqual( node.getXmin(), self.xmin )
        self.assertEqual( node.getXmax(), self.xmax )

        self.assertEqual( node.getYmin(), self.ymin )
        self.assertEqual( node.getYmax(), self.ymax )

    def test_twoD3(self):
        node = pycorgi.twoD.Node(self.Nx, self.Ny, self.Nz)
        self.assertEqual(node.getNx(), self.Nx)
        self.assertEqual(node.getNy(), self.Ny)
        self.assertEqual(node.getNz(), 1)

        #TODO setGridLims



if __name__ == '__main__':
    unittest.main()


