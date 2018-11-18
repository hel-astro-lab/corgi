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
        self.assertEqual(node.get_Nx(), self.Nx)

        node.set_grid_lims(self.xmin, self.xmax)
        self.assertEqual( node.get_xmin(), self.xmin )
        self.assertEqual( node.get_xmax(), self.xmax )


    def test_oneD2(self):
        node = pycorgi.oneD.Node(self.Nx, self.Ny)
        self.assertEqual(node.get_Nx(), self.Nx)
        self.assertEqual(node.get_Ny(), 1)

        node.set_grid_lims(self.xmin, self.xmax, self.ymin, self.ymax)
        self.assertEqual( node.get_xmin(), self.xmin )
        self.assertEqual( node.get_xmax(), self.xmax )

        self.assertEqual( node.get_ymin(), 0.0 )
        self.assertEqual( node.get_ymax(), 1.0 )

    def test_oneD3(self):
        node = pycorgi.oneD.Node(self.Nx, self.Ny, self.Nz)
        self.assertEqual(node.get_Nx(), self.Nx)
        self.assertEqual(node.get_Ny(), 1)
        self.assertEqual(node.get_Nz(), 1)

        #TODO set_grid_lims
    
    def test_twoD2(self):
        node = pycorgi.twoD.Node(self.Nx, self.Ny)
        self.assertEqual(node.get_Nx(), self.Nx)
        self.assertEqual(node.get_Ny(), self.Ny)

        node.set_grid_lims(self.xmin, self.xmax, self.ymin, self.ymax)
        self.assertEqual( node.get_xmin(), self.xmin )
        self.assertEqual( node.get_xmax(), self.xmax )

        self.assertEqual( node.get_ymin(), self.ymin )
        self.assertEqual( node.get_ymax(), self.ymax )

    def test_twoD3(self):
        node = pycorgi.twoD.Node(self.Nx, self.Ny, self.Nz)
        self.assertEqual(node.get_Nx(), self.Nx)
        self.assertEqual(node.get_Ny(), self.Ny)
        self.assertEqual(node.get_Nz(), 1)

        #TODO set_grid_lims



if __name__ == '__main__':
    unittest.main()


