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
        grid = pycorgi.oneD.Grid(self.Nx)
        self.assertEqual(grid.get_Nx(), self.Nx)

        grid.set_grid_lims(self.xmin, self.xmax)
        self.assertEqual( grid.get_xmin(), self.xmin )
        self.assertEqual( grid.get_xmax(), self.xmax )


    def test_oneD2(self):
        grid = pycorgi.oneD.Grid(self.Nx, self.Ny)
        self.assertEqual(grid.get_Nx(), self.Nx)
        self.assertEqual(grid.get_Ny(), 1)

        grid.set_grid_lims(self.xmin, self.xmax, self.ymin, self.ymax)
        self.assertEqual( grid.get_xmin(), self.xmin )
        self.assertEqual( grid.get_xmax(), self.xmax )

        self.assertEqual( grid.get_ymin(), 0.0 )
        self.assertEqual( grid.get_ymax(), 1.0 )

    def test_oneD3(self):
        grid = pycorgi.oneD.Grid(self.Nx, self.Ny, self.Nz)
        self.assertEqual(grid.get_Nx(), self.Nx)
        self.assertEqual(grid.get_Ny(), 1)
        self.assertEqual(grid.get_Nz(), 1)

        #TODO set_grid_lims
    
    def test_twoD2(self):
        grid = pycorgi.twoD.Grid(self.Nx, self.Ny)
        self.assertEqual(grid.get_Nx(), self.Nx)
        self.assertEqual(grid.get_Ny(), self.Ny)

        grid.set_grid_lims(self.xmin, self.xmax, self.ymin, self.ymax)
        self.assertEqual( grid.get_xmin(), self.xmin )
        self.assertEqual( grid.get_xmax(), self.xmax )

        self.assertEqual( grid.get_ymin(), self.ymin )
        self.assertEqual( grid.get_ymax(), self.ymax )

    def test_twoD3(self):
        grid = pycorgi.twoD.Grid(self.Nx, self.Ny, self.Nz)
        self.assertEqual(grid.get_Nx(), self.Nx)
        self.assertEqual(grid.get_Ny(), self.Ny)
        self.assertEqual(grid.get_Nz(), 1)

        #TODO set_grid_lims



if __name__ == '__main__':
    unittest.main()


