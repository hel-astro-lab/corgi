from mpi4py import MPI

import unittest
import numpy as np

import pycorgi



class Neighboords(unittest.TestCase):
    Nx = 3
    Ny = 3
    Nz = 3

    xmin = 0.0
    xmax = 1.0
    ymin = 2.0
    ymax = 3.0


    def test_nhood_1d(self):
        grid = pycorgi.oneD.Grid(self.Nx)

        if grid.master():
            for i in range(grid.get_Nx()):
                c = pycorgi.oneD.Tile()
                grid.add_tile(c, (i, ) ) 

        ref_nhood = [
        [(2,), (1,)],
        [(0,), (2,)],
        [(1,), (0,)],
        ]

        for i in range(grid.get_Nx()):
            cid = grid.id(i)
            c = grid.get_tile(cid)
            self.assertCountEqual( c.nhood(), ref_nhood[i] )
            



    def test_nhood_2d(self):
        grid = pycorgi.twoD.Grid(self.Nx, self.Ny)

        if grid.master():
            for i in range(grid.get_Nx()):
                for j in range(grid.get_Ny()):
                    c = pycorgi.twoD.Tile()
                    grid.add_tile(c, (i,j) ) 

        ref_nhood = [
            [(2, 2), (2, 0), (2, 1), (0, 2), (0, 1), (1, 2), (1, 0), (1, 1)],
            [(2, 0), (2, 1), (2, 2), (0, 0), (0, 2), (1, 0), (1, 1), (1, 2)],
            [(2, 1), (2, 2), (2, 0), (0, 1), (0, 0), (1, 1), (1, 2), (1, 0)],
            [(0, 2), (0, 0), (0, 1), (1, 2), (1, 1), (2, 2), (2, 0), (2, 1)],
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1), (2, 2)],
            [(0, 1), (0, 2), (0, 0), (1, 1), (1, 0), (2, 1), (2, 2), (2, 0)],
            [(1, 2), (1, 0), (1, 1), (2, 2), (2, 1), (0, 2), (0, 0), (0, 1)],
            [(1, 0), (1, 1), (1, 2), (2, 0), (2, 2), (0, 0), (0, 1), (0, 2)],
            [(1, 1), (1, 2), (1, 0), (2, 1), (2, 0), (0, 1), (0, 2), (0, 0)],
        ]


        q = 0
        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                #print("ij", i,j)
                cid = grid.id(i,j)
                c = grid.get_tile(cid)
                #print(c.nhood())
                self.assertCountEqual( c.nhood(), ref_nhood[q] )

                q += 1


    def skip_test_nhood_3d(self):
        grid = pycorgi.threeD.Grid(self.Nx, self.Ny, self.Nz)

        if grid.master():
            for i in range(grid.get_Nx()):
                for j in range(grid.get_Ny()):
                    for k in range(grid.get_Nz()):
                        c = pycorgi.threeD.Tile()
                        grid.add_tile(c, (i,j,k) ) 

        # TODO: add reference
        #ref_nhood = [
        #]

        q = 0
        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                for k in range(grid.get_Nz()):
                    print("ijk", i,j,k)
                    cid = grid.id(i,j,k)
                    c = grid.get_tile(cid)
                    print(c.nhood())
                    #self.assertCountEqual( c.nhood(), ref_nhood[q] )

                    q += 1







if __name__ == '__main__':
    unittest.main()
