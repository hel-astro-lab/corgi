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
        node = pycorgi.oneD.Node(self.Nx)

        if node.master():
            for i in range(node.getNx()):
                c = pycorgi.oneD.Tile()
                node.addTile(c, (i, ) ) 

        ref_nhood = [
        [(2,), (1,)],
        [(0,), (2,)],
        [(1,), (0,)],
        ]

        for i in range(node.getNx()):
            cid = node.id(i)
            c = node.getTile(cid)
            self.assertCountEqual( c.nhood(), ref_nhood[i] )
            



    def test_nhood_2d(self):
        node = pycorgi.twoD.Node(self.Nx, self.Ny)

        if node.master():
            for i in range(node.getNx()):
                for j in range(node.getNy()):
                    c = pycorgi.twoD.Tile()
                    node.addTile(c, (i,j) ) 

        ref_nhood = [
            [(2, 2), (2, 0), (2, 1), (0, 1), (1, 1), (0, 1), (1, 2), (0, 2)],
            [(2, 0), (2, 1), (2, 2), (0, 2), (1, 2), (0, 2), (1, 0), (0, 0)],
            [(2, 1), (2, 2), (2, 0), (0, 0), (1, 0), (0, 0), (1, 1), (0, 1)],
            [(0, 2), (0, 0), (0, 1), (1, 1), (2, 1), (1, 1), (2, 2), (1, 2)],
            [(0, 0), (0, 1), (0, 2), (1, 2), (2, 2), (1, 2), (2, 0), (1, 0)],
            [(0, 1), (0, 2), (0, 0), (1, 0), (2, 0), (1, 0), (2, 1), (1, 1)],
            [(1, 2), (1, 0), (1, 1), (2, 1), (0, 1), (2, 1), (0, 2), (2, 2)],
            [(1, 0), (1, 1), (1, 2), (2, 2), (0, 2), (2, 2), (0, 0), (2, 0)],
            [(1, 1), (1, 2), (1, 0), (2, 0), (0, 0), (2, 0), (0, 1), (2, 1)],
        ]

        q = 0
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                #print("ij", i,j)
                cid = node.id(i,j)
                c = node.getTile(cid)
                #print(c.nhood())
                self.assertCountEqual( c.nhood(), ref_nhood[q] )

                q += 1


    def skip_test_nhood_3d(self):
        node = pycorgi.threeD.Node(self.Nx, self.Ny, self.Nz)

        if node.master():
            for i in range(node.getNx()):
                for j in range(node.getNy()):
                    for k in range(node.getNz()):
                        c = pycorgi.threeD.Tile()
                        node.addTile(c, (i,j,k) ) 

        # TODO: add reference
        #ref_nhood = [
        #]

        q = 0
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                for k in range(node.getNz()):
                    print("ijk", i,j,k)
                    cid = node.id(i,j,k)
                    c = node.getTile(cid)
                    print(c.nhood())
                    #self.assertCountEqual( c.nhood(), ref_nhood[q] )

                    q += 1







if __name__ == '__main__':
    unittest.main()
