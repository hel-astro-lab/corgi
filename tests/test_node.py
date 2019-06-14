from mpi4py import MPI

import unittest

import numpy as np
import sys
import pycorgi.twoD as pycorgi

#sys.path.append('../lib')


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
        self.grid = pycorgi.Grid(self.Nx, self.Ny)

        self.grid.set_grid_lims(self.xmin, self.xmax,
                              self.ymin, self.ymax
                              )

    def test_size(self):
        nx = self.grid.get_Nx()
        ny = self.grid.get_Ny()

        self.assertEqual(self.grid.get_Nx(), self.Nx)
        self.assertEqual(self.grid.get_Ny(), self.Ny)

    def test_physicalSize(self):
        self.assertEqual( self.grid.get_xmin(), self.xmin )
        self.assertEqual( self.grid.get_xmax(), self.xmax )

        self.assertEqual( self.grid.get_ymin(), self.ymin )
        self.assertEqual( self.grid.get_ymax(), self.ymax )


def tile_id(i,j,Nx,Ny):
    return j*Nx + i


class Parallel(unittest.TestCase):
    
    Nx = 10
    Ny = 15

    xmin = 0.0
    xmax = 1.0
    ymin = 2.0
    ymax = 3.0


    def setUp(self):
        self.grid = pycorgi.Grid(self.Nx, self.Ny)
        self.grid.set_grid_lims(self.xmin, self.xmax, self.ymin, self.ymax)

    def test_mpiInitialization(self):

        self.refGrid = np.zeros((self.Nx, self.Ny), np.int)
        self.refGrid[0:5,   0:10] = 0
        self.refGrid[0:5,  10:15] = 1
        self.refGrid[5:10,  0:10] = 2
        self.refGrid[5:10, 10:15] = 3

        if self.grid.master():
            for j in range(self.grid.get_Ny()):
                for i in range(self.grid.get_Nx()):
                    val = self.refGrid[i,j]
                    self.grid.set_mpi_grid(i, j, val )
        self.grid.bcast_mpi_grid()

        for j in range(self.grid.get_Ny()):
            for i in range(self.grid.get_Nx()):
                val = self.grid.get_mpi_grid(i,j)
                self.assertEqual(val, self.refGrid[i,j])

    def test_cid(self):
        for j in range(self.grid.get_Ny()):
            for i in range(self.grid.get_Nx()):
                cid = self.grid.id(i, j)
                cidr = tile_id( i, j, self.grid.get_Nx(), self.grid.get_Ny() )
                self.assertEqual(cid, cidr)

    def test_loading(self):

        #load tiles
        k = 0
        for j in range(self.grid.get_Ny()):
            for i in range(self.grid.get_Nx()):
                c = pycorgi.Tile()
                self.grid.add_tile(c, (i,j) ) 
                k += 1
        self.assertEqual( k, self.Nx*self.Ny )

        cids = self.grid.get_tile_ids() 
        self.assertEqual( len(cids), self.Nx*self.Ny )

        #now try and get then back
        for cid in cids:
            c = self.grid.get_tile(cid)

            self.assertEqual(c.cid,   cid)
            self.assertEqual(c.communication.owner, self.grid.rank())
            #self.assertEqual(c.communication.local, True)

            #ci = c.i
            #cj = c.j
            #self.assertEqual(ci, ri)
            #self.assertEqual(cj, rj)


# advanced parallel tests
class Parallel2(unittest.TestCase):
    
    Nx = 10
    Ny = 2

    xmin = 0.0
    xmax = 1.0
    ymin = 2.0
    ymax = 3.0


    def test_send_recv(self):

        #set up
        grid = pycorgi.Grid(self.Nx, self.Ny)
        grid.set_grid_lims(self.xmin, self.xmax, self.ymin, self.ymax)

        # divide into upper and lower halfs
        #refGrid = np.zeros((self.Nx, self.Ny), np.int)
        #refGrid[:, 0] = 0
        #refGrid[:, 1] = 1

        #if grid.master():
        #    for j in range(grid.get_Ny()):
        #        for i in range(grid.get_Nx()):
        #            val = refGrid[i,j]
        #            grid.set_mpi_grid(i, j, val )
        #grid.bcast_mpi_grid()

        #load tiles
        if grid.rank() == 0:
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    #if grid.get_mpi_grid(i,j) == 0:
                    c = pycorgi.Tile()
                    grid.add_tile(c, (i,j) ) 

        #0 sends
        if grid.rank() == 0 and grid.size() > 1:
            #load cell with info
            cid = grid.id(2, 1)
            #print("0:  send............cid:", cid)
            c = grid.get_tile(cid)

            #communication object
            c.communication.top_virtual_owner           = 10
            c.communication.communications              = 11
            c.communication.number_of_virtual_neighbors = 12

            c.set_tile_mins([1.0,2.0])
            c.set_tile_maxs([1.1,2.1])



            #get the same cell and send
            grid.send_tile(cid, 1)


        #1 does nothing but receives
        if grid.rank() == 1 and grid.size() > 1:
            grid.recv_tile(0)
            #print("1 recv............")
        
        #grid.wait()

        #assert that we received the tile properly
        if grid.rank() == 1 and grid.size() > 1:
            cid = grid.id(2, 1)
            c = grid.get_tile(cid)
            #print("1:  cid=", cid)

            #cid
            #owner
            #top_virtual_owner
            #communications
            #number_of_virtual_neighbors

            self.assertEqual(c.communication.cid, cid)
            self.assertEqual(c.communication.owner, 0)
            self.assertEqual(c.communication.top_virtual_owner, 10)
            self.assertEqual(c.communication.communications,    11)
            self.assertEqual(c.communication.number_of_virtual_neighbors, 12)

            #indices
            #mins
            #maxs
            self.assertEqual(c.communication.indices, [2, 1, 0] )
            self.assertEqual(c.communication.mins,    [1.0, 2.0, 0.0] )
            self.assertEqual(c.communication.maxs,    [1.1, 2.1, 0.0] )

            #tile variables
            self.assertEqual(c.cid,     cid)
            self.assertEqual(c.mins,    [1.0, 2.0] )
            self.assertEqual(c.maxs,    [1.1, 2.1] )
            self.assertEqual(c.index,   (2,1) )



            #check that cell is reconstructed correctly from Communication obj
            self.assertEqual(c.cid, cid)


if __name__ == '__main__':
    unittest.main()


