from mpi4py import MPI

import unittest
import sys

import pycorgi.twoD
import pycorgitest


class Initialization(unittest.TestCase):

    i = 10
    j = 11
    o = 1

    Nx = 10
    Ny = 20

    def setUp(self):
        self.tile1 = pycorgitest.Welsh()
        self.tile2 = pycorgitest.Pembroke()

        self.tile1.index = (self.i, self.j)
        self.tile1.communication.owner = self.o

        self.tile2.index = (self.i, self.j)
        self.tile2.communication.owner = self.o

    #test that derived classes can inherit base class methods
    def test_inheritance(self):

        (i,j) = self.tile1.index
        self.assertEqual(i, self.i)
        self.assertEqual(j, self.j)

        (i,j) = self.tile2.index
        self.assertEqual(i, self.i)
        self.assertEqual(j, self.j)

    #tests that dserived classes can be extended
    def test_extending(self):

        self.assertEqual( self.tile1.bark(), "Woof!" )
        self.assertEqual( self.tile2.bark(), "Ruff!" )

        self.assertEqual( self.tile2.howl(), "Auuuuuu!" )


# Testing multiple inheritance bindings
class MultipleInheritance(unittest.TestCase):

    i = 10
    j = 11
    o = 1

    Nx = 10
    Ny = 20

    soc_num = 123456

    def setUp(self):
        self.swede  = pycorgitest.Swede(self.soc_num)
        self.viking = pycorgitest.Viking()

        self.vallhund = pycorgitest.Vallhund(self.soc_num)

    def test_extending(self):
        self.assertEqual( self.vallhund.bark(), "ruf ruf ruf" )

    def test_multiple_inheriting(self):

        self.assertEqual( self.vallhund.fika(), self.swede.fika() )
        self.assertEqual( self.vallhund.pray_for_odin(), self.viking.pray_for_odin() )

        self.assertEqual( self.vallhund.number, self.swede.number )



def tile_id(i,j,Nx,Ny):
    return j*Nx + i


class ParallelGrid(unittest.TestCase):
    
    Nx = 10
    Ny = 15

    xmin = 0.0
    xmax = 1.0
    ymin = 2.0
    ymax = 3.0


    def setUp(self):
        self.node = pycorgi.twoD.Node(self.Nx, self.Ny)
        self.node.set_grid_lims(self.xmin, self.xmax, self.ymin, self.ymax)


    #def test_extension(self):
    #    self.assertEqual( self.node.pet_shop(), "No Corgis for sale.")


    def test_cid(self):
        for j in range(self.node.get_Ny()):
            for i in range(self.node.get_Nx()):
                cid = self.node.id(i, j)
                cidr = tile_id( i, j, self.node.get_Nx(), self.node.get_Ny() )
                self.assertEqual(cid, cidr)

    def mpiInitialization(self):

        self.node.initMpi()

        self.refGrid = np.zeros((self.Nx, self.Ny), np.int)
        self.refGrid[0:5,   0:10] = 0
        self.refGrid[0:5,  10:15] = 1
        self.refGrid[5:10,  0:10] = 2
        self.refGrid[5:10, 10:15] = 3

        if self.node.master:
            for j in range(self.node.get_Ny()):
                for i in range(self.node.get_Nx()):
                    val = self.refGrid[i,j]
                    self.node.set_mpi_grid(i, j, val )
        self.node.bcast_mpi_grid()

        for j in range(self.node.get_Ny()):
            for i in range(self.node.get_Nx()):
                val = self.node.mpi_grid(i,j)
                self.assertEqual(val, self.refGrid[i,j])
        self.node.finalizeMpi()



    def test_loading(self):

        k = 0
        for j in range(self.node.get_Ny()):
            for i in range(self.node.get_Nx()):

                #initialize heterogeneous grid
                if (i % 2 == 0):
                    c = pycorgitest.Welsh()
                else:
                    c = pycorgitest.Pembroke()
                self.node.add_tile(c, (i,j) ) 
                k += 1

        self.assertEqual( k, self.Nx*self.Ny )

        cids = self.node.get_tile_ids() 
        self.assertEqual( len(cids), self.Nx*self.Ny )

        #now try and get them back
        for cid in cids:
            c = self.node.get_tile(cid)

            self.assertEqual(c.cid,   cid)
            self.assertEqual(c.communication.owner, self.node.rank())
            self.assertEqual(c.communication.local, True)

            # we need to be able to bark also after the getting.
            # This tests that operator slicing is not acting

            (i,j) = c.index
            if (i % 2 == 0): #Welsh
                self.assertEqual( c.bark(), "Woof!" )
            else:              #Pembroke
                self.assertEqual( c.bark(), "Ruff!" )





if __name__ == '__main__':
    unittest.main()

