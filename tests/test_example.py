import unittest
import sys

import pycorgi
import pycorgitest


class Initialization(unittest.TestCase):

    i = 10
    j = 11
    o = 1

    Nx = 10
    Ny = 20

    def setUp(self):
        self.cell1 = pycorgitest.Welsh(self.i, self.j, self.o, self.Nx, self.Ny)
        self.cell2 = pycorgitest.Pembroke(self.i, self.j, self.o, self.Nx, self.Ny)

    #test that derived classes can inherit base class methods
    def test_inheritance(self):

        (i,j) = self.cell1.index()
        self.assertEqual(i, self.i)
        self.assertEqual(j, self.j)

        (i,j) = self.cell2.index()
        self.assertEqual(i, self.i)
        self.assertEqual(j, self.j)

    #tests that dserived classes can be extended
    def test_extending(self):

        self.assertEqual( self.cell1.bark(), "Woof!" )
        self.assertEqual( self.cell2.bark(), "Ruff!" )

        self.assertEqual( self.cell2.howl(), "Auuuuuu!" )


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

        self.vallhund = pycorgitest.Vallhund(self.i, self.j, self.o, self.Nx, self.Ny, 
                                        self.soc_num)

    def test_extending(self):
        self.assertEqual( self.vallhund.bark(), "ruf ruf ruf" )

    def test_multiple_inheriting(self):

        self.assertEqual( self.vallhund.fika(), self.swede.fika() )
        self.assertEqual( self.vallhund.prayForOdin(), self.viking.prayForOdin() )

        self.assertEqual( self.vallhund.number, self.swede.number )



def cellID(i,j,Nx,Ny):
    return j*Nx + i


class ParallelGrid(unittest.TestCase):
    
    Nx = 10
    Ny = 15

    xmin = 0.0
    xmax = 1.0
    ymin = 2.0
    ymax = 3.0


    def setUp(self):
        self.node = pycorgitest.Grid(self.Nx, self.Ny)
        self.node.setGridLims(self.xmin, self.xmax, self.ymin, self.ymax)


    def test_extension(self):
        self.assertEqual( self.node.petShop(), "No Corgis for sale.")


    def test_cid(self):
        for j in range(self.node.getNy()):
            for i in range(self.node.getNx()):
                cid = self.node.cellId(i, j)
                cidr = cellID( i, j, self.node.getNx(), self.node.getNy() )
                self.assertEqual(cid, cidr)

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
                val = self.node.mpiGrid(i,j)
                self.assertEqual(val, self.refGrid[i,j])
        self.node.finalizeMpi()



    def test_loading(self):

        k = 0
        for j in range(self.node.getNy()):
            for i in range(self.node.getNx()):

                #initialize heterogeneous grid
                if (i%2 == 0):
                    c = pycorgitest.Welsh(i, j, 0, self.node.getNx(), self.node.getNy() )
                else:
                    c = pycorgitest.Pembroke(i, j, 0, self.node.getNx(), self.node.getNy() )
                self.node.addCell(c) 
                k += 1

        self.assertEqual( k, self.Nx*self.Ny )

        cids = self.node.getCellIds() 
        self.assertEqual( len(cids), self.Nx*self.Ny )

        #now try and get them back
        for cid in cids:
            c = self.node.getCellPtr(cid)

            self.assertEqual(c.cid,   cid)
            self.assertEqual(c.owner, self.node.rank)
            self.assertEqual(c.local, True)

            # we need to be able to bark also after the getting.
            # This tests that operator slicing is not acting

            if (c.i % 2 == 0): #Welsh
                self.assertEqual( c.bark(), "Woof!" )
            else:              #Pembroke
                self.assertEqual( c.bark(), "Ruff!" )





if __name__ == '__main__':
    unittest.main()

