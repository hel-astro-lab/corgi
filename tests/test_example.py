import unittest

import sys
sys.path.append('pycorgi')
sys.path.append('tests')


import example


class Initialization(unittest.TestCase):

    i = 10
    j = 11
    o = 1

    Nx = 10
    Ny = 20

    def setUp(self):
        self.cell1 = example.Welsh(self.i, self.j, self.o, self.Nx, self.Ny)
        self.cell2 = example.Pembroke(self.i, self.j, self.o, self.Nx, self.Ny)

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
         

