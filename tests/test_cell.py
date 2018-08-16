import unittest
import sys
import pycorgi 


class Params:
    mins = None
    maxs = None
    lens = None



class Initialization(unittest.TestCase):

    i = 10
    j = 11
    o = 1

    Nx = 10
    Ny = 20

    def setUp(self):
        self.cell = pycorgi.Cell(self.i, self.j, self.o, self.Nx, self.Ny)

    def test_indexing(self):
        (i,j) = self.cell.index()
        self.assertEqual(i, self.i)
        self.assertEqual(j, self.j)






if __name__ == '__main__':
    unittest.main()
