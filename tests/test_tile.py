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
        self.tile = pycorgi.Tile2D()
        self.tile.index = (self.i, self.j)
        self.tile.communication.owner = self.o

    def test_indexing(self):
        (i,j) = self.tile.index
        self.assertEqual(i, self.i)
        self.assertEqual(j, self.j)

        o = self.tile.communication.owner
        self.assertEqual(o, self.o)



if __name__ == '__main__':
    unittest.main()
