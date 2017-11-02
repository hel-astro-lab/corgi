import unittest

import sys
sys.path.append('pycorgi')

import corgi 


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
        self.node = corgi.Node(self.Nx, self.Ny)

        self.node.setGridLims(self.xmin, self.xmax,
                              self.ymin, self.ymax
                              )


    def test_size(self):
        nx = self.node.getNx()
        ny = self.node.getNy()

        self.assertEqual(self.node.getNx(), self.Nx)
        self.assertEqual(self.node.getNy(), self.Ny)

    def test_physicalSize(self):
        self.assertEqual( self.node.getXmin(), self.xmin )
        self.assertEqual( self.node.getXmax(), self.xmax )

        self.assertEqual( self.node.getYmin(), self.ymin )
        self.assertEqual( self.node.getYmax(), self.ymax )






if __name__ == '__main__':
    unittest.main()


