import unittest

import sys
sys.path.append('pycorgi')

import corgi 


class Params:
    mins = None
    maxs = None
    lens = None



class Initialization(unittest.TestCase):

    i = 10
    j = 11
    o = 1

    def setUp(self):
        self.cell = corgi.Cell(self.i, self.j, self.o)


    # test that we can inherit from the corgi::Cell base class
    def test_inheritance(self):

        (i,j) = self.cell.index()
        self.assertEqual(i, 10)
        self.assertEqual(j, 11)



if __name__ == '__main__':
    unittest.main()
