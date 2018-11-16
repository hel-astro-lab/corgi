from mpi4py import MPI

import unittest
import numpy as np
import os

import pycorgi

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import palettable as pal
    palette = pal.wesanderson.Moonrise1_5.mpl_colormap
    from visualize import plotNode
    from visualize import saveVisz

except:
    pass


class Conf:
    Nx = 3
    Ny = 3
    Nz = 1

    NxMesh = 5
    NyMesh = 5
    NzMesh = 5

    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0

    outdir = "out"

def readGrid(n, conf):
    tmp_grid = np.ones( (n.getNx(), n.getNy()) ) * -1.0

    for cid in n.getTileIds():
        c = n.getTile( cid )

        try:
            (i, j) = c.index
        except:
            (i,) = c.index
            j = 0
        tmp_grid[i,j] = c.communication.owner

    return tmp_grid



class Neighboords(unittest.TestCase):

    def test_boundary_construction(self):

        conf = Conf()

        # set up plotting and figure
        try:
            plt.fig = plt.figure(1, figsize=(3,3))
            plt.rc('font', family='serif', size=12)
            plt.rc('xtick')
            plt.rc('ytick')
            
            gs = plt.GridSpec(1, 1)
            gs.update(hspace = 0.5)
            
            axs = []
            for ai in range(1):
                axs.append( plt.subplot(gs[ai]) )

            if not os.path.exists( conf.outdir ):
                os.makedirs(conf.outdir)
        except:
            pass

        node = pycorgi.twoD.Node(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        # one tile surrounded by other rank
        refGrid = np.ones((conf.Nx, conf.Ny), np.int)
        refGrid[1,1] = 0

        #setup grid configuration
        if node.master():
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    val = refGrid[i,j]
                    node.setMpiGrid(i, j, val )
        node.bcastMpiGrid()

        # add tiles
        rank = node.rank()
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                if rank == refGrid[i,j]:
                    c = pycorgi.twoD.Tile()
                    node.addTile(c, (i,j) ) 

        try:
            plotNode(axs[0], node, conf)
            saveVisz(0, node, conf)
        except:
            pass

        if node.size() > 1:
            node.analyzeBoundaryTiles()
            #print(node.rank(), ":sq ", node.send_queue)
            #print(node.rank(), ":sqa", node.send_queue_address)

        if node.size() > 1:
            node.communicateSendTiles()
            node.communicateRecvTiles()

        plotNode(axs[0], node, conf)
        saveVisz(1, node, conf)

        cur = readGrid(node, conf)
        #print(cur)

        if node.size() > 1:
            for i in range(node.getNx()):
                for j in range(node.getNy()):
                    self.assertEqual(refGrid[i,j], cur[i,j])



if __name__ == '__main__':
    unittest.main()
