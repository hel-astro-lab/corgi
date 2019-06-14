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
    tmp_grid = np.ones( (n.get_Nx(), n.get_Ny()) ) * -1.0

    for cid in n.get_tile_ids():
        c = n.get_tile( cid )

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

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        # one tile surrounded by other rank
        refGrid = np.ones((conf.Nx, conf.Ny), np.int)
        refGrid[1,1] = 0

        #setup grid configuration
        if grid.master():
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    val = refGrid[i,j]
                    grid.set_mpi_grid(i, j, val )
        grid.bcast_mpi_grid()

        # add tiles
        rank = grid.rank()
        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                if rank == refGrid[i,j]:
                    c = pycorgi.twoD.Tile()
                    grid.add_tile(c, (i,j) ) 

        try:
            plotNode(axs[0], grid, conf)
            saveVisz(0, grid, conf)
        except:
            pass

        if grid.size() > 1:
            grid.analyze_boundaries()
            #print(grid.rank(), ":sq ", grid.send_queue)
            #print(grid.rank(), ":sqa", grid.send_queue_address)

        if grid.size() > 1:
            grid.send_tiles()
            grid.recv_tiles()

        try:
            plotNode(axs[0], grid, conf)
            saveVisz(1, grid, conf)
        except:
            pass

        cur = readGrid(grid, conf)
        #print(cur)

        if grid.size() > 1:
            for i in range(grid.get_Nx()):
                for j in range(grid.get_Ny()):
                    self.assertEqual(refGrid[i,j], cur[i,j])



if __name__ == '__main__':
    unittest.main()
