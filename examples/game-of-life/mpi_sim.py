from mpi4py import MPI

import numpy as np
import os

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import palettable as pal
    palette = pal.wesanderson.Moonrise1_5.mpl_colormap
except:
    pass

import pycorgi
import pyca

from simulation import imshow
from simulation import plotNode
from simulation import plotMesh
from simulation import plotMesh2
from simulation import saveVisz
from simulation import loadMpiRandomly
from simulation import loadMpiXStrides
from simulation import load_tiles
from simulation import randomInitialize

np.random.seed(0)


# make all tiles same type 
def initialize_virtuals(n, conf):

    for cid in n.get_virtual_tiles():
        c_orig = n.get_tile(cid)
        (i,j) = c_orig.index

        # new CA tile;
        # TODO: load_metainfo *HAS* to be after add_tile
        c = pyca.Tile()
        n.add_tile(c, (i,j)) 
        c_orig.communication.local = False;
        c.load_metainfo(c_orig.communication)
        print("{}: loading {} owned by {}".format(n.rank(), cid, c.communication.owner))

        mesh = pyca.Mesh( conf["NxMesh"], conf["NyMesh"] )
        c.add_data(mesh)
        c.add_data(mesh)



##################################################
##################################################

if __name__ == "__main__":

    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,4))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 2)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )
    
    
    
    #setup grid
    conf = {
            "Nx"     : 10,
            "Ny"     : 10,
            "NxMesh" : 10,
            "NyMesh" : 10,
            "dir"    : "out2",
            }
    
    grid = pycorgi.twoD.Grid( conf["Nx"], conf["Ny"] ) 
    grid.set_grid_lims(0.0, 1.0, 0.0, 1.0)
    
    loadMpiRandomly(grid)
    #loadMpiXStrides(grid)

    load_tiles(grid)
    randomInitialize(grid, conf)
    
    # Path to be created 
    if grid.master:
        if not os.path.exists( conf["dir"]):
            os.makedirs(conf["dir"])
    
    sol = pyca.Solver()
    
    #static setup; communicate neighbor info once
    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    initialize_virtuals(grid, conf)

    plotNode(axs[0], grid, conf)
    saveVisz(0, grid, conf)

    #grid.send_data(0)
    #grid.recv_data(0)


    for lap in range(1, 301):
        print("---lap: {}".format(lap))

        #send/recv boundaries
        grid.send_data(0)
        grid.recv_data(0)
        grid.wait_data(0)

        if (lap % 10 == 0):
            plotNode(axs[0], grid, conf)
            plotMesh(axs[1], grid, conf)
            saveVisz(lap, grid, conf)

        #update halo regions
        for cid in grid.get_local_tiles():
            c = grid.get_tile(cid)
            c.update_boundaries(grid)

        #progress one time step
        for cid in grid.get_local_tiles():
            c = grid.get_tile(cid)
            sol.solve(c)

        #cycle everybody in time
        for cid in grid.get_local_tiles():
            c = grid.get_tile(cid)
            c.cycle()

    
    
    
