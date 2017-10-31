from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import palettable as pal
cmap = pal.wesanderson.Moonrise1_5.mpl_colormap


import sys, os

import corgi


Nrank = 4

#make random starting order
def loadMpiRandomly(n):
    np.random.seed(4)
    if n.master:
        for i in range(corgi.Nx):
            for j in range(corgi.Ny):
                val = np.random.randint(n.Nrank)
                n.setMpiGrid(i, j, val)


#load nodes to be in stripe formation
def loadMpiStrides(n):
    if n.master: #only master initializes; then sends
        stride = np.zeros( (corgi.Ny), np.int64)
        dy = np.float(corgi.Ny) / np.float(n.Nrank) 
        for j in range(corgi.Ny):
            val = np.int( j/dy )
            stride[j] = val

        for i in range(corgi.Nx):
            for j in range(corgi.Ny):
                val = stride[j]
                n.setMpiGrid(i, j, val)
    n.bcastMpiGrid()


#load cells into each node
def loadCells(n):
    for i in range(corgi.Nx):
        for j in range(corgi.Ny):
            #print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.mpiGrid(i,j), ref[j,i]))
            if n.mpiGrid(i,j) == n.rank:
                c = corgi.Cell(i, j, n.rank)
                n.addLocalCell(c) #TODO load data to cell



##################################################
# plotting tools

# visualize matrix
def imshow(ax, grid):
    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(corgi.xmin, corgi.xmax)
    ax.set_ylim(corgi.ymin, corgi.ymax)

    extent = [corgi.xmin, corgi.xmax, corgi.ymin, corgi.ymax]

    mgrid = np.ma.masked_where(grid == -1.0, grid)
    ax.imshow(mgrid,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = cmap,
              vmin = 0.0,
              vmax = Nrank-1,
              aspect='auto',
              #vmax = Nrank,
              #alpha=0.5
              )


# Visualize current cell ownership on node
def plot_node(ax, n, lap):
    tmp_grid = np.ones( (corgi.Nx, corgi.Ny) ) * -1.0

    
    #for i in range(corgi.Nx):
    #    for j in range(corgi.Ny):
    #        cid = n.cell_id(i,j)
    #        if n.is_local(cid):
    #            tmp_grid[i,j] = 0.5


    for cid in n.getCells():
        c = n.getCell( cid )
        (i, j) = c.index()
        #check dublicates
        if tmp_grid[i,j] != -1.0:
            print("{}: ERROR in real cells at ({},{})".format(n.rank, i,j))
            sys.exit()
        tmp_grid[i,j] = c.owner


    for cid in n.getVirtuals():
        c = n.getCell( cid )
        (i,j) = c.index()
        if tmp_grid[i,j] != -1.0:
            print("{}: ERROR in virtual cells at ({},{})".format(n.rank, i,j))
            sys.exit()
        tmp_grid[i,j] = c.owner

    imshow(ax, tmp_grid)
    #imshow(ax, n.mpiGrid)


    # add text label about number of neighbors
    for cid in n.getCells():
        c = n.getCell( cid )
        (i, j) = c.index()
        ix = corgi.xmin + corgi.xmax*(i+0.5)/corgi.Nx
        jy = corgi.ymin + corgi.ymax*(j+0.5)/corgi.Ny

        #Nv = n.number_of_virtual_neighbors(c)
        Nv = c.number_of_virtual_neighbors
        label = str(Nv)
        #label = "{} ({},{})/{}".format(cid,i,j,Nv)
        #label = "({},{})".format(i,j)
        ax.text(jy, ix, label, ha='center',va='center', size=8)


    #for cid in n.getVirtuals():
    #    c = n.getCell( cid )
    #    (i,j) = c.index()
    #    ix = corgi.xmin + corgi.xmax*(i+0.5)/corgi.Nx
    #    jy = corgi.ymin + corgi.ymax*(j+0.5)/corgi.Ny
    #    label = "Vir"
    #    ax.text(jy, ix, label, ha='center',va='center')


    ax.set_title(str(len(n.getVirtuals() ))+"/"+str(len(n.getCells() )))


    #save
    slap = str(lap).rjust(4, '0')
    fname = fpath + '/node_{}_{}.png'.format(n.rank, slap)
    plt.savefig(fname)






if __name__ == "__main__":

    ################################################## 
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



    ################################################## 
    #init node
    node = corgi.Node()
    node.initMpi()
    loadMpiStrides(node)
    loadCells(node)


    # Path to be created 
    fpath = "out/"
    if node.master:
        if not os.path.exists(fpath):
            os.makedirs(fpath)


    ################################################## 
    #visualize as a test
    plot_node(axs[0], node, 0)


    ################################################## 
    # test step
    node.analyzeBoundaryCells()
    print("{}: send queue        : {}".format(node.rank, node.send_queue))
    print("{}: send queue address: {}".format(node.rank, node.send_queue_address))

    node.communicateSendCells()
    node.communicateRecvCells()
    plot_node(axs[0], node, 1)



    node.finalizeMpi()


