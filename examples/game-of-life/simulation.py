import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
import palettable as pal
palette = pal.wesanderson.Moonrise1_5.mpl_colormap


import pygol


# visualize matrix
def imshow(ax, 
           grid, xmin, xmax, ymin, ymax,
           cmap='plasma',
           vmin = 0.0,
           vmax = 1.0,
           clip = -1.0,
          ):

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    #ax.set_xlim(-3.0, 3.0)
    #ax.set_ylim(-3.0, 3.0)

    extent = [ xmin, xmax, ymin, ymax ]

    mgrid = np.ma.masked_where(grid <= clip, grid)
    
    mgrid = mgrid.T
    ax.imshow(mgrid,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = cmap,
              vmin = vmin,
              vmax = vmax,
              aspect='auto',
              #vmax = Nrank,
              #alpha=0.5
              )
 




# Visualize current cell ownership on node
def plotNode(ax, n, conf):
    tmp_grid = np.ones( (n.getNx(), n.getNy()) ) * -1.0
    
    #for i in range(n.getNx()):
    #    for j in range(n.getNy()):
    #        cid = n.cell_id(i,j)
    #        if n.is_local(cid):
    #            tmp_grid[i,j] = 0.5


    for cid in n.getCellIds():
        c = n.getCellPtr( cid )
        (i, j) = c.index()
        #check dublicates
        if tmp_grid[i,j] != -1.0:
            print("{}: ERROR in real cells at ({},{})".format(n.rank, i,j))
            sys.exit()
        tmp_grid[i,j] = c.owner

    #XXX add back
    #for cid in n.getVirtuals():
    #    c = n.getCell( cid )
    #    (i,j) = c.index()
    #    if tmp_grid[i,j] != -1.0:
    #        print("{}: ERROR in virtual cells at ({},{})".format(n.rank, i,j))
    #        sys.exit()
    #    tmp_grid[i,j] = c.owner

    imshow(ax, tmp_grid, 
            n.getXmin(), n.getXmax(), n.getYmin(), n.getYmax(),
            cmap = palette,
            vmin = 0.0,
            vmax = conf["Nrank"]-1
            )


    # add text label about number of neighbors
    for cid in n.getCellIds():
        c = n.getCellPtr( cid )
        (i, j) = c.index()
        dx = n.getXmax() - n.getXmin()
        dy = n.getYmax() - n.getYmin()

        ix = n.getXmin() + dx*(i+0.5)/n.getNx()
        jy = n.getYmin() + dy*(j+0.5)/n.getNy()

        #Nv = n.number_of_virtual_neighbors(c)
        Nv = c.number_of_virtual_neighbors
        label = str(Nv)
        #label = "{} ({},{})/{}".format(cid,i,j,Nv)
        #label = "({},{})".format(i,j)
        ax.text(ix, jy, label, ha='center',va='center', size=8)

    #for cid in n.getVirtuals():
    #    c = n.getCell( cid )
    #    (i,j) = c.index()
    #    ix = n.getXmin() + n.getXmax()*(i+0.5)/n.getNx()
    #    jy = n.getYmin() + n.getYmin()*(j+0.5)/n.getNy()
    #    label = "Vir"
    #    ax.text(jy, ix, label, ha='center',va='center')

    #XXX add back
    #ax.set_title(str(len(n.getVirtuals() ))+"/"+str(len(n.getCellIds() )))


def plotMesh(ax, n, conf):

    NxFull = conf["Nx"] * conf["NxMesh"]
    NyFull = conf["Ny"] * conf["NyMesh"]

    NxMesh = conf["NxMesh"]
    NyMesh = conf["NyMesh"]

    data = -1.0 * np.ones( (NxFull, NyFull) )

    for cid in n.getCellIds():
        c = n.getCellPtr( cid )
        (i, j) = c.index()

        mesh = c.getData()

        for k in range(NyMesh):
            for q in range(NxMesh):
                data[ i*NxMesh + q, j*NyMesh + k ] = mesh[q,k]


    imshow(ax, data,
            n.getXmin(), n.getXmax(), n.getYmin(), n.getYmax(),
            cmap = palette,
            vmin = 0.0,
            #vmax = data.max(),
            clip = 0,
            )

    #print "old mesh\n:"
    #print np.flipud(data.T)




def plotMesh2(ax, n, conf):

    NxFull = conf["Nx"] * conf["NxMesh"]
    NyFull = conf["Ny"] * conf["NyMesh"]

    NxMesh = conf["NxMesh"]
    NyMesh = conf["NyMesh"]

    data = -1.0 * np.ones( (NxFull, NyFull) )

    cid = n.cellId(1,1)
    c = n.getCellPtr( cid )
    (i, j) = c.index()

    mesh = c.getData()

    for k in range(-1, NyMesh+1, 1):
        for q in range(-1, NxMesh+1, 1):
            data[ i*NxMesh + q, j*NyMesh + k ] = mesh[q,k]


    imshow(ax, data,
            n.getXmin(), n.getXmax(), n.getYmin(), n.getYmax(),
            cmap = palette,
            vmin = 0.0,
            #vmax = data.max(),
            clip = 0,
            )

    #print "new mesh\n:"
    #print np.flipud(data.T)


def saveVisz(lap, n, conf):

    slap = str(lap).rjust(4, '0')
    fname = conf["dir"] + '/node_{}_{}.png'.format(n.rank, slap)
    plt.savefig(fname)




#make random starting order
def loadMpiRandomly(n):
    np.random.seed(4)
    if n.master:
        for i in range(n.getNx()):
            for j in range(n.getNy()):
                val = np.random.randint(n.Nrank)
                n.setMpiGrid(i, j, val)

#load nodes to be in stripe formation (splitted in X=horizontal direction)
def loadMpiXStrides(n):
    if n.master: #only master initializes; then sends
        stride = np.zeros( (n.getNx()), np.int64)
        dx = np.float(n.getNx()) / np.float(n.Nrank) 
        for i in range(n.getNx()):
            val = np.int( i/dx )
            stride[i] = val

        for j in range(n.getNy()):
            for i in range(n.getNx()):
                val = stride[i]
                n.setMpiGrid(i, j, val)
    n.bcastMpiGrid()


#load cells into each node
def loadCells(n):
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            #print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.getMpiGrid(i,j), ref[j,i]))

            if n.getMpiGrid(i,j) == n.rank:
                c = pygol.CellularAutomataCell(i, j, n.rank, n.getNx(), n.getNy())
                n.addCell(c) 



def randomInitialize(n, conf):

    val = 0
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            #if n.getMpiGrid(i,j) == n.rank:
            if True:
                cid = n.cellId(i,j)
                c = n.getCellPtr(cid) #get cell ptr

                mesh = pygol.Mesh( conf["NxMesh"], conf["NyMesh"] )

                # fill mesh
                if (i == 2) and (j == 2):
                    for q in range(conf["NxMesh"]):
                        for k in range(conf["NyMesh"]):
                            ref = np.random.randint(0,11)
                            if ref < 5:
                                mesh[q,k] = 1
                            else:
                                mesh[q,k] = 0



                        #mesh[q,k] = q + conf["NxMesh"]*k
                        #mesh[q,k] = val
                        val += 1

                #mesh[0,0] = 1

                c.addData(mesh)
                c.addData(mesh)



##################################################
##################################################

if __name__ == "__main__":

    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(4,8))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(2, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )
    
    
    
    #setup node
    conf = {
            "Nx"     : 5,
            "Ny"     : 5,
            "NxMesh" : 100,
            "NyMesh" : 100,
            "dir"    : "out",
            "Nrank"  : 1
            }
    
    node = pygol.Grid( conf["Nx"], conf["Ny"] ) 
    node.setGridLims(0.0, 1.0, 0.0, 1.0)
    
    loadCells(node)
    randomInitialize(node, conf)
    
    
    
    # Path to be created 
    if node.master:
        if not os.path.exists( conf["dir"]):
            os.makedirs(conf["dir"])
    
    sol = pygol.Solver()
    
    
    plotNode(axs[0], node, conf)
    plotMesh(axs[1], node, conf)
    saveVisz(0, node, conf)
    


    for lap in range(1, 1000):
        print "---lap: {}".format(lap)
    
        #update halo regions
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                c = node.getCellPtr(i,j) #get cell ptr
                c.updateBoundaries(node)

        #progress one time step
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                c = node.getCellPtr(i,j) #get cell ptr
                sol.solve(c)

        #cycle everybody in time
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                c = node.getCellPtr(i,j) #get cell ptr
                c.cycle()


        if (lap % 10 == 0):
            plotNode(axs[0], node, conf)
            plotMesh(axs[1], node, conf)
            saveVisz(lap, node, conf)
        

    #plotNode(axs[0], node, conf)
    #plotMesh2(axs[1], node, conf)
    #saveVisz(1, node, conf)
    
    
    
    
