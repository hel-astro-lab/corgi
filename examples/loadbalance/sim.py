from mpi4py import MPI

import numpy as np
import os, sys
import h5py

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import palettable as pal
    palette = pal.wesanderson.Moonrise1_5.mpl_colormap
except:
    pass

import pycorgi.twoD as pycorgi



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
              #alpha=0.5
              )
 

def saveVisz(lap, n, conf):
    slap = str(lap).rjust(4, '0')
    fname = conf.outdir + '/node_{}_{}.png'.format(n.rank(), slap)
    plt.savefig(fname)


#make random starting order
def loadMpiRandomly(n):
    np.random.seed(0)
    if n.master:
        for i in range(n.get_Nx()):
            for j in range(n.get_Ny()):
                val = np.random.randint( n.size() )
                n.set_mpi_grid(i, j, val)
    n.bcast_mpi_grid()

#load nodes to be in stripe formation (splitted in X=horizontal direction)
def loadMpiXStrides(n):
    if n.master: #only master initializes; then sends
        stride = np.zeros( (n.get_Nx()), np.int64)
        dx = np.float(n.get_Nx()) / np.float(n.size() ) 
        for i in range(n.get_Nx()):
            val = np.int( i/dx )
            stride[i] = val

        for j in range(n.get_Ny()):
            for i in range(n.get_Nx()):
                val = stride[i]
                n.set_mpi_grid(i, j, val)
    n.bcast_mpi_grid()

# Visualize current cell ownership on node
def plotNode(ax, n, conf, mpigrid=False):

    tmp_grid = np.ones( (n.get_Nx(), n.get_Ny()) ) * -1.0
    for cid in n.get_tile_ids():
        c = n.get_tile( cid )
        (i, j) = c.index

        #check dublicates
        if tmp_grid[i,j] != -1.0:
            print("{}: ERROR in real cells at ({},{})".format(n.rank(), i,j))
            sys.exit()
        tmp_grid[i,j] = c.communication.owner

    if not(mpigrid):
        imshow(ax, tmp_grid, 
               n.get_xmin(), n.get_xmax(), n.get_ymin(), n.get_ymax(),
               cmap = palette,
               vmin = 0.0,
               vmax = n.size(),
               )

    # internal mpi_grid
    if mpigrid:
        tmp_grid2 = np.ones( (n.get_Nx(), n.get_Ny()) ) * -1.0
        for i in range(n.get_Nx()):
            for j in range(n.get_Ny()):
                tmp_grid2[i,j] = node.get_mpi_grid(i,j)
                #print("{}: val={}".format(node.rank(), node.get_mpi_grid(i,j)))

        #for i in range(n.get_Nx()):
        #    for j in range(n.get_Ny()):
        #        if(tmp_grid2[i,j] != tmp_grid[i,j]):
        #            print(" mismatch {} vs {} at ({},{})".format(tmp_grid[i,j], tmp_grid2[i,j], i,j))

        imshow(ax, tmp_grid2, 
           n.get_xmin(), n.get_xmax(), n.get_ymin(), n.get_ymax(),
           cmap = palette,
           vmin = 0.0,
           vmax = n.size(),
           )


    # add text label about number of neighbors
    if True:
        if not(mpigrid):
            for cid in n.get_tile_ids():
                c = n.get_tile( cid )
                (i, j) = c.index
                dx = n.get_xmax() - n.get_xmin()
                dy = n.get_ymax() - n.get_ymin()

                ix = n.get_xmin() + dx*(i+0.5)/n.get_Nx()
                jy = n.get_ymin() + dy*(j+0.5)/n.get_Ny()

                Nv = c.communication.number_of_virtual_neighbors
                label = str(Nv)
                ax.text(ix, jy, label, ha='center',va='center', size=8)

    #mark boundaries with hatch
    if False:
        dx = n.get_xmax() - n.get_xmin()
        dy = n.get_ymax() - n.get_ymin()
        for cid in n.get_boundary_tiles():
            c = n.get_tile( cid )
            (i, j) = c.index

            ix0 = n.get_xmin() + dx*(i+0.0)/n.get_Nx()
            jy0 = n.get_ymin() + dy*(j+0.0)/n.get_Ny()

            ix1 = n.get_xmin() + dx*(i+1.0)/n.get_Nx()
            jy1 = n.get_ymin() + dy*(j+1.0)/n.get_Ny()

            #ax.fill_between([ix0,ix1], [jy0, jy1], hatch='///', alpha=0.0)

            ax.plot([ix0, ix0],[jy0, jy1], color='k', linestyle='dotted')
            ax.plot([ix1, ix1],[jy0, jy1], color='k', linestyle='dotted')
            ax.plot([ix0, ix1],[jy0, jy0], color='k', linestyle='dotted')
            ax.plot([ix0, ix1],[jy1, jy1], color='k', linestyle='dotted')

    if True:
        virs = n.get_virtual_tiles()
        boun = n.get_boundary_tiles()
        locs = n.get_local_tiles()
        ax.set_title(str(len(virs))+"/"+str(len(boun))+"/"+str(len(locs)))

def get_mpi_grid(n, conf):

    tmp = np.ones( (n.get_Nx(), n.get_Ny()) ) * -1.0
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            tmp[i,j] = node.get_mpi_grid(i,j)
    return tmp


def analyze(n, f5, lap, conf):

    virs = len(n.get_virtual_tiles()  )
    boun = len(n.get_boundary_tiles() )
    locs = len(n.get_local_tiles()    )

    f5['virtuals'][lap]   = virs
    f5['boundaries'][lap] = boun
    f5['locals'][lap]     = locs

    grid = get_mpi_grid(n, conf)
    f5['grid'][:,:,lap] = grid



#load tiles into each node
def load_tiles(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            #print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.get_mpi_grid(i,j), ref[j,i]))

            if n.get_mpi_grid(i,j) == n.rank():
                c = pycorgi.Tile()
                n.add_tile(c, (i,j)) 


#inject plasma into (individual) cells
def inject(node, ffunc, conf):

    #loop over all *local* cells
    for i in range(node.get_Nx()):
        for j in range(node.get_Ny()):
            if node.get_mpi_grid(i,j) == node.rank():
                #print("creating ({},{})".format(i,j))

                #get cell & its content
                cid    = node.id(i,j)
                c      = node.get_tile(cid) #get cell ptr

                # inject particles
                for ispcs in range(conf.Nspecies):
                    container = c.get_container(ispcs)

                    for n in range(conf.NzMesh):
                        for m in range(conf.NyMesh):
                            for l in range(conf.NxMesh):
                                #print(" sub mesh: ({},{},{})".format(l,m,n))
                                xloc = spatialLoc(node, (i,j), (l,m,n), conf)

                                for ip in range(conf.ppc):
                                    x0, u0 = ffunc(xloc, ispcs, conf)
                                    container.add_particle(x0, u0, 1.0)


def filler(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    #zz = xloc[2] + np.random.rand(1)
    zz = 0.0

    if ispcs == 0:
        vel = 0.1
    if ispcs == 1:
        vel = 0.5

    ang = 2.0*np.pi*np.random.rand(1)
    ux = vel*np.sin(ang)
    uy = vel*np.cos(ang)
    uz = vel*0.0

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0




# make all tiles same type 
def initialize_virtuals(n, conf):

    for cid in n.get_virtual_tiles():
        c_orig = n.get_tile(cid)
        (i,j) = c_orig.index

        # new prtcl tile;
        # TODO: load_metainfo *HAS* to be after add_tile
        c = pycorgi.Tile()
        n.replace_tile(c, (i,j)) 

        #c.communication.local = False;
        #c.load_metainfo(c_orig.communication)
        #print("{}: loading {} owned by {}".format(n.rank(), cid, c.communication.owner))
        
        #initialize_tile(c, i,j,n, conf)


class Conf:

    Nx  = 30
    Ny  = 30
    Nz  = 1

    NxMesh = 1
    NyMesh = 1
    NzMesh = 1

    outdir = "out"

    def update_bbox(self):
        self.xmin = 0.0
        self.xmax = self.Nx*self.NxMesh

        self.ymin = 0.0
        self.ymax = self.Ny*self.NyMesh

        self.zmin = 0.0
        self.zmax = self.Nz*self.NzMesh



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
    

    #setup node
    conf = Conf()
    conf.update_bbox()
    
    node = pycorgi.Node( conf.Nx, conf.Ny ) 
    node.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
    
    loadMpiRandomly(node)
    #loadMpiXStrides(node)

    load_tiles(node, conf)

    # Path to be created 
    if node.master:
        if not os.path.exists( conf.outdir):
            os.makedirs(conf.outdir)


    ################################################## 
    rank = str(node.rank())
    f5 = h5py.File(conf.outdir+"/run-"+rank+".h5", "w")

    Nsamples = 51
    f5.create_dataset("virtuals",   (Nsamples,), dtype='f')
    f5.create_dataset("locals",     (Nsamples,), dtype='f')
    f5.create_dataset("boundaries", (Nsamples,), dtype='f')
    f5.create_dataset("grid", (conf.Nx, conf.Ny, Nsamples), dtype='f')

    ################################################## 

    do_plots = True
    if do_plots:
        plotNode(axs[0], node, conf)
        plotNode(axs[1], node, conf, mpigrid=True)
        saveVisz(-1, node, conf)
    
    node.analyze_boundaries()
    node.send_tiles()
    node.recv_tiles()
    initialize_virtuals(node, conf)

    analyze(node, f5, 0, conf)

    if do_plots:
        plotNode(axs[0], node, conf)
        plotNode(axs[1], node, conf, mpigrid=True)
        saveVisz(0, node, conf)

    ##################################################
    for lap in range(1, Nsamples):
        print("---lap: {}".format(lap))

        if False:
            # corgi loadbalance 
            #print("adoption_council")
            node.adoption_council()
            #print("adopt")
            node.adopt()
            #print("communicate_adoptions")
            node.communicate_adoptions()
            #print("erase_virtuals")
            node.erase_virtuals()
        else:
            #print("adoption_council2")
            node.adoption_council2()
            #print("erase_virtuals")
            node.erase_virtuals()


        print("analyze_boundaries")
        node.analyze_boundaries()
        print("send_tiles")
        node.send_tiles()
        print("recv_tiles")
        node.recv_tiles()
        print("initialize")
        initialize_virtuals(node, conf)

        if (lap % 1 == 0):
            if do_plots:
                plotNode(axs[0], node, conf)
                plotNode(axs[1], node, conf, mpigrid=True)
                saveVisz(lap, node, conf)
    
        analyze(node, f5, lap, conf)

    f5.close()

    #plotNode(axs[0], node, conf)
    #plotNode(axs[1], node, conf, mpigrid=True)
    #saveVisz(, node, conf)


