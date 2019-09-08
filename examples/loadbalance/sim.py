from mpi4py import MPI

import numpy as np
import os, sys
import h5py
import argparse

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

# Visualize current cell ownership on grid
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
                tmp_grid2[i,j] = grid.get_mpi_grid(i,j)
                #print("{}: val={}".format(grid.rank(), grid.get_mpi_grid(i,j)))

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
            tmp[i,j] = grid.get_mpi_grid(i,j)
    return tmp

def get_work_grid(n, conf):

    tmp = np.ones( (n.get_Nx(), n.get_Ny()) ) * -1.0
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            tmp[i,j] = grid.get_work_grid(i,j)
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

    work = get_work_grid(n, conf)
    f5['work'][:,:,lap] = work


def plotWork(ax, n, conf):
    tmp_grid = np.ones( (n.get_Nx(), n.get_Ny()) ) * -1.0
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            tmp_grid[i,j] = grid.get_work_grid(i,j)

    #print("{}: min/max work {}/{}".format(n.rank(), np.min(tmp_grid), np.max(tmp_grid)))

    maxv = np.max(tmp_grid)

    imshow(ax, tmp_grid, 
           n.get_xmin(), n.get_xmax(), n.get_ymin(), n.get_ymax(),
           cmap = 'plasma',
           vmin = 0.0,
           vmax = maxv,
           )


#load tiles into each grid
def load_tiles(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            #print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.get_mpi_grid(i,j), ref[j,i]))

            if n.get_mpi_grid(i,j) == n.rank():
                c = pycorgi.Tile()
                n.add_tile(c, (i,j)) 

    n.update_work()


#inject plasma into (individual) cells
def inject(grid, ffunc, conf):

    #loop over all *local* cells
    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):
            if grid.get_mpi_grid(i,j) == grid.rank():
                #print("creating ({},{})".format(i,j))

                #get cell & its content
                cid    = grid.id(i,j)
                c      = grid.get_tile(cid) #get cell ptr

                # inject particles
                for ispcs in range(conf.Nspecies):
                    container = c.get_container(ispcs)

                    for n in range(conf.NzMesh):
                        for m in range(conf.NyMesh):
                            for l in range(conf.NxMesh):
                                #print(" sub mesh: ({},{},{})".format(l,m,n))
                                xloc = spatialLoc(grid, (i,j), (l,m,n), conf)

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


def add_virtual_work(n, lap, conf):

    ic = 0.0
    #ic = n.get_Nx()/2.0
    jc = n.get_Ny()/2.0

    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            oldw = grid.get_work_grid(i,j)

            #rvec2 = (i-ic)**2.0 + (j-jc)**2.0
            rvec2 = (i-ic)**2.0 #+ (j-jc)**2.0
            sig = 100.0
            nw = 4.0*np.exp(-rvec2/sig)

            grid.set_work_grid(i,j, nw)

            #grid.set_work_grid(i,j, 1.0)




class Conf:

    Nx  = 50
    Ny  = 50
    Nz  = 1

    NxMesh = 1
    NyMesh = 1
    NzMesh = 1

    def __init__(self):
        self.outdir = 'out2d_'+str(self.Nx)+'x'+str(self.Ny)+'n'+str(1)

    def __init__(self, Nx, Ny, Nz, Nr):
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz

        self.Nr = Nr

        self.outdir = 'out2d_'+str(Nx)+'x'+str(Ny)+'n'+str(Nr)


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


    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--Nx', dest='Nx', type=int, default=40)
    parser.add_argument('--Ny', dest='Ny', type=int, default=40)
    parser.add_argument('--Nz', dest='Nz', type=int, default=1)
    parser.add_argument('--Nt', dest='Nt', type=int, default=50)
    parser.add_argument('--Nr', dest='Nr', type=int, default=50)
    args = parser.parse_args()


    do_plots = True

    # set up plotting and figure
    if do_plots:
        plt.fig = plt.figure(1, figsize=(12,4))
        plt.rc('font', family='serif', size=12)
        plt.rc('xtick')
        plt.rc('ytick')
        
        gs = plt.GridSpec(1, 3)
        gs.update(hspace = 0.5)
        
        axs = []
        axs.append( plt.subplot(gs[0]) )
        axs.append( plt.subplot(gs[1]) )
        axs.append( plt.subplot(gs[2]) )
    

    #setup grid
    conf = Conf(args.Nx, args.Ny, args.Nz, args.Nr)
    conf.update_bbox()
    
    grid = pycorgi.Grid( conf.Nx, conf.Ny ) 
    grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
    
    loadMpiRandomly(grid)
    #loadMpiXStrides(grid)

    load_tiles(grid, conf)

    # Path to be created 
    if grid.master:
        if not os.path.exists( conf.outdir):
            os.makedirs(conf.outdir)
    MPI.COMM_WORLD.barrier()
    

    ################################################## 
    rank = str(grid.rank())
    f5 = h5py.File(conf.outdir+"/run-"+rank+".h5", "w")

    Nsamples = args.Nt + 1
    f5.create_dataset("virtuals",   (Nsamples,), dtype='f')
    f5.create_dataset("locals",     (Nsamples,), dtype='f')
    f5.create_dataset("boundaries", (Nsamples,), dtype='f')
    f5.create_dataset("grid", (conf.Nx, conf.Ny, Nsamples), dtype='f')
    f5.create_dataset("work", (conf.Nx, conf.Ny, Nsamples), dtype='f')

    ################################################## 

    if do_plots:
        plotNode(axs[0], grid, conf)
        plotNode(axs[1], grid, conf, mpigrid=True)
        plotWork(axs[2], grid, conf)
        saveVisz(-1, grid, conf)
    
    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    initialize_virtuals(grid, conf)
    grid.allgather_work_grid()

    #add_virtual_work(grid, 0, conf)

    analyze(grid, f5, 0, conf)
    if do_plots:
        plotNode(axs[0], grid, conf)
        plotNode(axs[1], grid, conf, mpigrid=True)
        plotWork(axs[2], grid, conf)
        saveVisz(0, grid, conf)


    ##################################################
    for lap in range(1, Nsamples):
        print("---lap: {}".format(lap))

        #add_virtual_work(grid, lap, conf)
        #grid.update_work()
        #grid.allgather_work_grid()
        MPI.COMM_WORLD.barrier()

        if False:
            # corgi loadbalance 
            #print("adoption_council")
            grid.adoption_council()
            #print("adopt")
            grid.adopt()
            #print("communicate_adoptions")
            grid.communicate_adoptions()
            #print("erase_virtuals")
            grid.erase_virtuals()
        else:
            print("adoption_council2")
            grid.adoption_council2()
            print("erase_virtuals")
            grid.erase_virtuals()


        print("analyze_boundaries")
        grid.analyze_boundaries()

        print("send_tiles")
        grid.send_tiles()
        print("recv_tiles")
        grid.recv_tiles()

        print("initialize")
        initialize_virtuals(grid, conf)

        analyze(grid, f5, lap, conf)

        #if (lap % 20 == 0):
        if True:
            if do_plots:
                print("visualizing...")
                plotNode(axs[0], grid, conf)
                plotNode(axs[1], grid, conf, mpigrid=True)
                plotWork(axs[2], grid, conf)
                saveVisz(lap, grid, conf)
            f5.flush()

    

    f5.close()

    #plotNode(axs[0], grid, conf)
    #plotNode(axs[1], grid, conf, mpigrid=True)
    #plotWork(axs[2], grid, conf)

    #saveVisz(, grid, conf)


