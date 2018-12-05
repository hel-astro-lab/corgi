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
import pyprtcls


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
def plotNode(ax, n, conf):
    tmp_grid = np.ones( (n.get_Nx(), n.get_Ny()) ) * -1.0

    for cid in n.get_tile_ids():
        c = n.get_tile( cid )
        (i, j) = c.index
        #check dublicates
        if tmp_grid[i,j] != -1.0:
            print("{}: ERROR in real cells at ({},{})".format(n.rank, i,j))
            sys.exit()
        tmp_grid[i,j] = c.communication.owner

    imshow(ax, tmp_grid, 
            n.get_xmin(), n.get_xmax(), n.get_ymin(), n.get_ymax(),
            cmap = palette,
            vmin = 0.0,
            vmax = n.size(),
            )

    # add text label about number of neighbors
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





class Particles:
    xs  = []
    ys  = []
    zs  = []

    uxs = []
    uys = []
    uzs = []

    wgt = []

    def clear(self):
        self.xs  = []
        self.ys  = []
        self.zs  = []

        self.uxs = []
        self.uys = []
        self.uzs = []

        self.wgt = []

def get_particles(node, conf, ip):
    prtcl = Particles()
    prtcl.clear()

    for i in range(conf.Nx):
        for j in range(conf.Ny):
            for k in range(conf.Nz):
                if node.get_mpi_grid(i,j) == node.rank():
                    cid = node.id(i,j)
                    c = node.get_tile(cid)

                    x, y, z, ux, uy, uz, wgt = get_particles_from_tile(c, ip)

                    prtcl.xs.extend(x)
                    prtcl.ys.extend(y)
                    prtcl.zs.extend(z)

                    prtcl.uxs.extend(ux)
                    prtcl.uys.extend(uy)
                    prtcl.uzs.extend(uz)

                    prtcl.wgt.extend(wgt)

    return prtcl


def get_particles_from_tile(tile, ispcs):
    container = tile.get_container(ispcs)
    x  = container.loc(0)
    y  = container.loc(1)
    z  = container.loc(2)

    ux = container.vel(0)
    uy = container.vel(1)
    uz = container.vel(2)

    wgt = container.wgt()

    return x, y, z, ux, uy, uz, wgt



def plotMesh(ax, n, conf, downsample=0):

    #ax.clear()
    #ax.cla()
    plotNode(ax, n, conf)
    #plotTileBoundaries(ax, n, conf)

    prtcl = get_particles(n, conf, 0)
    Np = len(prtcl.xs)
    print("particles to plot: {}".format(Np))

    if downsample > 0:
        rindxs = random.sample( range(0, Np-1), int(downsample*Np) )

        prtcl.xs = np.array( prtcl.xs )
        prtcl.ys = np.array( prtcl.ys ) 
        prtcl.zs = np.array( prtcl.zs )

        prtcl.xs = prtcl.xs[rindxs]
        prtcl.ys = prtcl.ys[rindxs]
        prtcl.zs = prtcl.zs[rindxs]

    ax.plot(prtcl.xs, prtcl.ys, ".", color='red')
    
    #second container
    prtcl1 = get_particles(n, conf, 1)
    ax.plot(prtcl1.xs, prtcl1.ys, ".", color='k')



def spatialLoc(node, Ncoords, Mcoords, conf):

    #node coordinates
    i, j    = Ncoords 
    Nx      = conf.Nx
    Ny      = conf.Ny

    #mesh coordinates
    l, m, n = Mcoords 
    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh
    NzMesh = conf.NzMesh

    #grid spacing
    xmin = node.get_xmin()
    ymin = node.get_ymin()

    dx = 1.0 #conf.dx
    dy = 1.0 #conf.dy
    dz = 1.0 #conf.dz


    #calculate coordinate extent
    x = xmin + i*(NxMesh)*dx + l*dx
    y = ymin + j*(NyMesh)*dy + m*dy
    z = 0.0                  + n*dz

    return [x, y, z]



def initialize_tile(c, i, j, n, conf):
    ppc = conf.ppc #/ conf.Nspecies
    
    # load particle containers
    for sps in range(conf.Nspecies):
        container = pyprtcls.ParticleBlock(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        
        #reserve memory for particles
        Nprtcls = conf.NxMesh*conf.NyMesh*conf.NzMesh*conf.ppc
        container.reserve(Nprtcls)

        c.set_container( container )
    
    #set bounding box of the tile 
    mins = spatialLoc(n, [i,j], [0,0,0], conf)
    maxs = spatialLoc(n, [i,j], [conf.NxMesh, conf.NyMesh, conf.NzMesh], conf)
    c.set_tile_mins(mins[0:2])
    c.set_tile_maxs(maxs[0:2])


#load tiles into each node
def load_tiles(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            #print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.get_mpi_grid(i,j), ref[j,i]))

            if n.get_mpi_grid(i,j) == n.rank():
                c = pyprtcls.Tile()

                initialize_tile(c, i, j, n, conf)

                #add it to the node
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

    for cid in n.get_virtuals():
        c_orig = n.get_tile(cid)
        (i,j) = c_orig.index

        # new prtcl tile;
        # TODO: load_metainfo *HAS* to be after add_tile
        c = pyprtcls.Tile()
        n.add_tile(c, (i,j)) 

        c_orig.communication.local = False;
        c.load_metainfo(c_orig.communication)
        print("{}: loading {} owned by {}".format(n.rank(), cid, c.communication.owner))
        
        initialize_tile(c, i,j,n, conf)


class Conf:

    Nx  = 10
    Ny  = 10
    Nz  = 1

    NxMesh = 1
    NyMesh = 1
    NzMesh = 1

    Nspecies = 2
    ppc = 1

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
    
    node = pycorgi.twoD.Node( conf.Nx, conf.Ny ) 
    node.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
    
    #loadMpiRandomly(node)
    loadMpiXStrides(node)

    load_tiles(node, conf)
    inject(node, filler, conf)
    
    # Path to be created 
    if node.master:
        if not os.path.exists( conf.outdir):
            os.makedirs(conf.outdir)
    
    pusher = pyprtcls.Pusher()
    
    #static setup; communicate neighbor info once
    node.analyze_boundaries()
    node.send_tiles()
    node.recv_tiles()
    initialize_virtuals(node, conf)

    plotNode(axs[0], node, conf)
    plotMesh(axs[1], node, conf)
    saveVisz(0, node, conf)


    node.send_data(0)
    node.recv_data(0)

    for lap in range(1, 101):
        print("---lap: {}".format(lap))

        #send/recv boundaries
        node.send_data(0)
        node.recv_data(0)

        if (lap % 1 == 0):
            plotNode(axs[0], node, conf)
            plotMesh(axs[1], node, conf)
            saveVisz(lap, node, conf)
    
        #move particles
        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            pusher.solve(tile)

        #communicate particles
        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            tile.check_outgoing_particles()

        #mpi communication

        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            tile.get_incoming_particles(node)

        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            tile.delete_transferred_particles()



