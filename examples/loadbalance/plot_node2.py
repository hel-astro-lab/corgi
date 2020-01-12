import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py
import sys, os
import matplotlib as mpl
import matplotlib.ticker as ticker
from scipy.stats import mstats
from scipy.optimize import curve_fit

from visualize import imshow

from configSetup import Configuration
from combine_files import get_file_list
from combine_files import combine_tiles

from scipy.ndimage.filters import gaussian_filter

from parser import parse_input
import argparse

import pycorgi.twoD as corgi
import initialize as init

import pyrunko.tools.twoD as pytools
import pyrunko.vlv.twoD as pyvlv
import pyrunko.pic.twoD as pypic
import pyrunko.fields.twoD as pyfld


# trick to make nice colorbars
# see http://joseph-long.com/writing/colorbars/
def colorbar(mappable, 
        loc="right", 
        orientation="vertical", 
        size="1%", 
        #pad=0.05, 
        pad=0.1, 
        ticklocation='right'):
        #loc="top", 
        #orientation="horizontal", 
        #size="8%", 
        #pad=0.03, 
        #ticklocation='top'):

    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(loc, size=size, pad=pad)
    return cax, fig.colorbar(mappable, cax=cax, orientation=orientation, ticklocation=ticklocation)


default_values = {
    #"cmap':"YlOrRd",
    #'cmap':"OrRd",
    'cmap':"viridis",
    'vmin': None,
    'vmax': None,
    'clip': None,
    'aspect':1,
    'vsymmetric':None,
    'winsorize_min':0.005,
    'winsorize_max':0.005,
    'title':'rank',
    'derived':False,
}

#load tiles into each grid
def loadTiles(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            ow = n.get_mpi_grid(i,j)
            c = corgi.Tile()
            c.communication.owner=ow
            n.add_tile(c, (i,j)) 
            n.set_mpi_grid(i,j,ow)
            #print("{} ({},{}) {} =? {}".format(n.rank(), i,j, n.get_mpi_grid(i,j),ow))

            #if n.get_mpi_grid(i,j) == n.rank():
                #c = pyrunko.vlv.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                #c = corgi.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

                #initialize_tile(c, i, j, n, conf)
                #add it to the grid

    for cid in n.get_tile_ids():
        c = n.get_tile( cid )
        (i, j) = c.index
        ow = n.get_mpi_grid(i,j)
        c.communication.owner=ow


def build_node_image(grid, conf):

    xmin = 0.0
    xmax = conf.Nx*conf.NxMesh #XXX scaled length
    ymin = 0.0
    ymax = conf.Ny*conf.NyMesh

    grid.set_grid_lims(xmin, xmax, ymin, ymax)

    init.loadMpi2D(grid, 200)
    loadTiles(grid, conf)
    grid.analyze_boundaries()

    #--------------------------------------------------
    tmp_grid = np.ones( (grid.get_Nx(), grid.get_Ny()) ) * -1.0
    #for cid in grid.get_tile_ids():
    #    c = grid.get_tile( cid )
    #    (i, j) = c.index

    #    #check dublicates
    #    if tmp_grid[i,j] != -1.0:
    #        print("{}: ERROR in real cells at ({},{})".format(grid.rank, i,j))
    #        sys.exit()
    #    tmp_grid[i,j] = c.communication.owner

    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):
            #val = igrid[i,j] 
            val = grid.get_mpi_grid(i,j)
            #n.set_mpi_grid(i, j, val)
            tmp_grid[i,j] = val

    return tmp_grid


def add_text_labels(ax, n, conf):

    # add text label about number of neighbors
    dx = n.get_xmax() - n.get_xmin()
    dy = n.get_ymax() - n.get_ymin()
    dx /= 10
    dy /= 10
    print(dx)
    print(dy)

    scalef = 2.0
    dx *= scalef
    dy *= scalef
    for cid in n.get_tile_ids():
        c = n.get_tile( cid )
        (i, j) = c.index

        ix = n.get_xmin() + dx*(i+0.5)/n.get_Nx()
        jy = n.get_ymin() + dy*(j+0.5)/n.get_Ny()
        #print(ix,jy,i,j)

        Nv = c.communication.number_of_virtual_neighbors
        if Nv == 0:
            continue
        if Nv == 8:
            continue

        label = str(Nv)
        ax.text(ix, jy, label, ha='center',va='center', size=5)




def plot2d_single(
        ax, 
        var,
        info, 
        conf,
        title= None,
        vmin = None,
        vmax = None,
        cmap = None,
        clip = None,
        ):

    #--------------------------------------------------
    # unpack incoming arguments that modify defaults
    args = {}

    #general defaults
    for key in default_values:
        args[key] = default_values[key]

    #TODO: read grid into val
    if True:
        grid = corgi.Grid(conf.Nx, conf.Ny, conf.Nz)
        val = build_node_image(grid, conf)
    else:
        outdir = "out500x500n1000"
        f5all = h5py.File(outdir+"/run-merged.h5", "r")
        val = f5all['grid'][:,:,:]
        val = val[:,:,180]
        nranks = 1000
        val = val / nranks
        val[0,0] = 1023

    nx, ny = np.shape(val)

    modv = 31
    args['vmin'] = np.min(val)
    args['vmax'] = np.max(val)
    #args['vmax'] = modv

    print(val)
    #val = np.mod(val, 20*np.ones(np.shape(val)))
    #val = np.mod(val, modv)   
    print("new")
    print(val)
    #val = np.reshape(np.random.shuffle(val),(nx,ny))

    xmin = 0.0
    ymin = 0.0
    #xmax = nx/info['skindepth']
    #ymax = ny/info['skindepth']
    xmax = nx*2
    ymax = ny*2

    xmax = ymax = 1024.0

    #cyc = plt.cycler(
    #        s=np.linspace(200, 50, 3),
    #        cmap=['viridis', 'magma', 'coolwarm'],
    #        alpha=[.25, .5, .75],
    #        lw=[0, .1, .5])

    #create cycling colormap array
    carr = []
    cmap_user = mpl.cm.get_cmap(args['cmap'])
    norm = mpl.colors.Normalize(vmin=0, vmax=modv)

    print(np.max(val))

    #for c in np.arange(0, np.max(val)+1):
    for c in np.arange(0, np.max(val)+1):
        v = np.mod(c, modv)
        v = norm(v)
        print(c, v, cmap_user(v))
        carr.append(cmap_user(v))
    cmap_cycl = mpl.colors.ListedColormap(carr)


    im = imshow(ax, val, xmin, xmax, ymin, ymax,
           #cmap = args['cmap'],
           cmap = cmap_cycl,
           vmin = args['vmin'],
           vmax = args['vmax'],
           clip = args['clip'],
           aspect=args['aspect'],
           )

    if False:
        ax.set_xlim((600.0, 900.0))
        ax.set_ylim((300.0, 600.0))

        #ax.set_xlim((0.0, 100.0))
        #ax.set_ylim((0.0, 100.0))

        #add_text_labels(ax, grid, conf)


    ax.set_xlabel(r"$x$ $(c/\omega_p)$")
    ax.set_ylabel(r"$y$ $(c/\omega_p)$")
    plt.subplots_adjust(left=0.15, bottom=0.10, right=0.87, top=0.97)


    wskip = 0.2
    pad = 0.01
    pos = ax.get_position()


    axleft   = pos.x0
    axbottom = pos.y0
    axright  = pos.x0 + pos.width
    axtop    = pos.y0 + pos.height

    cax = plt.fig.add_axes([axright+pad, axbottom+wskip, 0.01, axtop-axbottom-2*wskip])

    #norm = plt.Normalize()
    #colors = plt.cm.jet(norm(dz))

    #if False:
    #    #create cycling colormap array
    #    carr = []
    #    for c in np.arange(0, np.max(val), np.max(val)):
    #        carr.append(np.mod(c, modv))
    #    cmap = mpl.colors.ListedColormap(carr)
    #    cb = plt.fig.colorbar(
    #            im, 
    #            cax=cax, 
    #            cmap = cmap,
    #            orientation='vertical',
    #            ticklocation='right')

    #else:
    cb = plt.fig.colorbar(
            im, 
            cax=cax, 
            orientation='vertical',
            ticklocation='right')

    cax.text(1.0, 1.03, args['title'], transform=cax.transAxes)

    slap = "0"
    fname = var+'_{}.pdf'.format(slap)
    plt.savefig(fname)
    cb.remove()


if __name__ == "__main__":

    plt.fig = plt.figure(1, figsize=(4,3.5), dpi=200)
    plt.rc('font',  family='sans-serif')
    #plt.rc('text',  usetex=True)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)

    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.0)
    gs.update(wspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0,0]) )


    conf, fdir, args = parse_input()

    info = {}
    #info = quick_build_info(fdir, args.lap)
    info['skindepth'] = conf.c_omp/conf.stride

    fdir += '/'
    print("plotting {}".format(fdir))

    args.var = "rank"
    plot2d_single(axs[0], args.var, info, conf)







