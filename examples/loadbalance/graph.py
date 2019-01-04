import numpy as np
import h5py
import os

import palettable as pal
palette = pal.wesanderson.Moonrise1_5.mpl_colormap

import networkx as nx
from pylab import *
import matplotlib.pyplot as plt

class Conf:
    outdir = "out"


if __name__ == "__main__":

    # set up plotting and figure
    #plt.fig = plt.figure(1, figsize=(3.4,2.5))
    #fig = plt.figure(figsize=(3.54, 4.5)) #single column fig
    #fig = plt.figure(figsize=(7.48, 4.0))  #two column figure

    #plt.rc('font', family='serif', size=8)
    #plt.rc('xtick')
    #plt.rc('ytick')
    #
    #gs = plt.GridSpec(4, 1)
    #gs.update(hspace = 0.0)
    #
    #axs = []
    #axs.append( plt.subplot(gs[0]) )
    #axs.append( plt.subplot(gs[1]) )
    #axs.append( plt.subplot(gs[2]) )
    #axs.append( plt.subplot(gs[3]) )

    #for ax in axs:
    #    ax.minorticks_on()

    conf = Conf()


    ir = 0
    rank = str(ir)
    f5 = h5py.File(conf.outdir+"/run-"+rank+".h5", "r")

    virs = f5['virtuals']
    boun = f5['boundaries']
    locs = f5['locals']
        
    imgs = f5['grid'][:,:,:]

    Nx, Ny, Nt = np.shape(imgs)
    Nx, Ny = 10,10
    print("image size nx {} ny {} nt {}".format(Nx, Ny, Nt))
    #for t in range(Nt):
    for t in [0]:
        img = imgs[:,:,t]

        G = nx.grid_2d_graph(Nx,Ny, periodic=False) 
        
        pos = dict(zip(G.nodes(),G.nodes()))                                            
        ordering = [(y,Nx-1-x) for y in range(Ny) for x in range(Nx)]
        labels = dict(zip(ordering, range(len(ordering))))                              
        
        #nodes = G.nodes()
        #print(nodes)
        
        node_sizes = []
        node_cols = []
        labels = {}
        for (i,j) in G.nodes():
            #print(i,j)
            node_sizes.append(100*i+1)

            if i < 5:
                col = 'black'
            else:
                col = 'blue'
            node_cols.append(col)
        
        nx.draw_networkx(
                G, 
                #pos=pos, 
                pos = nx.spectral_layout(G),
                #pos = nx.circular_layout(G),
                #pos = nx.shell_layout(G),
                with_labels=False, 
                node_size = node_sizes,
                node_shape='s',
                node_color = node_cols,
                )
        
        #nx.draw_networkx_labels(
        #        G, 
        #        pos=pos, 
        #        labels=labels)                              
        
        
        plt.axis('off')                                                                 
        #plt.show()
        plt.savefig('graph.pdf')
