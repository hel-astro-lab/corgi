import numpy as np
import h5py
import os

import palettable as pal
palette = pal.wesanderson.Moonrise1_5.mpl_colormap

#import networkx as nx
#from networkx.drawing.nx_agraph import graphviz_layout
#import graph_tool.all as gt
from graph_tool.all import *

from pylab import *
import matplotlib.pyplot as plt
import matplotlib 

class Conf:
    outdir = "out2"


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

    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=4.0)

    virs = f5['virtuals']
    boun = f5['boundaries']
    locs = f5['locals']
        
    imgs  = f5['grid'][:,:,:]
    works = f5['work'][:,:,:]

    Nx, Ny, Nt = np.shape(imgs)
    #Nx, Ny = 10,10
    print("image size nx {} ny {} nt {}".format(Nx, Ny, Nt))
    #for t in range(Nt):
    #for t in range(0,101,1):
    for t in [0,1,2,3]:
        print("-------------------", t)

        slap = str(t).rjust(4, '0')
        fname = conf.outdir+'/gt'+slap+'.pdf'

        img  = imgs[:,:,t]
        work = works[:,:,t]
        g = lattice([Nx,Ny], periodic=False)

        ranks = g.new_vertex_property("int")
        loads = g.new_vertex_property("double")

        conns = g.new_edge_property("double")


        #pos = dict(zip(G.nodes(),G.nodes()))                                            
        #ordering = [(y,Nx-1-x) for y in range(Ny) for x in range(Nx)]
        #labels = dict(zip(ordering, range(len(ordering))))                              
        
        #nodes = G.nodes()
        #print(nodes)
        
        #node_sizes = []
        #node_cols = []
        #labels = {}
        #for (i,j) in g.nodes():

        #build vertexes
        j = 0
        for i,v in enumerate(g.vertices()):
            crank = img.flatten()[i]
            ranks[v] = crank

            load = works.flatten()[i]
            loads[v] = 10.0*load

            #crank = img[i,j]
            #col = palette(norm(crank))


        for e in g.edges():
            i = int(e.source() )
            j = int(e.target() )

            #w1 = work.flatten()[i]
            #w2 = work.flatten()[j]
            #d['weight'] = 2.0/(w1 + w2)

            r1 = img.flatten()[i]
            r2 = img.flatten()[j]
            if r1 == r2:
                v = 0.1
            else:
                v = 1.0
            conns[e] = v

        pos = sfdp_layout(g, 
                #eweight=conns,
                groups=ranks,
                cooling_step=0.95, 
                epsilon=1e-2)

        graph_draw(g, 
                pos=pos, 
                vertex_color=ranks,
                vertex_fill_color=ranks,
                vertex_size=loads,
                edge_pen_width=conns,
                output_size=(300,300), 
                output=fname)

        #nx.draw_networkx(
        #        G, 
        #        #pos=pos, 
        #        #pos = nx.spectral_layout(G, dim=2),
        #        #pos = nx.circular_layout(G),
        #        #pos = nx.shell_layout(G),
        #        #pos = nx.spring_layout(G,pos=pos, iterations=1, scale=10.0),
        #        #pos = graphviz_layout(G, prog='neato'),
        #        pos = nx.kamada_kawai_layout(G),
        #        with_labels=False, 
        #        node_size = node_sizes,
        #        node_shape='s',
        #        node_color = node_cols,
        #        )
        
        #nx.draw_networkx_labels(
        #        G, 
        #        pos=pos, 
        #        labels=labels)                              
        
        
        #plt.axis('off')                                                                 
        #plt.show()

        #slap = str(t).rjust(4, '0')
        #plt.savefig(conf.outdir+'/graph2_'+slap+'.png')


