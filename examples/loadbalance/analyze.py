import numpy as np
import os
import h5py

import matplotlib.pyplot as plt
from matplotlib import cm
import palettable as pal
palette = pal.wesanderson.Moonrise1_5.mpl_colormap


class Conf:
    outdir = "out"


if __name__ == "__main__":

    # set up plotting and figure
    #plt.fig = plt.figure(1, figsize=(3.4,2.5))
    fig = plt.figure(figsize=(3.54, 2.2)) #single column fig
    #fig = plt.figure(figsize=(7.48, 4.0))  #two column figure

    plt.rc('font', family='serif', size=8)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )

    for ax in axs:
        ax.minorticks_on()

    conf = Conf()

    ##################################################
    cols = ['k','b','r','g']
    for ir in range(2):
        rank = str(ir)
        f5 = h5py.File(conf.outdir+"/run-"+rank+".h5", "r")

        
        virs = f5['virtuals']
        boun = f5['boundaries']
        locs = f5['locals']
        
        axs[0].plot(locs, linestyle='solid' , color=cols[ir] )
        axs[0].plot(virs, linestyle='dashed', color=cols[ir] )
        #axs[0].plot(boun, "b-")


    axs[0].set_xlabel(r'Step $n$')
    #axs[0].set_yscale('log')

    plt.subplots_adjust(left=0.18, bottom=0.18, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
    plt.savefig('lbal.pdf')

