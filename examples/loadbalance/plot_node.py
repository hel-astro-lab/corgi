import numpy as np
import os
import h5py

import matplotlib.pyplot as plt
from matplotlib import cm
import palettable as pal
palette = pal.wesanderson.Moonrise1_5.mpl_colormap


class Conf:
    #outdir = "out200x200n10"
    outdir = "out200x200n100"
    #outdir = "out4_30x30"
    #outdir = "out4_100x100"
    #outdir = "out2_100x100"


def reduce_image(img, val):
    ret = np.empty_like(img)
    ind = np.where(img == val)

    ret[:] = 0
    ret[ind] = 255

    return ret


def imshow(ax, img):

    xmin, ymin = (0,0)
    xmax, ymax = np.shape(img)

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    extent = [ xmin, xmax, ymin, ymax ]

    mgrid = img
    mgrid = mgrid.T
    ax.imshow(mgrid,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = 'tab20c',
              vmin = 0,
              vmax = np.max(img),
              aspect='auto',
              )


def combine_ranks(fname):

    #read for base
    rank = str(0)
    fn = fname + rank + ".h5"
    f5 = h5py.File(fn, "r")
    imgs = f5['grid'][:,:,:]
    #imgs = np.where(imgs >= 0)
    imgs[imgs < 0] = 0

    #read whatever is left and combine
    ir = 1
    while(True):
        rank = str(ir)
        fn = fname + rank + ".h5"
        if not(os.path.isfile(fn)):
            break

        print("---reading rank: {} / file: {}".format(ir, fn))

        f5 = h5py.File(fn, "r")
        tmp = f5['grid'][:,:,:]
        tmp[tmp < 0] = 0

        imgs = imgs + tmp

        ir += 1

    return imgs


if __name__ == "__main__":

    # set up plotting and figure
    #plt.fig = plt.figure(1, figsize=(3.4,2.5))
    fig = plt.figure(figsize=(3.54, 3.54)) #single column fig
    #fig = plt.figure(figsize=(7.48, 4.0))  #two column figure

    plt.rc('font', family='serif', size=8)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )

    for ax in axs:
        ax.minorticks_on()

    conf = Conf()

    ##################################################
    nx, ny = 1,1
    cols = ['k','b','r','g']

    #ranks = range(1)
    #for ir in ranks:

    #fname = conf.outdir+"/run-"
    #imgs = combine_ranks(fname)

    nranks = 100
    f5all = h5py.File(conf.outdir+"/run-merged.h5", "r")
    imgs = f5all['grid'][:,:,:]
    imgs = imgs / nranks


    nx, ny, nt = np.shape(imgs)
    print("image size nx {} ny {} nt {}".format(nx, ny, nt))
    for t in range(nt):
    #for t in range(41):
        print("t={}".format(t))
        img = imgs[:,:,t]

        imshow(axs[0], img)

        slap = str(t).rjust(4, '0')
        plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
        plt.savefig(conf.outdir + "/node_" + slap + ".png")




    #axs[0].set_ylabel(r'$N_{\mathrm{vir}}/N_{\mathrm{loc}}$')
    #axs[0].set_xlabel(r'Step $n$')

    #plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
    #plt.savefig('nodes.pdf')
