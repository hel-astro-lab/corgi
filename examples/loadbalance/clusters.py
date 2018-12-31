import numpy as np
import os
import h5py
import cv2

import matplotlib.pyplot as plt
from matplotlib import cm
import palettable as pal
palette = pal.wesanderson.Moonrise1_5.mpl_colormap


class Conf:
    outdir = "out"


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
              cmap = palette,
              vmin = 0,
              vmax = 1,
              aspect='auto',
              )

def compute_clusters(img):

    #cv2.imwrite('out/im.png', img)
    #img2 = cv2.imread('out/im.png',0)
    #cv2.imwrite('img.png', img2)

    #imgray = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    img2 = cv2.convertScaleAbs(img)
    #print(img2.dtype)

    #ret,thresh = cv2.threshold(img2, 127,255,0)
    im2, contours, hierarchy = cv2.findContours(
            img2,
            cv2.RETR_TREE,
            cv2.CHAIN_APPROX_SIMPLE)

    cv2.drawContours(img2, contours, 1, (0,255,0), 3)

    #tresh_min = 0.5
    #tresh_max = 1.0
    #(thresh, imgray) = cv2.threshold(imgray, tresh_min, tresh_max, 0)
    #cv2.imwrite('bw.png', imgray)

    #contours, hierarchy = cv2.findContours(imgray, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    #cv2.drawContours(imgray, contours, -1, (0,255,0), 3)
    #cv2.imwrite('cnt.png', im2)

    #print("number of clusters: {}".format(len(contours)))

    return contours, hierarchy

def analyze_contours(cnts, hierarchy):
    res = {}

    N = len(cnts)
    res['N'] = N
    print("Number of clusters: {}".format(N))

    return res


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
        
        imgs = f5['grid'][:,:,:]

        nx, ny, nt = np.shape(imgs)
        print("image size nx {} ny {} nt {}".format(nx, ny, nt))
        for t in range(nt):
        #for t in [50]:
            print("t={}".format(t))
            img = imgs[:,:,t]
            img = reduce_image(img, ir)

            cnts, hier = compute_clusters(img)
            data = analyze_contours(cnts, hier)

            imshow(axs[0], img)



    #axs[0].set_xlabel(r'Step $n$')
    #axs[0].set_yscale('log')

    plt.subplots_adjust(left=0.18, bottom=0.18, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
    plt.savefig('clusters.pdf')

