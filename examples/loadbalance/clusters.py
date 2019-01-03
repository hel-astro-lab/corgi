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

    cv2.drawContours(img2, contours,-1, (0,255,0), 3)

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
    #print("Number of clusters: {}".format(N))

    res['areas'] = []
    res['arcs'] = []
    for cnt in cnts:
        res['areas'].append( cv2.contourArea(cnt) )
        res['arcs'].append(  cv2.arcLength(cnt, True) )

    res['area_mean'] = np.mean(res['areas'])
    res['area_min']  = np.min( res['areas'])
    res['area_max']  = np.max( res['areas'])
    #print("area mean={} min={} max={}".format(res['area_mean'], res['area_min'], res['area_max']))

    res['arc_mean'] = np.mean(res['arcs'])
    res['arc_min']  = np.min( res['arcs'])
    res['arc_max']  = np.max( res['arcs'])
    #print("arc mean={} min={} max={}".format(res['arc_mean'], res['arc_min'], res['arc_max']))

    return res


if __name__ == "__main__":

    # set up plotting and figure
    #plt.fig = plt.figure(1, figsize=(3.4,2.5))
    fig = plt.figure(figsize=(3.54, 4.5)) #single column fig
    #fig = plt.figure(figsize=(7.48, 4.0))  #two column figure

    plt.rc('font', family='serif', size=8)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(4, 1)
    gs.update(hspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )
    axs.append( plt.subplot(gs[2]) )
    axs.append( plt.subplot(gs[3]) )

    for ax in axs:
        ax.minorticks_on()

    conf = Conf()

    ##################################################
    nx, ny = 1,1
    cols = ['k','b','r','g']
    ranks = range(2)
    for ir in ranks:
        col = cols[ir]

        rank = str(ir)
        f5 = h5py.File(conf.outdir+"/run-"+rank+".h5", "r")

        virs = f5['virtuals']
        boun = f5['boundaries']
        locs = f5['locals']
        
        imgs = f5['grid'][:,:,:]

        N = []

        arc_mean = []
        arc_min  = []
        arc_max  = []

        area_mean = []
        area_min  = []
        area_max  = []

        nx, ny, nt = np.shape(imgs)
        print("image size nx {} ny {} nt {}".format(nx, ny, nt))
        for t in range(nt):
            #print("t={}".format(t))
            img = imgs[:,:,t]
            img = reduce_image(img, ir)

            cnts, hier = compute_clusters(img)
            data = analyze_contours(cnts, hier)

            #imshow(axs[0], img, lap)

            N.append(data['N'])

            arc_mean.append(data['arc_mean'])
            arc_min.append(data['arc_min'])
            arc_max.append(data['arc_max'])

            area_mean.append(data['area_mean'])
            area_min.append(data['area_min'])
            area_max.append(data['area_max'])

        #axs[0].plot(locs, linestyle='solid' , color=cols[ir] )
        #axs[0].plot(virs, linestyle='solid', color=cols[ir] )
        axs[0].plot(np.array(virs)/np.array(locs), linestyle='solid', color=cols[ir] )

        axs[1].plot(N, color=col)

        axs[2].plot(arc_mean, linestyle='solid', color=col)
        axs[2].fill_between(np.arange(len(arc_mean)), 
                arc_min, arc_max,
                color=col, alpha=0.5, edgecolor=None)

        axs[3].plot(area_mean, linestyle='solid', color=col)
        axs[3].fill_between(np.arange(len(area_mean)), 
                area_min, area_max,
                color=col, alpha=0.5, edgecolor=None)

    # some theoretical estimates for minimum curve
    area = nx*ny/len(ranks)
    print("nx={} ny={}".format(nx,ny))
    print("area per rank A_loc=",area)
    sidelen = 4*np.sqrt(area)
    print("arc length={}".format(sidelen))
    print("N_vir/N_loc",sidelen/area)


    axs[0].set_ylabel(r'$N_{\mathrm{vir}}/N_{\mathrm{loc}}$')

    axs[1].set_ylabel(r'$k$-clusters')
    axs[2].set_ylabel(r'Arc length')
    axs[3].set_ylabel(r'Area')

    #for ax in axs:
    #    ax.set_yscale('log')
    #    #ax.set_xscale('log')

    axs[1].set_yscale('log')
    axs[2].set_yscale('log')
    axs[3].set_yscale('log')

    axs[3].set_xlabel(r'Step $n$')


    plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
    plt.savefig('clusters.pdf')

