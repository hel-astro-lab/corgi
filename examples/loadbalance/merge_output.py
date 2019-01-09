import numpy as np
import os
import h5py


class Conf:


    outdir = "out200x200n10"
    #outdir = "out200x200n100"
    #outdir = "out4_30x30"
    #outdir = "out4_100x100"
    #outdir = "out2_100x100"



def combine_ranks(fname, f5all):

    #read for base
    rank = str(0)
    fn = fname + rank + ".h5"
    f5 = h5py.File(fn, "r")
    imgs = f5['grid'][:,:,:]
    #imgs = np.where(imgs >= 0)
    imgs[imgs < 0] = 0


    #virs/boun/locs to f5all
    virs = f5['virtuals']
    boun = f5['boundaries']
    locs = f5['locals']

    f5all['virtuals'][:,0]   = virs
    f5all['boundaries'][:,0] = boun
    f5all['locals'][:,0]     = locs

    f5.close()


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

        #virs/boun/locs to f5all
        virs = f5['virtuals']
        boun = f5['boundaries']
        locs = f5['locals']

        f5all['virtuals'][:,ir]   = virs
        f5all['boundaries'][:,ir] = boun
        f5all['locals'][:,ir]     = locs


        ir += 1
        f5.close()

    return imgs





if __name__ == "__main__":

    conf = Conf()
    fname = conf.outdir+"/run-"

    nranks = 10
    f5 = h5py.File(conf.outdir+"/run-merged.h5", "w")

    Nsamples = 201

    f5.create_dataset("virtuals",   (Nsamples,nranks), dtype='f')
    f5.create_dataset("locals",     (Nsamples,nranks), dtype='f')
    f5.create_dataset("boundaries", (Nsamples,nranks), dtype='f')




    ##################################################
    # img + virtuals/locals/boundaries
    imgs = combine_ranks(fname, f5)
    f5.create_dataset("grid", data=imgs)

    ##################################################
    # work
    rank = str(0)
    fn = fname + rank + ".h5"
    f5t = h5py.File(fn, "r")
    work = f5t['work'][:,:,:]
    f5.create_dataset("work", data=work)
    f5t.close()


    ##################################################
    #close
    f5.close()




