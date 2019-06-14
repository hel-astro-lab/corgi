import numpy as np
import math

#import scipy
#from pylab import *
#import palettable as pal
#from matplotlib import cm

#cmap = cm.get_cmap('inferno_r')
#cmap = pal.wesanderson.Moonrise1_5.mpl_colormap


import loadbalance as corgi



# setup configuration
corgi.xmin = corgi.ymin = 0.0
corgi.xmax = corgi.ymax = 1.0

corgi.Nrank = 4
corgi.Nx = 4
corgi.Ny = 4

corgi.grid = np.zeros( (corgi.Nx, corgi.Ny) )


corgi.grid[:2, :2] = 0
corgi.grid[:2, 2:] = 1

corgi.grid[2:, :2] = 2
corgi.grid[2:, 2:] = 3

print corgi.grid


# load nodes
nodes = []
for rank in range(corgi.Nrank):
    n = corgi.grid(rank)

    for i in range(corgi.Nx):
        for j in range(corgi.Ny):
            if corgi.grid[i,j] == rank:
                c = corgi.cell(i,j,rank)
                c.data = corgi.Nrank*1000 + i*100 + j*10

                print "inserting cell {} at ({},{}) d={}".format(rank, i, j, c.data)
                #n.cells = np.append(n.cells, c)
                n.cells.append( c ) 
    nodes.append(n)



print "\n\n\n"
print "Unit cell test"
print "--------------------------------------------------"
cell1 = nodes[0].cells[0]
print "0 0 data", cell1.data
cell1.data += 1

print "neighs"
print "0 0", cell1.neighs(0,0)
print "1 0", cell1.neighs(1,0)
print "0 1", cell1.neighs(0,1)
print "1 1", cell1.neighs(1,1)
print "-1 0", cell1.neighs(1,0)
print "0 -1", cell1.neighs(0,1)
print "-1 -1", cell1.neighs(1,1)

print "full neighborhood"
print cell1.full_neighborhood()


print "virtuals:"
print "--------------------------------------------------"
print nodes[0].get_all_virtuals()


print "returning cell and pointer-like behavior"
print "--------------------------------------------------"
print "0 1 index:", nodes[0].get_neighbor_index(cell1, 0, 1)

cell2 = nodes[0].get_neighbor_cell(cell1, 0, 1)
print "data value:", cell2.data
cell2.data += 1

cell2copy = nodes[0].get_neighbor_cell(cell1, 0, 1)
print "data value:", cell2copy.data
cell2.data += 1
print "data value:", cell2copy.data


print "virtual neighborhood:"
print "--------------------------------------------------"
corgi.communicate(nodes)
cell = nodes[0].cells[0]
print cell.owner
print "virtual neighborhood", nodes[0].virtual_neighborhood(cell)

print "communication"
print "--------------------------------------------------"
#nodes[0].pack_virtuals()
corgi.communicate(nodes)

print nodes[0].send_queue_address
print nodes[0].virtuals[0].data


print "--------------------------------------------------"
print "adoption"
corgi.adopt(nodes)






