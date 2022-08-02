import numpy as np 
import sys
import matplotlib.pyplot as plt 

with open(sys.argv[1],"r") as f:
    s = f.read().split()
    metavars={}
    for i in range(len(s)):
        if i%2==0:
            metavars[s[i]]=s[i+1]

Lx = int(metavars['Lx'])
Ly = int(metavars['Ly'])
Lz = int(metavars['Lz'])
Nh = int(int(metavars['number_head']))
Nt = int(int(metavars['number_tail']))

# data structures
type_matrix = np.zeros((Lx,Ly,Lz),dtype=str)

# fill all sites with waters first
cell_ndx = np.zeros((Lx*Ly*Lz,3),dtype=int)
i = 0
for x in range(Lx):
    for y in range(Ly):
        for z in range(Lz):
            type_matrix[x,y,z] = 'w'
            cell_ndx[i] = np.array([x,y,z])
            i += 1

Nh_alloc, Nt_alloc = 0, 0

ndx = np.arange(Lx*Ly*Lz)
heads_and_tails = np.random.choice(ndx, size=Nh+Nt, replace=False)
heads_ndx = cell_ndx[heads_and_tails[:Nh]]
tails_ndx = cell_ndx[heads_and_tails[Nh:]]

# allocate tails
for i in range(Nt):
    x = tails_ndx[i,0]
    y = tails_ndx[i,1]
    z = tails_ndx[i,2]
    type_matrix[x,y,z] = 't'
    Nt_alloc += 1

# allocate heads
for i in range(Nh):
    x = heads_ndx[i,0]
    y = heads_ndx[i,1]
    z = heads_ndx[i,2]
    type_matrix[x,y,z] = 'h'
    Nh_alloc += 1


"""
    filename = "squeezed_box.config"
    old_config = np.loadtxt(filename,dtype=str,skiprows=1)
    for i in range(old_config.shape[0]):
        t = old_config[i,3]
        x = int(old_config[i,0])
        y = int(old_config[i,1])+14
        z = int(old_config[i,2])
        if t == 'h':
            type_matrix[x,y,z] = 'h'
            Nh_alloc += 1
        elif t == 't':
            type_matrix[x,y,z] = 't'
            Nt_alloc += 1
"""
    
print("Requested {0:d} heads, {1:d} tails and {2:d} waters.".format(Nh,Nt,Lx*Ly*Lz-Nh-Nt))
print("Allocated {0:d} heads, {1:d} tails and {2:d} waters.".format(Nh_alloc,Nt_alloc,Lx*Ly*Lz-Nt_alloc-Nh_alloc))

with open(metavars['run_id']+"seq0.config","w") as f:
    f.write("{}\n".format("MC steps = 0"))
    # write charged particles first
    for y in range(Ly):
        for x in range(Lx):
            for z in range(Lz):
                if type_matrix[x,y,z] == 'h' or type_matrix[x,y,z] == 't':
                    f.write("{} {} {} {}\n".format(x,y,z,type_matrix[x,y,z]))
    # then write un-charged particles
    for y in range(Ly):
        for x in range(Lx):
            for z in range(Lz):
                if type_matrix[x,y,z] == 'o' or type_matrix[x,y,z] == 'v':
                    f.write("{} {} {} {}\n".format(x,y,z,type_matrix[x,y,z]))
    # then write solvent 
    for y in range(Ly):
        for x in range(Lx):
            for z in range(Lz):
                if type_matrix[x,y,z] == 'w':
                    f.write("{} {} {} {}\n".format(x,y,z,type_matrix[x,y,z]))
    
with open(metavars['run_id']+"seq0.trj","w") as f:
    f.write("{}\n".format("MC steps = 0"))
    for y in range(Ly):
        for x in range(Lx):
            for z in range(Lz):
                if type_matrix[x,y,z] != 'w': 
                    f.write("{} {} {} {}\n".format(x,y,z,type_matrix[x,y,z]))
