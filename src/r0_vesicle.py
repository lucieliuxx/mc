import numpy as np 
import sys
import matplotlib.pyplot as plt 

with open(sys.argv[1],"r") as f:
    s = f.read().split()
    metavars={}
    for i in range(len(s)):
        if i%2==0:
            metavars[s[i]]=s[i+1]

Lx, Ly, Lz = 30, 30, 30
d = 20

# data structures
type_matrix = np.zeros((Lx,Ly,Lz),dtype=str)

# fill all sites with waters first
for x in range(Lx):
    for y in range(Ly):
        for z in range(Lz):
            type_matrix[x,y,z] = 'w'
Nh_alloc, Nt_alloc, No_alloc = 0, 0, 0

# allocate tails and heads
x0, y0, z0 = Lx/2, Ly/2, Lz/2
for x in range(Lx):
    for y in range(Ly):
        for z in range(Lz):
            r = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
            if r > d/2 - 1 and r < d/2 + 1 and Nt_alloc < 2408 :
                type_matrix[x,y,z] = 't'
                Nt_alloc += 1
            elif r > d/2 + 1 and r < d/2 + 2:
                type_matrix[x,y,z] = 'h'
                Nh_alloc += 1
            elif r > d/2 - 2 and r < d/2 - 1:
                type_matrix[x,y,z] = 'h'
                Nh_alloc += 1

print("Allocated {0:d} heads, {1:d} tails, {2:d} oils and {3:d} waters.".format(Nh_alloc,Nt_alloc,No_alloc,Lx*Ly*Lz-No_alloc-Nt_alloc-Nh_alloc))

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
