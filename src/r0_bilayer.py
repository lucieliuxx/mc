import numpy as np 
import sys
import matplotlib.pyplot as plt 

make_new_config = True

with open(sys.argv[1],"r") as f:
    s = f.read().split()
    metavars={}
    for i in range(len(s)):
        if i%2==0:
            metavars[s[i]]=s[i+1]

n = 4 # number of bilayers to divide the system into
Lx = int(metavars['Lx'])
Ly = int(metavars['Ly'])
Lz = int(metavars['Lz'])
Nh = int(int(metavars['number_head']))
Nt = int(int(metavars['number_tail']))
No = int(int(metavars['number_oil']))

# data structures
type_matrix = np.zeros((Lx,Ly,Lz),dtype=str)
layer_sites = np.zeros((Lx*Lz,2))
for x in range(Lx):
    for z in range(Lz):
        layer_sites[x*Lz+z] = np.array([x,z])

# fill all sites with waters first
for x in range(Lx):
    for y in range(Ly):
        for z in range(Lz):
            type_matrix[x,y,z] = 'w'
Nh_alloc, Nt_alloc, No_alloc = 0, 0, 0

if make_new_config:
    # setting the desired layers
    N_tail_layers = 2*n
    N_head_layers = 2*n

    """
    top_tail_layers = np.arange(N_tail_layers/2,dtype=int)+int(oil_layers[-1]+1)
    bottom_tail_layers = np.arange(N_tail_layers/2,dtype=int)+int(oil_layers[0]-N_tail_layers/2)
    tail_layers = np.concatenate((bottom_tail_layers,top_tail_layers))
    top_head_layers = np.arange(N_head_layers/2,dtype=int)+int(top_tail_layers[-1]+1)
    bottom_head_layers = np.arange(N_head_layers/2,dtype=int)+int(bottom_tail_layers[0]-N_head_layers/2)
    head_layers = np.concatenate((bottom_head_layers,top_head_layers))
    """
    
    tail_layers = np.array([3,4,10,11,17,18,24,25])
    top_head_layers = tail_layers[::2] - 1
    bottom_head_layers = tail_layers[1::2] + 1
    head_layers = np.sort(np.concatenate((bottom_head_layers,top_head_layers)))

    print("tail", tail_layers, "head", head_layers)

    # allocate tails
    tails_per_layer = int(Nt/N_tail_layers)
    for m in range(len(tail_layers)):
        site_inds = np.random.choice(np.arange(Lx*Lz), size=tails_per_layer, replace=False)
        tail_sites_this_layer = layer_sites[site_inds]
        #y = (tail_layers[m]+int(curr_n*Ly/n)+Ly)%Ly
        y = tail_layers[m]
        for i in range(len(tail_sites_this_layer)):
            x = int(tail_sites_this_layer[i,0])
            z = int(tail_sites_this_layer[i,1])
            type_matrix[x,y,z] = 't'
            Nt_alloc += 1

    # allocate heads
    heads_per_layer = int(Nh/N_head_layers)
    for m in range(len(head_layers)):
        site_inds = np.random.choice(np.arange(Lx*Lz), size=heads_per_layer, replace=False)
        head_sites_this_layer = layer_sites[site_inds]
        #y = (head_layers[m]+int(curr_n*Ly/n)+Ly)%Ly
        y = head_layers[m]
        for i in range(len(head_sites_this_layer)):
            x = int(head_sites_this_layer[i,0])
            z = int(head_sites_this_layer[i,1])
            type_matrix[x,y,z] = 'h'
            Nh_alloc += 1

else:
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
        elif t == 'o':
            type_matrix[x,y,z] = 'o'
            No_alloc += 1

    
print("Requested {0:d} heads, {1:d} tails, {2:d} oils and {3:d} waters.".format(Nh,Nt,No,Lx*Ly*Lz-No-Nh-Nt))
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
