import numpy as np 
import sys
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d
plt.rcParams.update({'font.size': 13})
plt.rcParams.update({"errorbar.capsize": 2})

def main():
    with open(sys.argv[1],"r") as f:
        s = f.read().split()
        metavars={}
        for i in range(len(s)):
            if i%2==0:
                metavars[s[i]]=s[i+1]
    global Lx
    global Ly
    global Lz
    global Nsweeps
    global Nclusters
    global new_seq
    Lx, Ly, Lz = int(metavars['Lx']),int(metavars['Ly']),int(metavars['Lz'])
    Nsweeps = int(metavars['mc_sweeps'])
    Nclusters=1
    new_seq = int(metavars['old_seq'])+1

    logq_top_neat, loghq2_top_neat = parse_hq2(metavars,'top');
    logq_bot_neat, loghq2_bot_neat = parse_hq2(metavars,'bot');
    logq_med_neat, loghq2_med_neat = parse_hq2(metavars,'med');

    
    # make guide-lines for slope = -1, -2, -4 for the bilayer midpoint
    colors = ['tab:red','darkred']
    slopes = [-2,-4]
    labels = [r'$\sim 1/q^2$',r'$\sim 1/q^4$']
    for m in range(2):
        x0 = logq_med_neat[0]
        x1 = logq_med_neat[-1]
        y0 = loghq2_med_neat[0]
        y1 = y0 + slopes[m]*(x1-x0)
        plt.plot([x0,x1],[y0,y1],color=colors[m],label=labels[m])
     
    plt.plot(logq_top_neat, loghq2_top_neat, '2',color='tab:blue',markerfacecolor='none',label='upper interface')
    plt.plot(logq_bot_neat, loghq2_bot_neat, '1',color='tab:blue',markerfacecolor='none',label='lower interface')
    #plt.plot(logq_med_neat, loghq2_med_neat, '+',color='black',markerfacecolor='none',label='bilayer midpoint')
    
    loghq2_med_min, loghq2_med_max = errorbars_hq2(metavars,'med');
    loghq2_med_errs = np.vstack((-loghq2_med_min+loghq2_med_neat,loghq2_med_max-loghq2_med_neat))
    plt.errorbar(logq_med_neat, loghq2_med_neat, yerr=loghq2_med_errs, fmt='x',color='black',markerfacecolor='none',label='bilayer midpoint')


    # plot parameters
    ymin = min([np.min(loghq2_top_neat),np.min(loghq2_bot_neat),np.min(loghq2_med_neat)]) - 0.5
    ymax = max([np.max(loghq2_top_neat),np.max(loghq2_bot_neat),np.max(loghq2_med_neat)]) + 0.5
    plt.ylim([ymin,ymax])
    plt.xlabel(r'$\log{q}$')
    plt.ylabel(r'$\log{\langle |\hat{h}_{q}|^2 \rangle}$')
    plt.tight_layout()
    plt.legend(loc='best',fontsize=8)
    plt.savefig(metavars['run_id']+'seq{0:}'.format(new_seq)+'_loghq2.png', dpi=200)
     

def parse_hq2(metavars,side):
    fname = metavars['run_id']+'seq{0:}'.format(new_seq)+'_hq2_'+side+'.dat'
    nlines = (2*Lx-1)*Lz
    with open(fname,"r") as f:
        hq2_data = np.loadtxt(fname)

    hq2_data = np.reshape(hq2_data,(Nsweeps*Nclusters,nlines,3))
    #print(np.nonzero(np.isnan(hq2_data)))
    hq2_qxqy = np.mean(hq2_data,axis=0)

    # change coordinate from cartesian x-y to polar r
    hq2_qr = np.zeros((hq2_qxqy.shape[0],2))
    hq2_qr[:,0] = np.sqrt(hq2_qxqy[:,0]**2 + hq2_qxqy[:,1]**2) * 2*np.pi/float(Lx)
    hq2_qr[:,1] = hq2_qxqy[:,2]

    # sort and log data
    hq2_qr = hq2_qr[hq2_qr[:,0].argsort()]
    logq = np.log(hq2_qr[1:,:]).T[0]
    loghq2 = np.log(hq2_qr[1:,:]).T[1]

    # setting the upper and lower limit for q so we don't see finite-size effect
    expected_max_logq = np.log((Lx-1)*2*np.pi/Lx)
    expected_min_logq = np.log(2*np.pi/Lx)
    cond1 = logq > expected_min_logq - 0.0001
    cond2 = logq < expected_max_logq + 0.0001
    cond = cond1 & cond2
    logq = logq[cond] 
    loghq2 = loghq2[cond]
    #print(np.min(logq),np.max(logq))
    
    # average over degeneracies
    vals, idx_start, count = np.unique(logq, return_counts=True,return_index=True)
    global res
    res = np.split(np.arange(logq.shape[0]), idx_start[1:])
    loghq2_neat = []
    logq_neat = []
    for i in range(len(res)):
        logq_neat.append(np.mean(logq[res[i]]))
        loghq2_neat.append(np.mean(loghq2[res[i]]))
    loghq2_neat = np.array(loghq2_neat)
    logq_neat = np.array(logq_neat)
    return logq_neat, loghq2_neat

def errorbars_hq2(metavars,side):
    fname = metavars['run_id']+'seq{0:}'.format(new_seq)+'_hq2_'+side+'.dat'
    nlines = (2*Lx-1)*Lz
    with open(fname,"r") as f:
        hq2_data = np.loadtxt(fname)

    hq2_data = np.reshape(hq2_data,(Nsweeps*Nclusters,nlines,3))
    hq2_qxqy = np.mean(hq2_data,axis=0)
    dhq2_qxqy = stats.sem(hq2_data,axis=0)
    dhq2_qxqy[:,:2] = 0
    hq2_qxqy_max, hq2_qxqy_min = hq2_qxqy + dhq2_qxqy, hq2_qxqy - dhq2_qxqy

    hq2_qr_min, hq2_qr_max = np.zeros((hq2_qxqy.shape[0],2)),np.zeros((hq2_qxqy.shape[0],2))
    hq2_qr_min[:,0] = np.sqrt(hq2_qxqy[:,0]**2 + hq2_qxqy[:,1]**2) * 2*np.pi/float(Lx)
    hq2_qr_max[:,0] = np.sqrt(hq2_qxqy[:,0]**2 + hq2_qxqy[:,1]**2) * 2*np.pi/float(Lx)
    hq2_qr_min[:,1], hq2_qr_max[:,1] = hq2_qxqy_min[:,2], hq2_qxqy_max[:,2]

    hq2_qr_min = hq2_qr_min[hq2_qr_min[:,0].argsort()]
    logq_min = np.log(hq2_qr_min[1:,:]).T[0]
    loghq2_min = np.log(hq2_qr_min[1:,:]).T[1]
    hq2_qr_max = hq2_qr_max[hq2_qr_max[:,0].argsort()]
    logq_max = np.log(hq2_qr_max[1:,:]).T[0]
    loghq2_max = np.log(hq2_qr_max[1:,:]).T[1]
    
    # average over degeneracies
    #vals, idx_start, count = np.unique(logq_min, return_counts=True,return_index=True)
    #res = np.split(np.arange(logq_min.shape[0]), idx_start[1:])
    loghq2_min_neat = []
    loghq2_max_neat = []
    for i in range(len(res)):
        loghq2_min_neat.append(np.mean(loghq2_min[res[i]]))
        loghq2_max_neat.append(np.mean(loghq2_max[res[i]]))
    loghq2_min_neat = np.array(loghq2_min_neat)
    loghq2_max_neat = np.array(loghq2_max_neat)
    return loghq2_min_neat, loghq2_max_neat

main()

