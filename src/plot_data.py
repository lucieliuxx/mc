import numpy as np 
import sys
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt 
plt.rcParams.update({'font.size': 13})

def main():
    with open(sys.argv[1],"r") as f:
        s = f.read().split()
        metavars={}
        for i in range(len(s)):
            if i%2==0:
                metavars[s[i]]=s[i+1]
    global L
    L = int(metavars['boxside'])
    global Lf
    Lf = float(L)

    global qh
    qh = float(metavars['q_head'])
    global qt
    qt = float(metavars['q_tail'])

    global Nt
    global Nh
    Nt = float(metavars['number_tail'])
    Nh = float(metavars['number_head'])

    global rho_bulk
    rho_bulk = (Nh)/Lf**3

    global lattice_spacing # in nanometers
    lattice_spacing = 0.2

    global lam_debye
    # here I'm using the fact that only half of the box has solvent
    # the density used for ionic strength is the density in the aqueous phase
    # not the overall density N/L**3
    ionic_strength = (2*Nh/Lf**3)*qh**2 + (2*Nt/Lf**3)*qt**2
    lam_debye = 1./np.sqrt(4. * np.pi * ionic_strength)

    # can discard the early MC sweeps here to get better equilibrated data
    # note the trajectory might already have been trimmed so there might be no need for this
    equilibration_time = 0
    equilibration_time = int(equilibration_time/int(metavars['write_interval']))

    Average_Size_and_Max_Size(metavars)
    #Size_Distribution(metavars)
    #Center_of_Mass(metavars,equilibration_time)
    #Density_Profile(metavars,equilibration_time)       
    #Distance_Distribution(metavars,equilibration_time)       
    #Theoretical_Pair_Correlation(metavars,equilibration_time)       
    #Single_Particle_Diffusion(metavars)

def Average_Size_and_Max_Size(metavars):
    number_of_seqs = int(metavars['old_seq'])+1 
    ftj = np.zeros((1,int(Nh+Nt+1)),dtype=int)
    for new_seq in range(1,number_of_seqs+1):
        filename = metavars['run_id']+'seq{0:}'.format(new_seq)+'_ftj.dat'
        ftj_seq = np.loadtxt(filename,dtype=int)
        ftj = np.concatenate((ftj,ftj_seq),axis=0)
    ftj = ftj[1:]
    jftj = np.arange(ftj.shape[1])*ftj
    avg_jt = np.sum(jftj,axis=1)/np.sum(ftj,axis=1)
    max_jt = np.zeros(ftj.shape[0],dtype=int)
    for t in range(ftj.shape[0]):
        max_jt[t] = np.max(np.where(ftj[t]))
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('t / MC sweeps')
    ax1.set_ylabel(r'$\bar{j}(t)$',color=color)
    ax1.plot(avg_jt,color=color)
    #ax1.set_ylim([0,100])
    ax1.tick_params(axis='y', labelcolor=color)

    color = 'tab:blue'
    ax2 = ax1.twinx()
    ax2.set_ylabel(r'$j_{max}(t)$', color=color)
    ax2.plot(max_jt,color=color)
    #ax2.set_ylim([0,Nh+Nt])
    ax2.tick_params(axis='y', labelcolor=color)

    plt.tight_layout()
    new_seq = int(metavars['old_seq'])+1
    plt.savefig(metavars['run_id']+'seq{0:}'.format(new_seq)+'_jstats.png',dpi=200)

def Size_Distribution(metavars):
    new_seq = int(metavars['old_seq'])+1
    filename = metavars['run_id']+'seq{0:}'.format(new_seq)+'_ftj.dat'
    ftj = np.loadtxt(filename,dtype=int)
    jftj = np.arange(ftj.shape[1])*ftj
    plt.figure(0)
    for t in [0,5,10,49]:
        plt.plot(jftj[t],label=r'$jf({},j)$'.format(t))
    #plt.ylim([0])
    plt.legend(loc='best')
    plt.xlabel('cluster size j')
    plt.tight_layout()
    plt.savefig(metavars['run_id']+'seq{0:}'.format(new_seq)+'_ftj.png',dpi=200)

def pd(dx):
    if dx>L/2:
        dx-=L
    elif dx<-L/2:
        dx+=L
    return dx

def Center_of_Mass(metavars,equilibration_time):
    cmdat = np.loadtxt(metavars['run_id']+'_cm.dat', delimiter=' ')
    cmdat = cmdat[::3,:]
    plt.figure(0)
    plt.plot(cmdat[:,0],'k',label=r'$x_{CoM}$')
    plt.plot(cmdat[:,1],'gray',label=r'$y_{CoM}$')
    plt.plot(cmdat[:,2],'gray',linestyle='--',label=r'$z_{CoM}$')
    plt.ylim([0,L/2])
    plt.legend(loc='best')
    plt.xlabel('time / MC sweeps')
    plt.tight_layout()
    plt.savefig(metavars['run_id']+'_cm.png',dpi=200)

def Density_Profile(metavars,equilibration_time):
    # Parse data out of run_idseq*_rho.dat
    if int(metavars['read_config_switch']) == 1:
        filename = metavars['run_id']+'seq{0:}'.format(int(metavars['old_seq'])+1)+"_rho.dat"
    else:
        filename = metavars['run_id']+'seq{0:}'.format(int(metavars['old_seq']))+"_rho.dat"

    rhodata = np.genfromtxt(filename, delimiter=' ',dtype=int)
    
    volume_per_histogram_slab = float(L * L * 1) 
    rhodata = rhodata/(volume_per_histogram_slab)
    rhodata = rhodata[5*equilibration_time:,:]
    Nsweeps = rhodata.shape[0]/5 
    sample_n = 100
    Nsamples = int(Nsweeps/sample_n)

    s = 10000 #scale 
    
    # Parsing data from rhodata, and calculating errors
    # HEADS
    rho_h_sim = rhodata[0::5,1:]
    rho_h_sim_mu = np.mean(rho_h_sim,axis=0)
    sample_mu = np.zeros((Nsamples,L),dtype=float)
    for j in range(Nsamples):
        sample_mu[j,:] = np.copy(np.mean(rho_h_sim[j*sample_n:(j+1)*sample_n:,:],axis=0))
    rho_h_err = np.std(sample_mu, axis=0)/np.sqrt(sample_n)

    # TAILS
    rho_t_sim = rhodata[1::5,1:]
    rho_t_sim_mu = np.mean(rho_t_sim,axis=0)
    sample_mu = np.zeros((Nsamples,L),dtype=float)
    for j in range(Nsamples):
        sample_mu[j,:] = np.copy(np.mean(rho_t_sim[j*sample_n:(j+1)*sample_n:,:],axis=0))
    rho_t_err = np.std(sample_mu, axis=0)/np.sqrt(sample_n)
    
    # OIL AND WATER
    rho_o_sim = rhodata[2::5,1:]
    rho_o_sim_mean = np.mean(rho_o_sim,axis=0)
    rho_w_sim = rhodata[3::5,1:]
    rho_w_sim_mean = np.mean(rho_w_sim,axis=0)
    rho_v_sim = rhodata[4::5,1:]
    rho_v_sim_mean = np.mean(rho_v_sim,axis=0)

    # Theoretical predictions for net density
    x_lattice = np.arange(0,L)+0.5
    """
    prefac, prefac_err = fit_sim_data(x_lattice, rho_h_sim_mu - rho_t_sim_mu)
    print(prefac, prefac_err)
    print("error is {:.3f}%".format(prefac_err/prefac * 100))
    """

    # Plots of simulation data
    plt.figure(1)
    #plt.errorbar(x_lattice, s*rho_h_sim_mu, yerr = s*rho_h_err, fmt='none', ecolor='darkblue', elinewidth=1, capsize=2, label='head')
    #plt.errorbar(x_lattice, s*rho_t_sim_mu, yerr = s*rho_t_err, fmt='none', ecolor='orange', elinewidth=1, capsize=2, label='tail')
    #plt.errorbar(x_lattice, s*(rho_h_sim_mu - rho_t_sim_mu), yerr = s*(rho_t_err+rho_h_err), fmt='none', ecolor='black', elinewidth=1, capsize=2,label='net head')
    net_rho = s*(qh*rho_h_sim_mu + (qt/8)*rho_t_sim_mu);
    #head_charge_rho = s*qh*rho_h_sim_mu
    head_charge_rho = s*qh*rho_h_sim_mu
    tail_charge_rho = s*(qt/8.)*rho_t_sim_mu
    plt.plot(x_lattice*lattice_spacing, (net_rho+net_rho[::-1])/2,'k-',markerfacecolor='None',label='net positive')
    plt.fill_between(x_lattice*lattice_spacing, y1=(head_charge_rho+head_charge_rho[::-1])/2,y2=0,color='darkblue',alpha=0.5,label='head charge')
    plt.fill_between(x_lattice*lattice_spacing, y1=(tail_charge_rho+tail_charge_rho[::-1])/2,y2=0,color='orange',alpha=0.5,label='tail charge')

    #plt.plot(np.arange(0,L)+0.5,rho_w_sim_mean,'--',color='lightblue',label='{0:}'.format('water'))
    #plt.plot(np.arange(0,L)+0.5,rho_v_sim_mean,'--',color='darkgrey',label='{0:}'.format('empty'))
    #plt.plot(np.arange(0,L)+0.5,blow_up_multiple*rho_avg,'x',markerfacecolor='None',color='orange',linewidth=3)
    #plt.plot(np.arange(0,L)+0.5,blow_up_multiple*rho_avg,'x',markerfacecolor='None',color='darkblue',linewidth=3)

    # Plots of fit
    #x = np.linspace(0,L,1000)
    #nx = net_periodic(x, prefac)
    #plt.plot(x, s*nx, 'k-', label='fit')
    
    # Plot settings and output
    #plt.xlim([0,L])
    plt.xlim([0,L/2*lattice_spacing])
    plt.xlabel(r"$z$ / nm")
    plt.ylabel(r"$\rho_q(z) \times 10^{4}$")
    title_string = r"$|q_h\bar{\rho}_{h}| = |q_{t}\bar{\rho}_{t}| =$"+"${:.1f}$".format(Nh/Lf**3*qh * s) +r"$ \times 10^{-4}$"+ "\n"
    title_string += r"$q_{h}=$"+"{:.2f}".format(qh)+r", $q_{t}=$"+"{:.2f}".format(qt)
    title_string += r", $\lambda = ${:.1f}".format(lam_debye)
    
    #title_string = "charge-free heads and tails\n"+r"mean density $3.125\times 10^{-4}$"
    plt.title(title_string)
    plt.tight_layout()
    plt.legend(loc='upper left')
    if int(metavars['read_config_switch']) == 1:
        plt.savefig(metavars['run_id']+'seq{0:}'.format(int(metavars['old_seq'])+1)+'_rho.png',dpi=200)
    else:
        plt.savefig(metavars['run_id']+'seq{0:}'.format(int(metavars['old_seq']))+'_rho.png',dpi=200)

def fit_sim_data(sim_data_x, sim_data_y):
    popt, pcov = curve_fit(net_periodic, sim_data_x, sim_data_y)
    param = popt[0]
    perr = np.sqrt(np.diag(pcov))[0]
    return param, perr

def rho_head(x, prefactor):
    y = prefactor * (np.exp(-(x-L/2)/lam_debye) + np.exp((x-L)/lam_debye))
    y[x<L/2] = 0.
    return y

def rho_tail(x, prefactor):
    y = prefactor * (np.exp(-x/lam_debye) + np.exp((x-L/2.)/lam_debye))
    y[x>L/2] = 0.
    return y

def net_periodic(x, prefactor):
    net_p = net(x, prefactor) 
    for i in range(1,5):
        net_p += net(x + float(i)*L, prefactor)
        net_p += net(x - float(i)*L, prefactor)
    return net_p 

def net(x, prefactor):
    n_interface1 = prefactor * zigzag(x - L/4)
    n_interface2 = prefactor * zigzag(- (x - 3*L/4))
    return n_interface1 + n_interface2

def zigzag(x):
    n_zigzag = - np.exp(-np.absolute(x)/lam_debye)
    n_zigzag[x<0] = -n_zigzag[x<0]
    return n_zigzag

def Distance_Distribution(metavars,equilibration_time):
    hh = np.genfromtxt(metavars['run_id']+'_hh.dat', delimiter=' ',dtype=float)
    ht = np.genfromtxt(metavars['run_id']+'_ht.dat', delimiter=' ',dtype=float)
    tt = np.genfromtxt(metavars['run_id']+'_tt.dat', delimiter=' ',dtype=float)

    hh = np.hstack(hh)
    ht = np.hstack(ht)
    tt = np.hstack(tt)
    bins = np.linspace(0,int(L),21,endpoint=True)
    PRhh, PRhhx = np.histogram(hh,bins=bins)
    PRht, PRhtx = np.histogram(ht,bins=bins)
    PRtt, PRttx = np.histogram(tt,bins=bins)
    R = (bins[:-1]+bins[1:])/2
    dR = R[1] - R[0]
    PRhh = PRhh/np.sum(PRhh*dR)
    PRht = PRht/np.sum(PRht*dR)
    PRtt = PRtt/np.sum(PRtt*dR)

    plt.figure(2)
    plt.plot(R,PRhh,'-o',color='darkblue',linewidth=3,label='head-head')
    plt.plot(R,PRht,'-o',color='green',linewidth=3,label='head-tail')
    plt.plot(R,PRtt,'-o',color='orange',linewidth=3,label='tail-tail')
    plt.plot([L/2,L/2],[0,1],'--',color='grey')
    plt.xlabel(r'$R$')
    plt.ylabel(r'$P_{ij}(R)$')
    plt.ylim([0,0.1])
    plt.legend(loc='best')
    plt.title(r"Overall $\phi_t = 1.25\times 10^{-3}$, $\phi_h = 1.25\times 10^{-3}$"+"\n"+r"Charges $q_t=-1$, $q_h=+1$")
    plt.tight_layout()
    plt.savefig(metavars['run_id']+'_P(R).png',dpi=200)

    gRhh = PRhh*(Nh-1)/(4.*np.pi*R**2*Nh/L**3)
    gRht = PRht/(4.*np.pi*R**2/L**3)
    gRtt = PRtt*(Nt-1)/(4.*np.pi*R**2*Nt/L**3)
    plt.figure(3)
    plt.plot(R,gRhh,'-o',color='darkblue',linewidth=3,label='head-head')
    plt.plot(R,gRht,'-o',color='green',linewidth=3,label='head-tail')
    plt.plot(R,gRtt,'-o',color='orange',linewidth=3,label='tail-tail')
    plt.plot([L/2,L/2],[0,10],'--',color='grey')
    plt.xlabel(r'$R$')
    plt.ylabel(r'$g_{ij}(R)$')
    plt.ylim([0,4.0])
    plt.legend(loc='best')
    plt.title(r"Overall $\phi_t = 1.25\times 10^{-3}$, $\phi_h = 1.25\times 10^{-3}$"+"\n"+r"Charges $q_t=-1$, $q_h=+1$")
    plt.tight_layout()
    plt.savefig(metavars['run_id']+'_g(R).png',dpi=200)

def Theoretical_Pair_Correlation(metavars,equilibration_time):
    # load meta data 
    L = float(metavars['boxside'])
    Nt = float(metavars['number_tail'])
    Nh = float(metavars['number_head'])
    qt = float(metavars['q_tail'])
    qh = float(metavars['q_head'])

    # load pair distances from <analysis.c>
    hh = np.genfromtxt(metavars['run_id']+'_hh.dat', delimiter=' ',dtype=float)
    ht = np.genfromtxt(metavars['run_id']+'_ht.dat', delimiter=' ',dtype=float)
    tt = np.genfromtxt(metavars['run_id']+'_tt.dat', delimiter=' ',dtype=float)

    # parse the pair distances
    hh = np.hstack(hh)
    ht = np.hstack(ht)
    tt = np.hstack(tt)

    # Histogram and make center bins
    bins = np.linspace(0,int(L),21,endpoint=True)
    PRhh, PRhhx = np.histogram(hh,bins=bins)
    PRht, PRhtx = np.histogram(ht,bins=bins)
    PRtt, PRttx = np.histogram(tt,bins=bins)
    R = (bins[:-1]+bins[1:])/2
    dR = R[1] - R[0]
    PRhh = PRhh/np.sum(PRhh*dR)
    PRht = PRht/np.sum(PRht*dR)
    PRtt = PRtt/np.sum(PRtt*dR)

    # Convert P(R) into g(R)
    gRhh = PRhh*(Nh-1)/(4.*np.pi*R**2*Nh/L**3)
    gRht = PRht/(4.*np.pi*R**2/L**3)
    gRtt = PRtt*(Nt-1)/(4.*np.pi*R**2*Nt/L**3)

    # Calculate theoretical values for very dilute system
    r = np.linspace(0,int(L),100)[1:]
    grhh = np.exp(-qh*qh/r)
    grht = np.exp(-qh*qt/r)
    grtt = np.exp(-qt*qt/r)

    # Plot and save figure
    plt.figure(3)
    plt.plot([L/2,L/2],[0,10],'--',color='grey')

    plt.plot(R,gRhh,'-o',color='darkblue',linewidth=3,label='head-head')
    plt.plot(R,gRht,'-o',color='green',linewidth=3,label='head-tail')
    plt.plot(R,gRtt,'-o',color='orange',linewidth=3,label='tail-tail')
    
    plt.plot(r,grhh,'--',color='darkblue',linewidth=2)
    plt.plot(r,grht,'--',color='green',linewidth=2)
    plt.plot(r,grtt,'--',color='orange',linewidth=2)

    plt.xlabel(r'$R$')
    plt.ylabel(r'$g_{ij}(R)$')
    plt.ylim([0,4.0])
    plt.legend(loc='best')
    plt.title(r"Overall $\phi_t = 1.25\times 10^{-3}$, $\phi_h = 1.25\times 10^{-3}$"+"\n"+r"Charges $q_t=-1$, $q_h=+1$")
    plt.tight_layout()
    plt.savefig(metavars['run_id']+'_g(R).png',dpi=200)

def Single_Particle_Diffusion(metavars):
    #Nw = 7992
    Nw = 0
    No = 1

    # Parse data from .uw file
    uwtraj = np.genfromtxt(metavars['run_id']+'seq0.uw',dtype=np.str)
    Nframes = int(uwtraj.shape[0]/(Nw+No+1))
    MCstep_rows = uwtraj[:,0]=='MC'
    water_rows = uwtraj[:,3]=='w'
    oil_rows = uwtraj[:,3]=='o'

    # Parse water data
    water_all_frames = np.copy(uwtraj[water_rows,:3:]).astype(float)
    Delta_r2_w = np.zeros(Nframes)
    for frame in range(Nframes):
        water_this_frame = water_all_frames[frame*Nw:(frame+1)*Nw,:]
        if frame == 0:
            water_first_frame = water_this_frame
        Delta_r2_w[frame] = np.mean(np.sum((water_this_frame - water_first_frame)**2,axis=1))

    # Parse oil data
    oil_all_frames = np.copy(uwtraj[oil_rows,:3:]).astype(float)
    Delta_r2_o = np.zeros(Nframes)
    for frame in range(Nframes):
        oil_this_frame = oil_all_frames[frame*No:(frame+1)*No,:]
        if frame == 0:
            oil_first_frame = oil_this_frame
        Delta_r2_o[frame] = np.mean(np.sum((oil_this_frame - oil_first_frame)**2,axis=1))

    # Make plots of <dr2>(t)
    dr2_plot = False
    if dr2_plot:
        plt.figure(1)
        plt.plot(np.arange(Nframes), Delta_r2_w, label='water')
        plt.plot(np.arange(Nframes), Delta_r2_o, label='oil')
        plt.tight_layout()
        plt.xlabel('t / MC sweeps')
        plt.savefig('dr2.png',dpi=100)

    # Make plots of individual trajectories
    traj_plot = True
    if traj_plot:   
        plt.figure(2)
        plt.plot(np.arange(Nframes), oil_all_frames[:,0], label='oil x')
        plt.plot(np.arange(Nframes), oil_all_frames[:,1], label='oil y')
        plt.plot(np.arange(Nframes), oil_all_frames[:,2], label='oil z')
        plt.xlabel('t / MC sweeps')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig('diffusion_trajectories.png',dpi=100)

main()

