import numpy as np
import sys
import argparse
parser = argparse.ArgumentParser(description='This is an input script generator for an amphiphile lattice model.')
parser.add_argument('-nh', dest="Nh", type=int, default=0)
parser.add_argument('-nt', dest="Nt", type=int, default=0)
parser.add_argument('-j', dest="J", type=float, default=1)
parser.add_argument('-q', dest="q", type=float, default=1)
parser.add_argument('-deffnm', dest="run_id", type=str, default="unnamed")
args = parser.parse_args()

m = 4
J = args.J
Joo, Jot, Jtt = J, J, J
Lx, Ly, Lz = 20, 20, 28
Nh = args.Nh
Nt = args.Nt
N = Nh+Nt
No = 0*Lx*Lz
Nw = Lx*Ly*Lz - Nh - Nt - No
qh, qt = args.q, -args.q
Nsweeps = 100
rho = float(Nh)/float(Lx*Ly*Lz)

x = 13
istrength = (Nh*qh**2 + Nt*qt**2)/(Lx*Ly*Lz)
lambda_D = 1/np.sqrt(4*np.pi*istrength)
#Delta = np.sqrt(3./(4.*np.pi*rho*q2))
Lmin = x*lambda_D
print("Deybe length = {0:.1f}".format(lambda_D), "; the inter-layer separation is {0:.1f} times the debye length.".format(Ly/lambda_D/m))
#print("Delta = ", Delta)
if Ly/m < Lmin:
    #raise Exception("Box dimension too small compared to debye length.")
    print("Box dimension too small compared to debye length.")

s = 2.0
eps = np.exp(- s**2)/s**2
tau_r = 3.77
tau_f = 1.00
sqrta_max = (tau_r * np.pi**3 * N / tau_f / Lx**6.)**(1./6.)
sqrta_min = (tau_r * np.pi**3 * N / tau_f / Ly**6.)**(1./6.)
sqrta_avg = (tau_r * np.pi**3 * N / tau_f / (Lx**2 * Ly**2 * Lz**2))**(1./6.)
kc_conservative, rc_conservative = 2 * s * sqrta_max, s / sqrta_min
kc_avg, rc_avg = 2 * s * sqrta_avg, s / sqrta_avg

print("conservative choice: kc = {0:.3f}, rc = {1:.0f}".format(kc_conservative, rc_conservative))
print("average choice: kc = {0:.3f}, rc = {1:.0f}".format(kc_avg, rc_avg))

fprep_name = "{}.prep".format(args.run_id)
fin_name = "{}.in".format(args.run_id)

with open(fprep_name, "w") as fprep:
    #lines = "run_id\teps{0:.2f}_zh{1:.2f}_sig{2:.1f}\n".format(J, qh, -(qt/qh))
    lines = "run_id\t{}\n".format(args.run_id)
    lines += "old_seq\t0\n"
    lines += "temperature\t1\n"
    lines += "Lx\t{0:.0f}\n".format(Lx)
    lines += "Ly\t{0:.0f}\n".format(Lz)
    lines += "Lz\t{0:.0f}\n".format(Ly)
    lines += "mc_sweeps\t3\n"
    lines += "write_interval\t1\n"
    lines += "number_oil\t{0:.0f}\n".format(No)
    lines += "number_tail\t{0:.0f}\n".format(Nt)
    lines += "number_head\t{0:.0f}\n".format(Nh)
    lines += "number_water\t{0:.0f}\n".format(Nw)
    lines += "number_vacuum\t0\n"
    lines += "eps_ww\t{0:.4f}\n".format(-J)
    lines += "eps_ow\t{0:.4f}\n".format(J)
    lines += "eps_oo\t{0:.4f}\n".format(-Joo)
    lines += "eps_ot\t{0:.4f}\n".format(-Jot)
    lines += "eps_tt\t{0:.4f}\n".format(-Jtt)
    lines += "q_tail\t{0:.4f}\n".format(qt)
    lines += "q_head\t{0:.4f}\n".format(qh)
    lines += "cluster_move_freq\t100\n"
    lines += "random_seed\t8888\n"
    lines += "real_cutoff_switch\t1\n"
    lines += "energy_check_switch\t1\n"
    lines += "ewald_converge_switch\t1\n"
    lines += "cluster_move_switch\t0\n"
    lines += "read_config_switch\t1\n"
    lines += "solvent_in_trj_switch\t0\n"
    lines += "solvent_in_xyz_switch\t0\n"
    lines += "local_swap_switch\t0\n"
    lines += "surf_only_swap_switch\t0\n"
    fprep.write(lines)

with open(fin_name, "w") as fin:
    #lines = "run_id\teps{0:.2f}_zh{1:.2f}_sig{2:.1f}\n".format(J, qh, -(qt/qh))
    lines = "run_id\t{}\n".format(args.run_id)
    lines += "old_seq\t0\n"
    lines += "temperature\t1\n"
    lines += "Lx\t{0:.0f}\n".format(Lx)
    lines += "Ly\t{0:.0f}\n".format(Lz)
    lines += "Lz\t{0:.0f}\n".format(Ly)
    lines += "mc_sweeps\t{0:d}\n".format(Nsweeps)
    lines += "write_interval\t1\n"
    lines += "number_oil\t{0:.0f}\n".format(No)
    lines += "number_tail\t{0:.0f}\n".format(Nt)
    lines += "number_head\t{0:.0f}\n".format(Nh)
    lines += "number_water\t{0:.0f}\n".format(Nw)
    lines += "number_vacuum\t0\n"
    lines += "eps_ww\t{0:.4f}\n".format(-J)
    lines += "eps_ow\t{0:.4f}\n".format(J)
    lines += "eps_oo\t{0:.4f}\n".format(-Joo)
    lines += "eps_ot\t{0:.4f}\n".format(-Jot)
    lines += "eps_tt\t{0:.4f}\n".format(-Jtt)
    lines += "q_tail\t{0:.4f}\n".format(qt)
    lines += "q_head\t{0:.4f}\n".format(qh)
    lines += "cluster_move_freq\t100\n"
    lines += "random_seed\t8888\n"
    lines += "real_cutoff_switch\t1\n"
    lines += "energy_check_switch\t0\n"
    lines += "ewald_converge_switch\t0\n"
    lines += "cluster_move_switch\t0\n"
    lines += "read_config_switch\t1\n"
    lines += "solvent_in_trj_switch\t0\n"
    lines += "solvent_in_xyz_switch\t0\n"
    lines += "local_swap_switch\t0\n"
    lines += "surf_only_swap_switch\t0\n"
    fin.write(lines)

