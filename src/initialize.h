#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <complex.h>

typedef struct{
	int x,y,z;	 // coords
	char t;		 // type
	double q;	 // charge
} PARTICLE;

typedef struct{
	
	// fileI/O and run settings
	char run_id[30], traj_file[50], unwrapped_traj_file[50], output_file[50], new_config_file[50], old_config_file[50];
	FILE *trajptr, *unwraptr, *outptr, *oldconfigptr, *newconfigptr;
	int random_seed,real_cutoff_switch,energy_check_switch,ewald_converge_switch,cluster_move_switch,read_config_switch,solvent_in_trj_switch,solvent_in_xyz_switch,local_swap_switch,surf_only_swap_switch,old_seq;
	long int MCSTEPS,WRITE;
	
	// physical parameters
	double acc_ratio,T;						
	double eww, eow, eoo, eot, ett, a, E_k, E_k2, qt, qh;
	int Lx, Ly, Lz, Nh, Nt, No, Nw, Nv, N, first_step, rcut, cluster_move_freq;
	
	// data structures for all MC moves
	int ***id; 
	PARTICLE *Elem, *Elem2;
    int Elem2_inds_length, N_ids_moved;
    int *Elem2_inds;
    int * ids_moved;

    // data structures for MC cluster moves
    int **cl, *clsize, Nclusters, **neigh, *Nneigh;

    // data structures for MC dipole moves
    int **surr, Nsurr;
	
    // data structures for ewald summation
    double complex *rho_k, *rho_k2;
    double ** ksphere, kc;
    int total_vecs_in_ksphere;
    int nkx, nky, nkz;

    // data structures for dynamics
    PARTICLE *Unwrap;

	// files and data structures for analysis
	PARTICLE *Elemd;
	char cm_file[50], rho_file[50], cluster_file[50];
	FILE *cmptr, *rhoptr, *clusterptr;
	double xcm, ycm, zcm;

} SYSTEM;

SYSTEM * initialize_system(char* inputscript);

void read_config(SYSTEM *S);

void initialize_slab_config(SYSTEM *S);

void initialize_ewald(SYSTEM *S);

void initialize_clusters(SYSTEM *S);

void initialize_dipoles(SYSTEM *S);

void set_site_oil(int x, int y, int z, int m, SYSTEM *S);

void set_site_water(int x, int y, int z, int m, SYSTEM *S);

void set_site_vac(int x, int y, int z, int m, SYSTEM *S);

void set_site_tail(int x, int y, int z, int m, SYSTEM *S);

void set_site_head(int x, int y, int z, int m, SYSTEM *S);

#endif
