#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h> 
#include "utilities.h"
#include "initialize.h"
#include "energycalc.h"
#include "montecarlo.h"
#include "dipole.h"

int mc_dipole_move(SYSTEM *S){
    // dipole moves are inherently non-local, violates local dynamics
    if (S->local_swap_switch==1) return 0;

    double * U;
    U = (double *)calloc(2,sizeof(double));
    PARTICLE charges_old[4], charges_new[4];
    int * ntype1, *neigh1_w, *neigh1_h, *neigh1_t;
    int * ntype2, *neigh2_w, *neigh2_h, *neigh2_t;
    ntype1 = (int *)calloc(3,sizeof(int));
    neigh1_w = (int *)malloc(18*sizeof(int));
    neigh1_h = (int *)malloc(18*sizeof(int));
    neigh1_t = (int *)malloc(18*sizeof(int));
    ntype2 = (int *)calloc(3,sizeof(int));
    neigh2_w = (int *)malloc(18*sizeof(int));
    neigh2_h = (int *)malloc(18*sizeof(int));
    neigh2_t = (int *)malloc(18*sizeof(int));

	int n1,n2,n1e,n2e,nh,nt,nw,pi1,pj1,pi2,pj2,k;
	double p_accept, w2w1;
	
	// select first charged spin i1
    pi1 = rand()%(S->Nt + S->Nh);
	if (S->Elem[pi1].t=='w') return 0;

    ntype1 = count_types_around_particle(pi1,neigh1_w,neigh1_h,neigh1_t,S);
    nw = ntype1[0];
    nh = ntype1[1];
    nt = ntype1[2];
	n1e = nw + 1;
	
    // early rejections if there are no oppositely charged neighbors to i1; else, select cell j1
	if (S->Elem[pi1].t=='h'){
		if (nt==0) return 0;
		n1 = nt;
		pj1 = neigh1_t[rand()%nt];
	}
	if (S->Elem[pi1].t=='t'){
		if (nh==0) return 0;
		n1 = nh;
		pj1=neigh1_h[rand()%nh];
	}

	// Select empty positions to which the dipole should move
    pi2 = rand()%(S->Nw) + S->Nh + S->Nt;

    ntype2 = count_types_around_particle(pi2,neigh2_w,neigh2_h,neigh2_t,S);
    nw = ntype2[0];
    nh = ntype2[1];
    nt = ntype2[2];
	n2e = nw;

	if (n2e==0) return 0;
    pj2 = neigh2_w[rand()%nw];

	// calculate U(0) and note down the old configuration
    S->Elem2_inds_length = 0;
    int moved_list[4] = {pi1,pj1,pi2,pj2};
    for (int i=0; i<4; i++){
        U[0] += e1p(moved_list[i],S);
        for (int j=i+1; j<4; j++){ 
            U[0] -= e2p(moved_list[i],moved_list[j],S);
        }
        copy_particle_into_Elem2(moved_list[i],S);
        charges_old[i] = S->Elem[moved_list[i]];
    }
    U[0] += S->E_k;

    // Make the move
    int pid_moved_to_i2;
	if (drand48()>0.5){
        pid_moved_to_i2=pi1;
        generate_site_site_swap(pi1,pi2,S);
        generate_site_site_swap(pj1,pj2,S);
	}else{
        pid_moved_to_i2=pj1;
        generate_site_site_swap(pj1,pi2,S);
        generate_site_site_swap(pi1,pj2,S);
	}

    if (S->Elem[pid_moved_to_i2].t=='h'){
        n2 = nt + 1;
    }else if(S->Elem[pid_moved_to_i2].t=='t'){
        n2 = nh + 1;
    }

	// calculate U(2)
    for (int i=0; i<4; i++){
        U[1] += e1p(moved_list[i],S);
        for (int j=i+1; j<4; j++){ 
            U[1] -= e2p(moved_list[i],moved_list[j],S);
        }
        charges_new[i] = S->Elem[moved_list[i]];
    }
    S->E_k2 = e_recip(charges_old,charges_new,4,S);
    U[1] += S->E_k2;

	// Monte Carlo rule
	w2w1= (((double)(n1*n2e))/((double)(n2*n1e)))*exp(-(U[1]-U[0])/S->T);
	if (w2w1 < 1){
		p_accept = w2w1;
	} else{
		p_accept = 1;
	}
	if (drand48() < p_accept){
        k = accept_move(S);
	} else{
        k = reject_move(S);
	}

    return k;
}

int * count_types_around_particle(int p1, int * neigh_w, int * neigh_h, int * neigh_t, SYSTEM *S){
    // Count the neighbors according to their identity (charge) 
    int *ntype;
    ntype = (int *)calloc(3,sizeof(int));

    int x1,y1,z1,x2,y2,z2,p2;
    x1 = S->Elem[p1].x;
    y1 = S->Elem[p1].y;
    z1 = S->Elem[p1].z;
	for (int i=0; i<S->Nsurr; i++){
        x2 = p(x1+S->surr[i][0],S->Lx);
        y2 = p(y1+S->surr[i][1],S->Ly);
        z2 = p(z1+S->surr[i][2],S->Lz);
        p2 = S->id[x2][y2][z2];
		if (S->Elem[p2].t=='w'){
			neigh_w[ntype[0]]=p2;
			ntype[0]++;
		}else if(S->Elem[p2].t=='h'){
			neigh_h[ntype[1]]=p2;
			ntype[1]++;
		}else if(S->Elem[p2].t=='t'){
			neigh_t[ntype[2]]=p2;
			ntype[2]++;
		}
	}
    return ntype;
}
