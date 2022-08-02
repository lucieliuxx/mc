#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h> 
#include "utilities.h"
#include "initialize.h"
#include "energycalc.h"
#include "montecarlo.h"

int mc_neighbor_swap(SYSTEM *S){
	double * U;
	U = (double *)calloc(2,sizeof(double));
	int p1,p2,k,pid;
	int DiffusionMeasurement = 0;

	// 1: select a random particle and a random nearest-neighbor particle
    if (S->surf_only_swap_switch==1){
        p1 = rand()%(S->Nt+S->Nh);
    }else{
        p1 = rand()%(S->N);
    }

    if (S->local_swap_switch==1) p2 = get_neighbor(p1,S);
    else p2 = get_random(p1,S);     	

    // 2: do swaps 
	if ((S->Elem)[p1].t==(S->Elem)[p2].t){
		if (DiffusionMeasurement==1){
            accept_swap_of_identical_types_for_diffusion(p1,p2,S);
            apply_swap_without_pbc(p1,p2,S);
        }
        k = 1;
	}else{
		// 3: store old coordinates and calculate energy change
        S->Elem2_inds_length = 0;
        U = site_swap_get_energy(p1,p2,S);

		// 4: if accept, keep changes; if reject, reverse the swaps
		double p_accept = exp(-(U[1]-U[0])/S->T);
		if(drand48()<p_accept){
            //apply_swap_without_pbc(p1,p2,S);
            //fprintf(S->outptr,"accept\n");
			k = accept_move(S);
		}else{
            //fprintf(S->outptr,"reject\n");
			k = reject_move(S);
		}
	} 
    return k;
}

double * site_swap_get_energy(int p1, int p2, SYSTEM *S){
    PARTICLE charges_old[2], charges_new[2];
    int charge_flag = 0;
	double * U;
	U = (double *)calloc(2,sizeof(double));

    // look at the two particle coordinates and types
	if ((S->Elem)[p1].t=='h'||(S->Elem)[p1].t=='t'||(S->Elem)[p2].t=='h'||(S->Elem)[p2].t=='t') charge_flag = 1;

	// calculate old total energy of the relevant sites
	U[0] += e1p(p1,S) + e1p(p2,S) - e2p(p1,p2,S);
    if (charge_flag == 1){ 
		charges_old[0] = (S->Elem)[p1];
		charges_old[1] = (S->Elem)[p2];
        U[0] += S->E_k;
    }

    // note down the old configuration, and do the swap
    copy_particle_into_Elem2(p1,S);
    copy_particle_into_Elem2(p2,S);
	generate_site_site_swap(p1,p2,S);

    // calculate new total energy of the relevant sites
	U[1] += e1p(p1,S) + e1p(p2,S) - e2p(p1,p2,S);
    if (charge_flag == 1){ 
		charges_new[0] = (S->Elem)[p1];
		charges_new[1] = (S->Elem)[p2];
		S->E_k2 = e_recip(charges_old,charges_new,2,S);
        U[1] += S->E_k2;
    }

	return U;
}

void generate_site_site_swap(int p1, int p2, SYSTEM *S){
	// copy ids
	// generate new coordinates
    int x1_old = S->Elem[p1].x;
    int y1_old = S->Elem[p1].y;
    int z1_old = S->Elem[p1].z;
	S->Elem[p1].x = S->Elem[p2].x;
	S->Elem[p1].y = S->Elem[p2].y;
	S->Elem[p1].z = S->Elem[p2].z;
	S->Elem[p2].x = x1_old;
	S->Elem[p2].y = y1_old;
	S->Elem[p2].z = z1_old;
	// update neighbor lists 
	S->id[S->Elem[p1].x][S->Elem[p1].y][S->Elem[p1].z]=p1;
	S->id[S->Elem[p2].x][S->Elem[p2].y][S->Elem[p2].z]=p2;
}


int accept_move(SYSTEM *S){ 
    for (int i=0; i<S->total_vecs_in_ksphere; i++){
        S->rho_k[i] = S->rho_k2[i];
    }
	S->E_k = S->E_k2; 
    S->Elem2_inds_length = 0;
	return 1;
}

int reject_move(SYSTEM *S){ 
    int pid, x, y, z, i;
    for (i=0; i<S->Elem2_inds_length;i++){
        pid = (S->Elem2_inds)[i];
        (S->Elem)[pid] = (S->Elem2)[i];
        x = (S->Elem2)[i].x;
        y = (S->Elem2)[i].y;
        z = (S->Elem2)[i].z;
        S->id[x][y][z] = pid;
    }
    
    // There's actually no need to do anything with S->rho_k2 here because if will get overwritten next time
    ///*
    for (i=0; i<S->total_vecs_in_ksphere; i++){
        S->rho_k2[i] = S->rho_k[i];
    }
    //*/
	S->E_k2 = S->E_k;
    S->Elem2_inds_length = 0;
    return 0;
}

void copy_particle_into_Elem2(int p, SYSTEM *S){
    (S->Elem2_inds)[S->Elem2_inds_length] = p;
    (S->Elem2)[S->Elem2_inds_length] = (S->Elem)[p];
    S->Elem2_inds_length++;
}

int get_neighbor(int i, SYSTEM *S){
	int j;
    int dr[3] = {0};
    dr[rand()%3] = 2*(rand()%2)-1;
    j = S->id[p(S->Elem[i].x+dr[0],S->Lx)][p(S->Elem[i].y+dr[1],S->Ly)][p(S->Elem[i].z+dr[2],S->Lz)];
	return j;
}

int get_random(int i, SYSTEM *S){
	int j = i;
	while(j==i){
		j = rand()%(S->N);
	}
	return j;
}

/*************************FOR DIFFUSION MEASUREMENTS*****************************/
void move_block_of_sites(int *h, int dx, int dy, int dz, SYSTEM *S){
    int pt, s, sx, sy, sz;
    pt = S->id[h[0]][h[1]][h[2]];
    printf("%c\n", (S->Elem)[pt].t);
for (sx=0; sx<=1; sx++){
    for (sy=0; sy<=1; sy++){
            for (sz=0; sz<=1; sz++){
                s = S->id[p(h[0]+dx,S->Lx)][p(h[1]+dy,S->Ly)][p(h[2]+dz,S->Lz)];
                (S->Unwrap)[s].x += dx;
                (S->Unwrap)[s].y += dy;
                (S->Unwrap)[s].z += dz;
            }
        }
    }
}

void accept_swap_of_identical_types_for_diffusion(int p1, int p2, SYSTEM *S){
    // when dynamics are not important (during equilibration)
    // identical type swaps are ignored
    // but when dynamics is needed these swaps are added on
	int xs,ys,zs,dx,dy,dz;
    S->id[S->Elem[p1].x][S->Elem[p1].y][S->Elem[p1].z]=p2;
    S->id[S->Elem[p2].x][S->Elem[p2].y][S->Elem[p2].z]=p1;
    xs = (S->Elem)[p1].x; 
    ys = (S->Elem)[p1].y;
    zs = (S->Elem)[p1].z;
    (S->Elem)[p1].x = (S->Elem)[p2].x; 
    (S->Elem)[p1].y = (S->Elem)[p2].y;
    (S->Elem)[p1].z = (S->Elem)[p2].z;
    (S->Elem)[p2].x = xs; 
    (S->Elem)[p2].y = ys;
    (S->Elem)[p2].z = zs;
}

void apply_swap_without_pbc(int p1, int p2, SYSTEM *S){
    // Takes the change in coordinate of each particle, and apply them to the S.Unwrap[] coordinates without applying pbc
    int dx, dy, dz; // r12 - vector pointing from particle 1 to particle 2
    int *h1, *h2;
    dx = (S->Elem)[p2].x - (S->Elem)[p1].x;
    dy = (S->Elem)[p2].y - (S->Elem)[p1].y;
    dz = (S->Elem)[p2].z - (S->Elem)[p1].z;
    //printf("dx, dy, dz = %d, %d, %d\n",dx,dy,dz);     
    if (dx>2) dx -= S->Lx;
    else if (dx<-2) dx += S->Lx;
    if (dy>2) dy -= S->Ly;
    else if (dy<-2) dy += S->Ly;
    if (dz>2) dz -= S->Lz;
    else if (dz<-2) dz += S->Lz;
    //printf("dx, dy, dz = %d, %d, %d\n",dx,dy,dz);     
    (S->Unwrap)[p1].x += dx;
    (S->Unwrap)[p1].y += dy;
    (S->Unwrap)[p1].z += dz;
    (S->Unwrap)[p2].x -= dx;
    (S->Unwrap)[p2].y -= dy;
    (S->Unwrap)[p2].z -= dz;

    fprintf(S->unwraptr,"%d %d %d\n", dx,dy,dz);
}
