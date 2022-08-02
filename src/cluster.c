#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h> 
#include "utilities.h"
#include "initialize.h"
#include "energycalc.h"
#include "montecarlo.h"
#include "cluster.h"

int mc_cluster_move(SYSTEM *S){
	double * U;
	U = (double *)calloc(2,sizeof(double));
    reset_cluster_arrays(S);

	// step 1: generate old clusters
	S->Nclusters = generate_clusters(S);
	int i,j,pt,k,nS=0;
	int *dr = (int *)calloc(3, sizeof(int));
	int move_axis = rand()%3;
	dr[move_axis]=(rand()%2)*2-1;
    
	int select_cl;
    int attempts=0;
	do{
		select_cl = rand()%(S->Nclusters);
        attempts++;
        if (attempts>100) break;
	}while(S->clsize[select_cl]==1); // can do weighted selection here to immitate diffusion!

    int * Uold;
    Uold = S->cl[select_cl];
    int nUold = S->clsize[select_cl];

    // step 2: The old coordinates are stored into S->Elem2
    nS = S->clsize[select_cl];
    S->Elem2_inds_length = 0;

    // solvent and solute accounting to prepare for move
	int * VU = (int *)malloc(nS * sizeof(int)); 
	int ** UVr = (int **)malloc(nS * sizeof(int *)); 
	int ** VUr = (int **)malloc(nS * sizeof(int *)); 
    int nVU = compile_solvents_moved(dr, select_cl, VU, nS, UVr, VUr, S);
	if (nVU == -1){
        return reject_move(S);
	}

    // old energy is calculated
    compile_ids_moved(Uold, nUold, VU, nVU, nS, S);
    U[0] = e_real_cluster_and_solvent(S);
    PARTICLE * charges_old = make_list_of_cluster_charges(Uold, nUold, nS,S);
	U[0] += S->E_k;

	// step 3: move the solvents {VU}, and move the cluster {Uold} = {UU}+{UV} 
    displace_solvents(nVU, UVr, VUr, S);

	for (i=0; i<S->clsize[select_cl]; i++){
		pt = S->cl[select_cl][i];
		(S->Elem)[pt].x = p((S->Elem)[pt].x + dr[0], S->Lx);
		(S->Elem)[pt].y = p((S->Elem)[pt].y + dr[1], S->Ly);
		(S->Elem)[pt].z = p((S->Elem)[pt].z + dr[2], S->Lz);
        S->id[S->Elem[pt].x][S->Elem[pt].y][S->Elem[pt].z] = pt;
	}

        // no need to waste time calculating energy if clusters merged
	if (generate_clusters(S) != S->Nclusters){
        return reject_move(S);
	}

	// step 4: compile list of charge changes for ewald and calculate energy
    U[1] = e_real_cluster_and_solvent(S);
    PARTICLE * charges_new = make_list_of_cluster_charges(Uold, nUold, nS,S);
	S->E_k2 = e_recip(charges_old,charges_new,nS,S);
	U[1] += S->E_k2;

	    // release allocated memory
	free(charges_old); charges_old=NULL;
	free(charges_new); charges_new=NULL;
	free(VU); VU=NULL;
	free(VUr); VUr=NULL;
	free(UVr); UVr=NULL;
    free(S->ids_moved); S->ids_moved=NULL;
    
	// step 5: test acceptance
	double p_accept = exp(-(U[1]-U[0])/S->T);
	if (drand48()<p_accept){
        return accept_move(S);
	}else{
        return reject_move(S);
	}
}

void compile_ids_moved(int *Uold, int nUold, int *VU, int nVU, int nS, SYSTEM *S){
	int i, dx, dy, dz;
    S->N_ids_moved = 0;
    S->ids_moved = (int *)malloc((nS+nVU)*sizeof(int));

	for (i=0; i<nUold; i++){
        S->ids_moved[S->N_ids_moved] = Uold[i];
        S->N_ids_moved++;
	}

	for (i=0; i<nVU; i++){
		S->ids_moved[S->N_ids_moved] = VU[i];
		S->N_ids_moved++;
	}
}

double e_real_cluster_and_solvent(SYSTEM *S){
    double U;
    int xi,yi,zi,xj,yj,zj,i,j;
	for (int i=0; i<S->N_ids_moved; i++){
		U += e1p(S->ids_moved[i], S);
		for (int j=i+1; j<S->N_ids_moved; j++){
			U -= e2p(S->ids_moved[i], S->ids_moved[j], S);
		}
	}
	return U;
}

PARTICLE * make_list_of_cluster_charges(int *Uold, int nUold, int nS, SYSTEM *S){
	PARTICLE * charges = (PARTICLE *)malloc(nS*sizeof(PARTICLE));
	int i, j, pt;
	j = 0;
	for (i=0; i<nUold; i++){
		pt = Uold[i];
        charges[j] = S->Elem[pt];
        j++;
	}
    return charges;
}

int compile_solvents_moved(int *dr, int select_cl, int *VU, int nS, int ** UVr, int **VUr, SYSTEM * S){

	int nVU=0,nUV=0,nUU=0,nUnew=0,nUold=0;
	int i, j;
	int ** Uoldr = (int **)malloc(nS * sizeof(int *)); 
	int ** Unewr = (int **)malloc(nS * sizeof(int *)); 
	int ** UUr = (int **)malloc(nS * sizeof(int *)); 
	int * UV = (int *)malloc(nS * sizeof(int)); 
	
	// Find the ids on the new cluster sites {Uold}
	int pt;
	for (i=0; i<S->clsize[select_cl]; i++){
		pt = S->cl[select_cl][i];
        copy_particle_into_Elem2(pt,S);
        Uoldr[nUold] = (int *)malloc(4*sizeof(int));
        Uoldr[nUold][0] = S->id[S->Elem[pt].x][S->Elem[pt].y][S->Elem[pt].z];
        Uoldr[nUold][1] = p((S->Elem)[pt].x,S->Lx);
        Uoldr[nUold][2] = p((S->Elem)[pt].y,S->Ly);
        Uoldr[nUold][3] = p((S->Elem)[pt].z,S->Lz);
        nUold++;
	}

	// Find the ids on the new cluster sites {Unew}
	for (i=0; i<S->clsize[select_cl]; i++){
		pt = S->cl[select_cl][i];
        Unewr[nUnew] = (int *)malloc(4*sizeof(int));
        Unewr[nUnew][0] = S->id[p(S->Elem[pt].x+dr[0],S->Lx)][p(S->Elem[pt].y+dr[1],S->Ly)][p(S->Elem[pt].z+dr[2],S->Lz)];
        Unewr[nUnew][1] = p((S->Elem)[pt].x+dr[0],S->Lx);
        Unewr[nUnew][2] = p((S->Elem)[pt].y+dr[1],S->Ly);
        Unewr[nUnew][3] = p((S->Elem)[pt].z+dr[2],S->Lz);
        nUnew++;
	}

	// Reject right away if the new sites contain parts or all of a different cluster
	int InCurrCluster;
	for (i=0; i<nS; i++){
		// if Unew[i] is a tail or head
		if ((S->Elem)[Unewr[i][0]].t=='t' || (S->Elem)[Unewr[i][0]].t=='h'){
			// and if Unew[i] is not in the current cluster
			InCurrCluster=0;
			for (j=0; j<nS; j++){
				if (Uoldr[j][0]==Unewr[i][0]) InCurrCluster=1;
			}
			if (InCurrCluster==0){
				// then it must be in another cluster - reject
				return -1;
			}
		}
	}

	// Split {Uoldr} into {UUr} and {UVr}
	// {UUr} is the intersection of {Unewr} and {Uoldr}
	int InUU;
	for (i=0; i<nS; i++){
		InUU = 0;
		for (j=0; j<nS; j++){
			if (Unewr[j][1]==Uoldr[i][1] && Unewr[j][2]==Uoldr[i][2] && Unewr[j][3]==Uoldr[i][3]){
				InUU = 1;
				UUr[nUU] = (int *)malloc(4*sizeof(int));
				memcpy(UUr[nUU],Uoldr[i],4*sizeof(int));
				nUU++;
			} 
		}
		if (InUU==0) {
			UVr[nUV] = (int *)malloc(4*sizeof(int));
			memcpy(UVr[nUV],Uoldr[i],4*sizeof(int));
			nUV++;
		}
	}

	if (dr[0]!=0) qsort(UVr, nUV, sizeof(UVr[0]), compare_dx);
	else if (dr[1]!=0) qsort(UVr, nUV, sizeof(UVr[0]), compare_dy);
	else if (dr[2]!=0) qsort(UVr, nUV, sizeof(UVr[0]), compare_dz);

	// Split set {Unewr} = {VUr} + {UUr}, also compile VU
	for (i=0; i<nS; i++){
		InUU = 0;
		for (j=0; j<nUU; j++){
			if (Unewr[i][1]==UUr[j][1] && Unewr[i][2]==UUr[j][2] && Unewr[i][3]==UUr[j][3]){
				InUU = 1;
			} 
		}
		if (InUU==0) {
			VU[nVU] = Unewr[i][0];
			VUr[nVU] = (int *)malloc(4*sizeof(int));
			memcpy(VUr[nVU],Unewr[i],4*sizeof(int));
            copy_particle_into_Elem2(Unewr[i][0], S);
			nVU++;
		}
	}

	// int repeat, nV=0;
	// for (i=0; i<nVU; i++){
	// 	repeat=0;
	// 	for (j=i+1; j<nVU; j++){
	// 		if (VUr[i][0] == VUr[j][0]) repeat = 1;
	// 	}
	// 	if (repeat==0){
	// 		VU[nV]=VUr[i][0];
	// 		nV++;
	// 	}
	// }

	if (dr[0]!=0) qsort(VUr, nVU, sizeof(VUr[0]), compare_dx);
	else if (dr[1]!=0) qsort(VUr, nVU, sizeof(VUr[0]), compare_dy);
	else if (dr[2]!=0) qsort(VUr, nVU, sizeof(VUr[0]), compare_dz);


	// {VU} are solvents to be replaced by cluster
	// sites of {UV} are gap to be left by cluster
	// move the solvents and update to S->Elem and S->id 

	// printf("------------------------------------------------\n");

	// printf("dr = {%d, %d, %d}\n", dr[0], dr[1], dr[2]);

	// printf("\nnVU = %d sites\n",nVU);
	// for (i=0; i<nVU; i++){
	// 	printf("%d(%c): %d %d %d\n",VUr[i][0],(S->Elem)[VUr[i][0]].t,VUr[i][1],VUr[i][2],VUr[i][3]);
	// }

	// printf("\nnUV = %d sites\n",nUV);
	// for (i=0; i<nUV; i++){
	// 	printf("%d(%c): %d %d %d\n",UVr[i][0],(S->Elem)[UVr[i][0]].t,UVr[i][1],UVr[i][2],UVr[i][3]);
	// }

	if (nVU!=nUV){
		printf("Solvent gap (nUV = %d) != Solvent crushed (nVU = %d). Check code.\n", nUV, nVU);
		exit(1);
	}

	free(UUr); UUr=NULL;
	free(UV); UV=NULL;
	free(Uoldr); Uoldr=NULL;
	free(Unewr); Unewr=NULL;

	return nVU;
}
 
void displace_solvents(int nVU, int ** UVr, int **VUr, SYSTEM * S){
	for (int i=0; i<nVU; i++){
		(S->Elem)[VUr[i][0]].x = UVr[i][1];
		(S->Elem)[VUr[i][0]].y = UVr[i][2];
		(S->Elem)[VUr[i][0]].z = UVr[i][3];
		S->id[UVr[i][1]][UVr[i][2]][UVr[i][3]]=VUr[i][0];
	}
}

int generate_clusters(SYSTEM *S){
	int i,j,k;

	//construct merged neighbor lists
    reset_neighbor_arrays(S);
	get_neigh_lists(S);
	merge_lists(S);

	//fill in cluster lists using neighbor lists
	j=0;
	for(i=0; i<S->Nh+S->No+S->Nt; i++){
		if (S->Nneigh[i]>0){
			for(k=0; k<S->Nneigh[i]; k++){
				S->cl[j][S->clsize[j]] = S->neigh[i][k];
				S->clsize[j]++;
			}
		j++;
		}
	}

	return j;
}

void reset_cluster_arrays(SYSTEM *S){
    for (int i=0; i<S->Nclusters; i++){
        S->clsize[i] = 0;
    }
    S->Nclusters = 0;
}

void reset_neighbor_arrays(SYSTEM *S){
    for (int i=0; i<S->Nh+S->Nt+S->No; i++){
        S->Nneigh[i] = 0;
    }
}

void get_neigh_lists(SYSTEM * S){
	//construct neighbor lists for non-solvent particles
	int i,j,k,p1,p2,x1,y1,z1,d,dx,dy,dz;
	int AlreadyInList;
	int dr[3]={0};
	for(p1=0; p1<S->Nt+S->Nh+S->No; p1++){
		x1 = S->Elem[p1].x;
		y1 = S->Elem[p1].y;
		z1 = S->Elem[p1].z; 
		S->neigh[p1][0]=p1; 
		S->Nneigh[p1]=1;	
        for (j=0;j<3;j++){
            for (d=-1;d<=1;d+=2){
                dr[j]=d;
                p2 = S->id[p(x1+dr[0],S->Lx)][p(y1+dr[1],S->Ly)][p(z1+dr[2],S->Lz)];
                if (S->Elem[p2].t=='h' || S->Elem[p2].t=='t' || S->Elem[p2].t=='o'){
                    AlreadyInList=0;
                    for (k=0; k<S->Nneigh[p1]; k++){
                        if (S->neigh[p1][k]==p2) AlreadyInList=1;
                    }
                    if (AlreadyInList==0){
                        S->neigh[p1][S->Nneigh[p1]]=p2;
                        S->Nneigh[p1]++;	
                    }
                }
            }
            dr[j]=0;
        }
    }
}

void merge_lists(SYSTEM * S){
	int i,j,p1,p2;
	//repeated check if any pair has identical elements
	while(check_overlap(S)==1){
		for(p1=0; p1<S->Nh+S->Nt+S->No; p1++){
			for(p2=p1+1; p2<S->Nh+S->Nt+S->No; p2++){
				if ((S->Nneigh[p1] > 0)&&(S->Nneigh[p2] > 0)){
					merge_pair(p1,p2,S);
				}
			}
		}
	}
}

void merge_pair(int p1, int p2, SYSTEM * S){
	// merge two neighbor lists if they have common elements
	int m,n;
	int overlap=0,repeat;

	// check if they have common elements
	for(m=0;m<S->Nneigh[p1];m++){
		for(n=0;n<S->Nneigh[p2];n++){
			if(S->neigh[p1][m]==S->neigh[p2][n]) overlap=1;
		}
	}

	// if yes, merge two neighbor lists into one
	if(overlap==1){
		for(n=0;n<S->Nneigh[p2];n++){
			repeat=0;
			// check if n-th element exists in i's neighbor list
			for(m=0;m<S->Nneigh[p1];m++){
				if(S->neigh[p2][n]==S->neigh[p1][m]) repeat=1;
			}
			// if does not exist, add to i's neighbor list
			if(repeat==0){
				S->neigh[p1][S->Nneigh[p1]]=S->neigh[p2][n];
				S->Nneigh[p1]++;
			}
		}

		// after merging, kill j's neighbor list
		for(n=0;n<S->Nneigh[p2];n++){
			S->neigh[p2][n]=0;
		}
		S->Nneigh[p2]=0;
	}
}

int check_overlap(SYSTEM * S){
	// returns 1 if any pair of particles in the system have common elements in their neighbor lists
	int m,n,p1,p2;
	for(p1=0; p1<S->Nh+S->Nt+S->No; p1++){
		for(p2=p1+1; p2<S->Nh+S->Nt+S->No; p2++){
			for(m=0; m<S->Nneigh[p1]; m++){
				for(n=0; n<S->Nneigh[p2]; n++){
					if (S->neigh[p1][m]==S->neigh[p2][n]) return 1;
				}
			}
		}
	}

	return 0;
}


