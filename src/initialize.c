#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <complex.h>
#include <string.h>
#include <omp.h>
#include "utilities.h"
#include "initialize.h"
#include "fileio.h"

SYSTEM* initialize_system(char* inputscript){
	// Read from run_id.in (master input)
	SYSTEM* S = (SYSTEM *)malloc(sizeof(SYSTEM));
	infile_param_read(inputscript, S);

	// Calculate some parameters from the input
	srand(S->random_seed);
	srand48(S->random_seed);
	S->acc_ratio=0;
	S->N = S->Nv + S->No + S->Nw + S->Nh + S->Nt;

	if (S->real_cutoff_switch==1){
		S->rcut = (int)(2./sqrt(S->a)+1);
		if (S->rcut > S->Lx/2) S->rcut = S->Lx/2;
	}else{
		S->rcut = S->Lx/2;
	}
	// check charge neutrality
	double qtot = S->Nt * S->qt + S->Nh * S->qh;
	if (fabs(qtot)>0.00001){
		printf("System is not charge neutral.");
		exit(1);
	}

	// Initialize particle coordinates according to master input
	if(S->read_config_switch == 1){
		read_config(S);
	}else{
		printf("Input file did not specify whether to make or read configuration.\n");
		exit(1);
	}
    
    S->Elem2 = (PARTICLE *)malloc(S->N * sizeof(PARTICLE));
    S->Elem2_inds = (int *)malloc(S->N * sizeof(int));

    // data structures for tracking dynamics
    //S->Unwrap = (PARTICLE *)malloc(S->N * sizeof(PARTICLE));
    //memcpy(S->Unwrap, S->Elem, S->N*sizeof(PARTICLE));

	initialize_ewald(S);
    fprintf(S->outptr, "Ewald params: nc = (%d,%d,%d), rc = %d, alpha = %.3f\n", S->nkx, S->nky, S->nkz, S->rcut, S->a);

    initialize_clusters(S);
    initialize_dipoles(S);

	return S;
}

void read_config(SYSTEM* S){
	int i, j, k, outfile_exists;
    sprintf(S->traj_file, "%sseq%d.trj",S->run_id,S->old_seq+1);
    sprintf(S->unwrapped_traj_file, "%sseq%d.uw",S->run_id,S->old_seq+1);
    sprintf(S->old_config_file, "%sseq%d.config",S->run_id,S->old_seq);
    sprintf(S->new_config_file, "%sseq%d.config",S->run_id,S->old_seq+1);
    sprintf(S->output_file, "%s.out", S->run_id);


    if(access(S->output_file, R_OK)!=-1){
        // file exists
        outfile_exists = 1;
    }else{
        // file doesn't exist
        outfile_exists = 0;
    }

	S->trajptr = fopen(S->traj_file,"w");
	//S->unwraptr = fopen(S->unwrapped_traj_file,"w");
    S->outptr = fopen(S->output_file,"a");
	S->oldconfigptr = fopen(S->old_config_file,"r");
    
	S->id = (int ***)malloc((S->Lx)*sizeof(int **));
    if (S->id == NULL){
        fprintf(S->outptr,"Error -- out of memory.\n");
        exit(1);
    }
    for (i=0; i<S->Lx; i++){
        S->id[i] = (int **)malloc((S->Ly)*sizeof(int *));
        if (S->id[i] == NULL){
            fprintf(S->outptr,"Error -- out of memory.\n");
            exit(1);
        }
        for (j=0; j<S->Ly; j++){
            S->id[i][j] = (int *)malloc((S->Lz)*sizeof(int)); 
            if (S->id[i][j] == NULL){
                fprintf(S->outptr,"Error -- out of memory.\n");
                exit(1);
            }
            for (k=0; k<S->Lz; k++){
                S->id[i][j][k] == -1;
            }
        }
    }
	S->Elem = calloc(S->N, sizeof(PARTICLE));
	if (S->Elem == NULL){
	    fprintf(S->outptr,"Error -- out of memory.\n");
	    exit(1);
	}

	int oil_count=0, vac_count=0, water_count=0, head_count=0, tail_count=0;
	long int nstep;
	fscanf(S->oldconfigptr, "MC steps = %li\n", &nstep); // This is the "MCSTEPS = n" line
	int x,y,z;
	char t;
	for (i=0; i<S->N; i++){
		fscanf(S->oldconfigptr, "%d %d %d %c\n", &x, &y, &z, &t);
		if (t=='o'){
			set_site_oil(x,y,z,i,S);
			oil_count++;
		}else if (t=='t'){
			set_site_tail(x,y,z,i,S);
			tail_count++;
		}else if (t=='h'){
			set_site_head(x,y,z,i,S);
			head_count++;
		}else if (t=='w'){
			set_site_water(x,y,z,i,S);
			water_count++;
		}else if (t=='v'){
			set_site_vac(x,y,z,i,S);
			vac_count++;
		}
	} 

	if(outfile_exists == 0){
        outfile_param_write(S);
    }
    if(S->old_seq == 0){
        write_data(0,S);    
    }

	fprintf(S->outptr, "\nRead from file:\t%d total = %d vac + %d oil + %d water + %d head + %d tail\n\n",S->N,vac_count,oil_count,water_count,head_count,tail_count);

	fclose(S->trajptr);
	fclose(S->oldconfigptr);
	fclose(S->outptr);
	S->trajptr = fopen(S->traj_file,"a");
	S->outptr = fopen(S->output_file,"a");
	//fclose(S->unwraptr);
	//S->unwraptr = fopen(S->unwrapped_traj_file,"a");

	S->first_step = nstep+1;	
}

void initialize_clusters(SYSTEM *S){
    int i;

    S->Nclusters = 0;
    S->clsize = (int *)calloc(S->Nt+S->Nh+S->No,sizeof(int));
    S->cl = (int **)malloc((S->Nt+S->Nh+S->No)*sizeof(int *));
    for (i=0; i<S->Nt+S->Nh+S->No; i++){
        S->cl[i] = (int *)malloc((S->Nt+S->No+S->Nh)*sizeof(int));
    }

    S->Nneigh = (int *)calloc(S->Nh+S->Nt+S->No,sizeof(int));
    S->neigh = (int **)malloc((S->Nh+S->Nt+S->No)*sizeof(int *));
    for (i=0; i<S->Nh+S->Nt+S->No; i++){
        S->neigh[i] = (int *)malloc((S->Nt+S->Nh+S->No)*sizeof(int));
    }
    
    if (S->cl == NULL || S->clsize == NULL || S->neigh == NULL || S->Nneigh == NULL){
        printf("Error -- out of memory.\n");
        exit(1);
    }
}

void initialize_dipoles(SYSTEM *S){
    S->Nsurr = 18;

    S->surr = (int **)malloc(S->Nsurr * sizeof(int*));
    for (int i=0; i<S->Nsurr; i++){
        S->surr[i] = (int *)calloc(3, sizeof(int));
    }

    int n = 0;
    for (int i=0; i<3; i++){
        for (int d=-1; d<=1; d+=2){
            S->surr[n][i] = d;
            n++;
        }
    }

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            if (j!=i){
                for (int d=-1; d<=1; d+=2){
                    S->surr[n][j] = d;
                    n++;
                }
            }
        }
    }

   if (n!=S->Nsurr) printf("error in computing the initial surrounding displacements!\n");
}

void initialize_ewald(SYSTEM *S){
	double dkx = 2*M_PI/(double)(S->Lx);
	double dky = 2*M_PI/(double)(S->Ly);
	double dkz = 2*M_PI/(double)(S->Lz);
	double kx, ky, kz, k2, kr, degen, U_k=0;
	int i,j,k,m,pt;

    S->a = pow((S->Nh+S->Nt)*M_PI*M_PI*M_PI*3.77/(double)(S->Lx*S->Lx*S->Ly*S->Ly*S->Lz*S->Lz),1./3.);
    S->kc = 2 * 2.0 * sqrt(S->a);
    S->rcut = (int)(2.0/sqrt(S->a));

    S->nkx = (int)(S->kc / dkx)+1;
    S->nky = (int)(S->kc / dky)+1;
    S->nkz = (int)(S->kc / dkz)+1;
    // Step 1: counting how many terms we need to sum inside a k-space sphere with radius nkvecs
    S->total_vecs_in_ksphere = 0;
	for (i=S->nkx; i<2*S->nkx+1; i++){
		for (j=0; j<2*S->nky+1; j++){
			for (k=0; k<2*S->nkz+1; k++){
				if ((i==S->nkx)&&(j==S->nky)&&(k==S->nkz)) continue;
				kx = dkx*(double)(i-S->nkx);
				ky = dky*(double)(j-S->nky);
				kz = dkz*(double)(k-S->nkz);
				k2 = kx*kx + ky*ky + kz*kz;
                if (sqrt(k2) < S->kc) S->total_vecs_in_ksphere++;
            }
        }
    }

    // Step 2: making the data structures S->rho_k and S->ksphere for ewald sums
    S->ksphere = (double **)malloc(S->total_vecs_in_ksphere * sizeof(double *));
    if (S->ksphere == NULL){
        fprintf(S->outptr,"Error -- out of memory.\n");
        exit(1);
    }
    for (m = 0; m<S->total_vecs_in_ksphere; m++){
        S->ksphere[m] = (double *)malloc(5*sizeof(double));
        if (S->ksphere[m] == NULL){
            fprintf(S->outptr,"Error -- out of memory.\n");
            exit(1);
        }
    }

    m = 0;
	for (i=S->nkx; i<2*S->nkx+1; i++){
		for (j=0; j<2*S->nky+1; j++){
			for (k=0; k<2*S->nkz+1; k++){
				if ((i==S->nkx)&&(j==S->nky)&&(k==S->nkz)) continue;
				kx = dkx*(double)(i-S->nkx);
				ky = dky*(double)(j-S->nky);
				kz = dkz*(double)(k-S->nkz);
				k2 = kx*kx + ky*ky + kz*kz;
                if (i == S->nkx) degen = 1;
                else degen = 2;
                if (sqrt(k2) < S->kc){
                    S->ksphere[m][0] = kx;
                    S->ksphere[m][1] = ky;
                    S->ksphere[m][2] = kz;
                    S->ksphere[m][3] = k2;
                    S->ksphere[m][4] = degen;
                    m++;
                }
            }
        }
    }

    //for (i=0; i<S->total_vecs_in_ksphere; i++) printf("(%d) kx = %.2lf, ky = %.2lf, kz = %.2lf, degen = %.0f\n", i, S->ksphere[i][0],S->ksphere[i][1],S->ksphere[i][2],S->ksphere[i][4]);



    if (m!=S->total_vecs_in_ksphere) printf("Error generating k-sphere\n");

	S->rho_k = (double complex *)malloc(S->total_vecs_in_ksphere*sizeof(double complex));
	S->rho_k2 = (double complex *)malloc(S->total_vecs_in_ksphere*sizeof(double complex));

	if (S->rho_k2 == NULL || S->rho_k == NULL){
	    fprintf(S->outptr,"Error -- out of memory.\n");
	    exit(1);
	}

	// Step 3: calculating Uewald for the first time
    for (m=0; m<S->total_vecs_in_ksphere; m++){
        kx = S->ksphere[m][0];
        ky = S->ksphere[m][1];
        kz = S->ksphere[m][2];
        k2 = S->ksphere[m][3];
        degen = S->ksphere[m][4];

        // add up all delta_S(k) from individual atoms
        S->rho_k[m]=0.0+0.0*I;
        for (pt=0; pt<S->Nt+S->Nh; pt++){
            kr = kx*(double)(S->Elem)[pt].x+ky*(double)(S->Elem)[pt].y+kz*(double)(S->Elem)[pt].z;
            S->rho_k[m] += cexp(-I*kr)*(S->Elem)[pt].q;
        }
        
        U_k += degen*0.5/((double)(S->Lx*S->Ly*S->Lz)) * (4*M_PI/k2) * pow(cabs(S->rho_k[m]), 2) * exp(- k2/(4*S->a));
	}
	S->E_k = U_k;
	S->E_k2 = U_k;
}

void set_site_head(int x, int y, int z, int m, SYSTEM *S){
	(S->Elem)[m].x = x;
	(S->Elem)[m].y = y;
	(S->Elem)[m].z = z;
	(S->Elem)[m].t = 'h';
	(S->Elem)[m].q = S->qh;
	S->id[x][y][z]=m;
}

void set_site_tail(int x, int y, int z, int m, SYSTEM *S){
	(S->Elem)[m].x = x;
	(S->Elem)[m].y = y;
	(S->Elem)[m].z = z;
	(S->Elem)[m].t = 't';
	(S->Elem)[m].q = S->qt; 
    S->id[x][y][z]=m;
}

// DONT MESS WITH THE CHARGES BLOW THIS LINE!

void set_site_oil(int x, int y, int z, int m, SYSTEM *S){
	(S->Elem)[m].x = x;
	(S->Elem)[m].y = y;
	(S->Elem)[m].z = z;
	(S->Elem)[m].t = 'o';
	(S->Elem)[m].q = 0.; // ZERO
    S->id[x][y][z]=m;
}

void set_site_water(int x, int y, int z, int m, SYSTEM *S){
	(S->Elem)[m].x = x;
	(S->Elem)[m].y = y;
	(S->Elem)[m].z = z;
	(S->Elem)[m].t = 'w';
	(S->Elem)[m].q = 0.; //ZERO
	S->id[x][y][z]=m;
}

void set_site_vac(int x, int y, int z, int m, SYSTEM *S){
	(S->Elem)[m].x = x;
	(S->Elem)[m].y = y;
	(S->Elem)[m].z = z;
	(S->Elem)[m].t = 'v';
	(S->Elem)[m].q = 0.; //ZERO
	S->id[x][y][z]=m;
}

