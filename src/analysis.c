#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <complex.h>
#include <string.h> 
#include "utilities.h"
#include "initialize.h"
#include "fileio.h"
#include "energycalc.h"
#include "montecarlo.h"
#include "cluster.h"
#include "analysis.h"

int main(int argc, char *argv[]){
	int Equilibration = 0;
	int DensityProfile = 0;
	int DistanceDist = 0;
	int SizeProfile = 0;
    int HeightSpectrum = 1;
	int i, MCSTEPS, MCSTEPS_EQM=200;
	int j,k,l;

    FILE *ptr;
    char filename[50];
	SYSTEM* S = (SYSTEM *)malloc(sizeof(SYSTEM));
	infile_param_read(argv[1], S);
    
    // Overwrite previous calculation
    if (DensityProfile==1){
        if (S->read_config_switch == 1) sprintf(S->rho_file, "%sseq%d_rho.dat",S->run_id,S->old_seq+1);
        else sprintf(S->rho_file, "%sseq%d_rho.dat",S->run_id,S->old_seq);
        S->rhoptr = fopen(S->rho_file, "w");
    }
    if (DistanceDist==1){
        sprintf(filename,"%sseq%d_hh.dat",S->run_id, S->old_seq+1);
        ptr = fopen(filename,"w");
        fclose(ptr);
        sprintf(filename,"%sseq%d_ht.dat",S->run_id, S->old_seq+1);
        ptr = fopen(filename,"w");
        fclose(ptr);
        sprintf(filename,"%sseq%d_tt.dat",S->run_id, S->old_seq+1);
        ptr = fopen(filename,"w");
        fclose(ptr);
    }
    if (SizeProfile==1){
        initialize_clusters(S);
        char ftj_filename[50];
        sprintf(ftj_filename,"%sseq%d_ftj.dat",S->run_id, S->old_seq+1);
        FILE * ftj_file = fopen(ftj_filename, "w");
        fclose(ftj_file);
    }
    if (HeightSpectrum==1){
        initialize_clusters(S);
        char hq2_top_filename[50];
        char hq2_bot_filename[50];
        char hq2_med_filename[50];
        sprintf(hq2_top_filename,"%sseq%d_hq2_top.dat",S->run_id, S->old_seq+1);
        sprintf(hq2_bot_filename,"%sseq%d_hq2_bot.dat",S->run_id, S->old_seq+1);
        sprintf(hq2_med_filename,"%sseq%d_hq2_med.dat",S->run_id, S->old_seq+1);
        FILE * hq2_top_file = fopen(hq2_top_filename, "w");
        FILE * hq2_bot_file = fopen(hq2_bot_filename, "w");
        FILE * hq2_med_file = fopen(hq2_med_filename, "w");
        fclose(hq2_top_file);
        fclose(hq2_bot_file);
        fclose(hq2_med_file);
    }

    // load trajectory
    if (S->read_config_switch == 1) sprintf(S->traj_file, "%sseq%d.trj", S->run_id, S->old_seq+1);
    else sprintf(S->traj_file, "%sseq%d.trj", S->run_id, S->old_seq);
	S->trajptr = fopen(S->traj_file,"r");

    char CoMfile[50];
    if (S->read_config_switch == 1) sprintf(CoMfile,"%sseq%d_CoM.dat",S->run_id, S->old_seq+1);
    else sprintf(CoMfile,"%sseq%d_CoM.dat",S->run_id, S->old_seq);

    if (S->solvent_in_trj_switch == 0) S->cmptr = fopen(CoMfile, "r");
    else if (S->solvent_in_trj_switch == 1) S->cmptr = fopen(CoMfile, "w");

    // Switch between these two depending on if the trajectory contains solvents
	if (S->solvent_in_trj_switch == 1) S->N = S->Nv + S->No + S->Nw + S->Nh + S->Nt;
    else if (S->solvent_in_trj_switch == 0) S->N = S->Nh + S->Nt;

	/***** 			measurements 			*****/
	char type_of_interest;
	char type[5] = {'h','t','o','w','v'};

	// Discard the equilibration sweeps
	if (Equilibration==1){
		for (i=0; i<MCSTEPS_EQM; i++){
			fscanf(S->trajptr,"MC steps = %*i\n");
			printf("Equilibration step %d\n", i);
			for (j=0;j<S->N;j++){
				fscanf(S->trajptr,"%*d %*d %*d %*c\n");
			}
		}
	}else{
		MCSTEPS_EQM = 0;
	}

	for (i = MCSTEPS_EQM; i < S->MCSTEPS; i += S->WRITE){
		/******************************************************************/
		// density profile
		/******************************************************************/
		if(DensityProfile==1){
            read_next_frame_from_trj(S);

            // STEP 1: calculate or read CoM of the system 
            if (S->solvent_in_trj_switch == 0) read_CoM(S); 
            else if (S->solvent_in_trj_switch == 1) calculate_CoM(S); 

            for (j=0; j<5; j++){
                type_of_interest = type[j];
                measure_density_profile_along_x(type_of_interest,S);
            }
        }

		/******************************************************************/
		// Distance distribution by particle types
		/******************************************************************/
        if (DistanceDist==1){
            read_next_frame_from_trj(S);
            measure_distance_distribution('h','t',S);
            measure_distance_distribution('t','t',S);
            measure_distance_distribution('h','h',S);
        }

		/******************************************************************/
		// cluster size distribution
		/******************************************************************/
		if (SizeProfile==1){
            // read from the trajectory, and if solvent-free, fill in everything with water
            read_next_frame_from_trj_with_added_solvents(S);
            
            // parse the cluster size data into ftj array
            int * ftj = (int *)malloc((S->Nt+S->No+S->Nh+1)*sizeof(int));
            for (j=0; j<S->Nt+S->No+S->Nh+1; j++) ftj[j] = 0;
			S->Nclusters = generate_clusters(S);
            for (k=0; k<S->Nclusters; k++) ftj[S->clsize[k]]++;

            // print the ftj array to file
            char ftj_filename[50];
            sprintf(ftj_filename,"%sseq%d_ftj.dat",S->run_id, S->old_seq+1);
			FILE * ftj_file = fopen(ftj_filename, "a");
            for (j=0; j<S->Nt+S->No+S->Nh+1; j++) fprintf(ftj_file,"%d ", ftj[j]);
			fprintf(ftj_file,"\n");
			fclose(ftj_file);
            free(ftj);
            reset_cluster_arrays(S);
		}

		/******************************************************************/
		// Particle coordinates inside the largest cluster
		/******************************************************************/
		if (HeightSpectrum==1){
            read_next_frame_from_trj_with_added_solvents(S);
            measure_height_spectrum(S);
        }

		/******************************************************************/
		// clean up memory before reading the next frame
		/******************************************************************/
		free(S->Elem);
		free(S->id);
	}

	fclose(S->trajptr);
    fclose(S->cmptr);
    if (DensityProfile==1){
        fclose(S->rhoptr);
    }
    printf("\n");
}

void measure_height_spectrum(SYSTEM *S){
    int j, k, x, y, z, pt;

    // find the largest cluster(s)
    S->Nclusters = generate_clusters(S);
    struct data_index_duplex * all_clusters;
    all_clusters = (struct data_index_duplex *) malloc(S->Nclusters*sizeof(struct data_index_duplex));
    for (k=0; k<S->Nclusters; k++){
        all_clusters[k].data = S->clsize[k];
        all_clusters[k].index = k;
    }
    qsort(all_clusters, S->Nclusters, sizeof(struct data_index_duplex), compare_duplex);

    FILE *fvis;
    int n_largest_clusters = 1;
    int curr_max_cluster_size=0, curr_max_cluster_index;
    for (k=0; k<n_largest_clusters; k++){
        curr_max_cluster_size = all_clusters[k].data;
        curr_max_cluster_index = all_clusters[k].index;

        // determine if this slab is divided across periodic boundaries
        int divided_flag = slab_divided_across_periodic_box(curr_max_cluster_index,S);

        // identify this large cluster's center of mass for good interface location
        int cluster_CoM_y = 0;
        for (j=0; j<curr_max_cluster_size; j++){
            pt = S->cl[curr_max_cluster_index][j];
            cluster_CoM_y += S->Elem[pt].y;
        }
        cluster_CoM_y /= curr_max_cluster_size;
        if (divided_flag == 1) cluster_CoM_y = p(cluster_CoM_y + S->Ly/2, S->Ly);
        
        // Parse data to get h_{ij}
        double ** h_top = (double **)malloc(S->Lx*sizeof(double*));
        double ** h_bot = (double **)malloc(S->Lx*sizeof(double*));
        double ** h_med = (double **)malloc(S->Lx*sizeof(double*));
        for (x=0; x<S->Lx; x++){
            h_top[x] = (double *)malloc(S->Lz*sizeof(double));
            h_bot[x] = (double *)malloc(S->Lz*sizeof(double));
            h_med[x] = (double *)malloc(S->Lz*sizeof(double));
            for (y=0; y<S->Lz; y++){
                h_top[x][y] = 0;
                h_bot[x][y] = (double)(S->Ly)-1;
            }
        }
        /*
        if (k==2) {
            fvis = fopen("cluster_visual.xyz","w");
            fprintf(fvis, "%d\n\n", curr_max_cluster_size);
        }
        */

        // find the interface locations
        for (j=0; j<curr_max_cluster_size; j++){
            pt = S->cl[curr_max_cluster_index][j];
            x = S->Elem[pt].x;
            y = S->Elem[pt].z;
            z = p(S->Elem[pt].y - cluster_CoM_y, S->Ly);
            z = p(z + S->Ly/2, S->Ly);
            //if (k==2) fprintf(fvis,"5 %d %d %d\n", x, z, y);
            if ((double)z > h_top[x][y]) h_top[x][y] = (double)z;
            if ((double)z < h_bot[x][y]) h_bot[x][y] = (double)z;
        }

        // find the midpoint positions
        for (x=0; x<S->Lx; x++){
            for (y=0; y<S->Lz; y++){
                //if (h_top[x][y] == 0 || h_bot[x][y] == S->Ly-1) printf("Bilayer has hole at (%d, %d)\n",x,y);
                h_med[x][y] = 0.5*(h_top[x][y] + h_bot[x][y]);
            }
        }
        

        move_h_to_zero_mean(h_top, S);
        move_h_to_zero_mean(h_bot, S);
        move_h_to_zero_mean(h_med, S);
        
        double ** hq2_top = (double **)malloc((2*S->Lx-1)*sizeof(double *));
        double ** hq2_bot = (double **)malloc((2*S->Lx-1)*sizeof(double *));
        double ** hq2_med = (double **)malloc((2*S->Lx-1)*sizeof(double *));
        get_fourier_spectrum(h_top, hq2_top, S);
        get_fourier_spectrum(h_bot, hq2_bot, S);
        get_fourier_spectrum(h_med, hq2_med, S);

        // print to file
        char hq2_top_filename[50];
        char hq2_bot_filename[50];
        char hq2_med_filename[50];
        sprintf(hq2_top_filename,"%sseq%d_hq2_top.dat",S->run_id, S->old_seq+1);
        sprintf(hq2_bot_filename,"%sseq%d_hq2_bot.dat",S->run_id, S->old_seq+1);
        sprintf(hq2_med_filename,"%sseq%d_hq2_med.dat",S->run_id, S->old_seq+1);
        FILE * hq2_top_file = fopen(hq2_top_filename, "a");
        FILE * hq2_bot_file = fopen(hq2_bot_filename, "a");
        FILE * hq2_med_file = fopen(hq2_med_filename, "a");

        int nx, ny;
        for (nx=-S->Lx+1; nx<S->Lx; nx++){
            for (ny=0; ny<S->Lz; ny++){
                fprintf(hq2_top_file,"%d %d %lf\n", nx, ny, hq2_top[nx+S->Lx-1][ny]);
                fprintf(hq2_bot_file,"%d %d %lf\n", nx, ny, hq2_bot[nx+S->Lx-1][ny]);
                fprintf(hq2_med_file,"%d %d %lf\n", nx, ny, hq2_med[nx+S->Lx-1][ny]);
            }
        }
        fclose(hq2_top_file);
        fclose(hq2_bot_file);
        fclose(hq2_med_file);
        free(h_top);
        free(h_bot);
        free(h_med);
        free(hq2_top);
        free(hq2_bot);
        free(hq2_med);
    }
    
    free(all_clusters);
    reset_cluster_arrays(S);
}

void move_h_to_zero_mean(double **hij, SYSTEM *S){
    double sum_hij;
    int x,y;
    for (x=0; x<S->Lx; x++){
        for (y=0; y<S->Lz; y++){
            sum_hij += hij[x][y];
        }
    }
    for (x=0; x<S->Lx; x++){
        for (y=0; y<S->Lz; y++){
            hij[x][y] -= sum_hij/(S->Lx*S->Lz);
        }
    }
    return;
}

void get_fourier_spectrum(double **hij, double **hq2, SYSTEM *S){
    double ** Rehq = (double **)malloc((2*S->Lx-1)*sizeof(double *));
    double ** Imhq = (double **)malloc((2*S->Lx-1)*sizeof(double *));
    int nx, ny;
    for (nx=-S->Lx+1; nx<S->Lx; nx++){
        Rehq[nx+S->Lx-1] = (double *)malloc(S->Lz*sizeof(double));
        Imhq[nx+S->Lx-1] = (double *)malloc(S->Lz*sizeof(double));
        hq2[nx+S->Lx-1] = (double *)malloc(S->Lz*sizeof(double));
        for (ny=0; ny<S->Lz; ny++){
            Rehq[nx+S->Lx-1][ny] = 0;
            Imhq[nx+S->Lx-1][ny] = 0;
            for (int x=0; x<S->Lx; x++){
                for (int y=0; y<S->Lz; y++){
                    Rehq[nx+S->Lx-1][ny] += hij[x][y]*cos(2*M_PI/S->Lx*x*nx+2*M_PI/S->Lz*y*ny);
                    Imhq[nx+S->Lx-1][ny] += hij[x][y]*sin(2*M_PI/S->Lx*x*nx+2*M_PI/S->Lz*y*ny);
                }
            }
            hq2[nx+S->Lx-1][ny] = Rehq[nx+S->Lx-1][ny]*Rehq[nx+S->Lx-1][ny] + Imhq[nx+S->Lx-1][ny]*Imhq[nx+S->Lx-1][ny];
        }
    }
    free(Rehq);
    free(Imhq);
    return;
}

int slab_divided_across_periodic_box(int clid, SYSTEM *S){
    int i, pt, pty;

    // build an array of unique y values of the particles in the cluster
    int *unique_ys, N_unique_ys=0;
    unique_ys = (int *)malloc(S->Ly * sizeof(int));
    for (i=0; i<S->clsize[clid]; i++){
        pt = S->cl[clid][i];
        pty = S->Elem[pt].y;
        if (int_already_in_array(pty, unique_ys, N_unique_ys) == 0){
            unique_ys[N_unique_ys] = pty;
            N_unique_ys++;
        }
    }
    
    // sort the unique_ys in ascending order
    qsort(unique_ys, N_unique_ys, sizeof(int), compare_int_size);

    // iterate across unique_ys to determine if slab is split
    int split_flag = 0;
    for (i=0; i<N_unique_ys-1; i++){
        if (unique_ys[i+1] - unique_ys[i] > 1) split_flag = 1;
    }

    free(unique_ys);
    return split_flag;
}

int int_already_in_array(int elem, int * array, int len){
    int i,flag=0;
    for (i=0; i<len; i++){
        if (elem == array[i]) flag=1;
    }
    return flag;
}

void measure_density_profile_along_x(char type_of_interest, SYSTEM *S){
	int j;

	// STEP 2: build a frequency histogram for how many particles are at a certain x, y or z height
	double x_CoMframe; // system is setup to have a slab along x direction
	int *x_histogram = (int *)calloc(S->Ly,sizeof(int));
    double c;

	for (j=0;j<(S->Lx)*(S->Ly)*(S->Lz);j++){
		if( (S->Elem)[j].t == type_of_interest ){
			// calculate the coordinates in CoM frame
            c = S->xcm;
			x_CoMframe = (double)((S->Elem)[j].y) + (double)(S->Ly)/2 - c;

			if (x_CoMframe < 0) x_CoMframe += (double)S->Ly;
			else if (x_CoMframe >= (double)(S->Ly)) x_CoMframe -= (double)S->Ly;
            
			// place this particle into a histogram according to its x coordinate
			x_histogram[(int)(x_CoMframe)]++;
		}
	}

	fprintf(S->rhoptr, "%c ", type_of_interest);
	for (j = 0; j < S->Ly; j++){
		fprintf(S->rhoptr, "%d ", x_histogram[j]);
	}
	fprintf(S->rhoptr, "\n");

    free(x_histogram);
}

void measure_distance_distribution(char t1, char t2, SYSTEM *S){
    int i,j,dx,dy,dz;
    //double *R;
    double R;
    FILE *ptr;
    char filename[50];

    if ((t1=='h'&&t2=='t')||(t1=='t'&&t2=='h')){
        sprintf(filename,"%s_ht.dat",S->run_id);
        ptr = fopen(filename,"a");
        for (i=0; i<S->Nt; i++){
            for (j=S->Nt; j<S->Nh+S->Nt; j++){
                R = distance(i,j,S);
                fprintf(ptr,"%lf ", R);
                //fprintf(ptr,"(%.1f,%.1f,%.1f) ", R[0],R[1],R[2]);
            }
        }
    }else if(t1=='h'&&t2=='h'){
        sprintf(filename,"%s_hh.dat",S->run_id);
        ptr = fopen(filename,"a");
        for (i=S->Nt;i<S->Nh+S->Nt;i++){
            for (j=i+1; j<S->Nh+S->Nt; j++){
                R = distance(i,j,S);
                fprintf(ptr,"%lf ", R);
                //fprintf(ptr,"(%.1f,%.1f,%.1f) ", R[0],R[1],R[2]);
            }
        }
    }else if(t1=='t'&&t2=='t'){
        sprintf(filename,"%s_tt.dat",S->run_id);
        ptr = fopen(filename,"a");
        for (i=0; i<S->Nt; i++){
            for (j=i+1; j<S->Nt; j++){
                R = distance(i,j,S);
                fprintf(ptr,"%lf ", R);
                //fprintf(ptr,"(%.1f,%.1f,%.1f) ", R[0],R[1],R[2]);
            }
        }
    }else{
        printf("Type error!\n");
    }
    fprintf(ptr,"\n");
    fclose(ptr);
}

double distance(int i, int j, SYSTEM *S){
    double x1,x2,y1,y2,z1,z2,dx,dy,dz;
    double R;

    x1 = (double)((S->Elem)[i].x);
    y1 = (double)((S->Elem)[i].y);
    z1 = (double)((S->Elem)[i].z);
    x2 = (double)((S->Elem)[j].x);
    y2 = (double)((S->Elem)[j].y);
    z2 = (double)((S->Elem)[j].z);

    dx = x1 - x2;
    dy = y1 - y2;
    dz = z1 - z2;
    if (dx>0.5*(double)(S->Lx)) dx-=(double)(S->Lx);    
    if (dx<-0.5*(double)(S->Lx)) dx+=(double)(S->Lx);    
    if (dy>0.5*(double)(S->Ly)) dy-=(double)(S->Ly);    
    if (dy<-0.5*(double)(S->Ly)) dy+=(double)(S->Ly);    
    if (dz>0.5*(double)(S->Lz)) dz-=(double)(S->Lz);    
    if (dz<-0.5*(double)(S->Lz)) dz+=(double)(S->Lz);    

    R = sqrt(dx*dx + dy*dy + dz*dz);
    return R;
}

void read_CoM(SYSTEM *S){
    int i;
    fscanf(S->cmptr,"%d %lf %lf %lf\n",&i,&(S->xcm),&(S->ycm),&(S->zcm));
}

void calculate_CoM(SYSTEM *S){
   	S->xcm = 0;
	S->ycm = 0;
	S->zcm = 0;
    int Ntotal_for_CoM=0;
	for (int j=0;j<(S->Lx)*(S->Ly)*(S->Lz);j++){
		if( (S->Elem)[j].t == 'w' ){
			S->xcm += (double)((S->Elem)[j].x); 
			S->ycm += (double)((S->Elem)[j].y);
			S->zcm += (double)((S->Elem)[j].z);
			Ntotal_for_CoM++;
		}
	}
	S->xcm /= (double)(Ntotal_for_CoM);
	S->ycm /= (double)(Ntotal_for_CoM);
	S->zcm /= (double)(Ntotal_for_CoM);

	fprintf(S->cmptr,"%.3f %.3f %.3f\n", S->xcm, S->ycm, S->zcm);
}

void read_next_frame_from_trj(SYSTEM *S){
    int i, j, k;
	fscanf(S->trajptr,"MC steps = %d\n", &i);
	printf("\rMeasuring MCsweep %d...", i); 
	S->Elem = malloc((S->N)*sizeof(PARTICLE));
    S->id = (int ***)malloc((S->Lx)*sizeof(int **));
    for (i=0; i<S->Lx; i++){
        S->id[i] = (int **)malloc((S->Ly)*sizeof(int *));
        for (j=0; j<S->Ly; j++){
            S->id[i][j] = (int *)malloc((S->Lz)*sizeof(int));
            for (k=0; k<S->Lz; k++){
                S->id[i][j][k] = -1;
            }
        }
    }

	int x,y,z;
	char t;
	int site_counter = 0;
	for (int j = 0; j < S->N; j++){
		fscanf(S->trajptr,"%d %d %d %c\n",&x,&y,&z,&t);
		if (t=='o'){
		   	set_site_oil(x,y,z,site_counter,S);
			site_counter++;
		}else if (t=='t'){
		   	set_site_tail(x,y,z,site_counter,S);
			site_counter++;
		}else if (t=='w'){
		   	set_site_water(x,y,z,site_counter,S);
			site_counter++;
		}else if (t=='h'){
		   	set_site_head(x,y,z,site_counter,S);
			site_counter++;
		}else if (t=='v'){
		   	set_site_vac(x,y,z,site_counter,S);
			site_counter++;
		}
	}
}

void read_next_frame_from_trj_with_added_solvents(SYSTEM *S){
    // had to manually add in all the solvents for cluster generation to work correctly!
    int i, j, k;
	fscanf(S->trajptr,"MC steps = %d\n", &i);
	printf("Measuring MCsweep %d...\n", i); 
	S->Elem = malloc((S->Lx)*(S->Ly)*(S->Lz)*sizeof(PARTICLE));
    S->id = (int ***)malloc((S->Lx)*sizeof(int **));
    for (i=0; i<S->Lx; i++){
        S->id[i] = (int **)malloc((S->Ly)*sizeof(int *));
        for (j=0; j<S->Ly; j++){
            S->id[i][j] = (int *)malloc((S->Lz)*sizeof(int));
            for (k=0; k<S->Lz; k++){
                S->id[i][j][k] = -1;
            }
        }
    }

	int x,y,z;
	char t;
	int site_counter = 0;
	for (int j = 0; j < S->N; j++){
		fscanf(S->trajptr,"%d %d %d %c\n",&x,&y,&z,&t);
		if (t=='o'){
		   	set_site_oil(x,y,z,site_counter,S);
		}else if (t=='t'){
		   	set_site_tail(x,y,z,site_counter,S);
		}else if (t=='w'){
		   	set_site_water(x,y,z,site_counter,S);
		}else if (t=='h'){
		   	set_site_head(x,y,z,site_counter,S);
		}else if (t=='v'){
		   	set_site_vac(x,y,z,site_counter,S);
		}
        site_counter++;
	}
    int m = S->N;
    for (int i=0; i<S->Lx; i++){
        for (int j=0; j<S->Ly; j++){
            for (int k=0; k<S->Lz; k++){
                if (S->id[i][j][k]==-1){
                    set_site_water(i, j, k, m, S);
                    m++;
                }
            }
        }
    }
}


