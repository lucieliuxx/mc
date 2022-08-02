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
#include "dipole.h"
#include <time.h>

int main(int argc, char *argv[]){
	/***** 			initialize 				*****/
	SYSTEM* S = initialize_system(argv[1]);
	if (S->energy_check_switch==0) fprintf(S->outptr, "Initialization done.\n");
	// fprintf(S->outptr, "%d clusters.\n",generate_clusters(S->Elem,S->id,S));

	if (S->energy_check_switch==1){
		error_check(S); 
	}
	if (S->ewald_converge_switch==1) ewald_converge(S);

	/***** 		generate trajectory 		*****/
	long int i, AttemptsTotal=0;
    int accept_flag, j, k;
    int swap_freq;
    if (S->surf_only_swap_switch==1) swap_freq = S->Nh + S->Nt;
    else swap_freq = S->N;

    //swap_freq = 1; // for debugging

	int N_neighbor_swaps_per_cluster_move = swap_freq / S->cluster_move_freq;

	time_t now, start, finish;
	time(&start);

	for (i = S->first_step; i< S->first_step + S->MCSTEPS; i++){
        time(&now);
        fprintf(S->outptr, "Running MCSweep %li at %s", i,ctime(&now));

		if (S->cluster_move_switch == 0){
            for (j=0;j<swap_freq;j++){
                // spin swaps
                accept_flag = mc_neighbor_swap(S);
                AttemptsTotal++;
                S->acc_ratio += (double)accept_flag;
                // dipole moves
                accept_flag = mc_dipole_move(S);
                AttemptsTotal++;
                S->acc_ratio += (double)accept_flag;
            }
		} else if(S->cluster_move_switch == 1){
            for (j=0; j<S->cluster_move_freq; j++){
                for (k=0; k<N_neighbor_swaps_per_cluster_move; k++){
                    // spin swaps
                    accept_flag = mc_neighbor_swap(S);
                    AttemptsTotal++;
                    S->acc_ratio += (double)accept_flag;
                    // dipole moves
                    accept_flag = mc_dipole_move(S);
                    AttemptsTotal++;
                    S->acc_ratio += (double)accept_flag;
                }
                // cluster moves
                accept_flag = mc_cluster_move(S);
                AttemptsTotal++;
                S->acc_ratio += (double)accept_flag;
            }
		}
		
		// Write to .trj file
		if (i%S->WRITE==0){
			if (S->energy_check_switch==1){
				error_check(S); 
			}
			write_data(i,S);
            if (S->ewald_converge_switch==1) ewald_converge(S);
			//write_unwrapped_trajectory(i,S);
		}

		// Every 1000 WRITEs to .trj, write to .config so simulation can be restarted from this frame
		if (i%(S->WRITE*50)==0){
			write_final_config(i,S);
		}
	}
	time(&finish);
	fprintf(S->outptr,"Total running time = %.0f seconds, Average time per sweep = %.2f seconds.\n", difftime(finish,start), difftime(finish,start)/(float)S->MCSTEPS);

    //measure_tau(S);

	/***** 			acceptance stats 		*****/
	S->acc_ratio=S->acc_ratio/(double)(AttemptsTotal);
	fprintf(S->outptr, "Acceptance = %lf \n", S->acc_ratio);

	/*****write out last frame for continuation*****/
	write_final_config(i-1,S);
}
