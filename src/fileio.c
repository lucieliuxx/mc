#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utilities.h"
#include "initialize.h"
#include "fileio.h"

void write_data(long int nstep, SYSTEM *S){
	fprintf(S->trajptr, "MC steps = %li\n", nstep);
    
    if (S->solvent_in_trj_switch == 0){
        int Nsurf = S->Nh + S->Nt;
        for (int i=0; i<Nsurf + S->No; i++){
            fprintf(S->trajptr, "%d %d %d %c\n", (S->Elem)[i].x,(S->Elem)[i].y,(S->Elem)[i].z,(S->Elem)[i].t);
        }
        //write_CoM(nstep, S);
    }else if (S->solvent_in_trj_switch == 1){
        for (int i=0; i<(S->N); i++){
            fprintf(S->trajptr, "%d %d %d %c\n", (S->Elem)[i].x,(S->Elem)[i].y,(S->Elem)[i].z,(S->Elem)[i].t);
        }
    }
	
    fclose(S->trajptr);
	fclose(S->outptr);
	S->trajptr = fopen(S->traj_file,"a");
	S->outptr = fopen(S->output_file,"a");
}

void write_CoM(long int nstep, SYSTEM *S){
    char CoMfile[50];
    if (S->read_config_switch == 1){
        sprintf(CoMfile, "%sseq%d_CoM.dat",S->run_id,S->old_seq+1);
    }else{
        sprintf(CoMfile, "%sseq%d_CoM.dat",S->run_id,S->old_seq);
    }
    FILE *comptr = fopen(CoMfile,"a");
    
    double xc=0, yc=0, zc=0;
    for (int i=0; i<(S->N); i++){
        if ((S->Elem)[i].t == 'w'){
            xc += (double)((S->Elem)[i].x);
            yc += (double)((S->Elem)[i].y);
            zc += (double)((S->Elem)[i].z);
        }else{
            //do not count if it's vacuum
        }
    }
    xc /= (double)(S->Nw);
    yc /= (double)(S->Nw);
    zc /= (double)(S->Nw);

    fprintf(comptr, "%li %lf %lf %lf\n", nstep, xc, yc, zc);
    fclose(comptr);
}

void write_unwrapped_trajectory(long int nstep, SYSTEM *S){
	fprintf(S->unwraptr, "MC steps = %li\n", nstep);
	for (int i=0; i<(S->N); i++){
        if ((S->Elem)[i].t=='o') fprintf(S->unwraptr, "%d %d %d %c\n", (S->Unwrap)[i].x,(S->Unwrap)[i].y,(S->Unwrap)[i].z,(S->Unwrap)[i].t);
	}
	fclose(S->unwraptr);
	S->unwraptr = fopen(S->unwrapped_traj_file,"a");
}

void write_final_config(long int nstep, SYSTEM *S){
	S->newconfigptr = fopen(S->new_config_file,"w");
	fprintf(S->newconfigptr, "MC steps = %li\n", nstep);
	for (int i=0; i<(S->N); i++){
		fprintf(S->newconfigptr, "%d %d %d %c\n", (S->Elem)[i].x, (S->Elem)[i].y, (S->Elem)[i].z, (S->Elem)[i].t);
	}
	fclose(S->newconfigptr);
}
				
void outfile_param_write(SYSTEM *S){
	fprintf(S->outptr,"####OUTPUT SCRIPT. DO NOT SET PARAMETERS HERE.####\n\n");
	fprintf(S->outptr,"%-24s%s\n", "run_id",S->run_id);
	fprintf(S->outptr,"%-24s%d\n", "old_seq",S->old_seq);
	fprintf(S->outptr,"%-24s%lf\n", "temperature", S->T);
	fprintf(S->outptr,"%-24s%d\n", "Lx", S->Lx);
	fprintf(S->outptr,"%-24s%d\n", "Ly", S->Ly);
	fprintf(S->outptr,"%-24s%d\n", "Lz", S->Lz);
	fprintf(S->outptr,"%-24s%li\n", "mc_sweeps", S->MCSTEPS);
	fprintf(S->outptr,"%-24s%li\n", "write_interval", S->WRITE);
	fprintf(S->outptr,"%-24s%d\n%-24s%d\n%-24s%d\n%-24s%d\n%-24s%d\n", "number_oil", S->No,"number_tail", S->Nt,"number_head", S->Nh,"number_water", S->Nw,"number_vacuum", S->Nv);
	fprintf(S->outptr,"%-24s%lf\n%-24s%lf\n%-24s%lf\n%-24s%lf\n%-24s%lf\n", "eps_ww", S->eww,"eps_ow", S->eow,"eps_oo", S->eoo, "eps_ot", S->eot, "eps_tt", S->ett);
	fprintf(S->outptr,"%-24s%lf\n%-24s%lf\n", "q_tail", S->qt,"q_head", S->qh);
	fprintf(S->outptr,"%-24s%d\n", "cluster_move_freq", S->cluster_move_freq);
	fprintf(S->outptr,"%-24s%d\n", "random_seed", S->random_seed);
	fprintf(S->outptr,"%-24s%d\n", "real_cutoff_switch", S->real_cutoff_switch);
	fprintf(S->outptr,"%-24s%d\n", "energy_check_switch", S->energy_check_switch);
	fprintf(S->outptr,"%-24s%d\n", "ewald_converge_switch", S->ewald_converge_switch);
	fprintf(S->outptr,"%-24s%d\n", "cluster_move_switch", S->cluster_move_switch);
	fprintf(S->outptr,"%-24s%d\n", "read_config_switch", S->read_config_switch);
	fprintf(S->outptr,"%-24s%d\n", "solvent_in_trj_switch", S->solvent_in_trj_switch);
	fprintf(S->outptr,"%-24s%d\n", "solvent_in_xyz_switch", S->solvent_in_xyz_switch);
	fprintf(S->outptr,"%-24s%d\n", "local_swap_switch", S->local_swap_switch);
	fprintf(S->outptr,"%-24s%d\n", "surf_only_swap_switch", S->surf_only_swap_switch);
}

void infile_param_read(char *inputscript, SYSTEM *S){
	FILE* inptr = fopen(inputscript,"r");
	fscanf(inptr,"run_id\t%s\n", S->run_id);
	fscanf(inptr,"old_seq\t%d\n", &S->old_seq);
	fscanf(inptr,"temperature %lf\n", &S->T);
	fscanf(inptr,"Lx %d\n", &S->Lx);
	fscanf(inptr,"Ly %d\n", &S->Ly);
	fscanf(inptr,"Lz %d\n", &S->Lz);
	fscanf(inptr,"mc_sweeps %li\n", &S->MCSTEPS);
	fscanf(inptr,"write_interval %li\n", &S->WRITE);
	fscanf(inptr,"number_oil %d\nnumber_tail %d\nnumber_head %d\nnumber_water %d\nnumber_vacuum %d\n", &S->No,&S->Nt,&S->Nh,&S->Nw,&S->Nv);
	fscanf(inptr,"eps_ww %lf\neps_ow %lf\neps_oo %lf\neps_ot %lf\neps_tt %lf\n", &S->eww,&S->eow,&S->eoo,&S->eot,&S->ett);
	fscanf(inptr,"q_tail %lf\nq_head %lf\n",&S->qt,&S->qh);
	fscanf(inptr,"cluster_move_freq %d\n",&S->cluster_move_freq);
	fscanf(inptr,"random_seed %d\n",&S->random_seed);
	fscanf(inptr,"real_cutoff_switch %d\n",&S->real_cutoff_switch);
	fscanf(inptr,"energy_check_switch %d\n",&S->energy_check_switch);
	fscanf(inptr,"ewald_converge_switch %d\n",&S->ewald_converge_switch);
	fscanf(inptr,"cluster_move_switch %d\n",&S->cluster_move_switch);
	fscanf(inptr,"read_config_switch %d\n",&S->read_config_switch);
	fscanf(inptr,"solvent_in_trj_switch %d\n",&S->solvent_in_trj_switch);
	fscanf(inptr,"solvent_in_xyz_switch %d\n",&S->solvent_in_xyz_switch);
	fscanf(inptr,"local_swap_switch %d\n",&S->local_swap_switch);
	fscanf(inptr,"surf_only_swap_switch %d\n",&S->surf_only_swap_switch);
	fclose(inptr);
}
