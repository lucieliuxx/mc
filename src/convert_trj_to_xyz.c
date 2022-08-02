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

int main(int argc, char *argv[]){
	long int i, FirstFrame, LastFrame;
	int j,k,FrameInterval;
    int x,y,z,type;
    char t;
    char xyz_file[50];

	SYSTEM* S = (SYSTEM *)malloc(sizeof(SYSTEM));
	infile_param_read(argv[1], S);
    sprintf(S->traj_file,"%s",argv[2]);
    sprintf(xyz_file,"%s",argv[3]);
    FirstFrame = atoi(argv[4]);
    LastFrame = atoi(argv[5]);
    FrameInterval = atoi(argv[6]);

	S->trajptr = fopen(S->traj_file,"r");
	FILE * xyzptr = fopen(xyz_file,"a");

    if (S->solvent_in_trj_switch==1){
        S->N = S->Nv + S->No + S->Nw + S->Nh + S->Nt;
    }else{
        S->N = S->Nh + S->Nt + S->No;
    }

    printf("\n");
	while(!feof(S->trajptr)){
        fscanf(S->trajptr,"MC steps = %li\n", &i);
        if (i>=FirstFrame && i<LastFrame){
            if (i%FrameInterval==0){
                printf("\rReading MCsweep %li...", i); 
                if (S->solvent_in_xyz_switch == 1) fprintf(xyzptr,"%d\n\n",S->N);
                else fprintf(xyzptr,"%d\n\n",S->Nt+S->Nh+S->No);
                for (j = 0; j < S->N; j++){
                    fscanf(S->trajptr,"%d %d %d %c\n",&x,&y,&z,&t);
                    if (t=='o'){
                        type = 2;
                        //if (S->solvent_in_xyz_switch == 1) fprintf(xyzptr,"%d %d %d %d\n",type,x,y,z);
                        fprintf(xyzptr,"%d %d %d %d\n",type,x,y,z);
                    }else if (t=='t'){
                        type = 5;
                        fprintf(xyzptr,"%d %d %d %d\n",type,x,y,z);
                    }else if (t=='w'){
                        type = 1;
                        if (S->solvent_in_xyz_switch == 1) fprintf(xyzptr,"%d %d %d %d\n",type,x,y,z);
                    }else if (t=='h'){
                        type = 4;
                        fprintf(xyzptr,"%d %d %d %d\n",type,x,y,z);
                    }else if (t=='v'){
                        type = 3;
                        if (S->solvent_in_xyz_switch==1) fprintf(xyzptr,"%d %d %d %d\n",type,x,y,z);
                    }
                }
            }else{
                for (j = 0; j < S->N; j++) fscanf(S->trajptr,"%d %d %d %c\n",&x,&y,&z,&t);
            }
        }else{
            if (i>=LastFrame) break;
            for (j = 0; j < S->N; j++) fscanf(S->trajptr,"%d %d %d %c\n",&x,&y,&z,&t);
        }
	}
    printf("\n");

	fclose(S->trajptr);
	fclose(xyzptr);
}
