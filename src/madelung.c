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

// This program checks the Madelung constant using current the complete ewald calculation.
// If the Madelung constant calculated here is the correct answer, you will still need to
// 1) Check agreement between the complete ewald to the optimized ewald in simulation
// 2) Tune the ewald summation parameters for simulation

int main(int argc, char *argv[]){
	int i,j,k,m,p1,p2;
	int NaCl = 0;
	int PtQs = 1;

	SYSTEM* S = (SYSTEM *)malloc(sizeof(SYSTEM));
	S->T = 1;
	S->Lx = 2;
	S->Ly = 2;
	S->Lz = 2;
	S->N = S->Lx*S->Ly*S->Lz;
	S->rcut = S->Lx/2;
	S->eww = 0.; 
	S->eow = 0.;
	S->eoo = 0.;
	S->nkx=12;
	S->nky=12;
	S->nkz=12;
    S->kc = (12+1)*2*M_PI/(double)(S->Lx);
	S->a = 36;
    S->Nh = 8;
    S->Nt = 0;
	double U_self=0;

	S->Elem = malloc(S->N*sizeof(PARTICLE));
    S->id = (int ***)malloc((S->Lx)*sizeof(int **));
    for (i=0; i<S->Lx; i++){
        S->id[i] = (int **)malloc((S->Ly)*sizeof(int *));
        for (j=0; j<S->Ly; j++){
            S->id[i][j] = (int *)malloc((S->Lz)*sizeof(int));
        }
    }
	// m = 0;
	// for (i=0;i<S->Lx;i++){
	// 	for (j=0;j<S->Ly;j++){
	// 		for (k=0;k<S->Lz;k++){
	// 			set_site_vac(i,j,k,m,S);
	// 			m++;
	// 		}
	// 	}
	// }
	
	// for (int d=1;d<=S->Lx/2;d++){

	// 	p1=(S->id)[c(0,0,0,S->Lx)];
	// 	(S->Elem)[p1].q = 1.;
	// 	(S->Elem)[p1].t = 'h';
	// 	p2=(S->id)[c(0,0,d,S->Lx)];
	// 	(S->Elem)[p2].q = -1.;
	// 	(S->Elem)[p2].t = 'h';

	// 	for (i=0;i<S->N;i++) U_self += sqrt(S->a/M_PI)*(S->Elem)[i].q*(S->Elem)[i].q;

	// 	// printf("p1 = %d; p2 = %d\n",p1,p2);
	// 	printf("d/L = %.2f, (Er+Ek+Es)/N = %.3f\n",(double)d/(double)S->Lx, (e_real_uncut(S)+e_recip_abinitio(S)-U_self)/2.);

	// 	(S->Elem)[p1].q = 0.;
	// 	(S->Elem)[p1].t = 'v';
	// 	(S->Elem)[p2].q = 0.;
	// 	(S->Elem)[p2].t = 'v';
	// }

	set_site_head(0,0,0,0,S);
	S->Elem[0].q = 1.;
	set_site_head(1,0,0,1,S);
	S->Elem[1].q = -1.;
	set_site_head(0,1,0,2,S);
	S->Elem[2].q = -1.;
	set_site_head(0,0,1,3,S);
	S->Elem[3].q = -1.;
	set_site_head(0,1,1,4,S);
	S->Elem[4].q = 1.;
	set_site_head(1,0,1,5,S);
	S->Elem[5].q = 1.;
	set_site_head(1,1,0,6,S);
	S->Elem[6].q = 1.;
	set_site_head(1,1,1,7,S);
	S->Elem[7].q = -1.;

	for (i=0;i<S->N;i++) U_self += sqrt(S->a/M_PI)*S->Elem[i].q*S->Elem[i].q;

	printf("-2*(Er+Ek+Es)/ N = %.3f\n",-2*(e_real_uncut(S)+e_recip_abinitio(S)-U_self)/8.);
	printf("Madelung Constant = 1.748\n");
}

