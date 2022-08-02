#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <time.h>
#include "utilities.h"
#include "initialize.h"
#include "energycalc.h"

double e1p_ising(int pt, SYSTEM *S){
	// ising intn between site i and 6 nearest neighbors
	int m,d,d1,d2;
	double en=0; 
    
    // nearest neighbors
	for (d=-1;d<=1;d+=2){
		m = S->id[p(S->Elem[pt].x+d,S->Lx)][S->Elem[pt].y][S->Elem[pt].z];
		en += etype(S->Elem[pt].t,S->Elem[m].t,S);
		m = S->id[S->Elem[pt].x][p(S->Elem[pt].y+d,S->Ly)][S->Elem[pt].z];
		en += etype(S->Elem[pt].t,S->Elem[m].t,S);
		m = S->id[S->Elem[pt].x][S->Elem[pt].y][p(S->Elem[pt].z+d,S->Lz)];
		en += etype(S->Elem[pt].t,S->Elem[m].t,S);
	}
    // next-nearest neighbors
    /*    
    for (d1=-1;d1<=1;d1+=2){
        for (d2=-1;d2<=1;d2+=2){
            m = S->id[p(S->Elem[pt].x+d1,S->Lx)][p(S->Elem[pt].y+d2,S->Ly)][S->Elem[pt].z];
            en += (1/sqrt(2))*etype(S->Elem[pt].t,S->Elem[m].t,S);
            m = S->id[p(S->Elem[pt].x+d1,S->Lx)][S->Elem[pt].y][p(S->Elem[pt].z+d2,S->Lz)];
            en += (1/sqrt(2))*etype(S->Elem[pt].t,S->Elem[m].t,S);
            m = S->id[S->Elem[pt].x][p(S->Elem[pt].y+d1,S->Ly)][p(S->Elem[pt].z+d2,S->Lz)];
            en += (1/sqrt(2))*etype(S->Elem[pt].t,S->Elem[m].t,S);
        }
    }
   */
	return en;
}

double e1p_elec(int p1, SYSTEM *S){
	// screened electrostatics to all other charges
	double q1, en=0; 
    q1 = S->Elem[p1].q;

	if (S->rcut == S->Lx/2){
		if (fabs(q1)>0.0001){ 
            #pragma omp parallel for
			for (int dx=-S->rcut; dx<S->rcut; dx++){
				for (int dy=-S->rcut; dy<S->rcut; dy++){
					for (int dz=-S->rcut; dz<S->rcut; dz++){
						if (dx==0&&dy==0&&dz==0) continue;
						double d = sqrt((double)(dx*dx + dy*dy + dz*dz));
						if (d>(double)S->rcut) continue;

						int p2 = S->id[p(S->Elem[p1].x+dx,S->Lx)][p(S->Elem[p1].y+dy,S->Ly)][p(S->Elem[p1].z+dz,S->Lz)];
                        double q2 = S->Elem[p2].q;
                        char t2 = S->Elem[p2].t; 

						if (fabs(q2)<0.0001) continue;

                        #pragma omp atomic 
						en += q1*q2*erfc(sqrt(S->a)*d)/d;
					}
				}
			}
		}
	} else {
		if (fabs(q1)>0.0001){ 
            #pragma omp parallel for
			for (int dx=-S->rcut;dx<=S->rcut;dx++){
				for (int dy=-S->rcut;dy<=S->rcut;dy++){
					for (int dz=-S->rcut;dz<=S->rcut;dz++){
						if (dx==0&&dy==0&&dz==0) continue;

						double d = sqrt((double)(dx*dx + dy*dy + dz*dz));
						if (d>(double)S->rcut) continue;

						int p2 = S->id[p(S->Elem[p1].x+dx,S->Lx)][p(S->Elem[p1].y+dy,S->Ly)][p(S->Elem[p1].z+dz,S->Lz)];
                        double q2 = S->Elem[p2].q;
                        char t2 = S->Elem[p2].t; 

						if (fabs(q2)<0.0001) continue;

                        #pragma omp atomic 
						en += q1*q2*erfc(sqrt(S->a)*d)/d;
					}
				}
			}
		}
	}

	return en;
}

double e1p(int pt, SYSTEM *S){
	// combined ising and electrostatics 1 particle intns
	return e1p_ising(pt,S) + e1p_elec(pt,S);
}

double e2p(int p1, int p2, SYSTEM *S){
	// intn between site i and site j, depending on the distance r_ij
	double q1,q2,r,en=0;
	int dx = pdist(S->Elem[p1].x,S->Elem[p2].x,S->Lx);
	int dy = pdist(S->Elem[p1].y,S->Elem[p2].y,S->Ly);
	int dz = pdist(S->Elem[p1].z,S->Elem[p2].z,S->Lz);

	r = sqrt((double)(dx*dx + dy*dy + dz*dz));
	if (r<1.1) en += etype(S->Elem[p1].t,S->Elem[p2].t, S);
    //if (r<1.42 && r>1.40) en += (1/sqrt(2))*etype(S->Elem[p1].t,S->Elem[p2].t, S);

	if (fabs(S->Elem[p1].q*S->Elem[p2].q)>0.0001){
        q1 = S->Elem[p1].q;
        q2 = S->Elem[p2].q;
		en += q1*q2*erfc(sqrt(S->a)*r)/r;
	}

	return en;
}

double etype(char t1, char t2, SYSTEM *S){ 
	double en;
	if (t1=='h'||t1=='w'){
		if (t2=='h'||t2=='w') en = S->eww;
		else if(t2=='t'||t2=='o') en = S->eow;
		else en = 0;
	}else if(t1=='t'){
		if (t2=='h'||t2=='w') en = S->eow;
		else if(t2=='t') en = S->ett;
		else if(t2=='o') en = S->eot;
		else en = 0;
	}else if(t1=='o'){
		if (t2=='h'||t2=='w') en = S->eow;
		else if(t2=='t') en = S->eot;
		else if(t2=='o') en = S->eoo;
		else en = 0;
	}else{
		en = 0;
	}
	return en;
}

double e_recip(PARTICLE * charges_old, PARTICLE * charges_new, int ncharges, SYSTEM *S){
    double U_k=0;
	//double kx, ky, kz, k2, degen; 
	//double kr_new, kr_old;

    #pragma omp parallel for
    for (int i=0; i<S->total_vecs_in_ksphere; i++){
        double kx = S->ksphere[i][0];
        double ky = S->ksphere[i][1];
        double kz = S->ksphere[i][2];
        double k2 = S->ksphere[i][3];
        double degen = S->ksphere[i][4];

        S->rho_k2[i] = S->rho_k[i];
        for (int m=0; m<ncharges; m++){
            double kr_new = kx*(double)charges_new[m].x+ky*(double)charges_new[m].y+kz*(double)charges_new[m].z;
            double kr_old = kx*(double)charges_old[m].x+ky*(double)charges_old[m].y+kz*(double)charges_old[m].z;
            S->rho_k2[i] += charges_old[m].q*(cexp(-I*kr_new)-cexp(-I*kr_old));
        }
        
        double dU_k = degen*0.5/((double)(S->Lx*S->Ly*S->Lz)) * (4*M_PI/k2) * pow(cabs(S->rho_k2[i]), 2) * exp(- k2/(4*S->a));

        #pragma omp atomic
        U_k += dU_k;
	}

    /*
    for (int n=0; n<ncharges; n++){
        printf("n = %d: %d %d %d %.3lf\n", n, charges_old[n].x, charges_old[n].y, charges_old[n].z, charges_old[n].q);
    }
    printf("\n");
    for (int n=0; n<ncharges; n++){
        printf("n = %d: %d %d %d %.3lf\n", n, charges_new[n].x, charges_new[n].y, charges_new[n].z, charges_new[n].q);
    }
    */
    /*
    for (int i=0; i<S->total_vecs_in_ksphere; i++){
        printf("k = (%.2lf, %.2lf, %.2lf) rho(k) = %.2f + %.2f i\n", S->ksphere[i][0], S->ksphere[i][1], S->ksphere[i][2], creal(S->rho_k2[i]), cimag(S->rho_k2[i]));
    }
    */

	return U_k;
}

/***********************************/
//Functions for checking the energy calculation

double e_recip_abinitio(SYSTEM *S){
	double dkx = 2*M_PI/(double)(S->Lx);
	double dky = 2*M_PI/(double)(S->Ly);
	double dkz = 2*M_PI/(double)(S->Lz);
	double kx, ky, kz, k2, U_k=0, kr;
	double complex rho_k_elem;
	int i,j,k,m,dx,dy,dz;

	for (i=0; i<2*S->nkx+1; i++){
		for (j=0; j<2*S->nky+1; j++){
			for (k=0; k<2*S->nkz+1; k++){
				if ((i==S->nkx)&&(j==S->nky)&&(k==S->nkz)) continue;
				kx = dkx*(double)(i-S->nkx);
				ky = dky*(double)(j-S->nky);
				kz = dkz*(double)(k-S->nkz);
				k2 = kx*kx + ky*ky + kz*kz;
				// add up all delta_S(k) from individual atoms
                if (sqrt(k2) < S->kc){
                    rho_k_elem=0.0+0.0*I;
                    for (m=0;m<S->Nh+S->Nt;m++){
                        kr = kx*(double)(S->Elem)[m].x+ky*(double)(S->Elem)[m].y+kz*(double)(S->Elem)[m].z;
                        rho_k_elem+=cexp(-I*kr)*(S->Elem)[m].q;
                    }
                    U_k += 0.5/((double)(S->Lx*S->Ly*S->Lz)) * (4*M_PI/k2) * pow(cabs(rho_k_elem),2) * exp(- k2/(4*S->a));
                    //printf("k = (%.2lf, %.2lf, %.2lf) rho(k) = %.2f + %.2f i\n", kx, ky, kz, creal(rho_k_elem), cimag(rho_k_elem));
                }
			}
		}
	}

    //printf("id = %d: %d %d %d %.3lf\n",0, S->Elem[0].x, S->Elem[0].y, S->Elem[0].z, S->Elem[0].q);

	return U_k;
}

double e_real_uncut(SYSTEM *S){
	int temp_rcut = S->rcut;
	double U_r = 0;

	S->rcut = S->Lx/2;
	for (int pt=0; pt<S->N; pt++) U_r += 0.5*e1p(pt,S);
	S->rcut = temp_rcut;

	return U_r;
}

void error_check(SYSTEM * S){
	double U_r=0;

	for (int pt=0; pt<S->N; pt++) U_r += 0.5*e1p(pt,S);

	//double U_r_uncut = 0;
	double U_r_uncut = e_real_uncut(S);
	double U_k_abinitio = e_recip_abinitio(S);
    //fprintf(S->outptr, "--------Energy Calculation Checks--------\n");
    //fprintf(S->outptr, "%-15s%-15s%-15s%-15s%-15s%-15s\n","E_r(cut)","E_r(full)","RealError","E_k(opt)","E_k(abi)","RecipError");
	fprintf(S->outptr, "%-15.3lf%-15.3lf%-15.3lf%-15.3lf%-15.3lf%-15.3lf\n", U_r, U_r_uncut, fabs(U_r - U_r_uncut),S->E_k, U_k_abinitio, fabs(S->E_k - U_k_abinitio));
}

void ewald_converge(SYSTEM * S){
	int i, dx, dy, dz;
	double U_k, U_r, U_s;
	double temp_a = S->a, temp_kc = S->kc;
    int temp_rcut = S->rcut;
    double curr_sqrta, prev_E;

    fprintf(S->outptr, "--------Ewald Sum Convergence (nc = %d %d %d)--------\n", S->nkx, S->nky, S->nkz);
    //fprintf(S->outptr,"#sqrta, U_k, U_r, U_coulomb, (dU/da)*da\n");
    fprintf(S->outptr,"%-10s %-10s %-10s %-10s %-10s\n","sqrta","U_k","U_r","U_coulomb","(dU/da)*da");
    prev_E = 0;
    for (curr_sqrta=sqrt(temp_a)-0.02;curr_sqrta<=sqrt(temp_a)+0.015;curr_sqrta+=0.01){
        S->a = curr_sqrta*curr_sqrta;
        S->kc = 2*2.0 * sqrt(S->a);
        S->rcut = (int)(2.0/sqrt(S->a));
        
        U_r = 0;
        for (int pt=0; pt<S->N; pt++) U_r += 0.5*e1p_elec(pt,S);

        U_k = e_recip_abinitio(S);

        U_s = 0;
        for (i=0;i<S->N;i++) U_s += sqrt(S->a/M_PI)*(S->Elem)[i].q*(S->Elem)[i].q;
        

        fprintf(S->outptr,"%-10.3f %-10.3f %-10.3f %-10.3f %-10.3f\n", curr_sqrta, U_k, U_r, U_k + U_r - U_s, U_k + U_r - U_s - prev_E);
        prev_E = U_k + U_r - U_s;
        fclose(S->outptr);
        S->outptr = fopen(S->output_file,"a");
    }
    fprintf(S->outptr, "\n\n");

	S->a = temp_a;
	S->kc = temp_kc;
	S->rcut = temp_rcut;
}


void measure_tau(SYSTEM *S){
    S->rcut = S->Lx/2;
    double en=0;
    time_t real_ti, real_tf, recip_ti, recip_tf;
    double tau_R, tau_F;
    double NR = 0.5*(double)((S->Nh+S->Nt)*(S->Nh+S->Nt-1));
    double NF = (double)(S->total_vecs_in_ksphere * (S->Nt+S->Nh));

    time(&real_ti);
    for (int i=0; i<100;i++){
        #pragma omp parallel for
        for (int p1=0; p1<S->Nh+S->Nt; p1++){
            for (int p2=p1+1; p2<S->Nh+S->Nt; p2++){
                #pragma omp atomic
                en += e2p(p1, p2, S);
            }    
        }
    }
    time(&real_tf);
    tau_R = difftime(real_tf,real_ti)/NR;

    time(&recip_ti);
    for (int i=0; i<100;i++){
        #pragma omp parallel for
        for (int i=0; i<S->total_vecs_in_ksphere; i++){
            double kx = S->ksphere[i][0];
            double ky = S->ksphere[i][1];
            double kz = S->ksphere[i][2];
            double k2 = S->ksphere[i][3];
            double degen = S->ksphere[i][4];

            double complex rhok_elem=0.0+0.0I;
            for (int m=0; m<S->Nt+S->Nh; m++){
                double kr_new = kx*(double)S->Elem[m].x+ky*(double)S->Elem[m].y+kz*(double)S->Elem[m].z;
                double kr_old = kx*(double)S->Elem[m].x+ky*(double)S->Elem[m].y+kz*(double)S->Elem[m].z;
                rhok_elem += S->Elem[m].q*(cexp(-I*kr_new)+cexp(-I*kr_old));
                // note this is the wrong formula. using it to ensure non-zero rhok_elem, for time calculation only!
            }
            
            double dU_k = degen*0.5/((double)(S->Lx*S->Ly*S->Lz)) * (4*M_PI/k2) * pow(cabs(rhok_elem), 2) * exp(- k2/(4*S->a));

            #pragma omp atomic
            en += dU_k;
        }
    }
    time(&recip_tf);
    tau_F = difftime(recip_tf,recip_ti)/NF;

    fprintf(S->outptr, "tau_R = %lf, tau_F = %lf\ntau_R/tau_f = %lf\n",tau_R, tau_F, tau_R/tau_F);
}
