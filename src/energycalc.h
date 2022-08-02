#ifndef ENERGYCALC_H
#define ENERGYCALC_H

double e1p_ising(int pt, SYSTEM *S);

double e1p_elec(int pt, SYSTEM *S);

double e1p(int pt, SYSTEM *S);

double e2p(int p1, int p2, SYSTEM *S);

double etype(char s1, char s2, SYSTEM *S);

double e_recip(PARTICLE * charges_old, PARTICLE * charges_new, int ncharges, SYSTEM *S);
// list of old coords of moved charges (x1, y1, z1, q1) ... (xn, yn, zn, qn) 
// list of new coords of moved charges (x1', y1', z1', q1) ... (xn', yn', zn', qn) 
// length of these lists
// system params

/***********************************/
//Functions for checking the energy calculation

void error_check(SYSTEM * S);

double e_real_uncut(SYSTEM *S);

double e_recip_abinitio(SYSTEM *S);

void ewald_converge(SYSTEM * S);

void measure_tau(SYSTEM *S);

#endif
