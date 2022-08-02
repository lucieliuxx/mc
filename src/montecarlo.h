#ifndef MONTECARLO_H
#define MONTECARLO_H

int mc_neighbor_swap(SYSTEM *S);

void accept_swap_of_identical_types_for_diffusion(int p1, int p2, SYSTEM *S);

void apply_swap_without_pbc(int p1, int p2, SYSTEM *S);

int get_neighbor(int i, SYSTEM *S);

int get_random(int i, SYSTEM *S);

double * site_swap_get_energy(int p1, int p2, SYSTEM *S);

void generate_site_site_swap(int p1, int p2, SYSTEM *S);

int accept_move(SYSTEM *S);

int reject_move(SYSTEM *S);

void copy_particle_into_Elem2(int p, SYSTEM *S);

#endif
