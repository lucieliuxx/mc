#ifndef CLUSTER_H
#define CLUSTER_H

int mc_cluster_move(SYSTEM *S);

void compile_ids_moved(int *Uold, int nUold, int *VU, int nVU, int nS, SYSTEM *S);

double e_real_cluster_and_solvent(SYSTEM *S);

PARTICLE * make_list_of_cluster_charges(int *Uold, int nUold, int nS, SYSTEM *S);

void reset_cluster_arrays(SYSTEM *S);

void reset_neighbor_arrays(SYSTEM *S);

int compile_solvents_moved(int *dr, int select_cl, int *VU, int nS, int ** UVr, int **VUr, SYSTEM * S);

void displace_solvents(int nVU, int ** UVr, int **VUr, SYSTEM * S);

int generate_clusters(SYSTEM *S);

void get_neigh_lists(SYSTEM * S);

void merge_lists(SYSTEM * S);

void merge_pair(int p1, int p2, SYSTEM * S);

int check_overlap(SYSTEM * S);


#endif
