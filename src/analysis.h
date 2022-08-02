#ifndef ANALYSIS_H
#define ANALYSIS_H

void read_next_frame_from_trj(SYSTEM *S);

void read_next_frame_from_trj_with_added_solvents(SYSTEM *S);

void read_CoM(SYSTEM *S);

void calculate_CoM(SYSTEM *S);

void measure_density_profile_along_x(char type_of_interest, SYSTEM *S);

void measure_distance_distribution(char t1, char t2, SYSTEM *S);

void measure_height_spectrum(SYSTEM *S);

void move_h_to_zero_mean(double **hij, SYSTEM *S);

void get_fourier_spectrum(double **hij, double **hq2, SYSTEM *S);

int slab_divided_across_periodic_box(int clid, SYSTEM *S);

int int_already_in_array(int elem, int * array, int len);

double distance(int i, int j, SYSTEM *S);

#endif
