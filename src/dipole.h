#ifndef DIPOLE_H
#define DIPOLE_H

int mc_dipole_move(SYSTEM *S);

int * count_types_around_particle(int p1, int * neigh_w, int * neigh_h, int * neigh_t, SYSTEM *S);

#endif
