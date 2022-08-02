#ifndef UTILITIES_H
#define UTILITIES_H

int p(int x, int L);
// imposes periodic boundary onto a single coordinate

int pdist(int x1, int x2, int L);
// calculate one-dimensional displacement in the periodic box

int c(int x, int y, int z, int L);
// imposes peridic boundary and unravels x,y,z into x + y*L + z*L*L

int * make_vector(int x, int y, int z);

struct data_index_duplex {
    int data;
    int index;
};

int compare_dx(const void * pa, const void * pb);

int compare_dy(const void * pa, const void * pb);

int compare_dz(const void * pa, const void * pb);

int compare_duplex(const void * a, const void * b);

int compare_int_size(const void * a, const void * b);

#endif
