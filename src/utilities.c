#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "utilities.h"

int p(int x, int L){ 
    // folds displaced particle back into the box !!on the fine grid only!!
    /*
	if (x > L-1) x-=L;
	if (x < 0) x+=L;
	return x;
    */
    return (x+L)%L;
}

int pdist(int x1, int x2, int L){
	// maximum distance in the system is L/2 along each dimension
	int dx = x1 - x2;
	// minimum image distance criterion
	if (L%2 == 0){  // L = even
		if (dx > L/2) dx -= L;
		else if (dx < -(L/2)) dx += L; 
	} else {	// L = odd
		if (dx > (L-1)/2) dx -= L;
		else if (dx < -(L-1)/2) dx += L;
	}
	//return dx - round((dx - sign(dx)*0.5)/S->L)*S->L;
    return dx;
}

int c(int x, int y, int z, int L){
	return p(x,L)*L*L + p(y,L)*L + p(z,L);
}

int * make_vector(int x, int y, int z){
	int * v = (int *)malloc(3*sizeof(int));
	v[0]= x;
	v[1]= y;
	v[2]= z;
	return v;
}

int compare_dx(const void * pa, const void * pb) {
	const int *a = *(const int **)pa;
	const int *b = *(const int **)pb;

	if (a[2]==b[2]){
		if (a[3]==b[3]){
			return a[1] - b[1];
		}else{
			return a[3] - b[3];
		}
	}else{
		return a[2] - b[2];
	}

}

int compare_dy(const void * pa, const void * pb) {
	const int *a = *(const int **)pa;
	const int *b = *(const int **)pb;

	if (a[1]==b[1]){
		if (a[3]==b[3]){
			return a[2] - b[2];
		}else{
			return a[3] - b[3];
		}
	}else{
		return a[1] - b[1];
	}

}

int compare_dz(const void * pa, const void * pb) {
	const int *a = *(const int **)pa;
	const int *b = *(const int **)pb;

	if (a[1]==b[1]){
		if (a[2]==b[2]){
			return a[3] - b[3];
		}else{
			return a[2] - b[2];
		}
	}else{
		return a[1] - b[1];
	}

}

int compare_duplex (const void * a, const void * b)
{
   return ( ((struct data_index_duplex*)b)->data - ((struct data_index_duplex*)a)->data );
}

int compare_int_size (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}
