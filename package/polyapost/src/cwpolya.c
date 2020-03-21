
#include <R.h>
#include <R_ext/Utils.h>
#include "polyapost.h"

void cwpolya(double *x, double *w, int *nin, int *Nin)
{
	int i, j, k;
	int n = nin[0];
	int N = Nin[0];
	double a;

	GetRNGstate();

	for (i=0; i<N; i++) {
		a = w[n-1] * unif_rand();
		j = 0;
		while(a > w[j])
			j++;
		for(k=j; k <n; k++)
			w[k] = w[k] + 1;
		x[n+i] = x[j];
		R_CheckUserInterrupt();
	}

	PutRNGstate();
}
