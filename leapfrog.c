#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <omp.h>

#define pi 3.14159265359
#define pow3(x) ((x)*(x)*(x))

int main(){
	double       *x, *y, *vx, *vy;
	const double T = 365.256*86400; //[s]
	const double r = 149.6e+09;     //[m]
	const double M = 1.989e+30;   //[kg]
	const float  G = 6.67408e-11;
	const double v = 2*pi*r/T;
	int32_t      n,N;
	double       tstep;
	double       x0, y0, vx0, vy0;

	FILE *f;
	f = fopen("leapfrog.dat","w");

	N     = 48;
	tstep = (5*T)/(double)N;

	x  = (double *)calloc(N,sizeof(double));
	y  = (double *)calloc(N,sizeof(double));
	vx = (double *)calloc(N,sizeof(double));
	vy = (double *)calloc(N,sizeof(double));

	x0  = r;
	y0  = 0.0;
	vx0 = 0.0;
	vy0 = v;

	x[0]  = x0;
	y[0]  = y0;
	vx[0] = vx0-(tstep/2.0)*(G*M*x0/pow3(r)); //jumpstart
	vy[0] = vy0-(tstep/2.0)*(G*M*y0/pow3(r)); //jumpstart
	fprintf(f,"%f,%f\n",x[0],y[0]);

	for (n=1;n<N;n++){
		x[n]  = x[n-1]+tstep*vx[n-1];
		y[n]  = y[n-1]+tstep*vy[n-1];
		vx[n] = vx[n-1]-tstep*(G*M*x[n]/pow3(r));
		vy[n] = vy[n-1]-tstep*(G*M*y[n]/pow3(r));
		//fprintf(stderr,"v = %f\n",vx[n]);
		fprintf(f,"%f,%f\n",x[n],y[n]);
	}
	fclose(f);
	
	return 0;
}
