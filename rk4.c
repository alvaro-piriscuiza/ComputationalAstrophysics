#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <omp.h>

#define pi 3.14159265359
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

const double T = 365.256*86400; //[s]
const double r = 149.6e+09;     //[m]
const double M = 1.989e+30;     //[kg]
const float  G = 6.67408e-11;
const double v = 2*pi*r/T;

double accx(float,float);
double accy(float,float);

int main() {
	double       *x, *y, *vx, *vy;
	int32_t      n,N;
	double       tstep;
	double       x0, y0, vx0, vy0;
	double	     k1vx,k2vx,k3vx,k4vx,k1x,k2x,k3x,k4x,k1vy,k2vy,k3vy,k4vy,k1y,k2y,k3y,k4y;

	FILE *f;
	f = fopen("rk4.dat","w");

	N     = 10000;
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
	vx[0] = vx0;
	vy[0] = vy0;
	fprintf(f,"%f,%f\n",x[0],y[0]);

	for (n=1;n<N;n++){
		k1x   = vx[n-1];
		k1y   = vy[n-1];

		k1vx  = accx(x[n-1],y[n-1]);
		k1vy  = accy(x[n-1],y[n-1]);

		k2x   = vx[n-1]+(tstep/2.0)*k1vx;
		k2y   = vy[n-1]+(tstep/2.0)*k1vy;

		k2vx  = accx(x[n-1]+(tstep/2.0)*k1x,y[n-1]+(tstep/2.0)*k1y);
		k2vy  = accy(x[n-1]+(tstep/2.0)*k1x,y[n-1]+(tstep/2.0)*k1y);
		
		k3x   = vx[n-1]+(tstep/2.0)*k2vx;
		k3y   = vy[n-1]+(tstep/2.0)*k2vy;
	
		k3vx  = accx(x[n-1]+(tstep/2.0)*k2x,y[n-1]+(tstep/2.0)*k2y);
		k3vy  = accy(x[n-1]+(tstep/2.0)*k2x,y[n-1]+(tstep/2.0)*k2y);

		k4x   = vx[n-1]+tstep*k3vx;
		k4y   = vx[n-1]+tstep*k3vy;

		k4vx  = accx(x[n-1]+tstep*k3x,y[n-1]+tstep*k3y);
		k4vy  = accy(x[n-1]+tstep*k3x,y[n-1]+tstep*k3y);

		vx[n] = vx[n-1]+(tstep/6.0)*(k1vx+2*k2vx+2*k3vx+k4vx);
		vy[n] = vy[n-1]+(tstep/6.0)*(k1vy+2*k2vy+2*k3vy+k4vy);

		x[n]  = x[n-1]+(tstep/6.0)*(k1x+2*k2x+2*k3x+k4x);
		y[n]  = y[n-1]+(tstep/6.0)*(k1y+2*k2y+2*k3y+k4y);

		fprintf(f,"%f,%f\n",x[n],y[n]);
	}
	fclose(f);

	return 0;
}

double accx(float x, float y){
	float R = sqrt(pow2(x)+pow2(y));
	return -(G*M/pow3(R))*x;
}

double accy(float x, float y){
	float R = sqrt(pow2(x)+pow2(y));
	return -(G*M/pow3(R))*y;
}
