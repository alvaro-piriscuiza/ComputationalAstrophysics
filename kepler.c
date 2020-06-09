//INCLUDES
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <inttypes.h>
#include <libgen.h>

//DEFINES
#define pi 3.14159265359
#define G 6.67408e-11
#define pow3(x) ((x)*(x)*(x))
#define pow2(x) ((x)*(x))
#define MAXSTRING 2048

double **alloc_2D(double **, int, int);
double deg_to_radian(double);
double to_degree(double, double, double);
double get_period(double, double, double);
double get_coord(double, double);
double get_r(double, double, double, double, double, double);
double get_v(double,double);

struct Body{
	double m;
	double r;
	double i;
	double T;
	double v;
};

int main(){
	//Define variables and arrays. These include counters (n,p,r), integration steps (NSTEPS)
	int    n,NSTEPS,NPLANETS;
	int    p,r;
	double vx,vy,vz;
	double tstep;
	double *m_p;
	double *r_p;
	double *i_p;
	double **posX,**posY,**posZ,**vX,**vY,**vZ;
	double d_em,i_em;	

	//Set integration steps and number of planets
	NSTEPS   = 100000; //DUMMY
	NPLANETS = 11;
	tstep    = 86400.0; //DUMMY

	//Allocate memory for arrays
	m_p  = (double *)calloc(NPLANETS,sizeof(double));
	r_p  = (double *)calloc(NPLANETS,sizeof(double));
	i_p  = (double *)calloc(NPLANETS,sizeof(double));
	posX = alloc_2D(posX,NSTEPS,NPLANETS);
	posY = alloc_2D(posY,NSTEPS,NPLANETS);
	posZ = alloc_2D(posZ,NSTEPS,NPLANETS);
	vX   = alloc_2D(vX,NSTEPS,NPLANETS);
	vY   = alloc_2D(vY,NSTEPS,NPLANETS);
	vZ   = alloc_2D(vZ,NSTEPS,NPLANETS);

	//Open results file
	FILE *f;
	f = fopen("kepler_2.dat","wb");

	//Set initial conditions for the solar system Bodies, in the future these will be contained in a datafile or directory.
	const double msun = 1.989e+30;
	m_p[0] = msun;
	m_p[1] = 3.3e23;
	m_p[2] = 4.87e24;
	m_p[3] = 5.9742e24;
	m_p[4] = 6.42e23;
	m_p[5] = 1.90e27;
	m_p[6] = 5.69e26;
	m_p[7] = 8.68e25;
	m_p[8] = 1.03e26;
	m_p[9] = 1.46e22;

	r_p[0] = 0.0;
	r_p[1] = 5.79e10;
	r_p[2] = 1.082e11;
	r_p[3] = 1.496e11;
	r_p[4] = 2.279e11;
	r_p[5] = 7.786e11;
	r_p[6] = 1.433e12;
	r_p[7] = 2.873e12;
	r_p[8] = 4.495e12;
	r_p[9] = 5.906e12;

	i_p[0] = 0.0;
	i_p[1] = deg_to_radian(to_degree(7.0,0.0,16.0));
	i_p[2] = deg_to_radian(to_degree(3.0,23.0,40.0));
	i_p[3] = deg_to_radian(to_degree(0.0,0.0,0.0));
	i_p[4] = deg_to_radian(to_degree(1.0,50.0,59.0));
	i_p[5] = deg_to_radian(to_degree(1.0,18.0,0.0));
	i_p[6] = deg_to_radian(to_degree(2.0,29.0,0.0));
	i_p[7] = deg_to_radian(to_degree(0.77,0.0,0.0));
	i_p[8] = deg_to_radian(to_degree(1.77,0.0,0.0));
	i_p[9] = deg_to_radian(to_degree(17.1,0.0,0.0));

	//Set parameters for the Moon
	d_em    = 3.84e8;
	i_em    = deg_to_radian(to_degree(5.15,0.0,0.0));
	m_p[10] = 7.34767309e22;
	i_p[10] = atan((d_em*sin(i_em))/(r_p[3]+d_em*cos(i_em)));
	r_p[10] = d_em*sin(i_em)/sin(i_p[10]);

	//Create Solar System array
	struct Body b[NPLANETS];
	for(p=0;p<NPLANETS;p++){
		b[p].m = m_p[p];
		b[p].r = r_p[p];
		b[p].i = i_p[p];
		b[p].T = get_period(b[p].m,b[0].m,b[p].r);
		b[p].v = get_v(b[p].r,b[p].T);
	}
	//Correct for NaN velocity value for the Sun
	b[0].v=0.0;

	//Correct values for Moon's obirt around EARTH
	b[10].T=get_period(b[10].m,b[3].m,d_em);
	b[10].v=get_v(d_em,b[10].T)+b[3].v;

	//Set initial positions
	for(p=0;p<NPLANETS;p++){
		posX[0][p] = b[p].r*cos(b[p].i);
		posY[0][p] = 0.0;
		posZ[0][p] = b[p].r*sin(b[p].i);
		if(p==NPLANETS-1){
			fprintf(f,"%f,%f,%f",posX[0][p],posY[0][p],posZ[0][p]);
		}
		else{
			fprintf(f,"%f,%f,%f,",posX[0][p],posY[0][p],posZ[0][p]);		
		}
	}
	fprintf(f,"\n");

	//Set jumpstart in velocity
	for(p=0;p<NPLANETS;p++){
		vx=0.0;
		vy=0.0;
		vz=0.0;
		for(r=0;r<NPLANETS;r++){
			if(p==r){
				vx+=0.0;
				vy+=0.0;
				vz+=0.0;
			}
			else{
				vx+=G*(b[r].m)*get_coord(posX[0][p],posX[0][r])/pow3(get_r(posX[0][p],posX[0][r],posY[0][p],posY[0][r],posZ[0][p],posZ[0][r]));
				vy+=G*(b[r].m)*get_coord(posY[0][p],posY[0][r])/pow3(get_r(posX[0][p],posX[0][r],posY[0][p],posY[0][r],posZ[0][p],posZ[0][r]));
				vz+=G*(b[r].m)*get_coord(posZ[0][p],posZ[0][r])/pow3(get_r(posX[0][p],posX[0][r],posY[0][p],posY[0][r],posZ[0][p],posZ[0][r]));
			}

		}
		vX[0][p] = 0.0-(tstep/2.0)*vx;
		vY[0][p] = b[p].v-(tstep/2.0)*vy;
		vZ[0][p] = 0.0-(tstep/2.0)*vz;
	}

	//Leapfrog Algorithm
	for(n=1;n<NSTEPS;n++){
		//We first compute the positions at each integration step.
		for(p=0;p<NPLANETS;p++){
			posX[n][p] = posX[n-1][p]+tstep*vX[n-1][p];
			posY[n][p] = posY[n-1][p]+tstep*vY[n-1][p];
			posZ[n][p] = posZ[n-1][p]+tstep*vZ[n-1][p];
		}

		//Using the position allows us to find the velocities.
		for(p=0;p<NPLANETS;p++){
			vx=0.0;
			vy=0.0;
			vz=0.0;
			#pragma omp parallel for private (p,r) shared(NPLANETS) default(none) schedule(static) 
			for(r=0;r<NPLANETS;r++){
				if(p==r){
					vx+=0.0;
					vy+=0.0;
					vz+=0.0;
				}
				else{
					vx+=G*(b[r].m)*get_coord(posX[n][p],posX[n][r])/pow3(get_r(posX[n][p],posX[n][r],posY[n][p],posY[n][r],posZ[n][p],posZ[n][r]));
					vy+=G*(b[r].m)*get_coord(posY[n][p],posY[n][r])/pow3(get_r(posX[n][p],posX[n][r],posY[n][p],posY[n][r],posZ[n][p],posZ[n][r]));
					vz+=G*(b[r].m)*get_coord(posZ[n][p],posZ[n][r])/pow3(get_r(posX[n][p],posX[n][r],posY[n][p],posY[n][r],posZ[n][p],posZ[n][r]));
				}
			}
			vX[n][p] = vX[n-1][p]-tstep*vx;
			vY[n][p] = vY[n-1][p]-tstep*vy;
			vZ[n][p] = vZ[n-1][p]-tstep*vz;
		}
	}

	for(n=1;n<NSTEPS;n++){
		for(p=0;p<NPLANETS;p++){
			if (p==NPLANETS-1){
				fprintf(f,"%f,%f,%f",posX[n][p],posY[n][p],posZ[n][p]);
			}
			else{
				fprintf(f,"%f,%f,%f,",posX[n][p],posY[n][p],posZ[n][p]);
			}
		}
		fprintf(f,"\n");
	}

	return 0;
}

double **alloc_2D(double **array, int N1, int N2){
	int i;
	array = (double **)calloc(N1,sizeof(double *));
	for(i=0;i<N1;i++){
		array[i] = (double *)calloc(N2,sizeof(double *));
	}
	return array;
}

double deg_to_radian(double degrees){
	return degrees * pi / 180.0;
}

double to_degree(double degrees, double minutes, double seconds){
	return degrees + (minutes/60.0) + (seconds/3600.0);
}

double get_period(double mass1, double mass2, double a){
	return sqrt(4*pow2(pi)*pow3(a)/(G*(mass1+mass2)));
}

double get_v(double r,double T){
	return 2*pi*r/T;
}

double get_coord(double i, double j){
	return i-j;
}

double get_r(double xi, double xj, double yi, double yj, double zi, double zj){
	double xij, yij,zij;	
	xij = xi-xj;
	yij = yi-yj;
	zij = zi-zj;
	return sqrt(pow2(xij)+pow2(yij)+pow2(zij));
}
