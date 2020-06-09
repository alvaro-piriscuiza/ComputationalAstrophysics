#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <omp.h>

int main(){

	double 	*c_real, *c_img;
	int *notdiverged;
	int32_t n,i,j,N;
	int32_t niter,Niter;
	double 	re_min,re_max,img_min,img_max;
	double 	dre,dimg;
	double Zr,Zi,Zr2,Zi2, Zr_tmp, Zi_tmp;
	double Div,Div2;

	FILE *f;
	f = fopen("mandelbrot.dat","w");

	N     = 1000;
	Niter = 100;

	re_min  = -2.0;
	re_max  = 0.5;
	img_min = -1.1;
	img_max = 1.1;

	dre  = (re_max-re_min)/(double)N;
	dimg = (img_max-img_min)/(double)N;


	c_real      = (double *)calloc(N,sizeof(double));
	c_img       = (double *)calloc(N,sizeof(double));
	notdiverged = (int *)calloc(N*N,sizeof(int));
	for (n=0; n<N; n++){
		c_real[n] = re_min + n*dre;
		c_img[n]  = img_min + n*dimg;
	};

	Div = 2.0;
	Div2 = Div*Div;
	#pragma omp parallel for default(none) private(i,j,niter,Zr,Zi,Zr_tmp,Zi_tmp,Zr2,Zi2) shared(N,Niter,Div2,notdiverged,c_img,c_real,f) schedule(static)
	for (i=0; i<N; i++){
		for (j=0; j<N ;j++){
			Zr  = 0.0;
			Zi  = 0.0;
			Zr2 = Zr*Zr;
			Zi2 = Zi*Zi;
			for (niter=0;niter<Niter && ((Zr2+Zi2)<Div2);niter++){
        			Zr_tmp  = Zr2-Zi2 +c_real[i];
        			Zi_tmp  = 2*Zr*Zi +c_img[j];
        			Zr = Zr_tmp;
        			Zi = Zi_tmp;
        
        			Zr2 = Zr*Zr;
        			Zi2 = Zi*Zi;
        
			};
      
      			if (niter==Niter){
				notdiverged[i+(j*N)] = 1;
			}
			else {
				notdiverged[i+(j*N)] = 0;
			};
        		fprintf(f,"%f,%f,%d\n",c_real[i],c_img[j],notdiverged[i+(j*N)]);
		};
	};
	fclose(f);

	return 0;
}
