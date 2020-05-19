/*Fourier transform of exp(-x*x) using FFTW in C*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define xmax 20.0
#define xmin -20.0
#define n_pts 64
#define REAL 0
#define IMAG 1


void fx_data (fftw_complex* fxarr, double dx) /*the array which stores value of the function*/
{
  int i;
  double x;
  
  for (i=0; i<n_pts; i++)
  {
    x = xmin + i*dx;
    fxarr[i][REAL] = exp(-x*x);
    fxarr[i][IMAG] = 0;
  }
}

double fk (double k)					/*analytic fourier transform*/
{
  double x;

    return((1/sqrt(2))*exp(-0.25*k*k));
    
}
int main()									/*main prog starts*/
{
	double k_pts, fk_a;
	int i;
	double dx = (xmax-xmin)/(n_pts-1);		/*spacing*/
	double aft_real, aft_imag;
	
	FILE *mptr;
    mptr = fopen("fftw_2.dat", "w");		/*this file will store fourier transform of fx*/
		
	fftw_complex *fxarr;					/*stores fx values*/
	fftw_complex *fkarr;					/*stores Fourier transform of fx*/
	
	fxarr  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_pts);
	fkarr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_pts);
	
  	fftw_plan p;
	
	p = fftw_plan_dft_1d(n_pts, fxarr, fkarr, FFTW_FORWARD, FFTW_ESTIMATE);/*FFTW algorithm*/

	fx_data(fxarr ,dx);

  	fftw_execute(p); 						/*executing the fourier transform*/

	for (i=0; i<n_pts; i++)
  	{
	   if (i<=n_pts/2 - 1)
	   {
		k_pts = (2*M_PI/(n_pts*dx))*(i);	/*defining the k points*/
		fk_a = fk(k_pts);					/*analytic Fourier transform getting evaluated at given k*/
	   }
	   else
	   {
		k_pts = (2*M_PI/(n_pts*dx))*(i-n_pts);
		fk_a = fk(k_pts);
	   }
	  
	   aft_real = dx*sqrt(1/(2*M_PI))*(cos(k_pts*xmin)*fkarr[i][REAL] + sin(k_pts*xmin)*fkarr[i][IMAG]);/*fft according to our convention*/
	   aft_imag = dx*sqrt(1/(2*M_PI))*(cos(k_pts*xmin)*fkarr[i][IMAG] - sin(k_pts*xmin)*fkarr[i][REAL]);
	   
	   fprintf(mptr,"%e\t %e\t %e\t %e\n",k_pts, aft_real, aft_imag, fk_a);										/*printing the fft of fx into the file*/ 
	  																									/*the required fft is the real part i.e. aft_real*/
 	}
	
	fftw_destroy_plan(p);
	fftw_free(fxarr); 
	fftw_free(fkarr);
	fclose(mptr);
	return (0);
}


