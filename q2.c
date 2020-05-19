/*Fourier transform of sinx/x using FFTW*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define xmax 50.0
#define xmin -50.0
#define n_pts 256
#define REAL 0
#define IMAG 1


void fx_data (fftw_complex* fxarr, double dx) /*putting the value of the function which is to Fourier transformed into an array*/
{
  int i;
  double x;
  
  for (i=0; i<n_pts; i++)		
  {
    x = xmin + i*dx;
    fxarr[i][REAL] = sin(x)/x;
    fxarr[i][IMAG] = 0;
  }
}

int main()									/*main program starts*/
{
	double k_pts;
	int i;
	double dx = (xmax-xmin)/(n_pts-1);  	/*spacing*/
	double aft_real, aft_imag;
	
	FILE *mptr;
    mptr = fopen("fftw_1.dat", "w");		/*file in which the fourier transform will be stored*/
		
	fftw_complex *fxarr;
	fftw_complex *fkarr;
	
	fxarr  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_pts);	/*This array stores the value of the function*/
	fkarr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_pts);	/*This array stores the value of the dft of function*/
	
  	fftw_plan p;
	
	p = fftw_plan_dft_1d(n_pts, fxarr, fkarr, FFTW_FORWARD, FFTW_ESTIMATE);/*FFTW algorithm*/

	fx_data(fxarr ,dx);			

  	fftw_execute(p); /*fft is computed here*/

	for (i=0; i<n_pts; i++)
  	{
	   if (i<=n_pts/2 - 1)
	   {
		k_pts = (2*M_PI/(n_pts*dx))*(i);	/*defining the k points*/
	   }
	   else
	   {
		k_pts = (2*M_PI/(n_pts*dx))*(i-n_pts);
	   }
	  
	   aft_real = dx*sqrt(1/(2*M_PI))*(cos(k_pts*xmin)*fkarr[i][REAL] + sin(k_pts*xmin)*fkarr[i][IMAG]);/*fft according to our convention*/
	   aft_imag = dx*sqrt(1/(2*M_PI))*(cos(k_pts*xmin)*fkarr[i][IMAG] - sin(k_pts*xmin)*fkarr[i][REAL]);
	   
	   fprintf(mptr,"%e\t %e\t %e\n",k_pts, aft_real, aft_imag);	/*printing the fft of fx into the file*/ 
	  																/*the required fft is the real part i.e. aft_real*/
 	}
	
	fftw_destroy_plan(p);
	fftw_free(fxarr); 
	fftw_free(fkarr);
	fclose(mptr);
	return (0);
}


