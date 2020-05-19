/*Fourier transform of sinc function using gsl */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define x_max 50.0
#define x_min -50.0
#define n_pts 256
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

double fx (int i, double dx)						/*function which is to be Fourier transformed*/
{
  double x;

    x = x_min + i*dx;
    if(x!=0)
    {
        return(sin(x)/x);
    }
    else
        return (1.0);
}

int main ()
{
  FILE *mptr;
  mptr = fopen("fftgsl.dat","w");					/*stores the fourier transform of fx*/

  int i; 
  double k_pts, dx;
  double fx_data[2*n_pts];
  double aft_real, aft_imag;
  
  dx = (x_max-x_min)/(n_pts-1);						/*spacing*/
  	
  for (i = 0; i < n_pts; i++)
    {
       REAL(fx_data, i) = fx(i, dx); 				/*array containing value of function*/
       IMAG(fx_data, i) = 0.0;
    }

  gsl_fft_complex_radix2_forward (fx_data, 1, n_pts);/*Fourier transform using gsl fn*/

 for (i=0; i<n_pts; i++)
  	{
	   if (i <=n_pts/2 - 1)
	   {
		k_pts = (2*M_PI/(n_pts*dx))*(i);			/*defining the k points*/
	   }
	   else
	   {
		k_pts = (2*M_PI/(n_pts*dx))*(i-n_pts);
	   }
	   
	   
	   aft_real= dx*sqrt(1/(2*M_PI))*(cos(k_pts*x_min)*REAL(fx_data,i)+sin(k_pts*x_min)*IMAG(fx_data, i));/*fourier transform according to our convention*/
	   aft_imag= dx*sqrt(1/(2*M_PI))*(cos(k_pts*x_min)*IMAG(fx_data,i)-sin(k_pts*x_min)*REAL(fx_data, i));
	   
	   fprintf(mptr,"%e\t%e\t%e\n",k_pts ,aft_real ,aft_imag); /*printing fourier transform into the file*/
	  														/*the required fft is the real part i.e. aft_real*/
 	}
 	
  fclose(mptr);
  return 0;
}
