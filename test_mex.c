#include <math.h>
#include "mex.h"
#include <complex.h>

#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
//#define IS_REAL_2D_FULL_DOUBLE(P) (mxIsComplex(P)   && !mxIsSparse(P) && mxIsDouble(P))


typedef struct
{
	double real;
	double imaginary;

} signal;



int  power_of_two( int x)
{
    int y=x,next;
    while (((y % 2) == 0) && y > 1)  /* While x is even and > 1 */
    {
        y /= 2;
    }

    if(y != 1)
    {
	//		printf("y!=1");
		    next = pow(2, ceil(log(x)/log(2)));
  //          printf("next%d\n",next);

        return next;
    }
    else
		    return x;

}


signal *multiply(signal *c, signal *d)
{
	signal *a=( signal*)malloc( sizeof(signal));
	a->real=(c->real*d->real)-(c->imaginary*d->imaginary);
	a->imaginary=(c->real*d->imaginary)+(c->imaginary*d->real);
	return a;


}


signal *add(signal *c, signal *d)
{
	signal *a=( signal*)malloc( sizeof(signal));
	a->real= c->real + d->real;
	a->imaginary= c->imaginary + d->imaginary ;
	return a;
}



signal *subtract(signal *c, signal *d)
{
	signal *a=( signal*)malloc( sizeof(signal));
	a->real= c->real - d->real;
	a->imaginary= c->imaginary - d->imaginary ;
	return a;
}


signal  *fast_fourier_transform( signal *s, int signal_size, int step_size)
{
 signal *temp= ( signal*)malloc( sizeof(signal));
 signal *temp2= ( signal*)malloc( sizeof(signal));
 signal *temp3= ( signal*)malloc( sizeof(signal));
 signal *signal_array = ( signal*)malloc( sizeof(signal) * signal_size );
 signal *odd_indices, *even_indices;
 if ( signal_size == 1 )
   {
        signal_array[0] = s[0];
        return signal_array;
   }

  even_indices = fast_fourier_transform( s,signal_size/2, step_size * 2 );
  odd_indices = fast_fourier_transform( s + step_size, signal_size/2, step_size * 2 );

    for ( int k = 0; k < signal_size / 2; k++ ) {

 //       odd_indices[k].real = (odd_indices[k].real)*cos(2 * M_PI * k / signal_size)    -   (odd_indices[k].imaginary )* sin(2 * M_PI * k / signal_size);
 //       odd_indices[k].imaginary = (odd_indices[k].real)*sin(2 * M_PI * k / signal_size) +   (odd_indices[k].imaginary )* cos(2 * M_PI * k / signal_size);
              temp2->real=cos(2 * M_PI * k / signal_size);
              temp2->imaginary=sin(2 * M_PI * k / signal_size);

              temp3->real=odd_indices[k].real;
              temp3->imaginary=odd_indices[k].imaginary;
              temp=multiply(temp2,temp3);

              odd_indices[k].real=temp->real;
              odd_indices[k].imaginary=temp->imaginary;
          }


    for ( int k = 0; k < signal_size / 2; k++ ) {

            //even temp
              temp2->real=even_indices[k].real;
              temp2->imaginary=even_indices[k].imaginary;

              //odd temp
              temp3->real=odd_indices[k].real;
              temp3->imaginary=odd_indices[k].imaginary;

              temp=add(temp2,temp3);

             signal_array[k].real=temp->real;
             signal_array[k].imaginary=temp->imaginary;


    //     signal_array[k].real = even_indices[k].real + odd_indices[k].real;
    //    signal_array[k].imaginary = even_indices[k].imaginary + odd_indices[k].imaginary;


              temp=subtract(temp2,temp3);

             signal_array[k+ signal_size/2].real=temp->real;
             signal_array[k+ signal_size/2].imaginary=temp->imaginary;

             //signal_array[k+ signal_size/2].real = even_indices[k].real - odd_indices[k].real;
            //signal_array[k+ signal_size/2].imaginary = even_indices[k].imaginary - odd_indices[k].imaginary;


    }


    free(temp);
    free(temp2);
    free(temp3);

    free( odd_indices );
    free( even_indices );

    return signal_array;
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

#define A_IN prhs[0]
#define MAX_SIZE prhs[1]// first passed data

#define SIZE nrhs       //number of passed data array
double *output_signal, *A, *p;
int M, N, m, n,i;



//if(!IS_REAL_2D_FULL_DOUBLE(A_IN)) /* Check A */
//  mexErrMsgTxt("A must be a real 2D full double complex array.");
//else

//p = mxGetScalar(P_IN); /* Get p */

M = mxGetM(A_IN); /* Get the dimensions of A */
N = mxGetN(A_IN);

double *Ar=mxGetPr(A_IN); //real part
double *Ai=mxGetPi(A_IN); //imaginary part

double *max_size_real=mxGetPr(MAX_SIZE); //real part

///////////////////////////////////////////////////////////////
double max_size=*max_size_real;
printf("max_size %f\n",max_size);

////////////////////////////////////////
double temp1,temp2;

plhs[0] = mxCreateDoubleMatrix(M, 1, mxCOMPLEX); /* Create the output matrix */
plhs[1] = mxCreateDoubleMatrix(M, 1, mxCOMPLEX); /* Create the output matrix */

double *out_real = mxGetPr(plhs[0]); /* Get the pointer to the data of B */
double *out_imaginary = mxGetPi(plhs[0]); /* Get the pointer to the data of B */
double *inverse_real = mxGetPr(plhs[1]); /* Get the pointer to the data of B */
double *inverse_imaginary = mxGetPi(plhs[1]); /* Get the pointer to the data of B */

printf("m>> %d\n",M);

int x=power_of_two(M);
printf("x>> %d\n",x);
signal *s = (signal*)malloc(sizeof(signal)* x);
signal *output = (signal*)malloc(sizeof(signal)* x);
signal *ifft = (signal*)malloc(sizeof(signal)* x);





//C:\Users\Markeb\Desktop\dsp2


for (int i = 0; i < M; i++)
	{


		s[i].real =Ar[i] ;
		s[i].imaginary = Ai[i];

	}




   if(M!=x)
    {

			  for(int i=M; i<x; i++)
        {
            	s[i].real =0 ;
              s[i].imaginary =0;

		    }

    }


output= fast_fourier_transform( s,  x, 1)	;

printf("signal fft\n");


for (int i = 0; i < 5; i++)
	{
    printf("%.7f ", output[i].imaginary);
    printf("+ %.7fi\n", output[i].real);

	}

	for(int i=0;i<x;i++){
	    temp1=output[i].real;
	    temp2=output[i].imaginary;
	    output[i].imaginary= temp1;
	    output[i].real=temp2;
	}

	////////////////////////////////////////////////////////////////////////////

	for(i=0;i<x;i++){
	out_real[i]=output[i].real;
	out_imaginary[i]=output[i].imaginary;
	}

/////////////////////////////////////////////////////////////////////////
printf("zzzzzzz");
ifft= fast_fourier_transform( output,  x, 1)	;



for(int i=0;i<x;i++){
    temp1=ifft[i].real;
    temp2=ifft[i].imaginary;

    ifft[i].imaginary= temp1/x;
    ifft[i].real=temp2/x;
}

printf("signal ifft\n");

for (int i = 0; i < 5; i++)
	{

        printf("%.7f ", ifft[i].real);
		printf("%.7fi\n", ifft[i].imaginary);


	}


	for(i=0;i<x;i++){
	inverse_real[i]=ifft[i].real;
	inverse_imaginary[i]=ifft[i].imaginary;
}
return;
}
