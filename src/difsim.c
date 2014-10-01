#include <R.h>
#include <Rmath.h>
#include <math.h>
void difsim( double *v, double *a, double *z, double *s, double *var_v, double *range_z, double *Ter, double *range_Ter, int *N, int *err, double *t_top, double *t_bot)
{
/*==========================================================================
# generate data from model by random walk approximation (Ratcliff &
# Tuerlinckx, 2002, pp. 441-442)
# Note. It is more efficient to program loops such as those occurring in this
# function in C, and then link the C code to R. ==> this is that code
#
# Raoul Grasman 2005
#==========================================================================*/
	GetRNGstate(); 					/* get random number seed 			*/

	*err = 0;
    if (*a < 0)
    {
        *err = -1;
        return;
    }
    if (*s < 0)
    {
        *err = -2;
	return;
    }

    double h           = 0.05 * (1.0/1000.0);          // gives stepsize h in ms
    double delta       = (*s) * sqrt(h);               // distance up or down
    double pstepdown   = 0.5*(1.0 - ((*v)*sqrt(h)/(*s)));  // N.B. Ratcliff and Tuerlinckx (2002) mistakenly give sqrt(h/s)!
    int i = 0;
    int j = 0;
    int k = 0;

    for(i=0;i < *N; i++) 
    { 
	t_top[i] = -1;
	t_bot[i] = -1;
    }
	
    for (i=0;i < (*N); i++)
    {
        double position = *z;
        double hit      =  0;
        double RT       =  0;
        while (hit == 0)
        {
            double u = runif(0.0,1.0);  	/*  use R function to draw uniform random variable 		*/
            if ( u < pstepdown)
            {
                position = position - delta;  // step down
            }
            else
            {
                position = position + delta; // step up
            }
            if (position >= *a)               // walker hits top
            {
                t_top[j] = RT;
                j        = j + 1;
                hit      = 1;
            }
			if (position <= 0) // walker hits bottom
            {
                t_bot[k] = RT;
                k        = k + 1;
                hit      = 1;
            }
            RT = RT + h;
        }
		R_CheckUserInterrupt();		/* allow user to interrupt computations */
    }
    for(i=0;i<k;i++) t_top[i] = t_top[i] + *Ter;
    for(i=0;i<j;i++) t_bot[i] = t_bot[i] + *Ter;
	PutRNGstate();						/* write back random number seed		*/
    return;
}
