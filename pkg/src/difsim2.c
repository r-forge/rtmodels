#include <R.h>
#include <Rmath.h>
#include <math.h>
void difsim2( double *v, double *a, double *z, double *s, 
  double *var_v, double *range_z, double *Ter, double *range_Ter, 
  int *N, int *err, double *stepsize, double *t_top, double *t_bot)
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
		error("Upper bound a should be greater than lower bound 0.");
    }
    if (*s <= 0)
    {
		error("s should be greater than 0.");
    }

    double h = 0.05 * (1.0/1000.0);          // gives stepsize h in ms
	if((*stepsize)>0){
		h = (*stepsize);
	}
    double delta       = (*s) * sqrt(h);               // distance up or down
    int i = 0, j = 0, k = 0;

	/*
	for(i=0;i < *N; i++) // zou weg kunnen omdat t_top en t_bot kennelijk op 0 geinitialiseerd worden
    { 
	t_top[i] = -1;
	t_bot[i] = -1;
    }
	*/
	
    for (i=0;i < (*N); i++)
    {
		double ksi			= rnorm(*v,sqrt(*var_v));
		double ter         = *Ter + runif( -(*range_Ter)/2.0, (*range_Ter)/2.0 );
        double position 	= *z + runif(-(*range_z)/2.0, (*range_z)/2.0);
		double pstepdown	= 0.5*(1 - (ksi*sqrt(h)/(*s)));  // N.B. Ratcliff and Tuerlinckx (2002) mistakenly give sqrt(h/s)!
//warning("ksi = %6.3f, ter = %6.3f, position = %6.3f, pstepdown = %6.3f.\n",ksi,ter,position,pstepdown);
//        double hit      =  0;
//        double RT       =  0;

		int m = 0;
		while (position < *a && position > 0) { 
			position += (runif(0.0,1.0)<pstepdown) ? -delta : delta; // use R function to draw uniform random variable
			m++;
//Rprintf("%9.5f %s",position,((m%10==0)?"\n":" "));
		}
		if(position >= *a){					// walker hits top
			t_top[j++] = ter + ((double) m) * h;
		}
		else{								// walker hits bottom
			t_bot[k++] = ter + m * h;
		}
		R_CheckUserInterrupt();			/* allow user to interrupt computations */
    }
	PutRNGstate();						/* write back random number seed		*/
    return;
}
