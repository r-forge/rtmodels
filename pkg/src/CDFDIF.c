#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>
#define ANSI_C_TARGET 1 /*ansi C (or better) required*/
#define IS_ODD(j)   ((j) & 1 )

/*The following definitions can be changed for compilation without R.h and Rmath.h */
#define powi(x,n) ((double)(R_pow_di((double)(x),(int)(n))))         /* returns:  x^n */
#define ipow(m,n) ((long)(R_pow_di((double)(m),(int)(n))))           /* returns:  m^n */

        /*PARAMETER definitions*/
#define DELTA   1.e-29 /* convergence value fo*/
#define EPSILON 1.e-7 /* value to check devi*/
#define MIN_RT  .001 /* realistic minimum tim*/
#define NR_NU   10 /* nr of Gauss-Hermite quadrature no*/
#define NR_Z    10 /* nr of Gaussian quadrature nodes fo*/
#ifndef PI
  #define PI  3.141592653589793238462643383279
#endif
#define S2  0.01 /* scale factor*/
#define TRUNC_  1.e-16 /* truncation value (se*/
#define V_MAX   5000 /* maximum number of terms in a pa*/
        /*end of PARAMETER definitions*/


void /*FUNCTION*/ cdfdif(double *t, double *x, double par[],  /* Input :: x, t, Par(7)*/
      double *cdf) /* Output :: CDF*/
{
    /*  CDFDif    Cumulative Distribution Function for the Diffusion model
     *   Subroutine that computes the cumulative distribution function of the
     *   diffusion model with random drift, random starting point and random non-
     *   decisional component ter.
     *   The input arguments are : t (the reaction time)
     *                             x (the choice response)
     *                             Par (the parameter vector)
     *                              Par(1) = a (the boundary separation)
     *                              Par(2) = Ter (the mean non-decisional
     *                                            component time)
     *                              Par(3) = eta (the standard deviation of
     *                                            the normal drift distribution)
     *                              Par(4) = z (the mean starting point)
     *                              Par(5) = sZ (the spread of the starting
     *                                           point distribution)
     *                              Par(6) = st (the spread of the non-decisional
     *                                           component time distribution)
     *                              Par(7) = nu (the mean drift rate)
     *
     *   The output argument is : CDF (the value of the cumulative distribution
     *                                 function evaluated at t and x for parameters
     *                                 Par, F_{X,T}(x,t) = Pr(X=x,T<=t) )
     *
     *   Author original FORTRAN 90 source: Francis Tuerlinckx
     *   Last change original FORTRAN 90 source: May 4, 2004
     *
     *
     *   Tuerlinckx, F. (2004a). The efficient computation of the cumulative
     *   distribution and probability density functions in the diffusion model.
     *   Behavior Research Methods, Instruments, and Computers, 36(4), 702-716.
     *
     *   Tuerlinckx, F. (2004b). Tuerlinckx-BRMIC-2004. Retrieved Nov 20, 2006 from
     *   Psychonomic Society Web Archive: http://www.psychonomic.org/ARCHIVE/.
     *
     *   Translation to C: Raoul Grasman
     *   Last change: Nov 23, 2006
     */

    /*  !!!!!!!!!!! defining the constants and variables*/
    long int i, i_, m, m_, v, v_;
    double a, a2, denom, eta, exdif, fact, fnew, gk[NR_NU], gz[NR_Z],
      low, lower_t, nu, p0, p1, prob0, ser, sifa, sl, st, su, sum_hist[3],
      sum_nu, sum_z, sz, ter, upp, upper_t, wk[NR_NU], wz[NR_Z], z,
      z_l, z_u, zzz;

    /*  !!!!!!!!!!! 10 Gauss-Hermite quadrature nodes and weights*/
    static double gr_gh[NR_NU]={-3.43615911883774,-2.53273167423279,
      -1.75668364929988,-1.03661082978951,-0.34290132722370,0.34290132722370,
      1.03661082978951,1.75668364929988,2.53273167423279,3.43615911883774};
    static double w_gh[NR_NU]={0.00000764043286,0.00134364574678,0.03387439445548,
      0.24013861108230,0.61086263373530,0.61086263373530,0.24013861108230,
      0.03387439445548,0.00134364574678,0.00000764043286};

    /*  !!!!!!!!!!! 10 Gaussian quadrature nodes and weights*/
    static double gr_g[NR_Z]={-0.973906528517172,-0.865063366688985,
      -0.679409568299024,-0.433395394129247,-0.148874338981631,0.148874338981631,
      0.433395394129247,0.679409568299024,0.865063366688985,0.973906528517172};
    static double w_g[NR_Z]={0.066671344308688,0.149451349150581,0.219086362515982,
      0.269266719309996,0.295524224714753,0.295524224714753,0.269266719309996,
      0.219086362515982,0.149451349150581,0.066671344308688};

    /*  initialization of the parameters */
    fnew = 0.0e0;

    a = par[0];
    a2 = powi(a,2);
    ter = par[1];
    eta = par[2];
    z = par[3];
    sz = par[4];
    st = par[5];
    nu = par[6];

    z_u = (1 - *x)*z + *x*(a - z) + sz/2.; /* upper bound starting point distribution */
    z_l = (1 - *x)*z + *x*(a - z) - sz/2.; /* lower bound starting point distribution */
    lower_t = ter - st/2.;                 /* lower bound non-decisional component distribution */

    /*   appropriate scaling of GH quadrature nodes and weights based on
     *   normal kernel
     */
    for( i=1; i <= NR_NU; i++ ){
        i_ = i - 1;
        gk[i_] = sqrt(2.0e0)*gr_gh[i_]*eta + nu;
        wk[i_] = w_gh[i_]/sqrt(PI);
        }

    /*   appropriate scaling of the gaussian quadrature nodes and weights */
    for( i=1; i <= NR_Z; i++ ){
        i_ = i - 1;
        gz[i_] = sz/2.*gr_g[i_] + z;
        wz[i_] = w_g[i_];
        }
    /*   is t larger than lower boundary of non-decisional component time
     *   distribution?
     */
    if( *t - lower_t > MIN_RT ){
        /*   compute integrated probability Pr(X=0) */
        sum_z = 0.;
        /*   numerical integration with respect to z0 */
        for( i=1; i <= NR_Z; i++ ){
            i_ = i - 1;
            sum_nu = 0.;
            for( m=1; m <= NR_NU; m++ ){ /* numerical integration with respect to xi */
                m_ = m - 1;
                if( fabs(gk[m_]) >= EPSILON ){ /* gk[m] close to zero? */
                    sum_nu += (exp(-2.*a*gk[m_]/S2) - exp(-2.*gz[i_]*
                      gk[m_]/S2))/(exp(-2.*a*gk[m_]/S2) - 1.)*wk[m_];
                    }
                else if( fabs(gk[m_]) < EPSILON ){
                    sum_nu += (a - gz[i_])/a*wk[m_];
                    }
                }
            sum_z += sum_nu*wz[i_]/2.;
            }
        prob0 = sum_z;
        /*   compute second part cumulative distribution function */
        upper_t = fmin(*t,ter+st/2.);
        /*   integrate probability with respect to t */
        p0 = prob0*(upper_t - lower_t)/st;
        p1 = (1 - prob0)*(upper_t - lower_t)/st;
        /*   is t larger than upper boundary non-decisional component time
         *   distribution?
         */
        if( *t > ter + st/2. ){
            sum_hist[0] = 0.0e0;
            sum_hist[1] = 0.0e0;
            sum_hist[2] = 0.0e0;
            for( v=1; v <= V_MAX; v++ ){ /* infinite series */
                v_ = v - 1;
                sum_hist[0] = sum_hist[1]; /* shift history for convergence check */
                sum_hist[1] = sum_hist[2];
                sum_nu = 0.;
                for( m=1; m <= NR_NU; m++ ){ /* numerical integration with respect to xi */
                    m_ = m - 1;
                    denom = powi(gk[m_],2)/S2 + (powi(PI,2))*(ipow(v,2))*
                      S2/a2;
                    upp = exp((2.**x-1)*z_u*gk[m_]/S2-3.*log(denom)+
                      log(wk[m_])+2.*log(S2));
                    low = exp((2.**x-1)*z_l*gk[m_]/S2-3.*log(denom)+
                      log(wk[m_])+2.*log(S2));
                    sifa = PI*v/a;
                    fact = upp*((2.**x - 1)*gk[m_]*sin(sifa*z_u)/S2 -
                      sifa*cos(sifa*z_u)) - low*((2.**x - 1)*gk[m_]*
                      sin(sifa*z_l)/S2 - sifa*cos(sifa*z_l));
                    exdif = exp((-0.5*denom*(*t-upper_t))+log(1-exp(-0.5*
                      denom*(upper_t-lower_t))));
                    sum_nu += exdif*fact;
                    }
                sum_hist[2] = sum_hist[1] + v*sum_nu;
                if( fabs(sum_hist[0]-sum_hist[1]) < DELTA ){
                    if( fabs(sum_hist[1]-sum_hist[2]) < DELTA ){
                        if( sum_hist[2] > 0. )
                            break;
                        }
                    }
                }
            /*   cumulative distribution function for t and x given Par */
            fnew = (p0*(1 - *x) + p1**x) - sum_hist[2]*4.*PI/(a2*sz*
              st);
            /*   is t less than upper boundary non-decisional component time
             *   distribution?
             */
            }
        else if( *t <= ter + st/2. ){
            sum_nu = 0.;
            for( m=1; m <= NR_NU; m++ ){
                m_ = m - 1;
                if( fabs(gk[m_]) > EPSILON ){
                    sum_z = 0.;
                    for( i=1; i <= NR_Z; i++ ){
                        i_ = i - 1;
                        zzz = (a - gz[i_])**x + gz[i_]*(1 - *x);
                        ser = -((powi(a,3))/((1 - 2**x)*gk[m_]*PI*
                          S2))*sinh(zzz*(1-2**x)*gk[m_]/S2)/(powi(sinh((1-
                          2**x)*gk[m_]*a/S2),2)) + (zzz*a2)/((1 -
                          2**x)*gk[m_]*PI*S2)*cosh((a-zzz)*(1-2**x)*
                          gk[m_]/S2)/sinh((1-2**x)*gk[m_]*a/S2);
                        sum_hist[0] = 0.0e0;
                        sum_hist[1] = 0.0e0;
                        sum_hist[2] = 0.0e0;
                        for( v=1; v <= V_MAX; v++ ){
                            v_ = v - 1;
                            sum_hist[0] = sum_hist[1]; /* shift history for convergence check */
                            sum_hist[1] = sum_hist[2];
                            sifa = PI*v/a;
                            denom = powi(gk[m_],2)/S2 + powi(PI*v,2)*
                              S2/a2;
                            sum_hist[2] = sum_hist[1] + v*sin(sifa*
                              zzz)*exp(-0.5*denom*(*t-lower_t)-2*log(denom));
                            if( ((fabs(sum_hist[0]-sum_hist[1]) <
                              DELTA) && (fabs(sum_hist[1]-sum_hist[2]) <
                              DELTA)) && (sum_hist[2] > 0.) )
                                break;
                            }
                        sum_z += 0.5*wz[i_]*(ser - 4.*sum_hist[2])*
                          (PI*S2)/(a2*st)*exp((2**x-1)*zzz*gk[m_]/
                          S2);
                        }
                    }
                else if( fabs(gk[m_]) <= EPSILON ){
                    sum_hist[0] = 0.0e0;
                    sum_hist[1] = 0.0e0;
                    sum_hist[2] = 0.0e0;
                    su = -(powi(z_u,2))/(12*a2) + (powi(z_u,3))/(12*
                      powi(a,3)) - (powi(z_u,4))/(48*powi(a,4));
                    sl = -(powi(z_l,2))/(12*a2) + (powi(z_l,3))/(12*
                      powi(a,3)) - (powi(z_l,4))/(48*powi(a,4));
                    for( v=1; v <= V_MAX; v++ ){
                        v_ = v - 1;
                        sifa = PI*v/a;
                        denom = powi(PI*v,2)*S2/a2;
                        sum_hist[2] = sum_hist[1] + powi(PI*v,-4)*
                          (cos(sifa*z_l) - cos(sifa*z_u))*exp(-0.5*
                          denom*(*t-lower_t));
                        if( fabs(sum_hist[0]-sum_hist[1]) < DELTA ){
                            if( fabs(sum_hist[1]-sum_hist[2]) < DELTA ){
                                if( sum_hist[2] > 0. )
                                    break;
                                }
                            }
                        }
                    sum_z = 4.*(powi(a,3))/(st*sz*S2)*(sl - su - sum_hist[2]);
                    }
                sum_nu += sum_z*wk[m_];
                }
            fnew = (p0*(1 - *x) + p1**x) - sum_nu;
            }
        /*   is t less than lower boundary of non-decisional component time
         *   distribution?
         */
        }
    else if( *t - lower_t <= MIN_RT ){
        fnew = 0.;
        }
    /*  values that are smaller than 1.E-16 will be set to zero*/
    if( fnew < TRUNC_ )
        fnew = 0.;
    *cdf = fnew;

    return;
} /*end of function*/


void /*FUNCTION*/ cdfdifn(double t[], int *n, double *x, double par[],  /* Input :: x, t(n), n, Par(7)*/
      double cdf[]) /* Output :: cdf(n) */
{
    int i;
    double cdfvalue;
    for(i=0;i<*n;i++){
        cdfdif(&t[i], x, par, &cdfvalue);
        cdf[i] = cdfvalue;
 		R_CheckUserInterrupt();		/* allow R-user to interrupt computations */
   }
    return;
} /*end of function*/
