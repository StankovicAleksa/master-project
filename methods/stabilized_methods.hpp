#ifndef STABILIZED_MODELS_HPP
#define STABILIZED_MODELS_HPP

#include "core_methods.hpp"
class RKC:public StepIntegrator
{
  private:
    const Real nu=0.95;
    static const int maxs=200; // maximal number of stages
    Real *T,*U;
    Real w0;
    Real w1;
    Real** K;  // holding values at different stages
    Real** F;  // holding values of different stages
    int s;
    int f_evals=0;
  public:
    RKC(Equation *e):StepIntegrator(e)
    {
      K = new Real*[maxs+1];
      F = new Real*[maxs+1];
      for ( int i=0;i<=maxs;i++) // probably just s?
      {
        K[i]=new Real[e->get_dimension()];
        F[i]=new Real[e->get_dimension()];
      }
      /* construct Chebyshev polynomial values...*/
      T=new Real[maxs+1];
      U=new Real[maxs+1];
    }
    ~RKC(){}
    int get_conv_order();
    int get_f_evals();
    void calc_number_of_stages(Real *x0, Real h);
    void prepare_chebychev_info();
    void operator()(Real *x0, Real *x1, Real h);
};

class ROCK2:public StepIntegrator
{
  private:
    static const int maxs=200; // maximal number of stages
    Real** K;  // holding values at different stages
    Real** F;  // holding values of different stages
    Real *K_error;
    Real *F_error;
    int f_evals=0;
    int s; // use safest value if no other is specified
    
    // static things
    static const int available_stages_cnt=46; // number of stages for which we have calculated coefficients
    static int ms[46];  // stage number
    // final stage coefficients
    static Real fp1[46];
    static Real fp2[46];
    
    // mu, nu, kappa
    static Real recf[4476];
    // mu, nu, kappa damped
    static Real recf2[184];
    static Real recalph[46];
    static bool static_coefs_constructed;
    

    // coefs to be used in ROCK2 (as in paper, unwrapping)
    static Real **nu;
    static Real **mu;
    static Real **kappa;
    static Real precomputed_alpha[maxs+1];
    static Real precomputed_sigma[maxs+1];
    static Real precomputed_tau[maxs+1];
    static void initialize_static_coefficients();
    void calc_number_of_stages(Real *x0, Real h);
    int get_f_evals();


  public:
    ROCK2(Equation *e):StepIntegrator(e)
    {
      // initialize coefficients if not initialized already
      initialize_static_coefficients();
      K = new Real*[maxs+1];
      F = new Real*[maxs+1];
      for ( int i=0;i<=maxs;i++) // probably just s?
      {
        K[i]=new Real[e->get_dimension()];
        F[i]=new Real[e->get_dimension()];
      }
      K_error = new Real[n];
      F_error = new Real[n];
    }
    ~ROCK2()
    {
      for ( int i=0;i<=maxs;i++) // probably just s?
      {
        delete []K[i];
        delete []F[i];
      }
      delete[] K;
      delete[] F;
      delete[] K_error;
      delete []F_error;
    }
    int get_conv_order();
    void operator()(Real *x0, Real *x1, Real h);
};

#endif
