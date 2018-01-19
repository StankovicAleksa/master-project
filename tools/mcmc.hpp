#ifndef MCMC_H
#define MCMC_H
#include "../config.h"
#include "../tools/random_variable.hpp"

class MCMC_Sampler
{
  private:
    RandomVariable* rand_var;
    Real gamma;
    int no_iter;
    Real *init_val;
    Real*x;
    Real **sampled_array;
    int n;
    int no_rejects,no_acc;
    NormalDistribution** gamma_normal;
    UniformDistribution U;
  public:
  MCMC_Sampler(RandomVariable* rand_var,Real gamma, int no_iter,Real*init_val):rand_var(rand_var),gamma(gamma),no_iter(no_iter),init_val(init_val),n(rand_var->dim()),U(0,1)
  {
    sampled_array=new Real*[no_iter];
    for ( int i=0;i<no_iter;i++)
    sampled_array[i]=new Real[rand_var->dim()];
    gamma_normal=new NormalDistribution*[n];
    for ( int i=0;i<n;i++)
    {
      gamma_normal[i]=new NormalDistribution(0,gamma*gamma);
    }

  }
  ~MCMC_Sampler()
  {
    for ( int i=0;i<no_iter;i++)
      delete[]sampled_array[i];
    for ( int i=0;i<n;i++)
      delete gamma_normal[i];
    delete[] gamma_normal;
    delete[] sampled_array;
  }
  void sample();
  Real **get_sampled_array(){return sampled_array;}
};


#endif
