#ifndef PCN_HPP
#define PCN_HPP
#include "../config.h"
#include "../tools/random_variable.hpp"
#include "../models/inverse_model.hpp"
class PCN_Sampler
{
  private:
    InverseProblem* rand_var;
    Real gamma;
    int no_iter;
    Real *init_val;
    Real*x;
    Real **sampled_array;
    int n;
    int no_rejects,no_acc;
    NormalDistribution** gamma_normal;
    UniformDistribution U;
    Real beta; 
  public:
    PCN_Sampler(InverseProblem* rand_var,Real gamma, int no_iter,Real*init_val):rand_var(rand_var),gamma(gamma),no_iter(no_iter),init_val(init_val),n(rand_var->dim()),U(0,1)
    {
      //beta = 0.64;
      beta = 0.1;
      sampled_array=new Real*[no_iter];
      for ( int i=0;i<no_iter;i++)
      sampled_array[i]=new Real[rand_var->dim()];
      gamma_normal=new NormalDistribution*[n];
      for ( int i=0;i<n;i++)
      {
        gamma_normal[i]=new NormalDistribution(0,gamma*gamma);
      }

    }
    void sample();
    Real **get_sampled_array(){return sampled_array;}
};


#endif

