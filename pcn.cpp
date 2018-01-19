#include <cstdio>
#include <cstdlib>
#include "models/inverse_model.hpp"
#include "mcmc/pcn.hpp"
#include "tools/init.hpp"
#include <iostream>

#include "models/inverse_model.hpp"
#include "methods/methods.hpp"
#include "config.h"
#include "smc/smc.hpp"
int main()
{
  Real gamma=2.0;
  Real prior_variance=100.0;
  Real noise_variance=0.02*0.02;
  int no_kl_elm=10;
  int no_params = no_kl_elm*2;

  InverseProblem* rand_var =BrusselatorInverse::get_instance(prior_variance,noise_variance,no_kl_elm);
  rand_var->set_include_prior(false);
  int no_samples=10000;
  //int no_samples=100;
  Real* x0=new Real[no_kl_elm*2];
  for ( int i=0;i<2*no_kl_elm;i++)
    x0[i]=0;
  //x0[0]=1;
  //x0[no_kl_elm]=3;

  //MCMC_Sampler* mcmc_sampler = new MCMC_Sampler(rand_var,gamma,no_samples,x0);
  PCN_Sampler* mcmc_sampler = new PCN_Sampler(rand_var,gamma,no_samples,x0);
  mcmc_sampler->sample();
  Real **arr= mcmc_sampler->get_sampled_array();  
  /*
  FILE *fout = fopen("pcn_res.txt","w");
  for ( int i=0;i<rand_var->get_prior_distribution().size();i++)
  {
    if ( i > 0 ) fprintf(fout,",");
    fprintf(fout, "x%d", i) ;
  }
  
  for ( int i=0;i<no_samples;i++)
  {
    for ( int k=0;k<no_params;k++)
    {
      if ( k!=0 ) fprintf(fout,",");
        fprintf(fout,"%lf",arr[i][k]);

    }
    fprintf(fout,"\n");
  }
  fclose(fout);
  */
  return 0;
}

