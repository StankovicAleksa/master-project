#include <cstdio>
#include <cstdlib>
#include "tools/mcmc.hpp"

int main()
{
  Real gamma=0.15;
  RandomVariable* rand_var =new SarkaDistribution();
  int no_samples=500000;
  Real* x0=new Real[2];
  x0[0]=1;
  x0[1]=-1;
  //MCMC_Sampler* mcmc_sampler = new MCMC_Sampler(rand_var,gamma,no_samples,x0);
  MCMC_Sampler* mcmc_sampler = new MCMC_Sampler(rand_var,gamma,no_samples,x0);
  mcmc_sampler->sample();
  Real **arr= mcmc_sampler->get_sampled_array();  
  FILE *fout = fopen("mcmc_out.txt","w");
  fprintf(fout,"x1,x2\n");
  for ( int i=no_samples/2;i<no_samples;i++)
  {
    fprintf(fout,"%lf,%lf\n",arr[i][0],arr[i][1]);
  }
  return 0;
}
