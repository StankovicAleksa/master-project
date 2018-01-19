#include <cstdio>
#include <cstdlib>
#include "tools/mcmc.hpp"
#include "tools/init.hpp"
#include "models/inverse_model.hpp"
#include "methods/methods.hpp"
#include "config.h"
#include "smc/smc.hpp"

int main()
{
  int no_steps=500;
  int no_particles=2000;
  int J=200;
  //InverseProblem *inverse_problem=TestInverse::get_instance();
  Real noise_variance = 0.05*0.05;

  int no_integration_steps = 500;
  std::string forward_op_type = "det";
  Real prior_variance = 0.01; 
  //InverseProblem* inverse_problem=FitzHughInverse::get_instance(prior_variance,noise_variance,no_integration_steps,forward_op_type);
  InverseProblem* inverse_problem=TestInverse::get_instance(prior_variance,noise_variance,no_integration_steps,forward_op_type);
  inverse_problem->init_x_obs();
  Stuart_SMC smc_solver(no_particles,inverse_problem,J,no_steps);
  std::vector<Particle> result=smc_solver.estimate_parameters();

  std::vector<Real> weights;
  std::vector<Real*> vals;
  int p_size =inverse_problem->get_prior_distribution().size();

  for (auto p = result.begin();p<result.end();p++)
  {
    weights.push_back(p->get_weight());
    Real *x=new Real[p_size];
    for ( int i=0;i<p_size ;i++)
      x[i]=p->get_params()[i];
    //printf("%lf\n",p->get_x()[0]);
    vals.push_back(x);
  }

  DiscreteDirac rand_var(weights,vals,inverse_problem->get_prior_distribution().size());
  for ( int i=0;i<p_size;i++)
  {
    if ( i!=0 ) printf(",");
    printf("x%d",i);
  }
  printf("\n");

  for ( int i=0;i<100000;i++)
  {
    std::vector<Real> sampled = rand_var.sample();
    for (auto it =sampled.begin();it< sampled.end();it++)
    {
      if ( it != sampled.begin() ) printf(",");
      printf("%lf",*it);
    }
    printf("\n");
  }
  return 0;
}
