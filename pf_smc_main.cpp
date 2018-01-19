#include <cstdio>
#include <cstdlib>
#include "tools/mcmc.hpp"
#include "tools/init.hpp"
#include "models/inverse_model.hpp"
#include "methods/methods.hpp"
#include "config.h"
#include "smc/smc.hpp"

int main(int argc, char *argv[])
{

  int no_steps=1000;
  int no_particles=10000;
  Real a=0.97;
  //Real a=1.0;
  //InverseProblem *inverse_problem=TestInverse::get_instance();
  Real noise_variance = 0.05*0.05;

  int no_integration_steps = 500;
  std::string forward_op_type = "det";
  Real prior_variance = 0.01; 
  //InverseProblem* inverse_problem=FitzHughInverse::get_instance(prior_variance,noise_variance,no_integration_steps,forward_op_type);
  //InverseProblem *inverse_problem=FitzHughInverse::get_instance();
  InverseProblem* inverse_problem=TestInverse::get_instance(prior_variance,noise_variance,no_integration_steps,forward_op_type);
  inverse_problem->init_x_obs();
  
  PF_SMC smc_solver(no_particles,inverse_problem,a,no_steps);
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
    vals.push_back(x);
  }
  DiscreteDirac rand_var(weights,vals,inverse_problem->get_prior_distribution().size());
  for ( int i=0;i<p_size;i++)
  {
    if ( i!=0 ) printf(",");
    printf("x%d",i);
  }
  printf("\n");
  for ( int i=0;i<10000;i++)
  {
    std::vector<Real> sampled = rand_var.sample();
    for (auto it =sampled.begin();it< sampled.end();it++)
    {
      if ( it != sampled.begin() ) printf(",");
      printf("%lf",*it);
    }
    printf("\n");
  }
  
  /*
  std::vector<Real*>&calc_vals  = inverse_problem->get_x_obs();

  std::vector<Real> params;
  double a1,b,c;
  a1=-0.055201;
  b=-0.998293;
  c=0.837150;
  printf("%lf,%lf,%lf\n",exp(a1),exp(b),exp(c));
  a1=log(0.2);
  b=log(0.2);
  c=log(3);
  params.push_back(a1);
  params.push_back(b);
  params.push_back(c);
  std::vector<Real*>calc_vals2 = inverse_problem->construct_x_obs(inverse_problem->get_model(params));
  for ( size_t i=0;i<vals.size();i++)
  {
    //printf("%lf,%lf",calc_vals[i][0]-calc_vals2[i][0],calc_vals[i][1]-calc_vals2[i][1]);
    //printf("\n");
  }

  for ( auto part = result.begin(); part < result.end();part++)
  {
    //if ( part->get_weight() > 0.001)
    //printf("%.16lf,%lf\n",part->get_weight(),part->get_params()[0]);
  }
  */

  return 0;
}
