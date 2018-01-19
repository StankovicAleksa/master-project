#include "smc.hpp"
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>

struct ShuffleSelect
{
  ShuffleSelect(double up_bound,int index):up_bound(up_bound),index(index){}
  double up_bound;
  int index;
};
bool operator <( ShuffleSelect s1, ShuffleSelect s2) { return s1.up_bound< s2.up_bound;}

void Stuart_SMC::init_particles()
{
  particles.clear();
  // init particles from prior distribution;
  for ( int i=0;i<no_particles;i++)
  {
    std::vector<Real> v;
    inverse_problem->get_prior_distribution();
    for ( auto it= inverse_problem->get_prior_distribution().begin();it<inverse_problem->get_prior_distribution().end();it++)
    {
      RandomVariable* a=*it;
      RandomVariable& rv=*a;
      v.push_back(rv());
    }
    std::vector<Real> x_tmp; 
    Real*x0 = true_model->get_x0();
    for ( int i=0;i<d;i++)
      x_tmp.push_back( x0[i]);
    particles.push_back(Particle(1.0/no_particles,v,x_tmp));
  }
}
std::vector<Particle> Stuart_SMC::estimate_parameters()
{
  /*********************************************************************************
   * 1. Initialization.
   *********************************************************************************/
  Real h = 1.0/J;
  init_particles();
  inverse_problem->set_include_prior(false);
  int p_size = inverse_problem->get_prior_distribution().size() ;
  int d = inverse_problem->get_model(particles[0].get_params())->get_dim();
  int n=inverse_problem->dim();
  Real *w=new Real[n];
  Real *y=new Real[n];
  Real *x=new Real[n];

  UniformDistribution U(0,1);
  // calculate mean and covariance of parameter random variablec

  Real gamma =0.01;
  NormalDistribution normal_distribution(0,gamma*gamma);
  /*********************************************************************************
                             *** Main loop ***
   *********************************************************************************/
  for ( int j=1;j<=J;j++)
  {
    fflush(stdout);
    printf("J: %d\n",j);
  /*********************************************************************************
   * 2. Apply P (Bayes theorem to calculate probabilities)
   *********************************************************************************/
    if ( j>=1)
    {
       for ( size_t ind=0;ind<particles.size();ind++)
       {
        for ( int i=0;i<n;i++) 
        {
          w[i]=normal_distribution();
        }
        for ( int i=0;i<n;i++) y[i]=(particles[ind].get_params()[i])+w[i];
        for ( int i=0;i<n;i++) x[i]=(particles[ind].get_params()[i]);
        Real py_pot = inverse_problem->potential(y);
        Real px_pot = inverse_problem->potential(x);
        Real py_prior = inverse_problem->prior_potential(y); 
        Real px_prior = inverse_problem->prior_potential(x); 
        Real alpha= rmin(1,exp( h*(j-1)* (-py_pot+px_pot) + (-py_prior+px_prior)  ) );
        if ( U() < alpha ) //accept the move
        {
          for ( int i=0;i<n;i++)
            particles[ind].get_params()[i]=y[i];
        }
       }
    }

  /*********************************************************************************
   * 3. Apply L (Bayes theorem to calculate probabilities)
   *********************************************************************************/

    std::vector<Real> fitness_weights;
    std::vector<Real> t_vec;
    // get observed x
    std::vector<Real*> obs_x=(inverse_problem->get_x_obs());
    Real w_sum =0;
    Real p_sum = 0;
    for ( size_t ind=0;ind<particles.size();ind++)
    {
      Real potential=inverse_problem->potential(particles[ind].get_params());
      Real prior_potential=inverse_problem->prior_potential(particles[ind].get_params());
      Real w = exp(-h*potential-prior_potential)*particles[ind].get_weight();
      p_sum += potential; 
      fitness_weights.push_back(w);
      w_sum +=w;
    }
    //printf("%lf\n",w_sum); 
    t_vec = fitness_weights;
    fitness_weights.clear();

    // normalize fitness weights
    for ( auto it =t_vec.begin();it<t_vec.end();it++)
    {
      fitness_weights.push_back(*it/w_sum);
    }

  /*********************************************************************************
   * 3. Shuffle
   *********************************************************************************/
    
    std::vector<ShuffleSelect> shufle_helper;
    Real sum=0;
    for ( size_t i=0;i<particles.size();i++)
    {
      sum+=fitness_weights[i];
      shufle_helper.push_back(ShuffleSelect(sum,i));
    }
    std::sort(shufle_helper.begin(),shufle_helper.end());

    if ( j<= J)
    {
      // Reshuffle particles
      UniformDistribution uniform_distribution(0,1);
      std::vector<Particle> sh_particles; 
      std::vector<Real*> sh_hat_x_j1; 
      for ( size_t i=0;i<particles.size();i++)
      {
        Real u=uniform_distribution();
        auto f_it = std::upper_bound(shufle_helper.begin(),shufle_helper.end(),ShuffleSelect(u,0));
        int index=f_it -> index;

        sh_particles.push_back(Particle(1.0/particles.size(),particles[index].get_params(),particles[index].get_x()));
      }
      particles=sh_particles;
    }
    else 
    {
      for ( size_t i=0;i<particles.size();i++)
        particles[i].set_weight(fitness_weights[i]);
    }
  }

  delete []x;
  delete []y;
  delete []w;
  return particles;
}
