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

  int p_size = inverse_problem->get_prior_distribution().size() ;
  int d = inverse_problem->get_model(particles[0].get_params())->get_dim();

  // calculate mean and covariance of parameter random variablec
  std::vector<Real> mean(p_size);
  Eigen::MatrixXd cov(p_size,p_size);

  for ( auto it=particles.begin();it!=particles.end();it++)
  {
    for ( int i=0;i<p_size;i++)
      mean[i]+=it->get_weight()*(it->get_params())[i];
  }

  for ( auto it=particles.begin();it!=particles.end();it++)
  {
    std::vector<Real> &x = it->get_params();
    for ( int i=0;i<p_size;i++)
    {
      for ( int j=0;j<p_size;j++)
      {
        cov(i,j)+=it->get_weight()*(x[i]-mean[i])*(x[j]-mean[j]);
      }
    }
  }

  
  /*********************************************************************************
                             *** Main loop ***
   *********************************************************************************/
  for ( int j=1;j<=J;j++)
  {
  fflush(stdout);
  /*********************************************************************************
   * 2. Apply L (Bayes theorem to calculate probabilities)
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
      Real w = exp(-h*potential)*particles[ind].get_weight();
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

    // draw new indices
    
    std::vector<ShuffleSelect> shufle_helper;
    Real sum=0;
    for ( size_t i=0;i<particles.size();i++)
    {
      sum+=fitness_weights[i];
      shufle_helper.push_back(ShuffleSelect(sum,i));
    }
    std::sort(shufle_helper.begin(),shufle_helper.end());

    // Reshuffle hat_x_j1 and particles
    UniformDistribution uniform_distribution(0,1);
    std::vector<Particle> sh_particles; 
    std::vector<Real*> sh_hat_x_j1; 

    for ( size_t i=0;i<particles.size();i++)
    {
      Real u=uniform_distribution();
      auto f_it = std::upper_bound(shufle_helper.begin(),shufle_helper.end(),ShuffleSelect(u,0));
      int index=particles.size()-1; 
      if ( f_it  != shufle_helper.end())
        index=f_it -> index;

      Real *x = new Real[d];
      for ( int j=0;j<d;j++)
        x[j]=particles[index].get_x()[j];
      sh_particles.push_back(Particle(1.0/particles.size(),particles[index].get_params(),particles[index].get_x()));
    }

    particles=sh_particles;
    
    /* update cov and mean matrices */
    if ( j <= J )
    {

      for ( int i=0;i<p_size;i++)
        mean[i]=0;
      for ( int i=0;i<p_size;i++)
        for( int j=0;j<p_size;j++)
          cov(i,j)=0;

      for ( auto it=particles.begin();it!=particles.end();it++)
      {
        for ( int i=0;i<p_size;i++)
          mean[i]+=it->get_weight()*(it->get_params())[i];
      }

      for ( auto it=particles.begin();it!=particles.end();it++)
      {
        std::vector<Real> &x = it->get_params();
        for ( int i=0;i<p_size;i++)
        {
          for ( int j=0;j<p_size;j++)
          {
            cov(i,j)+=it->get_weight()*(x[i]-mean[i])*(x[j]-mean[j]);
          }
        }
      }
      Real a=0.97;
      Real s2 = 1-a*a;
      EigenNormalRandomVariable par_est(s2*cov);

      std::vector<Particle> pert_particles;
      for ( size_t i=0;i<particles.size();i++)
      {
        Eigen::VectorXd param_noise =  par_est();
        std::vector<Real> new_params;
        for ( int j=0;j<p_size;j++)
        {
          new_params.push_back(particles[i].get_params()[j]+param_noise[j]);
        }
        pert_particles.push_back(Particle(particles[i].get_weight(),new_params,particles[i].get_x()));
      }

      particles=pert_particles;

      /* shrink params */
      std::vector<Particle> new_particles;

      for ( auto it = particles.begin();it!=particles.end();it++)
      {
        std::vector<Real>& prev_params=it->get_params();
        std::vector<Real> shrinked_params;
        for ( int i=0;i<p_size;i++)
        {
          shrinked_params.push_back(prev_params[i]*a+(1-a)*mean[i]);
        }
        new_particles.push_back(Particle(it->get_weight(),shrinked_params,it->get_x()));
      }
      particles=new_particles;
    }
  }


  return particles;
}
