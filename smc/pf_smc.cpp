#include "smc.hpp"
#include <vector>
#include <set>
#include <algorithm>
void PF_SMC::init_particles()
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
struct ShuffleSelect
{
  ShuffleSelect(double up_bound,int index):up_bound(up_bound),index(index){}
  double up_bound;
  int index;
};

bool operator <( ShuffleSelect s1, ShuffleSelect s2) { return s1.up_bound< s2.up_bound;}
std::vector<Particle> PF_SMC::estimate_parameters()
{
  /*********************************************************************************
   * 1. Initialization.
   *********************************************************************************/

  init_particles();

  int p_size = inverse_problem->get_prior_distribution().size() ;
  int d = inverse_problem->get_model(particles[0].get_params())->get_dim();

  // calculate mean and covariance of parameter random variablec
  std::vector<Real> mean(p_size);
  Eigen::MatrixXd cov(p_size,p_size);
  /* just in case, init cov and mean */
  for ( int i=0;i<p_size;i++)
  {
    for (int j=0;j<p_size;j++)
      cov(i,j)=0;
    mean[i]=0;
  }


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
   * 1*. Calculate costs for each step
   *********************************************************************************/

  Real prev_t = true_model->get_t0();
  std::vector<Real> costs;
  std::vector<std::pair<Real, Real*>> result;
  std::vector<Real> &t_obs = inverse_problem->get_t_obs();
  for ( auto t = t_obs.begin();t<t_obs.end();t++)
  {
    costs.push_back(*t-prev_t);
    prev_t=*t;
  }
  Real cost_sum = true_model->get_tf()-true_model->get_t0();

  prev_t = true_model->get_t0();
  
  int t_cnt =0;
  int km=0;
  for ( auto t_it = inverse_problem->get_t_obs().begin();t_it < inverse_problem->get_t_obs().end();t_it++,t_cnt++)
  {
    //printf("%d\n",km);
    km++;
    //calculate k, the number of steps for this round
    int k = (*t_it-prev_t)/cost_sum * no_steps; if ( k==0 ) k++;
    Real h0 = cost_sum/k;

    /*********************************************************************************
     * 2. Propagation.
     *********************************************************************************/
    //printf("Propagation ... \n"); 
    // shrink parameters
    std::vector<Particle> new_particles;

    for ( auto it = particles.begin();it!=particles.end();it++)
    {
      std::vector<Real>& prev_params=it->get_params();
      std::vector<Real> shrinked_params;
      for ( int i=0;i<p_size;i++)
        shrinked_params.push_back(prev_params[i]*a+(1-a)*mean[i]);
      new_particles.push_back(Particle(it->get_weight(),shrinked_params,it->get_x()));
    }

    particles=new_particles;

    // predict new states
    std::vector<Real*> hat_x_j1;
    for ( auto it = particles.begin();it!=particles.end();it++)
    {
      Real *x=new Real[d];
      Model* model=inverse_problem->get_model(it->get_params());
      StepIntegrator* s_integrator = new Heun(model->get_equation());
      ClassicalNumericalMethod* mm= new ClassicalNumericalMethod(*model,h0,*s_integrator);
      FixedTimeSteppingMethod *cnm=mm;
      cnm->integrate_part(it->get_x(),x,prev_t,*t_it,k);
      delete mm;
      delete s_integrator;
      /*
      FixedTimeSteppingMethod* cnm= new ClassicalNumericalMethod(*model,h0,*s_integrator);
      jnm->integrate_part(true_model->get_x0(),x,true_model->get_t0(),*t_it,(*t_it-true_model->get_t0())/cost_sum*no_steps+1); 
      */

      hat_x_j1.push_back(x);
    }
    /*********************************************************************************
     * 3. Survival of the fittest
     *********************************************************************************/
    //printf("Survival of the fittest... \n"); 
    // Compute the fitness weights

    std::vector<Real> fitness_weights;
    std::vector<Real> t_vec;
    // get observed x
    Real* obs_x=(inverse_problem->get_x_obs())[t_cnt];

    Real w_sum =0;
    for ( size_t ind=0;ind<particles.size();ind++)
    {
      Real prob=1.0;
      Real *temp_var = new Real[1];
      for( int i=0;i<d;i++)
      {
        NormalDistribution n_dist(obs_x[i],inverse_problem->get_noise_variance());
        temp_var[0]=hat_x_j1[ind][i];
        prob*=n_dist.density(temp_var);
      }
      delete[] temp_var;
      Real w = particles[ind].get_weight() * prob;
      fitness_weights.push_back(w);
      w_sum +=w;
    }
    t_vec = fitness_weights;
    fitness_weights.clear();
    // normalize fitness weights
    for ( auto it =t_vec.begin();it<t_vec.end();it++)
    {
      fitness_weights.push_back(*it/w_sum);
    }

    // draw new indices
    //printf("Drawing new indices\n"); 
    //std::set<ShuffleSelect> shufle_helper;
    std::vector<ShuffleSelect> shufle_helper;
    Real sum=0;
    for ( size_t i=0;i<particles.size();i++)
    {
      sum+=fitness_weights[i];
      //shufle_helper.insert(ShuffleSelect(sum,i));
      shufle_helper.push_back(ShuffleSelect(sum,i));
    }
    std::sort(shufle_helper.begin(),shufle_helper.end());

    //printf("Reshuffling sh_hat_x_j1\n"); 
    // Reshuffle hat_x_j1 and particles
    UniformDistribution uniform_distribution(0,1);
    std::vector<Particle> sh_particles; 
    std::vector<Real*> sh_hat_x_j1; 
    //printf("%d\n",shufle_helper.size());
    for ( size_t i=0;i<particles.size();i++)
    {
      //printf("move: %d\n",i);
      Real u=uniform_distribution();
      auto f_it = std::upper_bound(shufle_helper.begin(),shufle_helper.end(),ShuffleSelect(u,0));
      int index=particles.size()-1; 
      if ( f_it  != shufle_helper.end())
        index=f_it -> index;
      sh_particles.push_back(particles[index]);
      Real *x = new Real[d];
      for ( int j=0;j<d;j++)
        x[j]=hat_x_j1[index][j];
      sh_hat_x_j1.push_back(x); 
    }

    /*********************************************************************************
     * 4. Proliferation
     *********************************************************************************/
    //printf("Proliferation ... \n");
    Real s2 = 1-a*a;
    EigenNormalRandomVariable par_est(s2*cov);

    std::vector<Particle> proliferated_particles;
    for ( size_t i=0;i<sh_particles.size();i++)
    {
      Eigen::VectorXd param_noise =  par_est();
      std::vector<Real> new_params;
      for ( int j=0;j<p_size;j++)
      {
        new_params.push_back(sh_particles[i].get_params()[j]+param_noise[j]);
      }
      /*
      for ( int j=0;j<p_size;j++)
      {
        printf("%lf ",new_params[j]);
      }
      printf("\n");
      */
      // at this point, random time-stepping I guess 
      Model *m=inverse_problem->get_model(sh_particles[i].get_params()); 
      StepIntegrator *s_integrator = new Heun(m->get_equation());
      Real *next_points=new Real[d];
      FixedTimeSteppingMethod* ftsm = & AbdulleGaregnaniMethod::createWithUniformPerturbations(*m,k,*s_integrator,1);
      //FixedTimeSteppingMethod* ftsm = new ClassicalNumericalMethod(*m,(*t_it-prev_t)/k,*s_integrator);
      ftsm->integrate_part(sh_particles[i].get_x(),next_points,prev_t,*t_it,k); 

      /* 
      FixedTimeSteppingMethod* ftsm = new ClassicalNumericalMethod(*m,no_steps,*s_integrator);
      ftsm->integrate_part(true_model->get_x0(),next_points,true_model->get_t0(),*t_it,(*t_it-true_model->get_t0())/cost_sum*no_steps+1); 
      */
      
      delete s_integrator;
      delete ftsm;

      std::vector<Real> v_next_point;
      for ( int i=0;i<d;i++)
        v_next_point.push_back(next_points[i]);
      delete []next_points;
      proliferated_particles.push_back(Particle(sh_particles[i].get_weight(),new_params,v_next_point));
    }
    /*********************************************************************************j1
     * 5. Weight updating
     *********************************************************************************/
    //printf("Weight updating ... \n"); 


    Real sum_weights=0;
    for ( size_t ind=0;ind<proliferated_particles.size();ind++)
    {
      
      Real prob=1.0;
      Real *temp_var = new Real[1];
      for( int i=0;i<d;i++)
      {
        NormalDistribution n_dist(obs_x[i],inverse_problem->get_noise_variance());
        temp_var[0]=sh_hat_x_j1[ind][i];
        prob/=n_dist.density(temp_var);
        temp_var[0]=proliferated_particles[ind].get_x()[i];
        prob*=n_dist.density(temp_var);
      }
      delete[] temp_var;
      Real w = proliferated_particles[ind].get_weight() * prob;
      proliferated_particles[ind].set_weight(w);
      sum_weights+=w;
    }


    // normalize weights
    for ( size_t ind=0;ind<=proliferated_particles.size();ind++)
      proliferated_particles[ind].set_weight(proliferated_particles[ind].get_weight()/sum_weights);

    /*********************************************************************************
     * 6. Update mean/covariance of parameters
     *********************************************************************************/
    //printf("Updating mean and covariance of parameters ... \n");
    particles=proliferated_particles;

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

    prev_t=*t_it;

    /*********************************************************************************
     * ### Free memory. ###
     *********************************************************************************/
    for ( size_t i = 0;i<particles.size();i++)
    {
      delete [] sh_hat_x_j1[i];
      delete [] hat_x_j1[i];
    }
  }
  return particles;
}
