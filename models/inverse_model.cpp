#include "inverse_model.hpp"
#include "models.hpp"
//#include "../methods/methods.hpp"
#include <cmath>
#include <vector>

Real InverseProblem::potential(std::vector<Real> & params)
{
  /*
  TestEquation* eq = new TestEquation(1,x[0]);
  Real *x0 = new Real[1]{1};
  */
  Model* model=get_model(params);
  if ( x_obs.size() == 0) 
    init_x_obs();
  //Real *x1=new Real[model->get_dim()];

  //int N=integration_steps;
  int N=2000; 
  // consider only Heun's method for now... pfff... need to make it better...
  //StepIntegrator *s_integrator=new EulerExplicit(model->get_equation());
  //StepIntegrator *s_integrator=new Heun(model->get_equation());
  StepIntegrator*s_integrator=new ROCK2(model->get_equation());
  //StepIntegrator*s_integrator=new RKC(model->get_equation());

  Real h0 = (model->get_tf()-model->get_t0())/N;
  NumericalMethod *method=0;
  if ( forward_op_type == "det" )
    method= new ClassicalNumericalMethod(*model,h0,*s_integrator);
  else  if ( forward_op_type == "rts_rk")
    method= &AbdulleGaregnaniMethod::createWithUniformPerturbations(*model,N,*s_integrator,1);
  else if ( forward_op_type == "rts_rk_skewed")
    method= &AbdulleGaregnaniMethod::createWithQuasiUniformPerturbations(*model,N,*s_integrator,1);
  else
  { fprintf(stderr,"Type of forward operator(det,rts_rk,rts_rk_skewed) not recognized. Given operator is %s. Aborting...",forward_op_type.c_str()); abort();}
  std::vector<std::pair<Real,Real*>> res= method->integrate(t_obs);
  //method= &AbdulleGaregnaniMethod::createWithUniformPerturbations(*model,N,*s_integrator,1);

  // bayesian step...

  Real P_A_B=0;
  auto t_it = t_obs.begin();auto r_it = res.begin(); auto obs_it = x_obs.begin();
  for (;t_it < t_obs.end();t_it++,r_it++,obs_it++)
  {
    Real *observed_x=*obs_it;
    Real *suggested_x= r_it->second;
    //printf("obs/sug: %lf %lf\n",observed_x[0],suggested_x[0]);
    for (int i=0;i<model->get_dim();i++)
    {
      Real*x=new Real[1];
      NormalDistribution normal_distribution(observed_x[i],noise_variance);
      x[0]=suggested_x[i];
      //P_A_B += -log(normal_distribution.density(x));
      P_A_B += normal_distribution.potential(x);
      delete[] x;
    }
    delete[] suggested_x;
  }
  Real P_A=0;
  auto params_it = params.begin();
  for ( auto prior_it = prior_distribution.begin();prior_it < prior_distribution.end();prior_it++,params_it++)
  {
    RandomVariable * p_dist= *prior_it;
    Real * x = new Real[1];
    x[0]=*params_it;
    P_A += -log((*p_dist).density(x));
    delete [] x;
  }
  
  //printf("%lf %lf %lf\n",params[0],P_A,P_A_B);
  //P_A/=pow(prior_variance,model->get_dim()); //P_A=exp(-0.5*P_A)/sqrt(pow(2*M_PI*prior_variance,dim()));
  // watch out for memory issues, not fixed yet...
  
  
  delete method;
  delete s_integrator;
  delete model;
  if ( include_prior) 
    return P_A+P_A_B;
  else
    return P_A_B;
}
Real InverseProblem::prior_potential(std::vector<Real> & params)
{
  Real P_A=0;
  auto params_it = params.begin();
  for ( auto prior_it = prior_distribution.begin();prior_it < prior_distribution.end();prior_it++,params_it++)
  {
    RandomVariable * p_dist= *prior_it;
    Real * x = new Real[1];
    x[0]=*params_it;
    P_A += -log((*p_dist).density(x));
    delete [] x;
  }
  return P_A;
}

Real InverseProblem::density(std::vector<Real> & params)
{
  return exp(-potential(params));
}


Model* BrusselatorInverse::get_model(const std::vector<Real>& params)
{
  int no_nodes = 30;

  Real *x0 = new Real[no_nodes*2+1];
  // initialize
  for ( int i=0;i<no_nodes*2+1;i++) x0[i]=0;

  Real sqroot2 = sqrt(2);
  // note that we do take bc as known!
  // Apply trigonometric basis
  for ( int i=0;i<no_nodes;i++)
  {
    Real x = 1.0/(no_nodes-1)*i;
    int cut_at=(no_kl_elm+1)/2;
    for ( int j=0;j<cut_at;j++)
    {
      if ( j == 0 ) 
        //x0[i] += params[j]+1 ;
        x0[i] += 1 ;
      else
      {
        Real eigen_val = 1.0 / (2*M_PI*j) / (2*M_PI*j);
        x0[i] += 0*params[j] * sqroot2 * cos(2*M_PI*j*x) *  eigen_val ;
      }
    }
    for ( int j=cut_at;j<no_kl_elm;j++)
    {
      Real eigen_val = 1.0;
      eigen_val/= (2*M_PI*(j-cut_at+1)) * (2*M_PI*(j-cut_at+1));
      x0[i] += params[j] * sqroot2 * sin(2*M_PI*(j-cut_at+1)*x) * eigen_val;
    }
  }

  for ( int i=0;i<no_nodes;i++)
  {
    Real x = 1.0/(no_nodes-1)*i;
    int cut_at=(no_kl_elm+1)/2;
    for ( int j=0;j<cut_at;j++)
    {
      if ( j == 0 ) 
        //x0[i+no_nodes] +=params[j+no_kl_elm]+3 ;
        x0[i+no_nodes] +=3;
      else
      {
        Real eigen_val = 1.0 / (2*M_PI*j) / (2*M_PI*j);
        x0[no_nodes+i] +=0* params[j+no_kl_elm] * sqroot2 * cos(2*M_PI*j*x) *  eigen_val ;
      }
    }
    for ( int j=cut_at;j<no_kl_elm;j++)
    {
      Real eigen_val = 1.0;
      eigen_val /= (2*M_PI*(j-cut_at+1)) * (2*M_PI*(j-cut_at+1));
      x0[no_nodes+i] += params[j+no_kl_elm] * sqroot2 * sin(2*M_PI*(j-cut_at+1)*x) * eigen_val;
    }
  }

  /** apply translation **/

  x0[2*no_nodes]=0;

  Brusselator* eq = new Brusselator(no_nodes);
  Model* brusselator_model=new Model(x0,0,10,eq);
  return brusselator_model;
};
