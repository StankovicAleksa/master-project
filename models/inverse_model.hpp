#ifndef INVERSE_MODEL_HPP
#define INVERSE_MODEL_HPP

#include "../tools/random_variable.hpp" 
#include "../config.h"
#include "../methods/methods.hpp"
#include <boost/filesystem.hpp>
#include <cstdio>
#include <vector>


class InverseProblem:public RandomVariable
{
  /* Inverse problem. For now only multidimensional distributions with diagonal covariance matrix
   * are considered.
   * */
  private:
    std::vector<Real> Real_params;
    std::vector<Real> t_obs;
    std::vector<Real*> x_obs;
    std::vector<RandomVariable*> prior_distribution;
    Real noise_variance;
    Model *true_model=0;
    bool include_prior=true;
    int integration_steps;
    std::string forward_op_type;
  public:
    InverseProblem(const std::vector<Real> & parameters, const std::vector< Real>& t_obs ,
        std::vector<RandomVariable*>& prior_distribution, Real noise_variance = 10.0, 
        int integration_steps=200, std::string forward_op_type="det")  // probably I should make prior and posterior noises nicer..
      :Real_params(parameters),t_obs(t_obs),prior_distribution(prior_distribution),noise_variance(noise_variance),
      integration_steps(integration_steps),forward_op_type(forward_op_type) {true_model=0;}
    void init_x_obs()
    {
      Model *m=get_true_model(); // initializes true model in case it is not initialized
      x_obs= construct_x_obs(m);
    }
    void set_include_prior(bool flag){include_prior=flag;}
    std::vector<Real*> construct_x_obs(Model*m)
    {
      int N=integration_steps*2000;
      Real h0 = (m->get_tf()-m->get_t0())/N;
      std::vector<Real*> tmp_x_obs;
      StepIntegrator *s_integrator=new Heun(m->get_equation());
      //StepIntegrator *s_integrator=new ROCK2(m->get_equation());
      NumericalMethod *method= new ClassicalNumericalMethod(*m,h0,*s_integrator);
      std::vector<std::pair<Real,Real*>> res= method->integrate(t_obs);
      for ( auto it=res.begin();it<res.end();it++)
      {
        tmp_x_obs.push_back(it->second);
      } 
      delete method;
      delete s_integrator;
      return tmp_x_obs;
    }
    Real prior_potential(Real *x)
    {
      std::vector<Real> vc;
      for ( size_t i=0;i<Real_params.size();i++)
        vc.push_back(x[i]);
      return prior_potential(vc);
    }
    Real potential(Real *x)
    {
      std::vector<Real> vc;
      for ( size_t i=0;i<Real_params.size();i++)
        vc.push_back(x[i]);
      return potential(vc);
    }
    Real density(Real *x)
    {
      std::vector<Real> vc;
      for ( size_t i=0;i<Real_params.size();i++)
        vc.push_back(x[i]);
      return density(vc);
    }
    virtual ~InverseProblem()
    {
      for (auto it=prior_distribution.begin();it!=prior_distribution.end();it++)
        delete (*it);
      for ( auto it=x_obs.begin();it<x_obs.end();it++)
        delete[] (*it);
      if ( true_model != 0 )
      delete true_model;
    }
    virtual Real density(std::vector<Real> &params); // x is set of params, for which we calculate density using Bayes theorem
    Real potential(std::vector<Real> &params); // x is set of params, for which we calculate density using Bayes theorem
    Real prior_potential(std::vector<Real> &params); // x is set of params, for which we calculate density using Bayes theorem
    virtual Model* get_model(const std::vector<Real>& params)=0;
    int dim(){return Real_params.size() ;}
    std::vector<RandomVariable*>& get_prior_distribution(){return prior_distribution;}
    std::vector<Real>& get_Real_params(){return Real_params;}
    std::vector<Real>& get_t_obs(){return t_obs;}
    std::vector<Real*>& get_x_obs(){return x_obs;}
    Real get_noise_variance(){return noise_variance;}

    Model *get_true_model()
    {
      if ( true_model == 0 )
        true_model=get_model(Real_params);
      return true_model;
    }
};

class TestInverse:public InverseProblem
{
  private:
    Real *x_Real;
    Real *obs;
    TestInverse(const std::vector<Real> & parameters, const std::vector< Real>& t_obs , std::vector<RandomVariable*>& prior_distribution, Real noise_variance, int integration_steps,std::string forward_op_type)  // probably I should make prior and posterior noises nicer..
      :InverseProblem(parameters,t_obs,prior_distribution,noise_variance,integration_steps,forward_op_type){}
  public:
    static TestInverse* get_instance(Real prior_variance=0.5,Real noise_variance=1.0, int integration_steps=5,std::string forward_op_type="det")
    {
      // true parameters
      //std::vector<Real> params; params.push_back(1);
      std::vector<Real> params; params.push_back(1);

      // observation times
      std::vector<Real> t_obs;
      int no_obs=5;
      for ( int i=1;i<=no_obs;i++)
      t_obs.push_back(1.0*i/no_obs);

      // create prior
      std::vector<RandomVariable*> prior_distribution;
      //prior_distribution.push_back(new NormalDistribution(2,prior_variance));
      prior_distribution.push_back(new NormalDistribution(1,prior_variance));
      return new TestInverse(params,t_obs,prior_distribution,noise_variance,integration_steps,forward_op_type);
    }
    virtual Model* get_model(const std::vector<Real>& params)
    {
      Real lambda = params[0];
      Real* test_eq_x0=new Real[1]{1};
      TestEquation *test_eq = new TestEquation(1,lambda); 
      Model* test_equation=new Model(test_eq_x0,0,1,test_eq);
      delete[] test_eq_x0;
      return test_equation;
    };
};


class FitzHughInverse:public InverseProblem
{
  private:
    Real *x_Real;
    Real *obs;
    FitzHughInverse(const std::vector<Real> & parameters, const std::vector< Real>& t_obs ,
        std::vector<RandomVariable*>& prior_distribution, Real noise_variance,int integration_steps, std::string forward_op_type)  // probably I should make prior and posterior noises nicer..
      :InverseProblem(parameters,t_obs,prior_distribution,noise_variance,integration_steps,forward_op_type){}
  public:
    static FitzHughInverse* get_instance(Real prior_variance=0.5,Real noise_variance=1.0, int integration_steps=5,std::string forward_op_type="det")
    {
      // true parameters
      std::vector<Real> params; 
      params.push_back(log(0.2));
      params.push_back(log(0.2));
      params.push_back(log(3));

      // observation times
      std::vector<Real> t_obs;
      Real final_t= get_fitz()->get_tf();
      for ( int i=1;i<=100;i++)
      t_obs.push_back(final_t*i/100.0);

      // create prior
      std::vector<RandomVariable*> prior_distribution;
      prior_distribution.push_back(new NormalDistribution(log(0.2),prior_variance));
      prior_distribution.push_back(new NormalDistribution(log(0.2),prior_variance));
      prior_distribution.push_back(new NormalDistribution(log(3),prior_variance));
      return new FitzHughInverse(params,t_obs,prior_distribution,noise_variance,integration_steps,forward_op_type);
    }
    Model* get_model(const std::vector<Real>& params)
    {
      Real a = exp(params[0]);
      Real b = exp(params[1]);
      Real c = exp(params[2]);

      Real *x0 = new Real[2]{-1,1};
      FitzHughNagumo* eq = new FitzHughNagumo(a,b,c);
      Model* fitz_model=new Model(x0,0,20,eq);
      return fitz_model;
    };
};

class BrusselatorInverse:public InverseProblem
{
  private:
    Real *x_Real;
    Real *obs;
    int no_kl_elm;
    BrusselatorInverse(const std::vector<Real> & parameters, const std::vector< Real>& t_obs ,
        std::vector<RandomVariable*>& prior_distribution, Real noise_variance, int no_kl_elm=5):  // probably I should make prior and posterior noises nicer..
      InverseProblem(parameters,t_obs,prior_distribution,noise_variance),no_kl_elm(no_kl_elm){}
  public:
    static BrusselatorInverse* get_instance(Real prior_variance=0.5,Real noise_variance=1.0, int no_kl_elm=5, int integration_steps=5)
    {
      // true parameters
      // slightly tricky, in case the params are non-constant

      std::vector<Real> params; 
      params.push_back(1.0);
      for ( int i=1;i<no_kl_elm;i++)
        params.push_back(0.0);

      params[(no_kl_elm+1)/2]=1*2*M_PI*2*M_PI;

      params.push_back(3.0);
      for ( int i=1;i<no_kl_elm;i++)
        params.push_back(0.0);

      // observation times
      std::vector<Real> t_obs;
      Real final_t= get_fitz()->get_tf();
      for ( int i=1;i<=10;i++)
      t_obs.push_back(final_t*i/10.0);

      // create prior
      std::vector<RandomVariable*> prior_distribution;
     
      // set prior variance
      for ( int i=0;i<2*no_kl_elm;i++)
        prior_distribution.push_back(new NormalDistribution(0,prior_variance));  


      return new BrusselatorInverse(params,t_obs,prior_distribution,noise_variance,no_kl_elm);
    }

    // trigonometric basis
    Model* get_model(const std::vector<Real>& params);
};

#endif
