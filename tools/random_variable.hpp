#ifndef RAND_VAR_H
#define RAND_VAR_H
#include <iostream>
#include <algorithm>
#include <random>
#include "../config.h"
#include <cstdio>
#include <chrono>
#include <vector>
#include <set>
#include "../Eigen/Eigen"
class RandomVariable
{
  private:
  public:
  int virtual dim() = 0;
  Real virtual operator()(){fprintf(stderr,"Sampling from random variable not implemented");abort(); return 1.0;}
  Real virtual density(Real*x){fprintf(stderr,"Density calculation of random variable not implemented");abort();return 1.0;}
  std::vector<Real> sample(){fprintf(stderr,"Multidimensional sampling from random variable not implemented");abort(); return std::vector<Real>();}
  virtual Real potential(Real *x){fprintf(stderr,"Potential calculation not implemented"); abort(); return -1.0; };
};

class NormalDistribution: public RandomVariable
{
  private:
    Real mean,variance;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;
    long long seed; 
  public:
    NormalDistribution(Real mean, Real variance):seed(std::chrono::system_clock::now().time_since_epoch().count())
    {
      generator=std::default_random_engine(seed);
      this->mean=mean;this->variance=variance;
      this->distribution = std::normal_distribution<double>(mean,sqrt(variance));
      //this->distribution = std::normal_distribution<double>(mean,(variance));
    }
    Real operator()() 
    {
      return this->distribution(this->generator);
    }
    Real density(Real *x)
    {
      Real val =0;
      for ( int i=0;i<dim();i++)
      {
        val+=(mean-x[i])*(mean-x[i]);
      }

      val/=pow(variance,dim());
      val=exp(-0.5*val)/sqrt(pow(2*M_PI*variance,dim()));
      return val;
    }
    virtual Real potential(Real *x)
    {
      Real val =0;
      for ( int i=0;i<dim();i++)
      {
        val+=(mean-x[i])*(mean-x[i]);
      }

      val/=pow(variance,dim());
      //val=exp(-0.5*val)/sqrt(pow(2*M_PI*variance,dim()));
      val = -0.5*val;
      Real n_pot = val -log(sqrt(pow(2*M_PI*variance,dim())));
      return -n_pot;
    }
    int dim() { return 1;}
};

class UniformDistribution:public RandomVariable
{
  private:
    Real a,b;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
    long long seed; 
  public:
    UniformDistribution(Real a, Real b):generator(seed),seed(std::chrono::system_clock::now().time_since_epoch().count())
    {
      generator=std::default_random_engine(seed);
      this->a=a;this->b=b;
      //this->distribution = std::uniform_Real_distribution<Real>(a,b);
      this->distribution = std::uniform_real_distribution<double>(0,1);
    }
    Real operator()() 
    {
      return a+this->distribution(this->generator)*(b-a);
    }
    Real density(Real *x) {fprintf(stderr,"Density is not defined for UniformDistribution"); abort(); }
    int dim() { return 1;}
};

class AssymetricUniformDistribution:public RandomVariable
{
  private:
    Real a,b;
    int mult_factor;
    bool left_side;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
    double ratio;
    long long seed; 
  public:
    AssymetricUniformDistribution(Real a, Real b):generator(seed),seed(std::chrono::system_clock::now().time_since_epoch().count())
    {
      generator=std::default_random_engine(seed);
      this->a=a;this->b=b;
      ratio =a/b;
      if ( ratio < 0 ) ratio = -ratio;
      //this->distribution = std::uniform_Real_distribution<Real>(a,b);
      this->distribution = std::uniform_real_distribution<double>(0,1);
    }
    Real operator()() 
    {
      Real choice_sample = this->distribution(this->generator);
      if ( choice_sample < 1/(1+ratio) )
        return this->distribution(this->generator)*a;
      else
        return this->distribution(this->generator)*b;
    }
    Real density(Real *x) {fprintf(stderr,"Density is not defined for UniformDistribution"); abort(); }
    int dim() { return 1;}
};

class SarkaDistribution:public RandomVariable
{

  private:
  public:
    SarkaDistribution() { }
    Real operator()() 
    {
      fprintf(stderr,"Density is not defined for UniformDistribution"); 
      abort(); 
    }
    Real density(Real *x) { return exp(-10*pow(x[0]*x[0]-x[1],2)-pow(x[1]-1/4.0,4)); } // up to multiplicative constant, which is irrelevant for MCMC method..
    int dim() { return 2;}
};

void calc_mean_and_covariance(std::vector<std::vector<Real>> data,std::vector<Real> & mean, std::vector<std::vector<Real>>&covariance);
struct EigenNormalRandomVariable
{
    EigenNormalRandomVariable(Eigen::MatrixXd const& covar)
        : EigenNormalRandomVariable(Eigen::VectorXd::Zero(covar.rows()), covar)
    {}

    EigenNormalRandomVariable(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar)
        : mean(mean)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
        transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    Eigen::VectorXd mean;
    Eigen::MatrixXd transform;

    Eigen::VectorXd operator()() const
    {
        static std::mt19937 gen{ std::random_device{}() };
        static std::normal_distribution<> dist;

        return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](auto x) { return dist(gen); });
    }
};

class DirSample
{
  public:
    DirSample(int ind,Real up_bound):index(ind),up_bound(up_bound){}
    int index;
    Real up_bound;
};
inline bool operator< (DirSample ds1, DirSample ds2){return ds1.up_bound<ds2.up_bound;}
class DiscreteDirac:public RandomVariable
{
  std::vector<Real*> vals;
  std::set<DirSample> indexer;
  int dim_val;
  public:
    DiscreteDirac(std::vector<Real> weights, std::vector<Real*> vals,int dim_val):vals(vals),dim_val(dim_val)
    {
      Real w_sum=0;
      int i=0;
      for ( auto it=weights.begin();it<weights.end();it++,i++)
      {
        w_sum+=*it;
        indexer.insert(DirSample(i,w_sum));
      }
    }
    Real operator()() 
    {
      UniformDistribution ud(0,1);
      Real w=ud();
      auto s=std::upper_bound(indexer.begin(),indexer.end(),DirSample(0,w));
      return vals[s->index][0];
    }
    std::vector<Real> sample()
    {
      UniformDistribution ud(0,1);
      Real w=ud();
      auto s=std::upper_bound(indexer.begin(),indexer.end(),DirSample(0,w));
      std::vector<Real> v_s;
      for ( int i=0;i<dim_val;i++)
        v_s.push_back(vals[s->index][i]);
      return v_s;
    }
    int dim(){return dim_val;}
};
#endif
