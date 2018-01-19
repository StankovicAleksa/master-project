#ifndef CORE_METHODS_H
#define CORE_METHODS_H
#include "../models/models.hpp"
#include "../tools/random_variable.hpp"
#include "../config.h"
#include <vector>
#include <utility>

class StepIntegrator
{
  public:
    Equation *eq;
    int n;
    bool do_error_estimation;
    Real error;
  public:
    StepIntegrator(Equation *e, bool error_estimation_active = false): eq(e),do_error_estimation(error_estimation_active) {n=eq->get_dimension(); }
    void set_use_error_estimation(bool flag){ do_error_estimation=flag;}
    virtual void operator()(Real *x0, Real *x1,Real h)=0;
    Real get_error(){return error;}
    virtual int get_conv_order() = 0;
    virtual int get_f_evals() = 0;
    virtual ~StepIntegrator(){}
};

class EulerExplicit:public StepIntegrator 
{
  public:
  EulerExplicit(Equation *e):StepIntegrator(e){}
  void operator()(Real *x0, Real *x1, Real h);
  int get_conv_order(){return 1;}
  int get_f_evals();
};

class Heun:public StepIntegrator 
{
  private:
    Real *x_m,*x_t;
    Real *x_tmp;
  public:
  Heun(Equation *e):StepIntegrator(e)
  {
    x_m=new Real[e->get_dimension()];
    x_t=new Real[e->get_dimension()];
    x_tmp = new Real[e->get_dimension()];
  }
  ~Heun(){delete[] x_m;delete[] x_t;delete[] x_tmp;}
  void operator()(Real *x0, Real *x1, Real h);
  int get_f_evals();
  int get_conv_order(){return 2;}
};

class RK3:public StepIntegrator
{
  private:
    Real *k1,*k2,*k3;
    Real *tmp;
  public:
    RK3(Equation *e):StepIntegrator(e)
    {
      k1=new Real[e->get_dimension()];
      k2=new Real[e->get_dimension()];
      k3=new Real[e->get_dimension()];
      tmp = new Real[e->get_dimension()];
    }
    ~RK3(){delete[] k1;delete[] k2;delete[] k3;;delete[] tmp;}
    void operator()(Real *x0, Real *x1, Real h);
    int get_conv_order(){return 3;}
    int get_f_evals();
};

class RK4:public StepIntegrator
{
  private:
    Real *k1,*k2,*k3,*k4;
    Real *tmp;
  public:
    RK4(Equation *e):StepIntegrator(e)
    {
      k1=new Real[e->get_dimension()];
      k2=new Real[e->get_dimension()];
      k3=new Real[e->get_dimension()];
      k4=new Real[e->get_dimension()];
      tmp = new Real[e->get_dimension()];
    }
    ~RK4(){delete[] k1;delete[] k2;delete[] k3;delete[] k4;delete[] tmp;}
    void operator()(Real *x0, Real *x1, Real h);
    int get_conv_order(){return 4;}
    int get_f_evals();
};


class NumericalMethod
{
  protected:
    int p_count=0;
    Model& model;
    Real h0;
    Real cur_t;
    Real *cur_x;
    Real *next_x;
    Real *x_tmp;
    int f_evals=0;
    int n;
    StepIntegrator& base_integrator;
    FILE *f_print=0;
  public:
    NumericalMethod(Model& m, Real h0, StepIntegrator& base_integrator);
    virtual ~NumericalMethod(){delete[] cur_x; delete [] next_x; delete[] x_tmp;}
    static Real calc_h(Model &m, int N) { return (m.get_tf()-m.get_t0())/N; }
    virtual void get_final(Real *x1)=0;  // executes all steps, returns final solution
    virtual std::vector<std::pair<Real,Real*>> integrate(const std::vector<Real>& t_i){fprintf(stderr,"Integration for multiple t_i not implemented.");abort();};  // executes all steps, returns final solution
    Real get_current_t() {return cur_t;}
    virtual Real average_timestep(){return -1.0;}
    void print_pos();
    void print_pos(Real *data);
    void close_f(){ if (f_print!=0) fclose(f_print); }
    int get_f_evals(){return f_evals;}
    void set_f_print(std::string path);
};

class FixedTimeSteppingMethod: public NumericalMethod
{
  public:
    virtual void integrate_part(Real *x_start, Real *x_end,Real t_start,Real t_end, int N)=0;
    void integrate_part(std::vector<Real>x_start, Real *x_end,Real t_start,Real t_end, int N);
    FixedTimeSteppingMethod(Model&m, Real h0, StepIntegrator& base_integrator):NumericalMethod(m,h0,base_integrator){}
    void get_final(Real *x1);
    std::vector<std::pair<Real,Real*>> integrate(const std::vector<Real>& t_i);
    Real average_timestep(){ return h0;}
};

class ClassicalNumericalMethod:public FixedTimeSteppingMethod
{
  public:
    ClassicalNumericalMethod(Model&m, Real h0, StepIntegrator& base_integrator):FixedTimeSteppingMethod(m,h0,base_integrator){};
    void integrate_part(Real *x_start, Real *x_end,Real t_start,Real t_end, int N);
};

class ConradNumericalMethod:public FixedTimeSteppingMethod 
{
  private:
    RandomVariable &step_perturbator;
    ConradNumericalMethod(Model&m, Real h0, StepIntegrator& base_integrator, RandomVariable& step_perturbator):FixedTimeSteppingMethod(m,h0,base_integrator), step_perturbator(step_perturbator){};
  public:
    void integrate_part(Real *x_start, Real *x_end,Real t_start,Real t_end, int N);
    static ConradNumericalMethod& createWithNormalPerturbations(Model &m, int N, StepIntegrator& base_integrator);

    Real average_timestep(){ return h0;}
};


class AbdulleGaregnaniMethod:public FixedTimeSteppingMethod
{
  private:
    Real perturb_step(Real h);
    RandomVariable* step_distribution; 
    //UniformDistribution uniform_distribution; 
    int reduce_on,reduce_factor;
    int no_paths;
    AbdulleGaregnaniMethod(Model&m,Real h0, StepIntegrator& base_integrator,RandomVariable* step_distribution,int no_paths=1):FixedTimeSteppingMethod(m,h0,base_integrator),step_distribution(step_distribution)
    {
    }
  public:
    void integrate_part(Real *x_start, Real *x_end,Real t_start,Real t_end, int N);
    static AbdulleGaregnaniMethod& createWithUniformPerturbations(Model&m,int N, StepIntegrator& base_integrator,int no_paths);
    static AbdulleGaregnaniMethod& createWithQuasiUniformPerturbations(Model&m,int N, StepIntegrator& base_integrator,int no_paths);
    //static AbdulleGaregnaniMethod& createWithLogNormalPerturbations(Model&m,int N, StepIntegrator& base_integrator);
    Real average_timestep(){ return h0;}
    virtual ~AbdulleGaregnaniMethod(){delete step_distribution;}
};



#endif
