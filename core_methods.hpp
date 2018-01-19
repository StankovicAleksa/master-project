
#ifndef METHODS_H
#define METHODS_H
#include "../models/models.hpp"
#include "../tools/random_variable.hpp"

class StepIntegrator
{
  protected:
    Equation *eq;
    int n;
    bool do_error_estimation;
    double error;
  public:
    StepIntegrator(Equation *e, bool error_estimation_active = false): eq(e),do_error_estimation(error_estimation_active) {n=eq->get_dimension(); }
    virtual void operator()(double *x0, double *x1,double h)=0;
    double get_last_step_error(){return error;}
    virtual int get_conv_order() = 0;
};


class NumericalMethod
{
  protected:
    Model& model;
    double h0;
    double cur_t;
    double *cur_x;
    double *next_x;
    int n;
    StepIntegrator& base_integrator;
  public:
    NumericalMethod(Model& m, double h0, StepIntegrator& base_integrator);
    static double calc_h(Model &m, int N)
    {
      return (m.get_tf()-m.get_t0())/N;
    }
    virtual void get_final(double *x1)=0;  // executes all steps, returns final solution
    double get_current_t() {return cur_t;}
};



class ClassicalNumericalMethod:public NumericalMethod
{
  public:
    ClassicalNumericalMethod(Model&m, double h0, StepIntegrator& base_integrator):NumericalMethod(m,h0,base_integrator){};
    void get_final(double *x1);
};

class AdaptedClassicalNumericalMethod:public NumericalMethod
{
  private:
  AdaptedClassicalNumericalMethod(Model &m, double h0,StepIntegrator& base_integrator):NumericalMethod(m,h0,base_integrator){}
  public:
  static AdaptedClassicalNumericalMethod getInstance(Model&m, StepIntegrator & base_integrator);
  void get_final(double *x1)
  {
  };
};


class ConradNumericalMethod:public NumericalMethod
{
  private:
    RandomVariable &step_perturbator;
    ConradNumericalMethod(Model&m, double h0, StepIntegrator& base_integrator, RandomVariable& step_perturbator):NumericalMethod(m,h0,base_integrator), step_perturbator(step_perturbator){};
  public:
    static ConradNumericalMethod& createWithNormalPerturbations(Model &m, int N, StepIntegrator& base_integrator);
    void get_final(double *x1);
};


class AbdulleGaregnaniMethod:public NumericalMethod
{
  private:
    RandomVariable &step_perturbator; 
    AbdulleGaregnaniMethod(Model&m,double h0, StepIntegrator& base_integrator,RandomVariable& step_perturbator):NumericalMethod(m,h0,base_integrator),step_perturbator(step_perturbator){}
  public:
    static AbdulleGaregnaniMethod& createWithUniformPerturbations(Model&m,int N, StepIntegrator& base_integrator);
    //static AbdulleGaregnaniMethod& createWithLogNormalPerturbations(Model&m,int N, StepIntegrator& base_integrator);
    void get_final(double *x1);
};



class EulerExplicit:public StepIntegrator 
{
  public:
  EulerExplicit(Equation *e):StepIntegrator(e){}
  void operator()(double *x0, double *x1, double h);
  int get_conv_order(){return 1;}
};

class Heun:public StepIntegrator 
{
  private:
    double *x_m,*x_t;
    double *x_tmp;
  public:
  Heun(Equation *e):StepIntegrator(e)
  {
    x_m=new double[e->get_dimension()];
    x_t=new double[e->get_dimension()];
    x_tmp = new double[e->get_dimension()];
  }
  ~Heun(){delete[] x_m;delete[] x_t;delete[] x_tmp;}
  void operator()(double *x0, double *x1, double h);
  int get_conv_order(){return 2;}
};

class RK3:public StepIntegrator
{
  private:
    double *k1,*k2,*k3;
    double *tmp;
  public:
    RK3(Equation *e):StepIntegrator(e)
    {
      k1=new double[e->get_dimension()];
      k2=new double[e->get_dimension()];
      k3=new double[e->get_dimension()];
      tmp = new double[e->get_dimension()];
    }
    ~RK3(){delete[] k1;delete[] k2;delete[] k3;;delete[] tmp;}
    void operator()(double *x0, double *x1, double h);
    int get_conv_order(){return 3;}
};

class RK4:public StepIntegrator
{
  private:
    double *k1,*k2,*k3,*k4;
    double *tmp;
  public:
    RK4(Equation *e):StepIntegrator(e)
    {
      k1=new double[e->get_dimension()];
      k2=new double[e->get_dimension()];
      k3=new double[e->get_dimension()];
      k4=new double[e->get_dimension()];
      tmp = new double[e->get_dimension()];
    }
    ~RK4(){delete[] k1;delete[] k2;delete[] k3;delete[] k4;delete[] tmp;}
    void operator()(double *x0, double *x1, double h);
    int get_conv_order(){return 4;}
};

#endif
