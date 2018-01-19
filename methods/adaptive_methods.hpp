#ifndef ADAPTIVE_METHODS_HPP
#define ADAPTIVE_METHODS_HPP
#include "core_methods.hpp"

class ClassicalAdaptiveScheme:public NumericalMethod
{
  private:
    int f_rho_evals=0;
    Real tol;
    const Real fac = 0.8;
    double sum_h=0;
    int n_h=0;
  public:
    ClassicalAdaptiveScheme(Model&m, Real h0, StepIntegrator& base_integrator, Real tol):NumericalMethod(m,h0,base_integrator),tol(tol)
    {
      base_integrator.set_use_error_estimation(true);
    }
    void get_final(Real *x1);
    int get_f_evals(){return f_evals;}
    int get_f_rho_evalations(){ return f_rho_evals; }
    int get_total_f_evaluations(){ return f_evals+f_rho_evals;}
    

    // we can try different sigmas here!
    Real sigma(Real *x0, Real *x1,Real h);
    double average_timestep(){ return sum_h/n_h;}
};


class GiacomoRandomPathsAdaptiveScheme:public NumericalMethod
{
  private:
    int f_rho_evals=0;
    Real tol;
    const Real fac = 0.8;
    int no_paths;
    Real **cur_points, ** next_points;
    Real *h_path;
    double sum_h=0;
    int n_h=0;
  public:
    GiacomoRandomPathsAdaptiveScheme(Model&m, Real h0, StepIntegrator& base_integrator, Real tol, int no_paths):NumericalMethod(m,h0,base_integrator),tol(tol),no_paths(no_paths)
    {
      base_integrator.set_use_error_estimation(true);
      cur_points = new Real*[no_paths];
      next_points = new Real*[no_paths];
      for ( int i=0;i<no_paths;i++)
      {
        cur_points[i]=new Real[n];
        next_points[i]=new Real[n];
      }
      h_path = new Real[no_paths];
    }
    ~GiacomoRandomPathsAdaptiveScheme()
    {
      for ( int i=0;i<no_paths;i++)
      {
        delete[] cur_points[i];
        delete[] next_points[i];
      }
      delete [] cur_points;
      delete [] next_points;
      delete []h_path;
    }
    void get_final(Real *x1);
    int get_f_evals(){return f_evals;}
    int get_f_rho_evalations(){ return f_rho_evals; }
    int get_total_f_evaluations(){ return f_evals+f_rho_evals;}

    // we can try different sigmas here!
    Real sigma(Real *x0, Real *x1,Real h);
    double average_timestep(){ return sum_h/n_h;}
};

class GiacomoRandomPathsAdaptiveScheme2:public NumericalMethod
{
  private:
    int f_rho_evals=0;
    Real tol;
    const Real fac = 0.8;
    int no_paths;
    Real **cur_points, ** next_points;
    Real *h_path;
    double sum_h=0;
    int n_h=0;
  public:
    GiacomoRandomPathsAdaptiveScheme2(Model&m, Real h0, StepIntegrator& base_integrator, Real tol, int no_paths):NumericalMethod(m,h0,base_integrator),tol(tol),no_paths(no_paths)
    {
      base_integrator.set_use_error_estimation(true);
      cur_points = new Real*[no_paths];
      next_points = new Real*[no_paths];
      for ( int i=0;i<no_paths;i++)
      {
        cur_points[i]=new Real[n];
        next_points[i]=new Real[n];
      }
      h_path = new Real[no_paths];
    }
    ~GiacomoRandomPathsAdaptiveScheme2()
    {
      for ( int i=0;i<no_paths;i++)
      {
        delete[] cur_points[i];
        delete[] next_points[i];
      }
      delete [] cur_points;
      delete [] next_points;
      delete []h_path;
    }
    void get_final(Real *x1);
    int get_f_evals(){return f_evals;}
    int get_f_rho_evalations(){ return f_rho_evals; }
    int get_total_f_evaluations(){ return f_evals+f_rho_evals;}

    // we can try different sigmas here!
    Real sigma(Real *x0, Real *x1,Real h);
    double average_timestep(){ return sum_h/n_h;}
};

class AleksaRandomPathsAdaptiveScheme:public NumericalMethod
{
  private:
    int f_rho_evals=0;
    Real tol;
    const Real fac = 0.8;
    int no_paths;
    Real **cur_points, ** next_points;
    Real *h_path;
    double sum_h=0;
    int n_h=0;
  public:
    AleksaRandomPathsAdaptiveScheme(Model&m, Real h0, StepIntegrator& base_integrator, Real tol, int no_paths):NumericalMethod(m,h0,base_integrator),tol(tol),no_paths(no_paths)
    {
      base_integrator.set_use_error_estimation(true);
      cur_points = new Real*[no_paths];
      next_points = new Real*[no_paths];
      for ( int i=0;i<no_paths;i++)
      {
        cur_points[i]=new Real[n];
        next_points[i]=new Real[n];
      }
      h_path = new Real[no_paths];
    }
    ~AleksaRandomPathsAdaptiveScheme()
    {
      for ( int i=0;i<no_paths;i++)
      {
        delete[] cur_points[i];
        delete[] next_points[i];
      }
      delete [] cur_points;
      delete [] next_points;
      delete []h_path;
    }
    void get_final(Real *x1);
    int get_f_rho_evalations(){ return f_rho_evals; }
    int get_total_f_evaluations(){ return f_evals+f_rho_evals;}
    int get_f_evals(){return f_evals;}

    // we can try different sigmas here!
    Real sigma(Real *x0, Real *x1,Real h);
    double average_timestep(){ return sum_h/n_h;}
};
#endif
