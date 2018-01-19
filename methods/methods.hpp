#ifndef METHODS_HPP
#define METHODS_HPP
#include "core_methods.hpp"
#include "adaptive_methods.hpp"
#include "stabilized_methods.hpp"
#include "../models/models.hpp"

class NumericalMethodFactory
{
  private:
    std::string method_name;
  public:
    NumericalMethodFactory(std::string method_name ):method_name(method_name){}
    NumericalMethod* get_method(Model *model,int N,StepIntegrator*step_integrator,double tol)
    {
      Real h0 = (model->get_tf()-model->get_t0())/N;
      if ( method_name == "fixed_timestepping" )
      {
        return new ClassicalNumericalMethod(*model,h0,*step_integrator);
      }
      else if ( method_name == "abdulle_garegnani")
      {
        //return &AbdulleGaregnaniMethod::createWithUniformPerturbations(*model,N,*step_integrator);
        return &AbdulleGaregnaniMethod::createWithQuasiUniformPerturbations(*model,N,*step_integrator,100);
      }
      else if ( method_name == "conrad")
      {
        return &ConradNumericalMethod::createWithNormalPerturbations(*model,N,*step_integrator);
      }
      else if ( method_name == "classical_adaptive")
      {
        return new ClassicalAdaptiveScheme(*model,h0,*step_integrator,tol);
      }
      else if ( method_name == "r_adaptive")
      {
        return new GiacomoRandomPathsAdaptiveScheme(*model,h0,*step_integrator,tol,5);
      }
      else if ( method_name == "aleksa_adaptive")
      {
        return new AleksaRandomPathsAdaptiveScheme(*model,h0,*step_integrator,tol,1000);
      }
      else if ( method_name == "giacomo2")
      {
        return new GiacomoRandomPathsAdaptiveScheme2(*model,h0,*step_integrator,tol,1000);
      }
      else
      {
        fprintf(stderr,"Numerical integration scheme is not recognized");
        abort();
      }
    }
    std::string get_name(){ return method_name; }
};
#endif
