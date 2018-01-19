#include "tools/random_variable.hpp"
#include "models/models.hpp"
#include "methods/methods.hpp"
#include "methods/stabilized_methods.hpp"
#include "models/inverse_model.hpp"
#include "tools/init.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "config.h"
#include <boost/filesystem.hpp>


extern Model test_equation;
int main()
{
  /* Try out explicit Euler method */
  Model *model; 
  StepIntegrator* step_integrator;
 
  //model=get_brusselator(300);
  //
  Real prior_variance=0.5;
  Real noise_variance=1.0;
  int no_kl_elm=3;

  InverseProblem* rand_var =BrusselatorInverse::get_instance(prior_variance,noise_variance,no_kl_elm);

  model = rand_var->get_true_model();
  
  //step_integrator=new ROCK2(model->get_equation());
  step_integrator=new EulerExplicit(model->get_equation());

  int N = 200;
  Real h0 = (model->get_tf()-model->get_t0())/N;
  ClassicalNumericalMethod cnm(*model,h0,*step_integrator);
  cnm.set_f_print("test.txt");
  Real*x=new Real[model->get_dim()];
  cnm.get_final(x);
  return 0;
}
