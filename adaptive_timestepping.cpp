#include "tools/random_variable.hpp"
#include "models/models.hpp"
#include "methods/methods.hpp"
#include "methods/stabilized_methods.hpp"
#include "tools/init.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "config.h"
#include <boost/filesystem.hpp>


extern Model test_equation;
Real e = 2.7182818284;
void calc_exact(Model *model, StepIntegrator* s)
{
  // check if the folder exists...
  // if it does, assume exact solution is already calculated
  boost::filesystem::path dir(model->get_name().c_str());
  if ( boost::filesystem::is_directory (dir) )
    return;
  else
    boost::filesystem::create_directory(dir);
  int N = 800000;
  Real h0 = (model->get_tf()-model->get_t0())/N;
  ClassicalNumericalMethod cnm(*model,h0,*s);
  cnm.set_f_print(model->get_name()+"/exact.txt");
  Real*x=new Real[model->get_dim()];
  cnm.get_final(x);
}
int main(int argc, char *argv[])
{
  /* Try out explicit Euler method */
  Model *model; 
  StepIntegrator* step_integrator;
  NumericalMethodFactory *numerical_method_factory; 
 
  process_command_line(argc,argv,model,step_integrator,numerical_method_factory);
  int n = model->get_dim();
  Real* final_x =new Real[n];

  //for (int i=0;i<n;i++)
    //printf("%lf\n",final_x[i]);

  /* research over N */
 
  calc_exact(model,step_integrator); 
  // do  
  int N=2000;
  int maxi=1;
  Real tol = 1.0;
  //Real tol = 0.00001;
  std::ofstream f(model->get_name()+ "/stats_"+numerical_method_factory->get_name()+".txt");
  f << "i,f_evals"<< std::endl;
  for ( int i=1;i<=maxi;i++)
  {
    for ( int j=0;j<1;j++)
    {
      NumericalMethod* method=numerical_method_factory->get_method(model,N,step_integrator,tol);//(*model,h,step_integrator);
      std::stringstream ss;
      ss << model->get_name()+"/"+numerical_method_factory->get_name();
      ss<<i;
      ss<<"_" << j << "_";
      ss<<".txt";
      method->set_f_print(ss.str());
      method->get_final(final_x);
      if ( j==0)
      {
      std::cout << i << ","<< method->get_f_evals()  << std::endl;
      f << i << ","<< method->get_f_evals()  << std::endl;
      }
      //delete method;
    }
    N*=(1/0.8);
    //N*=2;
    tol/=2;
  }
  /*
  std::stringstream ss,ss2;
  ss << model->get_name()+"/"+numerical_method_factory->get_name();
  ss<<maxi;
  ss<<".txt";
  ss2 << model->get_name()+"/approximation";
  ss2<<".txt";
  if (boost::filesystem::exists (ss2.str()))
    boost::filesystem::remove (ss2.str());
  boost::filesystem::copy_file(ss.str(),ss2.str());
  */
  // copy last file to approximation
	return 0;
}
