#include <cstdio>
#include <cstdlib>
#include "tools/mcmc.hpp"
#include "tools/init.hpp"
#include "models/inverse_model.hpp"
#include "methods/methods.hpp"
#include "config.h"
#include "mcmc/phmc.hpp"
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

void study_distribution(InverseProblem *rand_var);
int main(int argc,char *argv[])
{
  int no_samples=1000000;
  //int no_samples=100; 
  bool study_dist = false;
  Real prior_variance = 1; 
  //Real noise_variance = (0.01)*0.01; 
  Real noise_variance = (0.03)*0.03;
  //Real noise_variance = 0.01;

  //std::string forward_op_type = "det";
  //std::string forward_op_type = "rts_rk";
  std::string forward_op_type = "rts_rk_skewed";
  
  bool just_dist = true;

  //for ( int no_integration_steps=50;no_integration_steps<=100;no_integration_steps+=2)
  for ( int no_integration_steps=5;no_integration_steps<=50;no_integration_steps+=2)
  {

    if ( just_dist ) no_integration_steps = 20;
    printf("no_integration_steps: %d\n",no_integration_steps);
    InverseProblem* rand_var =TestInverse::get_instance(prior_variance,noise_variance,no_integration_steps,forward_op_type);
   
    //InverseProblem* rand_var =FitzHughInverse::get_instance(prior_variance,noise_variance,no_integration_steps,forward_op_type);
    rand_var->set_include_prior(true);
    //InverseProblem* rand_var =FitzHughInverse::get_instance();
    rand_var->init_x_obs();

    if ( study_dist) study_distribution(rand_var);

    /*
    Real* x0=new Real[3];
    x0[0]=0;
    x0[1]=0;
    x0[2]=0;
    */
    Real* x0=new Real[1];
    x0[0]=2;
  
    //Real gamma=0.5;
    Real gamma=0.05;
    PHMC_Sampler * mcmc_sampler = new PHMC_Sampler(rand_var,gamma,no_samples,x0);
    //MCMC_Sampler * mcmc_sampler = new MCMC_Sampler(rand_var,gamma,no_samples,x0);
   
    mcmc_sampler->sample();
    Real **arr= mcmc_sampler->get_sampled_array();  

    std::string model_name=rand_var->get_true_model()->get_name();
    std::string model_inv_name = model_name + "_inverse";
    boost::filesystem::path dir(model_inv_name.c_str());
    if ( !boost::filesystem::is_directory (dir) )
      boost::filesystem::create_directory(dir);
    std::stringstream ss;
    ss<<model_inv_name+"/"+forward_op_type+"_dist_";
    ss<<no_integration_steps;
    ss<<".txt";

    std::ofstream fout;
    if ( just_dist )
    {
      fout = std::ofstream("distribution.txt");
    }
    else
    {
      fout = std::ofstream(ss.str());
    }
    //FILE *fout = fopen("distribution_.txt","w");
    //fprintf(fout,"x0,x1,x2\n");
    //fprintf(fout,"x0\n");
    for ( int i=0;i<rand_var->get_prior_distribution().size();i++)
    {
      if ( i>0) fout << ",";
      fout <<  "x" << i ;
    }

    fout<< std::endl;
    for ( int i=no_samples/10;i<no_samples;i++)
    {
      //fprintf(fout,"%lf,%lf,%lf\n",exp(arr[i][0]),exp(arr[i][1]),exp(arr[i][2]));
      //fprintf(fout,"%lf,%lf,%lf\n",(arr[i][0]),(arr[i][1]),(arr[i][2]));
      //fprintf(fout,"%lf\n",arr[i][0]);
      
      for ( int j=0;j<rand_var->get_prior_distribution().size();j++)
      {
        if ( j>0) fout << ",";
        fout << arr[i][j];
      }
      fout << std::endl;
    }
    printf("\n\n");
    delete[] x0;
    delete mcmc_sampler;
    delete rand_var;
    if ( just_dist )
      break;
  }
  return 0;
}


void study_distribution(InverseProblem *rand_var)
{
  Real *x=new Real[1];
  FILE *f_distribution=fopen("distribution.txt","w");
  fprintf(f_distribution,"x1,x2\n");
  int bound = 400;
  double min_pot=1000000000.0;
  for ( int i=0;i<bound;i++)
  {
    x[0]=i/300.0;
    min_pot=rmin(min_pot,rand_var->potential(x));
  }
  for ( int i=0;i<bound;i++)
  {
    x[0]=i/300.0;
    fprintf(f_distribution,"%.16lf,%.16lf\n",x[0],exp(-rand_var->potential(x)+min_pot));
  }
  delete[] x;
  fclose(f_distribution);
}
