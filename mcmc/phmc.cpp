#include "phmc.hpp"

void PHMC_Sampler::sample()
{
  int no_samples=1;
  no_acc=0;
  no_rejects=0;
  Real prev_p=0;
  Real prev_p_normalization;
  
  for ( int i=0;i<n;i++)
    sampled_array[0][i]=init_val[i];
  Real *w=new Real[n];
  Real *y=new Real[n];
  Real *x=new Real[n];
  for ( int i=0;i<n;i++)
    x[i]=sampled_array[0][i];

  prev_p_normalization = rand_var->potential(x);
  //prev_p=1;

  //prev_p=exp(-rand_var->potential(x));
  prev_p = rand_var->potential(x);

  for ( int k=1;k<no_samples;k++)
    prev_p+=rand_var->potential(x);
    //prev_p+=exp(-rand_var->potential(x)+prev_p_normalization);
    //prev_p+=exp(-rand_var->potential(x));
  prev_p/=no_samples;
  //prev_p-=log(no_samples);
  long long last=0;
  for ( long long  iter=1;iter<no_iter;iter++)
  {
    if (iter*100/no_iter > last  )
    {
      last = iter*100/no_iter;
      printf("Percentage of accepts: %lf\n",100.0*no_acc/(no_rejects+no_acc));
      printf("%d\%\n",last);
    }
    Real px=prev_p;
    //Real px_normalization = prev_p_normalization;

    for ( int i=0;i<n;i++) 
    {
      w[i]=(*gamma_normal[0])();
      //printf("%lf ",w[i]);
    }
    //printf("\n");
    for ( int i=0;i<n;i++) y[i]=x[i]+w[i];
    //Real py_normalization = rand_var->potential(y);

    Real py;
    //py=1;

    //py=py_normalization;
    //py = exp(-rand_var->potential(y));
    py = rand_var->potential(y); 

    for ( int k=1;k<no_samples;k++)
      py+=rand_var->potential(y);
      //py += exp(-rand_var->potential(y)+py_normalization);
      //py += exp(-rand_var->potential(y));
    py/=no_samples;
    //py -= log(no_samples);
    //Real py=rand_var->potential(y);
    //printf("%lf\n",py);
    Real alpha= rmin(1,exp(-py+px));
    //Real alpha= rmin(1,py/px);
    //Real alpha= rmin(1,py/px*exp(-py_normalization+px_normalization));
    //printf("%lf %.16lf\n",y[0],alpha);
    //printf("%lf %lf %lf %lf %.16lf\n",alpha,px_normalization,py_normalization,px,py);
    double u = U();
    if ( u < alpha) // accept
    {
      no_acc++;
      for ( int i=0;i<n;i++)
        sampled_array[iter][i]=x[i];
      for ( int i=0;i<n;i++)
        x[i]=y[i];
      prev_p=py;
      //prev_p_normalization = py_normalization;
    }
    else
    {
      no_rejects++;
      for ( int i=0;i<n;i++)
        sampled_array[iter][i]=x[i];
      prev_p=px;
    }
  }
  printf("Percentage of accepts: %lf\n",100.0*no_acc/(no_rejects+no_acc));
  delete []x;
  delete []y;
  delete []w;
}

PHMC_Sampler::~PHMC_Sampler()
{
  for ( int i=0;i<no_iter;i++)
    delete []sampled_array[i];
  delete []sampled_array;
  for ( int i=0;i<n;i++)
  {
    delete gamma_normal[i];
  }
  delete [] gamma_normal;
}
