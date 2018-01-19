#include "mala.hpp"

void MALA_Sampler::sample()
{
  int no_samples=10;
  no_acc=0;
  no_rejects=0;
  Real prev_p=0;

  for ( int i=0;i<n;i++)
    sampled_array[0][i]=init_val[i];
  Real *w=new Real[n];
  Real *y=new Real[n];
  Real *x=new Real[n];
  for ( int i=0;i<n;i++)
    x[i]=sampled_array[0][i];
  for ( int k=0;k<no_samples;k++)
  {
    prev_p+=rand_var->potential(x);
  }
  prev_p-=log(no_samples);
  for ( int iter=1;iter<no_iter;iter++)
  {
    printf("%d\n",iter);
    Real px=0;
    px=prev_p;
    for ( int i=0;i<n;i++) 
    {
      w[i]=(*gamma_normal[0])();
      //printf("%lf ",w[i]);
    }
    //printf("\n");
    for ( int i=0;i<n;i++) y[i]=x[i]+w[i];
    Real py=rand_var->potential(y);
    Real alpha= rmin(1,exp(-py+px));
    double u = U();
    if ( u < alpha) // accept
    {
      no_acc++;
      for ( int i=0;i<n;i++)
        sampled_array[iter][i]=x[i];
      for ( int i=0;i<n;i++)
        x[i]=y[i];
      prev_p=py;
    }
    else
    {
      no_rejects++;
      for ( int i=0;i<n;i++)
        sampled_array[iter][i]=x[i];
      prev_p=px;
    }
  }
  delete []x;
  delete []y;
  delete []w;
}

