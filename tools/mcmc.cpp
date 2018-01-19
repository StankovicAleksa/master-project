#include "mcmc.hpp"

void MCMC_Sampler::sample()
{
  no_acc=0;
  no_rejects=0;
  for ( int i=0;i<n;i++)
    sampled_array[0][i]=init_val[i];
  Real *w=new Real[n];
  Real *y=new Real[n];
  Real *x=new Real[n];
  Real *tx= new Real[n];

  for ( int i=0;i<n;i++)
    x[i]=sampled_array[0][i];
  for ( int iter=1;iter<no_iter;iter++)
  {
    if ( iter %100 == 0)
    printf("%d\n",iter);
    Real px=rand_var->potential(x);
    for ( int i=0;i<n;i++) 
    {
      w[i]=(*gamma_normal[0])();
      //printf("%lf ",w[i]);
    }
    //printf("\n");


    for ( int i=0;i<n;i++) y[i]=x[i]+w[i];
    Real py=rand_var->potential(y);
    
    
    for ( int i=0;i<n;i++)
      tx[i]=0;
    Real pt=rand_var->potential(tx);

    //Real alpha= rmin(1,py/px);
    Real alpha= rmin(1,exp(-py+px));
    double u = U();
    if ( u < alpha) // accept
    {
      no_acc++;
      for ( int i=0;i<n;i++)
        sampled_array[iter][i]=x[i];
      for ( int i=0;i<n;i++)
        x[i]=y[i];
    }
    else
    {
      no_rejects++;
      for ( int i=0;i<n;i++)
        sampled_array[iter][i]=x[i];
    }
  }
  delete []x;
  delete []y;
  delete []w;
  delete []tx;
}
