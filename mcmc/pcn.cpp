#include "pcn.hpp"

void PCN_Sampler::sample()
{
  rand_var->set_include_prior(false);
  no_acc=0;
  no_rejects=0;
  for ( int i=0;i<n;i++)
    sampled_array[0][i]=init_val[i];

  Real *w=new Real[n];
  Real *y=new Real[n];
  Real *x=new Real[n];
  Real *tx= new Real[n];

  int model_dim = rand_var->get_true_model()->get_dim();
  UniformDistribution un_dist(0,5);
 


  FILE *pcn_f_out = fopen("pcn_res.txt","w");
  for ( int i=0;i<model_dim;i++)
  {
    if ( i!=  0 ) fprintf(pcn_f_out,",");
    fprintf(pcn_f_out,"x%d",i);
  }

  fprintf(pcn_f_out,"\n");
  for ( int i=0;i<n;i++)
    x[i]=sampled_array[0][i];

  for ( int iter=1;iter<no_iter;iter++)
  {
    std::vector<Real> params ;
    for ( int i=0;i<n;i++)
      params.push_back(x[i]);
    Model * m = rand_var->get_model(params);
    for ( int i=0;i<m->get_dim();i++)
    {
      if ( i!=  0 ) fprintf(pcn_f_out,",");
      fprintf(pcn_f_out,"%lf",m->get_x0()[i]);
    }
    fprintf(pcn_f_out,"\n");
    Real px=rand_var->potential(x);
    Real ud_sample = un_dist();
    for ( int i=0;i<n;i++)   
    {
      //w[i]=ud_sample*(*gamma_normal[0])();
      w[i]=(*gamma_normal[0])();
    }
    //printf("\n");
    if ( iter%100==0)
    {
      printf("%d\n",iter);
      printf("No accepts: %d, no rejects: %d, acceptance ratio: %lf \%\n",no_acc,no_rejects,((double)no_acc)/(no_acc+no_rejects));
    }
    // Adaptive beta
    //beta = 1.0/sqrt(iter);
    //beta = 1.0/(iter);
    //if ( beta < 0.02 ) beta = 0.02;
    beta = 0.2;
    for ( int i=0;i<n;i++) 
    {
      y[i]=sqrt(1-beta*beta)*x[i]+beta*w[i];
    }
      //y[i]=x[i]+w[i];
    
    /*
    printf("\n");
    printf("u: ");
    for ( int i=0;i<n;i++)
    {
      if (2*i == n )
        printf("  v: ");
      else if ( i!=0 ) printf(",");
      printf("%lf",x[i]);
    }
    printf("\n\n");
    printf("u: ");
    for ( int i=0;i<n;i++)
    {
      if (2*i == n )
        printf("  v: ");
      else if ( i!=0 ) printf(",");
      printf("%lf",y[i]);
    }
    printf("\n\n");
    */
    Real py=rand_var->potential(y);

    for ( int i=0;i<n;i++)
      tx[i]=0;
    //Real pt=rand_var->potential(tx);
    //printf("Probability of  x: %lf .Probability of y: %lf. Probability of t: %lf\n\n",px,py,pt);
    Real alpha= rmin(1,exp(+px-py));
    double u = U();
    if ( u < alpha) // accept
    {
      no_acc++;
      for ( int i=0;i<n;i++)
        sampled_array[iter][i]=x[i];
      for ( int i=0;i<n;i++)
        x[i]=y[i];

      // output 
      
    }
    else
    {
      //printf("Rejected!");
      no_rejects++;
      for ( int i=0;i<n;i++)
        sampled_array[iter][i]=x[i];
    }
  }
  printf("\n");
  printf("u: ");
  for ( int i=0;i<n;i++)
  {
    Real eigen_val = 1.0;
    if (i >0 && i < 5  )
    {
       eigen_val /= (2*M_PI*i) * (2*M_PI*i);
       continue;
    }
    else if ( i >= 5 && i < 10)
       eigen_val /= (2*M_PI*(i-4)) * (2*M_PI*(i-4));

    if ( 2*i == n) break;
    if (2*i == n )
      printf("  v: ");
    else if ( i!=0 ) printf(",");
    //printf("%lf",x[i]);
    printf("%lf",x[i]*eigen_val);
  }
  printf("\n\n");
  delete []x;

  delete []y;
  delete []w;
}

