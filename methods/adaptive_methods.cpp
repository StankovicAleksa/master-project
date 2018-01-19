#include "adaptive_methods.hpp"
#include "../config.h" 
#include <cmath>

void ClassicalAdaptiveScheme::get_final(Real *x1)
{
  Real h=h0;
  Real prev_err_factor;
  Real prev_h;
  Real new_h;
  bool first_iteration= true;
  print_pos(); 
  double eps=1e-8;
  while (cur_t < model.get_tf() )
  {
    if ( cur_t+h > model.get_tf() )
    {
      h=model.get_tf()-cur_t;
    }
    this->base_integrator(cur_x,next_x,h);
    f_evals+=this->base_integrator.get_f_evals();
    Real sgma = sigma(cur_x,next_x,h);
    Real error = base_integrator.get_error();
    Real err_factor = error/sgma;
    if ( error < sgma ) // accept time-step
    {
      for ( int i=0;i<n;i++)
        x_tmp[i]=cur_x[i];
      print_pos();
      //printf("%lf %lf\n",cur_t,h);
      for ( int i=0;i<n;i++)
        cur_x[i]=next_x[i];
      cur_t+=h;
      // new time step is 
      if (!first_iteration )  
        new_h = fac*h*sqrt(1/err_factor)*h/prev_h*sqrt(prev_err_factor/err_factor);
      else
        new_h=h;
      sum_h +=h;
      n_h ++;
    }
    else 
    {
      //printf("Rejected %g %g\n",cur_t,error,h);
      new_h = h*sqrt(1/err_factor);
    }
    if ( new_h>0.1)
    {
      //printf("New h proposal: %lf\n",new_h);
      new_h=0.1; 
    }
    prev_h=h;
    h=new_h;
    first_iteration = false;
    prev_err_factor= err_factor;
    // we move iff the error is smaller than tol...
  }
  for ( int i=0;i<n;i++)
    x_tmp[i]=next_x[i];
  print_pos();
  for ( int i=0;i<n;i++)
    x1[i]=next_x[i];
}

Real ClassicalAdaptiveScheme::sigma(Real *x0,Real *x1,Real h)
{
  Real n0,n1; n0 = n1 = 0;
  for ( int i=0;i<n;i++) n0+=x0[i]*x0[i],n1+=x1[i]*x1[i];
  n0=sqrt(n0), n1=sqrt(n1);
  return (rmax(n0,n1))*tol+tol;
}



void GiacomoRandomPathsAdaptiveScheme::get_final(Real *x1)
{
  UniformDistribution  ud(-1.0,1.0);
  Real h=h0;
  Real prev_err_factor;
  Real prev_h;
  Real new_h;
  Real h_p;
  bool first_iteration= true;
  double emb_error=0;
  double err_factor;
  double rand_error=0;
  double abs_error =0 ;
  double E[n];
  double variance[n][n];
  double tmp[n];
  double prev_error=0;
  //double C=tol;
  double lya_dist = tol;
  Real max_eig;//=eq->largest_eigenval(x0)*safety;
  std::string f_err_path = model.get_name()+"/err.txt";
  FILE *f_errors=fopen(f_err_path.c_str(),"w");
  fprintf(f_errors,"t,var,emb_err,sqrt(var)\n");
  /** generate starting points **/
  for ( int i=0;i<no_paths;i++)
  {
    for ( int j=0;j<n;j++)
      cur_points[i][j]=cur_x[j];
    //h_path_sum[i]=0;
  }
  print_pos(cur_x); 
  while (cur_t < model.get_tf() )
  {
    if ( cur_t+h > model.get_tf() )
    {
      h=model.get_tf()-cur_t;
    }
    // random integration step, Assyr-Garegnani style 
    h_p=pow(h,1*(base_integrator.get_conv_order()+0.5));
    emb_error=0;
    max_eig=0;
    for ( int i=0;i<no_paths;i++)
    {
      if ( i==0 ) 
        h_path[i]=h;
      else
        h_path [i] = h+ ud()*h_p;
      this->base_integrator(cur_points[i],next_points[i],h_path[i]);
      f_evals+=this->base_integrator.get_f_evals();
      if ( i ==0 )
      {
        emb_error = base_integrator.get_error();
      }

      double t_eig=base_integrator.eq->directional_eigenval(cur_points[i],next_points[i]);
      //printf("%lf\n",t_eig);
      //double t_eig=base_integrator.eq->largest_eigenval(cur_points[i])*safety;
      if ( t_eig > max_eig) max_eig = t_eig;
      else 
        emb_error = rmax(emb_error,base_integrator.get_error());
    }
    for ( int i=0;i<n;i++)
      E[i]=0;
    for ( int i=0;i<no_paths;i++)
      for ( int j=0;j<n;j++)
        E[j]+=cur_points[i][j];
    for ( int i=0;i<n;i++)
      E[i]/=no_paths;
    for ( int i=0;i<no_paths;i++)
    {
      double tmp =0;
      for ( int j=0;j<n;j++)
        tmp+=(E[j]-cur_points[i][j])*(E[j]-cur_points[i][j]);
      tmp=sqrt(tmp);
      if ( tmp > abs_error )
        abs_error=tmp;
    }
    for ( int i=0;i<n;i++)
      for ( int j=0;j<n;j++)
        variance[i][j]=0;

    for ( int i=0;i<no_paths;i++)
    {
      for ( int j=0;j<n;j++)
      {
        tmp[j]=cur_points[i][j]-E[j];
      }
      // add up to variance 
      for ( int j=0;j<n;j++)
        for ( int k=0;k<n;k++)
        {
          variance[j][k]+=tmp[j]*tmp[k];
        }
    }
    for ( int i=0;i<n;i++)
      for ( int j=0;j<n;j++)
        variance[i][j]/=no_paths-1;
    rand_error=0; // variance
    for ( int i=0;i<n;i++)
    {
      rand_error += variance[i][i];
    }
    rand_error=sqrt(rand_error);
    //rand_error=sqrt(rand_error);
    double fact = 1;
    Real sigma_min=0;
    if (max_eig <=  0 ) max_eig=0.02;
    double err;
    //sigma_min = C*exp(lambda*(cur_t+h));
    double new_lya_dist = lya_dist*(1+max_eig*h);
    sigma_min=new_lya_dist;
    if  (true)
    {
      err=emb_error;
      //sigma_min = sigma(cur_points[0],next_points[0],h_path[0]);
    }
    else
    {
      err=rand_error;
      //err=abs_error;
      //fact = (cur_t-model.get_t0())/(model.get_tf()-model.get_t0());
      //sigma_min = prev_error;
    }
    fflush(stdout);
    
    sigma_min = sigma(cur_points[0],next_points[0],h_path[0]);
    for ( int i=0;i<no_paths;i++)
      sigma_min = rmin(sigma_min,sigma(cur_points[i],next_points[i],h_path[i]));
    err=emb_error;
    //printf("%g %g %.4g %.4g %g %g\n",cur_t,h,emb_error,abs_error,sigma_min,new_lya_dist);
    fprintf(f_errors,"%g,%g,%g,%g\n",cur_t,rand_error,emb_error,sqrt(rand_error));
    //printf("%g,%g,%g,%g\n",cur_t,rand_error,emb_error,sqrt(rand_error));
    //for ( int i=1;i<no_paths;i++)
      //sigma_min = rmin(sigma_min,sigma(cur_points[i],next_points[i],h_path[i])*fact);
    err_factor = err/sigma_min;
    if ( err < sigma_min  ) // accept time-step
    {
      if (!first_iteration )  
      {
        new_h = fac*h*sqrt(1/err_factor)*h/prev_h*sqrt(prev_err_factor/err_factor);
      }
      else
        new_h=h;
      print_pos(cur_points[0]);
      //printf("%lf %lf\n",cur_t,h);
      for ( int i=0;i<n;i++)
        for ( int j=0;j<no_paths;j++)
          cur_points[j][i]=next_points[j][i];
      for ( int i=0;i<n;i++)
        cur_x[i]=cur_points[0][i];
      cur_t+=h;
      sum_h +=h;
      n_h++;
      // new time step is 
      first_iteration = false;
      lya_dist = new_lya_dist;
    }
    else 
    {
      //printf("Rejected %g\n",err_factor);
      new_h = h*sqrt(1/err_factor);
    }
    /* fixed time-stepping */
    //new_h=h;
    //print_pos(); 


    if ( new_h > 0.1) new_h=0.1; 
    prev_h=h;
    h=new_h;
    prev_err_factor= err_factor;
    prev_error = rand_error;


    // we move iff the error is smaller than tol...
  }

  print_pos(cur_points[0]);
  for ( int i=0;i<n;i++)
  {
    x1[i]=0;
    for ( int j=0;j<no_paths;j++)
      x1[i]+=cur_points[j][i];
  }
  for ( int i=0;i<n;i++)
      x1[i]/=no_paths;
  fclose(f_errors);
  close_f();
}

Real GiacomoRandomPathsAdaptiveScheme::sigma(Real *x0,Real *x1,Real h)
{
  Real n0,n1; n0 = n1 = 0;
  for ( int i=0;i<n;i++) n0+=x0[i]*x0[i],n1+=x1[i]*x1[i];
  n0=sqrt(n0), n1=sqrt(n1);
  n1=0.0; // TODO -> change this
  return rmax(n0,n1)*tol+tol;
}


void AleksaRandomPathsAdaptiveScheme::get_final(Real *x1)
{
  UniformDistribution  ud(-1.0,1.0);
  Real h=h0;
  Real prev_err_factor;
  Real prev_h;
  Real new_h;
  Real h_p;
  bool first_iteration= true;
  double emb_error=0;
  double err_factor;
  double rand_error=0;
  double abs_error =0 ;
  double E[n];
  double variance[n][n];
  double tmp[n];
  double prev_error=0;
  //double C=tol;
  Real max_eig;//=eq->largest_eigenval(x0)*safety;
  std::string f_err_path = model.get_name()+"/err.txt";
  FILE *f_errors=fopen(f_err_path.c_str(),"w");
  fprintf(f_errors,"t,var,emb_err,sqrt(var)\n");
  /** generate starting points **/
  for ( int i=0;i<no_paths;i++)
  {
    for ( int j=0;j<n;j++)
      cur_points[i][j]=cur_x[j];
    //h_path_sum[i]=0;
  }
  print_pos(); 
  while (cur_t < model.get_tf() )
  {
    double min_h=1e-8;
    if ( h< min_h) h=min_h;
    if ( cur_t+h > model.get_tf() )
    {
      h=model.get_tf()-cur_t;
    }

    /******8
     * JUST SOME UNIMPORTANT CALCULATIONS
     *
     **/
    // random integration step, Assyr-Garegnani style 
    h_p=pow(h,1*(base_integrator.get_conv_order()+0.5));
    emb_error=0;
    max_eig=0;
    for ( int i=0;i<no_paths;i++)
    {
      if ( i==0 ) 
        h_path[i]=h;
      else
        h_path [i] = h+ ud()*h_p;
      this->base_integrator(cur_points[i],next_points[i],h_path[i]);
      f_evals+=this->base_integrator.get_f_evals();
      if ( i ==0 )
      {
        emb_error = base_integrator.get_error();
      }

      double t_eig=base_integrator.eq->directional_eigenval(cur_points[i],next_points[i]);
      //printf("%lf\n",t_eig);
      //double t_eig=base_integrator.eq->largest_eigenval(cur_points[i])*safety;
      if ( t_eig > max_eig) max_eig = t_eig;
      //else 
      //emb_error = rmax(emb_error,base_integrator.get_error());
    }
    for ( int i=0;i<n;i++)
      E[i]=0;
    for ( int i=0;i<no_paths;i++)
      for ( int j=0;j<n;j++)
        E[j]+=next_points[i][j];
    for ( int i=0;i<n;i++)
      E[i]/=no_paths;
    for ( int i=0;i<no_paths;i++)
    {
      double tmp =0;
      for ( int j=0;j<n;j++)
        tmp+=(E[j]-next_points[i][j])*(E[j]-next_points[i][j]);
      tmp=sqrt(tmp);
      if ( tmp > abs_error )
        abs_error=tmp;
    }
    for ( int i=0;i<n;i++)
      for ( int j=0;j<n;j++)
        variance[i][j]=0;

    for ( int i=0;i<no_paths;i++)
    {
      for ( int j=0;j<n;j++)
      {
        tmp[j]=next_points[i][j]-E[j];
      }
      // add up to variance 
      for ( int j=0;j<n;j++)
        for ( int k=0;k<n;k++)
        {
          variance[j][k]+=tmp[j]*tmp[k];
        }
    }
    for ( int i=0;i<n;i++)
      for ( int j=0;j<n;j++)
        variance[i][j]/=no_paths-1;
    rand_error=0; // variance
    for ( int i=0;i<n;i++)
    {
      rand_error += variance[i][i];
    }
    rand_error=sqrt(rand_error);
    //rand_error/=h;
    //printf("%lf\n",h);
    /******8
     * END OF JUST SOME UNIMPORTANT CALCULATIONS
     *
     **/
    
    Real sigma_min=0;
    double err;
    sigma_min = sigma(cur_points[0],next_points[0],h_path[0]);
    err=rand_error-prev_error;
    
    //err=sqrt(rand_error)-sqrt(prev_error);

    if ( err < 0 ) err = -err;
    //err=sqrt(err); 

    double min_err=1e-10*tol;
    //printf("%g\n",err);
    if ( err < min_err) 
    {
      //err=min_err;
      //printf("MIN_ERR %g\n",rand_error-prev_error);
    }
    //err=emb_error; 
    fprintf(f_errors,"%g,%g,%g,%g\n",cur_t,rand_error,emb_error,sqrt(rand_error));
    err_factor = err/sigma_min;
    printf("%lf\n",cur_t);
    //printf("%g %g %g\n" ,h,prev_error,rand_error);
    //if ( err < sigma_min ) // accept time-step
    if ( true )
    {
      if ( !first_iteration)
      {
        //printf(" Adapting %g!\n",fac*sqrt(1/err_factor)*h/prev_h*sqrt(prev_err_factor/err_factor));
        //new_h = fac*h*sqrt(1/err_factor)*h/prev_h*sqrt(prev_err_factor/err_factor);
        new_h = h*fac*sqrt(sqrt(1/err_factor))*h/prev_h*sqrt(sqrt(prev_err_factor/err_factor));
        //new_h=fact*h;
      }
      else
        new_h=h;

      //printf("%lf %lf\n",cur_t,h);
      for ( int i=0;i<n;i++)
        for ( int j=0;j<no_paths;j++)
          cur_points[j][i]=next_points[j][i];
      
      cur_t+=h;
      for ( int i=0;i<n;i++)
      {
        cur_x[i]=0;
        for ( int j=0;j<no_paths;j++)
          cur_x[i]+=cur_points[0][i];
        cur_x[i]/=no_paths;
      }
      print_pos(cur_points[0]);
      
      
      sum_h +=h;
      n_h++;
      // new time step is 
      first_iteration = false;
    }
    else 
    {
      new_h = h*sqrt(1/err_factor);
      //new_h = h/fact;
    }
    /* fixed time-stepping */
    //new_h=h;
    if ( new_h > 0.1) 
    {
      new_h=0.1; 
    }
    prev_h=h;
    h=rmax(0.8*prev_h,rmin(1.2*prev_h,new_h));
    //h=new_h ;
    //printf("%lf\n",h);
    h=prev_h;
    prev_err_factor= err_factor;
    prev_error = rand_error;

    // we move iff the error is smaller than tol...
  }
  for ( int i=0;i<n;i++)
  {
    x1[i]=0;
    for ( int j=0;j<no_paths;j++)
      x1[i]+=cur_points[j][i];
    x1[i]/=no_paths;
  }
  print_pos(x1);
  fclose(f_errors);
  close_f();
}

Real AleksaRandomPathsAdaptiveScheme::sigma(Real *x0,Real *x1,Real h)
{
  Real n0,n1; n0 = n1 = 0;
  for ( int i=0;i<n;i++) n0+=x0[i]*x0[i],n1+=x1[i]*x1[i];
  n0=sqrt(n0), n1=sqrt(n1);
  return rmax(n0,n1)*tol+tol;
}
