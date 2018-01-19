#include "methods.hpp"
#include "stabilized_methods.hpp"
#include <algorithm>
/********************************************************************************
 * RKC method for Stiff methods
 *******************************************************************************/

int RKC::get_conv_order(){ return 1;}
void RKC::operator()( Real *x0, Real *x1, Real h )
{
  f_evals=0;
  // initialize s
  //calc_number_of_stages(x0,h);
  s=100;
  prepare_chebychev_info();
  for ( int i=0;i<n;i++)
    K[0][i]=x0[i];
  eq->eval(K[0],F[0]);

  f_evals++;
  for ( int i=0;i<n;i++)
  {
    K[1][i]=K[0][i]+h*w1/w0*F[0][i];
    
    // no dumping 
    //K[1][i]=K[0][i]+h/s/s*F[0][i];
  }
  for ( int j=2;j<=s;j++)
  {
    eq->eval(K[j-1],F[j-1]);
    f_evals++;
    for ( int i=0;i<n;i++)
    {
      K[j][i]=2*h*w1*T[j-1]/T[j]*F[j-1][i]+2*w0*T[j-1]/T[j]*K[j-1][i]-T[j-2]/T[j]*K[j-2][i];
      // no dumping
      //K[j][i]=2*h/s/s*F[j-1][i]+2*K[j-1][i]-K[j-2][i];
    }
  }
  // final stage..
  for ( int i=0;i<n;i++)
    x1[i]=K[s][i];
  if ( do_error_estimation)
  {
    error =0 ;
    for ( int i=0;i<n;i++)
    {
      error += (x1[i]-x0[i])*(x1[i]-x0[i]);
    }
    error = sqrt(error);
  }
}

int RKC::get_f_evals()
{
  return f_evals;
}

void RKC::prepare_chebychev_info()
{
  w0= 1+nu/s/s;
  T[0]=1;
  T[1]=w0;
  U[0]=1;
  U[1]=2*w0;
  for ( int i=2;i<=s;i++)
  {
    T[i]=2*w0*T[i-1]-T[i-2];
    U[i]=2*w0*U[i-1]-U[i-2];
  }
  w1=T[s]/(s*U[s-1]);
}
void RKC::calc_number_of_stages(Real *x0, Real h)
{
  static const Real safety = 1.05;
  Real max_eig=eq->largest_eigenval(x0);
  max_eig*=safety;
  s=maxs;
  Real r_tmp = (sqrt((1.5+h*max_eig)/0.65));
#ifdef BOOST_FLOAT
  s = (int) r_tmp.convert_to<double>();
#else
  s = (int) r_tmp;
#endif
  s+=1;
  if ( s < 3 ) s= 3; // at least 3 stages
  if ( s > maxs ) s= maxs;
  //printf("%g %d\n",max_eig,s);
  // find first bigger or equal s for which we have info about coefficients..
}
Real** ROCK2::mu;
Real** ROCK2::nu;
Real** ROCK2::kappa;
int ROCK2::get_conv_order(){ return 2;}

Real ROCK2::precomputed_alpha[maxs+1];
Real ROCK2::precomputed_sigma[maxs+1];
Real ROCK2::precomputed_tau[maxs+1];
bool ROCK2::static_coefs_constructed = false;

void ROCK2::operator() ( Real *x0, Real *x1, Real h )
{
  f_evals=0;
  //calc_number_of_stages(x0,h);
  //printf ("s is %d\n",s);
  s=8;
  Real alpha = precomputed_alpha[s];
  Real sigma = precomputed_sigma[s];
  Real tau   = precomputed_tau[s];
  //printf("%lf %lf %lf\n",alpha,tau,sigma);
  for ( int i=0;i<n;i++)
    K[0][i]=x0[i];
  eq->eval(K[0],F[0]);
  f_evals++;
  for ( int i=0;i<n;i++)
    K[1][i]=K[0][i]+alpha*mu[s][1]*h*F[0][i];
  for ( int j=2;j<=s-2;j++)
  {
    eq->eval(K[j-1],F[j-1]); 
    f_evals++;
    for ( int i=0;i<n;i++)
      K[j][i]=alpha*mu[s][j]*h*F[j-1][i]-nu[s][j]*K[j-1][i]-kappa[s][j]*K[j-2][i];
  }
  //** final stages **//
  eq->eval(K[s-2],F[s-2]);
  f_evals++;

  //printf(" S is %d, %lf\n",s,alpha);
  Real sigma_a=(1-alpha)/2+alpha*sigma;
  Real tau_a = (alpha-1)*(alpha-1)/2.0+2*alpha*(1-alpha)*sigma+alpha*alpha*tau;
  for ( int i=0;i<n;i++)
    K[s-1][i]=K[s-2][i]+2*tau_a*h*F[s-2][i];
  eq->eval(K[s-1],F[s-1]);
  f_evals++;
  for ( int i=0;i<n;i++)
    x1[i]=K[s-2][i]+(2*sigma_a-0.5)*h*F[s-2][i]+0.5*h*F[s-1][i];

  if ( do_error_estimation )
  {
    // construct K_{s-1}^{**} 
    for ( int i=0;i<n;i++)
      K_error[i] = K[s-2][i]+2*tau_a*h*F[s-2][i];
    eq->eval(K_error,F_error);
    //f_evals++;
    error=0;
    for ( int i=0;i<n;i++)
    {
      Real v=(F[s-2][i]-F_error[i]); 
      error+= v*v;
    }
    error = sqrt(error);
    error *= 0.5*h*(1-sigma_a*sigma_a/tau_a);
  }
}

void ROCK2::calc_number_of_stages(Real *x0,Real h)
{
  static const Real safety = 1.05;
  Real max_eig=eq->largest_eigenval(x0)*safety;
  //printf("max eigenval is %lf\n",max_eig);
  Real r_tmp = (sqrt((1.5+h*max_eig)/0.811));
#ifdef BOOST_FLOAT
  s = (int) r_tmp.convert_to<double>();
#else
  s = (int) r_tmp;
#endif
  s-=2;
  if ( s >= maxs-2 )
  {
    s= maxs;
  }
  auto it = std::upper_bound(ms,ms+46,s);
  if ( it!= ms+46 )
    s=*it;
  else 
    s= maxs-2;
  s+=2;

  if ( s < 3 ) s= 3; // at least 3 stages
}
int ROCK2::get_f_evals()
{
  return f_evals;
}
