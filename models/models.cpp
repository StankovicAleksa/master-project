#include "models.hpp"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <limits>
#include "../config.h"


Model::Model (const Model &model)
{
  x0=new Real[model.dim];
  t0=model.t0;
  dim=model.dim;
  tf=model.tf;
}


Real eps=std::numeric_limits<Real>::epsilon();
Real sqrt_eps = sqrt(eps);
// not completely the same as in the paper!


/**** Power method *****/
Real Equation::largest_eigenval(const Real *y0)
{
  if ( eigvec_prev == 0  )
  {
    eigvec_prev = new Real[dim];
    eigvec_tmp1 = new Real[dim];
    eigvec_tmp2 = new Real[dim];
    eigvec_tmp3 = new Real[dim];
    // initialize eigval_prev as random vector
    for ( int i=0;i<dim;i++)
    {
      eigvec_prev[i] = (-RAND_MAX/2+rand()) /(Real) RAND_MAX;
      //printf("%lf ",eigvec_prev[i]); 
    }
    //printf("\n");
  }
  Real *delta=eigvec_prev;
  Real *z= eigvec_tmp1;
  Real *y=eigvec_tmp2;
  Real *f_y0 = eigvec_tmp3;
  Real norm_y0=0;
  Real norm_delta=0;
  static const int max_iter = 10000;
  static const Real tol = 1e-3;
  eval(y0,f_y0);

  // calculate norms
  for ( int i=0;i<dim;i++)
  {
    norm_y0+=y0[i]*y0[i];
    norm_delta+=delta[i]*delta[i];
  }
  norm_y0=sqrt(norm_y0);
  norm_delta=sqrt(norm_delta);
  
  Real h;
  Real coef; 
  if ( norm_y0 != 0.0 && norm_delta !=0.0  )
  {
    h=sqrt_eps*norm_y0;
    coef = 1.0/norm_delta;
  }
  else if ( norm_y0 != 0.0 )
  {
    for ( int i=0;i<dim;i++)
      delta[i]=y0[i];
    coef = 1;
    h=sqrt_eps;
  }
  else if ( norm_delta != 0.0)
  {
    coef = 1.0/ norm_delta;
    h=sqrt_eps;
  }
  else
  {
    h=eps*sqrt(dim);
    for ( int i=0;i<dim;i++)
      delta[i]+=sqrt_eps;  
    coef = 1.0;
  }
  for ( int i=0;i<dim;i++)
    z[i]=y0[i]+h*coef*delta[i];

  Real prev_eigmax;
  Real eigmax = 0;
  Real new_norm; 
  for ( int i=0;i<max_iter;i++)
  {
    eval(z,y);
    for ( int i=0;i<dim;i++)
      y[i]-=f_y0[i];
    new_norm = 0;
    for ( int i=0;i<dim;i++)
      new_norm+=y[i]*y[i];
    new_norm= sqrt(new_norm);
    prev_eigmax = eigmax;
    eigmax = new_norm/h;
    //printf("%g\n",new_norm);
    if ( i >= 2 && fabs(eigmax-prev_eigmax) <= eigmax * tol + tol )
    {
      for ( int i=0;i<dim;i++)
        eigvec_prev[i] = z[i]-y0[i];
      return fabs(eigmax);
    }
    if ( new_norm != 0.0)
    {
      for ( int i=0;i<dim;i++)
        z[i]=y0[i]+y[i]*h/new_norm;
    }
    else break;
  }
  for ( int i=0;i<dim;i++)
  {
    if ( i!=0 ) printf(",");
    printf("%lf",y0[i]);
  }
  printf("\n");
  printf("Searching for maximum eigenvalue failed\n");
  return -1.0;
  // calculate offset
}

Real Equation::directional_eigenval(const Real *y0,const Real* y1)
{
  if ( eigvec_prev == 0  )
  {
    eigvec_prev = new Real[dim];
    eigvec_tmp1 = new Real[dim];
    eigvec_tmp2 = new Real[dim];
    eigvec_tmp3 = new Real[dim];
    // initialize eigval_prev as random vector
    for ( int i=0;i<dim;i++)
    {
      eigvec_prev[i]=y1[i]-y0[i];
      //printf("%lf ",eigvec_prev[i]); 
    }
    //printf("\n");
  }
  Real *delta=eigvec_prev;
  Real *z= eigvec_tmp1;
  Real *y=eigvec_tmp2;
  Real *f_y0 = eigvec_tmp3;
  Real norm_y0=0;
  Real norm_delta=0;
  static const Real tol = 1e-3;
  eval(y0,f_y0);

  // calculate norms
  for ( int i=0;i<dim;i++)
  {
    norm_y0+=y0[i]*y0[i];
    norm_delta+=delta[i]*delta[i];
  }
  norm_y0=sqrt(norm_y0);
  norm_delta=sqrt(norm_delta);
  
  Real h;
  Real coef; 
  if ( norm_y0 != 0.0 && norm_delta !=0.0  )
  {
    h=sqrt_eps*norm_y0;
    coef = 1.0/norm_delta;
  }
  else if ( norm_y0 != 0.0 )
  {
    for ( int i=0;i<dim;i++)
      delta[i]=y0[i];
    coef = 1;
    h=sqrt_eps;
  }
  else if ( norm_delta != 0.0)
  {
    coef = 1.0/ norm_delta;
    h=sqrt_eps;
  }
  else
  {
    h=eps*sqrt(dim);
    for ( int i=0;i<dim;i++)
      delta[i]+=sqrt_eps;  
    coef = 1.0;
  }
  for ( int i=0;i<dim;i++)
    z[i]=y0[i]+h*coef*delta[i];
  
  eval(z,y);
  for ( int i=0;i<dim;i++)
    y[i]-=f_y0[i];
  Real new_norm = 0;
  for ( int i=0;i<dim;i++)
    new_norm+=y[i]*y[i];
  new_norm= sqrt(new_norm);
  Real eigmax = new_norm/h;
  norm_delta = 0;
  for ( int i=0;i<dim;i++)
    norm_delta+=delta[i]*delta[i];
  norm_delta=sqrt(norm_delta);
  for ( int i=0;i<dim;i++)
    delta[i]/=norm_delta;
  for ( int i=0;i<dim;i++)
    y[i]/=new_norm;
  double tmp=0;
  for ( int i=0;i<dim;i++)
    tmp+=y[i]*delta[i];
  // find sign...
  return eigmax*tmp;
  printf("Searching for directional eigenvalue failed\n");
  return -1.0;
  // calculate offset
}
void TestEquation::eval(const Real *y0, Real *y1 )
{
  for ( int i=0;i<dim;i++)
  {
    y1[i]=y0[i]*lambda;
  }
}

Real TestEquation::largest_eigenval(const Real *y0)
{
  return fabs(lambda);
}

Model* get_test_equation()
{
  Real* test_eq_x0=new Real[1]{1};
  TestEquation *test_eq = new TestEquation(1,4); 
  Model* test_equation=new Model(test_eq_x0,0,1,test_eq);
  return test_equation;
}


void PeroxideOxideEquation::eval(const Real *y0,Real *y1)
{
  y1[0]=k[7]*(A0-y0[0])-k[3]*y0[0]*y0[1]*y0[3];
  y1[1]=k[8]*B0-k[1]*y0[1]*y0[2]-k[3]*y0[0]*y0[1]*y0[3];
  y1[2]=k[1]*y0[1]*y0[2]-2*k[2]*y0[2]*y0[2]+3*k[3]*y0[0]*y0[1]*y0[3]-k[4]*y0[2]+k[6]*X0;
  y1[3]=2*k[2]*y0[2]*y0[2]-k[5]*y0[3]-k[3]*y0[0]*y0[1]*y0[3];
}

Model* get_po_equation()
{
  Real *k=new Real[9]{0,0.35,250,0.035,20,5.35,1e-5,0.1,0.825};
  PeroxideOxideEquation* eq = new PeroxideOxideEquation(8,1,1,k); 
  Real *x0 = new Real[4]{6,58,0,0};
  //Model* peroxide_oxide_model=new Model(x0,0,10,eq);
  Model* peroxide_oxide_model=new Model(x0,0,200,eq);
  return peroxide_oxide_model;
}


void Lorenz::eval(const Real *y0,Real *y1)
{
  y1[0]=sigma*(y0[1]-y0[0]);
  y1[1]=y0[0]*(rho-y0[2])-y0[1];
  y1[2]=y0[0]*y0[1]-beta*y0[2];
}

Model* get_lorenz_system()
{
  Lorenz* eq = new Lorenz(10,28,8.0/3); 
  Real *x0 = new Real[3]{-10,-1,40};
  Model* lorenz_model=new Model(x0,0,40,eq);
  //Model* peroxide_oxide_model=new Model(x0,0,4,eq);
  return lorenz_model;
}


void FitzHughNagumo::eval(const Real *y0,Real *y1)
{
  y1[0]=c*(y0[0]-y0[0]*y0[0]*y0[0]/3.0+y0[1]);
  y1[1]=-1/c*(y0[0]-a+b*y0[1]);
}

Model* get_fitz()
{
  FitzHughNagumo* eq = new FitzHughNagumo(0.2,0.2,3);
  Real *x0 = new Real[2]{-1,1};
  Model* fitz_model=new Model(x0,0,20,eq);
  //Model* peroxide_oxide_model=new Model(x0,0,4,eq);
  return fitz_model;
}
Brusselator::Brusselator(int disc_nodes) :Equation(2*(disc_nodes)+1),dn(disc_nodes){
  a=1;
  b=3;
  nu=1.0/50;
  dx=1.0/(disc_nodes-1);
}  // last node is time $t$
void Brusselator::eval(const Real *y0, Real *y1)
{
  /*********************************************************************************
   *
   * y[0,1,...,disc_nodes-1] -> u_n^i
   *  y[disc_nodes,disc_nodes+2,...,2*disc_nodes-1] -> v_n^i
   *  y[2disc_nodes] -> t
   *
   *********************************************************************************/

  /*********************************************************************************
   *
   * Apply boundary condition first. Since the bondary conditions are constant functions,
   * we take here their first derivative to be zero.
   * 
   *
   *********************************************************************************/
  y1[0]=0; y1[dn-1]=0;
  y1[dn]=0; y1[2*dn-1]=0;
 
  /* Apply one step */
  for ( int i=1;i<dn-1;i++)
    y1[i]=a+y0[i]*y0[i]*y0[dn+i]-(b+1)*y0[i]+nu/dx*(y0[i+1]-2*y0[i]+y0[i-1])/dx; 

  for ( int i=dn+1;i<2*dn-1;i++)
  {
    y1[i]=b*y0[i-dn]-y0[i-dn]*y0[i-dn]*y0[i]+nu/dx*(y0[i+1]-2*y0[i]+y0[i-1])/dx; 
  }
  //printf("%lf ",y0[20]);
  //printf("\n");
  
  /*********************************************************************************
   *
   * The derivate w.r.t time $t$ is 1.
   *
   *********************************************************************************/
  y1[2*dn]=1;

}

Model *get_brusselator(int no_nodes)
{

  Brusselator* eq = new Brusselator(no_nodes);
  /* create initial condition */
  Real *x0 = new Real[no_nodes*2+1];
  
  // Apply initial condition: u(x,0) = 1 + sin(2*pi*x), v(x,0) = 3 
  for ( int i=0;i<no_nodes;i++)
  {
    Real x = 1.0/(no_nodes-1)*i;
    x0[i] = 1+sin(2*M_PI*x);
    x0[no_nodes+i] = 3;
  }
  x0[2*no_nodes]=0;
  
  Model* brusselator_model=new Model(x0,0,10,eq);
  return brusselator_model;
}
