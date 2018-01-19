#include "core_methods.hpp"
#include <algorithm>
using namespace std;
/* Basic numerical method class */
NumericalMethod::NumericalMethod(Model& m, Real h0, StepIntegrator& base_integrator) :model(m),h0(h0),base_integrator(base_integrator)
{
  n=model.get_equation()->get_dimension();
  cur_t = model.get_t0();
  cur_x=new Real[n];
  next_x=new Real [n];
  x_tmp = new Real[n];
  // init cur_x
  for ( int i=0;i<n;i++)
  {
    cur_x[i]=model.get_x0()[i];
  }
}
void NumericalMethod::print_pos()
{
  if ( f_print == 0 )
    return;
  p_count++;
  //if (p_count % 10 != 1 ) return;
  //if ( f_print ==0 )
    //return;
  fprintf(f_print,"%.16g,",cur_t);
  //fprintf(f_print,"%g,",cur_t);
  for ( int i=0;i<n;i++)
  {
    fprintf(f_print,"%.16g",x_tmp[i]);
   // fprintf(f_print,"%g",x_tmp[i]);
    if ( i < n-1)
      fprintf(f_print,",");
  }
  fprintf(f_print,"\n");
}

void NumericalMethod::print_pos(Real *data){
  p_count++;
  //if (p_count % 10 != 1 ) return;
  //if ( f_print ==0 )
    //return;
  fprintf(f_print,"%.16g,",cur_t);
  //fprintf(f_print,"%g,",cur_t);
  for ( int i=0;i<n;i++)
  {
    fprintf(f_print,"%.16g",data[i]);
   // fprintf(f_print,"%g",data[i]);
    if ( i < n-1)
      fprintf(f_print,",");
  }
  fprintf(f_print,"\n");
}
void NumericalMethod::set_f_print(std::string path){ 
  f_print=fopen(path.c_str(),"w");
  fprintf(f_print,"t,");
  for ( int i=0;i<n;i++)
  {
    fprintf(f_print,"v%d",i+1);
    if ( i != n-1 )
      fprintf(f_print,",");
  }
  fprintf(f_print,"\n");
}

void FixedTimeSteppingMethod::get_final(Real *x1)
{
  vector<Real> t_vec;
  t_vec.push_back(model.get_tf());
  auto res= integrate(t_vec);
  Real *x_res = res.begin()->second;
  for ( int i=0;i<n;i++)
    x1[i]=x_res[0];
  delete x_res;
}


void FixedTimeSteppingMethod::integrate_part(std::vector<Real> x_start, Real *x_end,Real t_start,Real t_end, int N)
{
  Real* x_s=new Real[n];
  for ( int i=0;i<n;i++)
    x_s[i]=x_start[i];
  integrate_part(x_s,x_end,t_start,t_end,N);
  delete[] x_s;
}

std::vector<std::pair<Real,Real*>> FixedTimeSteppingMethod::integrate(const std::vector<Real>& t_i)
{
  Real prev_t=model.get_t0();
  std::vector<Real> costs;
  std::vector<std::pair<Real, Real*>> result;
  //std::sort(t_i.begin(),t_i.end());
  for ( auto t = t_i.begin();t<t_i.end();t++)
  {
    costs.push_back(*t-prev_t);
    prev_t=*t;
  }
  Real cost_sum = model.get_tf()-model.get_t0();
  
  int N = cost_sum / h0;
  if ( N == 0 ) N++;
  prev_t = model.get_t0();

  // set cur_x, cur_t
  for ( auto t = t_i.begin();t<t_i.end();t++)
  {
    int k = (*t-prev_t)/cost_sum * N;
    //printf("%d\n",k);
    if ( k ==0 ) k =1 ;
    integrate_part(cur_x,cur_x,prev_t,*t,k);
    Real *x_t = new Real[n];
    for ( int i=0;i<n;i++)
      x_t[i]=cur_x[i];
    result.push_back(std::pair<Real,Real*>(*t,x_t));
    prev_t=*t;
  }
  print_pos();
  close_f();
  return result;
}

const double eps=1e-8;
/* Typical numerical method */
void ClassicalNumericalMethod::integrate_part( Real *x_start, Real *x_end,Real t_start,Real t_end,int N)
{
  for ( int i=0;i<n;i++)
    x_tmp[i] = x_start[i];

  double h=(t_end-t_start)/N;
  for ( int i=0;i<N;i++)
  {
    print_pos();
    this->base_integrator(x_tmp,next_x,h);
    f_evals+=this->base_integrator.get_f_evals();
    for ( int i=0;i<n;i++)
      x_tmp[i]=next_x[i];
    cur_t+=h;
  }
  for ( int i=0;i<n;i++)
    x_end[i]=x_tmp[i];
}



/* Conrad method with perturbed drift */

ConradNumericalMethod& ConradNumericalMethod::createWithNormalPerturbations( Model &m, int N, StepIntegrator &base_integrator )
{
  // calculate h
  Real h0 =  calc_h(m,N);
  Real h_2p_1=pow(h0,base_integrator.get_conv_order());
  NormalDistribution *rand_var = new NormalDistribution(0,h_2p_1*h_2p_1);
  ConradNumericalMethod* res = new ConradNumericalMethod(m,h0,base_integrator,*rand_var);
  return *res;
}

void ConradNumericalMethod::integrate_part( Real *x_start, Real *x_end,Real t_start,Real t_end,int N)
{
  for ( int i=0;i<n;i++)
    x_tmp[i] = x_start[i];
  Real h=(t_end-t_start)/N;
  for ( int i=0;i<N;i++)
  {
    print_pos();
    this->base_integrator(x_tmp,next_x,h);
    for ( int i=0;i<n;i++)
    {
      x_tmp[i]=next_x[i]+step_perturbator();
    }
    f_evals+=this->base_integrator.get_f_evals();
    cur_t+=h;
  }
  for ( int i=0;i<n;i++)
    x_end[i]=x_tmp[i];
}

void AbdulleGaregnaniMethod::integrate_part( Real *x_start, Real *x_end,Real t_start,Real t_end,int N)
{
  for ( int i=0;i<n;i++)
    x_tmp[i] = x_start[i];
  Real h=(t_end-t_start)/N;
  for ( int i=0;i<N;i++)
  {
    print_pos();
    Real next_h=perturb_step(h);
    this->base_integrator(x_tmp,next_x,next_h);
    for ( int i=0;i<n;i++)
      x_tmp[i]=next_x[i];
    f_evals+=this->base_integrator.get_f_evals();
    cur_t+=h;
  }
  for ( int i=0;i<n;i++)
    x_end[i]=x_tmp[i];
}

AbdulleGaregnaniMethod& AbdulleGaregnaniMethod::createWithUniformPerturbations( Model &m, int N, StepIntegrator &base_integrator,int no_paths )
{
  Real h0 =  calc_h(m,N);
  Real h_p = pow(h0,base_integrator.get_conv_order()+0.5);
  //UniformDistribution* step_distribution= new UniformDistribution(h0-h_p,h0+h_p);
  double coeff=1;
  UniformDistribution* step_distribution= new UniformDistribution(-coeff,coeff);
  AbdulleGaregnaniMethod* res = new AbdulleGaregnaniMethod(m,h0,base_integrator,step_distribution,no_paths);
  return *res;
}

AbdulleGaregnaniMethod& AbdulleGaregnaniMethod::createWithQuasiUniformPerturbations( Model &m, int N, StepIntegrator &base_integrator, int no_paths )
{
  Real h0 =  calc_h(m,N);
  Real h_p = pow(h0,base_integrator.get_conv_order()+0.5);
  double left_fact = 0.1;
  double right_fact = 10;
  //double right_fact = 1;
  AssymetricUniformDistribution *step_distribution=new AssymetricUniformDistribution(-left_fact,right_fact);
  AbdulleGaregnaniMethod* res = new AbdulleGaregnaniMethod(m,h0,base_integrator,step_distribution,no_paths);
  return *res;
}

Real AbdulleGaregnaniMethod::perturb_step(Real h)
{
  Real h_p=pow(h,base_integrator.get_conv_order()+0.5);
  return h+h_p*(*step_distribution)();
  //return (*step_distribution)();
}

void EulerExplicit::operator()( Real *x0, Real *x1, Real h)
{
  eq->eval(x0,x1);
  //printf("%d\n",eq->get_dimension());
  for ( int i=0;i<eq->get_dimension();i++)
  {
    x1[i]*=h;
    x1[i]+=x0[i];
  }
}
int EulerExplicit::get_f_evals(){return 1;}

void Heun::operator()( Real *x0, Real *x1, Real h )
{
  eq->eval(x0,x_m);
  for ( int i=0;i<eq->get_dimension();i++)
    x_tmp[i]=x0[i]+x_m[i]*h;
  eq->eval(x_tmp,x_t);
  for ( int i=0;i<eq->get_dimension();i++)
    x1[i]=x0[i]+h/2*(x_m[i]+x_t[i]);
  if (do_error_estimation)
  {
    error = 0;
    for ( int i=0;i<n;i++)
      error +=h*h/4*(x_m[i]-x_t[i])*(x_m[i]-x_t[i]); 
    error=sqrt(error);
  }
}

int Heun::get_f_evals(){return 2;}

void RK3::operator()(Real *x0, Real *x1, Real h)
{
  eq->eval(x0,k1);
  for ( int i=0;i<n;i++)
    tmp[i]=x0[i]+h/2*k1[i];
  eq->eval(tmp,k2);
  for ( int i=0;i<n;i++)
    tmp[i]=x0[i]+h/2*k2[i];
  eq->eval(tmp,k3);
  for ( int i=0;i<n;i++)
    x1[i]=x0[i]+h/3*(2*k2[i]+1*k3[i]);
}
int RK3::get_f_evals(){return 3;}

void RK4::operator()(Real *x0, Real*x1, Real h)
{
  eq->eval(x0,k1);
  for ( int i=0;i<n;i++)
    tmp[i]=x0[i]+h/2*k1[i];
  eq->eval(tmp,k2);
  for ( int i=0;i<n;i++)
    tmp[i]=x0[i]+h/2*k2[i];
  eq->eval(tmp,k3);
  for ( int i=0;i<n;i++)
    tmp[i]=x0[i]+h*k3[i];
  eq->eval(tmp,k4); 
  for ( int i=0;i<n;i++)
    x1[i]=x0[i]+h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
}

int RK4::get_f_evals(){return 4;}
