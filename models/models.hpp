#ifndef MODELS_H
#define MODELS_H

#include <cstdio>
#include "../config.h"
#include <string>
class Equation
{
  protected:
    int dim;
    Real *eigvec_tmp1,*eigvec_tmp2,*eigvec_tmp3,*eigvec_prev;
  public:
    Equation(int dim):dim(dim),eigvec_prev(0)
    { 
      eigvec_tmp1=eigvec_tmp2=eigvec_tmp3=0; 
    }
    Equation(const Equation& equation)
    {
      dim=equation.dim;
      eigvec_tmp1=eigvec_tmp2=eigvec_tmp3=0; 
      eigvec_prev=0;
    }
    virtual ~Equation()
    {
      if ( eigvec_prev != 0 ) delete eigvec_prev;
      if ( eigvec_tmp1 != 0 ) delete eigvec_tmp1;
      if ( eigvec_tmp2 != 0 ) delete eigvec_tmp2;
      if ( eigvec_tmp3 != 0 ) delete eigvec_tmp3;
    }
    int get_dimension() {return dim;}
    virtual void eval(const Real *y0, Real *y1)=0;
    Real largest_eigenval(const Real *y0);
    Real directional_eigenval(const Real *y0,const Real *y1);
    virtual std::string get_name()=0;
};

class Model
{
  private:
    Real *x0;
    Real t0,tf;
    Equation* equation;
    int dim;
  public:
    Model(Real *x0,Real t0, Real tf,Equation *equation):t0(t0), tf(tf), equation(equation)
    {
      dim = equation->get_dimension();
      this->x0 = new Real[dim];
      for ( int i=0;i<dim;i++)
      {
        this->x0[i]=x0[i];
      }
    }
    Model (const Model &model );
    virtual ~Model()
    {
      delete equation;
      delete[] x0;
    }
    Real get_t0() const { return t0;} 
    Real get_tf() const { return tf;} 
    Real* get_x0(){ return x0;}
    int get_dim(){return dim;}
    Equation* get_equation() {return equation;}
    std::string get_name(){ return equation->get_name();}
};


class TestEquation:public Equation
{
  private:
    Real lambda;
  public:
    TestEquation(int n,Real lambda):Equation(n),lambda(lambda) {}
    void eval (const Real *y0, Real *y1);
    Real largest_eigenval(const Real *y0);
    std::string get_name(){ return "test" ; }
};


class PeroxideOxideEquation:public Equation
{
  private:
    Real A0,B0,X0;
    Real k[9];
  public:
    PeroxideOxideEquation(Real A0,Real B0,Real X0,Real*k):Equation(4),A0(A0),B0(B0),X0(X0)
    {
      for ( int i=1;i<=8;i++)
        this->k[i]=k[i];
    }
    void eval(const Real *y0, Real *y1);
    Real largest_eigenval(const Real *y0);
    std::string get_name(){ return "po" ; }
};

class Lorenz:public Equation
{
  private:
    Real sigma,rho,beta;
  public:
    Lorenz(Real sigma,Real rho,Real beta):Equation(3),sigma(sigma),rho(rho),beta(beta) { }
    void eval(const Real *y0, Real *y1);
    std::string get_name(){ return "lorenz" ; }
};

class FitzHughNagumo:public Equation
{
  private:
    Real a,b,c;
  public:
    FitzHughNagumo(Real a,Real b,Real c):Equation(2),a(a),b(b),c(c) { }
    void eval(const Real *y0, Real *y1);
    std::string get_name(){ return "fitz" ; }
};

class Brusselator:public Equation
{
  /*********************************************************************************
   *
   * Brusselator problem: 
   *                        du/dt = a + u*u*v-(b+1)*u+\nu*d^2u/dx^2,
   *                        dv/dt = b*u-u*u*v+\nu*d^2v/dx^2,
   * where x \in ]0,1[
   *
   *
   * For now, only with constant boundary condition:
   *
   *                        u(0,t)= 1, u(1,t) = 3, 
   * is implemented. Furthermore, at every point $t$, the solution (u,v)(t) is discretized 
   * on a grid with disc_nodes+2 points.
   *
   *
   *********************************************************************************/
  private:
    int dn;
    Real a,b;
    Real nu;
    Real dx;
  public:
    Brusselator(int disc_nodes=30);

    void eval(const Real *y0, Real *y1);
    std::string get_name(){ return "brusselator" ; }
};

Model* get_test_equation();
Model* get_po_equation();
Model* get_lorenz_system();
Model* get_fitz();
Model* get_brusselator(int no_nodes);

#endif
