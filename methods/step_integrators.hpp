class EulerExplicit:public StepIntegrator 
{
  public:
  EulerExplicit(Equation *e):StepIntegrator(e){}
  void operator()(Real *x0, Real *x1, Real h);
  int get_conv_order(){return 1;}
  int get_f_evals();
};

class Heun:public StepIntegrator 
{
  private:
    Real *x_m,*x_t;
    Real *x_tmp;
  public:
  Heun(Equation *e):StepIntegrator(e)
  {
    x_m=new Real[e->get_dimension()];
    x_t=new Real[e->get_dimension()];
    x_tmp = new Real[e->get_dimension()];
  }
  ~Heun(){delete[] x_m;delete[] x_t;delete[] x_tmp;}
  void operator()(Real *x0, Real *x1, Real h);
  int get_f_evals();
  int get_conv_order(){return 2;}
};

class RK3:public StepIntegrator
{
  private:
    Real *k1,*k2,*k3;
    Real *tmp;
  public:
    RK3(Equation *e):StepIntegrator(e)
    {
      k1=new Real[e->get_dimension()];
      k2=new Real[e->get_dimension()];
      k3=new Real[e->get_dimension()];
      tmp = new Real[e->get_dimension()];
    }
    ~RK3(){delete[] k1;delete[] k2;delete[] k3;;delete[] tmp;}
    void operator()(Real *x0, Real *x1, Real h);
    int get_conv_order(){return 3;}
    int get_f_evals();
};

class RK4:public StepIntegrator
{
  private:
    Real *k1,*k2,*k3,*k4;
    Real *tmp;
  public:
    RK4(Equation *e):StepIntegrator(e)
    {
      k1=new Real[e->get_dimension()];
      k2=new Real[e->get_dimension()];
      k3=new Real[e->get_dimension()];
      k4=new Real[e->get_dimension()];
      tmp = new Real[e->get_dimension()];
    }
    ~RK4(){delete[] k1;delete[] k2;delete[] k3;delete[] k4;delete[] tmp;}
    void operator()(Real *x0, Real *x1, Real h);
    int get_conv_order(){return 4;}
    int get_f_evals();
};


