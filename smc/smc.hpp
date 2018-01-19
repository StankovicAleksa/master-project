#ifndef SMC_H
#define SMC_H
#include "../models/models.hpp"
#include "../tools/mcmc.hpp"
#include "../models/inverse_model.hpp"
#include <vector>
#include <utility>
class Particle
{
  private:
    Real weight;
    std::vector<Real>params;
    std::vector<Real> x;
  public:
    Particle(Real w,std::vector<Real> params, std::vector<Real>x):weight(w),params(params),x(x) { }
    std::vector<Real>& get_params() {return params;}
    Real get_weight()const {return weight;}
    std::vector<Real>& get_x(){return x;}
    void set_weight(Real w){weight=w;}
};

class PF_SMC
{
  private:
    int no_particles;
    InverseProblem* inverse_problem;
    Model* true_model;
    Real a;
    std::vector<Particle> particles;
    int d;
    int no_steps;
  public:
    PF_SMC(int no_particles,InverseProblem* inverse_problem,Real a=0.97, int no_steps=50000):no_particles(no_particles),inverse_problem(inverse_problem),a(a),no_steps(no_steps)
    {
      true_model = inverse_problem->get_true_model();
      d=true_model->get_dim();
    }
    void init_particles();
    std::vector<Particle> estimate_parameters();
};

class Stuart_SMC
{
  private:
    int no_particles;
    InverseProblem* inverse_problem;
    Model* true_model;
    std::vector<Particle> particles;
    int d;
    int no_steps;
    int J;
  public:
    Stuart_SMC(int no_particles,InverseProblem* inverse_problem,int J=200,int no_steps=50000):no_particles(no_particles),inverse_problem(inverse_problem),no_steps(no_steps),J(J)
    {
      true_model = inverse_problem->get_true_model();
      d=true_model->get_dim();
    }
    void init_particles();
    std::vector<Particle> estimate_parameters();
};

#endif
