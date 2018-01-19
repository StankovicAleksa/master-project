#ifndef INIT_H
#define INIT_H
#include "../models/models.hpp"
#include <getopt.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "../methods/methods.hpp"
void process_command_line(int argc, char * argv[],Model*&model, StepIntegrator *&step_integrator, NumericalMethodFactory* &numerical_method_factory)
{
  int index;
  int c;

  model=0;
  step_integrator=0;
  numerical_method_factory=0;
  char *s_step_integrator=0;
  opterr = 0;
  while ((c = getopt (argc, argv, "m:s:f:")) != -1)
    switch (c)
      {
      case 'm':
      {
        if ( strcmp(optarg,"po")== 0 )
            model=get_po_equation();
        else if ( strcmp(optarg,"lorenz")== 0 )
            model=get_lorenz_system();
        else if ( strcmp(optarg,"test")== 0 )
            model=get_test_equation();
        else if ( strcmp(optarg,"fitz") == 0 )
          model=get_fitz();
        else if ( strcmp(optarg,"brusselator") == 0 )
          model=get_brusselator(30);
        break;
      }
      case 's':
      {
        s_step_integrator=optarg;
        break;
      }
      case 'f':
      {
        if ( strcmp(optarg,"fixed")==0)
          numerical_method_factory=new NumericalMethodFactory("fixed_timestepping");
        else if ( strcmp(optarg,"agm")==0 )
          numerical_method_factory=new NumericalMethodFactory("abdulle_garegnani");
        else if ( strcmp(optarg,"add")==0 )
          numerical_method_factory=new NumericalMethodFactory("conrad");
        else if ( strcmp(optarg,"adapt")==0 )
          numerical_method_factory=new NumericalMethodFactory("classical_adaptive");
        else if ( strcmp(optarg,"r_adapt")==0 )
          numerical_method_factory=new NumericalMethodFactory("r_adaptive");
        else if ( strcmp(optarg,"a_adapt")==0 )
          numerical_method_factory=new NumericalMethodFactory("aleksa_adaptive");
        else if ( strcmp(optarg,"giacomo2")==0 )
          numerical_method_factory=new NumericalMethodFactory("giacomo2");
        break;
      }
      case '?':
        if (optopt == 'm' || optopt == 's' || optopt=='f' )
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return;
      default:
        abort ();
      }
  
  for (index = optind; index < argc; index++)
    fprintf (stderr,"Non-option argument %s\n", argv[index]);
  if ( model ==0 || numerical_method_factory ==0 || s_step_integrator==0)
  {
    fprintf(stderr,"Bad command line arguments\n");
    abort();
  } 
  if( strcmp(s_step_integrator, "ee")==0 )
  {
    step_integrator=new EulerExplicit(model->get_equation());
  }
  else if( strcmp(s_step_integrator, "rk4")==0 )
  {
   step_integrator = new RK4(model->get_equation());
  }
  else if( strcmp(s_step_integrator, "rk3")== 0 )
  {
   step_integrator = new RK3(model->get_equation());
  }
  else if( strcmp(s_step_integrator, "heun")== 0)
  {
   step_integrator = new Heun(model->get_equation());
  }
  else if( strcmp(s_step_integrator, "rkc")==0 )
  {
   step_integrator = new RKC(model->get_equation());
  }
  else if( strcmp(s_step_integrator, "rock2")== 0 )
  {
   step_integrator =new  ROCK2(model->get_equation());
  }
  return;
}
#endif
