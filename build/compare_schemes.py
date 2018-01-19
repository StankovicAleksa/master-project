import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
if __name__ == "__main__":
    #types=["fixed_timestepping","r_adaptive","aleksa_adaptive","classical_adaptive"]
    types=["fixed_timestepping","classical_adaptive","aleksa_adaptive"]
    #types=["fixed_timestepping","classical_adaptive","r_adaptive"]
    #types=["abdulle_garegnani","fixed_timestepping"]
    #types=["abdulle_garegnani","conrad","fixed_timestepping"]
    if ( len(sys.argv) < 2 ):
        raise ValueError("You need to specify the name of a problem") 
    problem = sys.argv[1]

    files = os.listdir(problem)

    r_plots = []
    exact_sol = pd.read_csv("{}/exact.txt".format(problem)).iloc[-1]
    for tp in types:
        f_files =[]
        for f in files:
            if f.startswith(tp):
                f_files.append(f)
        f_stats=pd.read_csv("{}/stats_{}.txt".format(problem,tp))
        print f_stats
        for i in range(len(f_stats)):
            a =0 
            no_paths = 1
            for j in range(no_paths):
                df = pd.read_csv("{}/{}{}_{}_.txt".format(problem,tp,i+1,j))
                if ( i==0 ):
                    err = pd.DataFrame(columns=df.columns)
                a+=df.iloc[-1]
            a/=no_paths
            err=err.append(a)
        # calculate err
        print err
        print exact_sol
        err=err-exact_sol
        err["err"]=err.apply(lambda x: np.sqrt(x.dot(x)), axis=1)
        if tp.startswith("aleksa_adaptive"):
            f_stats["f_evals"]/=1000
        if tp.startswith("giacomo2"):
            pass
            #f_stats["f_evals"]/=1000
        print err["err"]
        r_plots.append(plt.loglog(err["err"],f_stats["f_evals"],"-o")[0])
    # we also need number of $f$ 
    #types[1]="adaptive algorithm"
    #types[0]="fixed time-stepping algorithm"
    #types[2]="Runge-Kutta"
    #types[0] = "RTS-RK"
    #types[1] = "AN-RK"
    #types[2]="Heun"
    types[0] = "fixed time-stepping"
    types[1] = "A-ADAPT"
    types[2]="VAR-ADAPT"
    plt.xlabel('Error')
    plt.ylabel('Number of f evaluations')
    plt.gca().invert_xaxis()
    plt.legend(r_plots,types)
    plt.show()
