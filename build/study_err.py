import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
# read problem from command line 
if __name__ == "__main__":

    if ( len(sys.argv) < 2 ):
        raise ValueError("You need to specify the name of a problem") 
    problem = sys.argv[1]

    r1 = pd.read_csv("{}/err.txt".format(problem))
    r_exact = pd.read_csv("{}/exact.txt".format(problem))
    r_approximated = pd.read_csv("{}/approximation.txt".format(problem))


    r_approximated.set_index("t",inplace=True)
    r_exact.set_index("t",inplace=True)
    r_approximated[r_approximated.index.duplicated()]
    r_exact.index.union(r_approximated.index)
    r_exact = r_exact.reindex(r_exact.index.union(r_approximated.index),method="ffill")



    #r_interp = r_exact.interpolate(method='piecewise_polynomial',order=1)
    r_interp = r_exact
    diff =r_approximated-r_interp
    diff.dropna(inplace=True)
    print len(diff)
    diff["err"]=diff.apply(lambda x: np.sqrt(x.dot(x)), axis=1)
    diff["err"][0] = 2e-5
    #print r_interp
    d = r1["var"]-r1["var"].shift(1)
    e = pd.ewma(r1["var"],span=2.3)
    
    r_plots=[] 
    legends=[]

    #r_plots.append(plt.plot( r1["t"], np.log(r1["var"]),"-")[0])
    #legends.append("$\sqrt{tr(Var(Y_n))}$")
    r_plots.append(plt.plot( r1["t"], np.log(np.clip(np.abs(r1["var"]-r1["var"].shift(1)),pow(10,-10),None)),"-")[0])
    legends.append("$\hat{E_n}$")

    r_plots.append(plt.plot( r1["t"], np.log(r1["emb_err"]),"-")[0])
    legends.append("embedded error")
    #r_plots.append(plt.plot( r1["t"], np.log(r1["sqrt(var)"]),"-")[0])
    #r_plots.append(plt.plot( r1["t"], np.log(r1["var"]),"-")[0])
    #legends.append(r"$tr(Var(Y_n))^{0.25}$")
    r_plots.append(plt.plot( diff.index, np.log(diff["err"]),"-")[0])
    legends.append("true error")

    plt.legend(r_plots,legends)
    plt.show()   
    # add time-steps for most precise simulation
    types=["fixed_timestepping","classical_adaptive","aleksa_adaptive"]
    #types=["aleksa_adaptive"]
    #types=[] 
    files = os.listdir(problem)
    fig = plt.figure()
    for tp in types:
        f_files =[]
        for f in files:
            if f.startswith(tp):
                f_files.append(f)
        f_stats=pd.read_csv("{}/stats_{}.txt".format(problem,tp))
        df = pd.read_csv("{}/{}{}_0_.txt".format(problem,tp,len(f_stats)))
        # calculate err
        r_plots.append(plt.plot(df["t"],-df["t"].shift(1)+df["t"],"-")[0])
        #legends.append("h_{}".format(tp))
        print df


    # adding legends
    plt.xlabel('t')
    plt.ylabel('h_n')
    plt.legend(r_plots,["fixed time-stepping", "A-ADAPT","VAR-ADAPT" ])
    #plt.legend( loc="upper right" )
    plt.show()
