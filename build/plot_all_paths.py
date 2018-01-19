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
    scheme = "abdulle_garegnani"
    r_plots=[]
    # add time-steps for most precise simulation
    f_stats=pd.read_csv("{}/stats_{}.txt".format(problem,scheme))
    t=1.0
    for i in range(len(f_stats)):
        if ( i>=10):
            break
        print i
        df = pd.read_csv("{}/{}{}.txt".format(problem,scheme,i+1))
        # calculate err
        if (i==0):
            r=2.0
        else:
            r=1.4
        col = (t/r, t/r, t/r)
        r_plots.append(plt.plot(df["t"],df["v1"],"-",color=col,linewidth=0.5)[0])
        #legends.append("ts_{}".format(tp))

    # adding legends
    plt.xlabel('t')
    plt.ylabel('x')    
    #plt.legend(r_plots,legends)
    
    plt.show()
