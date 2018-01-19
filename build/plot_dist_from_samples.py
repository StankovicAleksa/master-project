import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import seaborn as sns
import os
# read problem from command line 
if __name__ == "__main__":
    sns.set(color_codes=True)
    mpl.rcParams["savefig.directory"] = ("/home/aleksa/Dropbox/master/latex/thesis/Chapter3/Figs")
    r1 = pd.read_csv("distribution.txt")
    #r1 = pd.read_csv("samples.txt")
    #r1 = pd.read_csv("mcmc_out.txt")
    #r1 = pd.read_csv("det_out.txt")
    #r1 = pd.read_csv("pmmh_out.txt")
    #r1 = pd.read_csv("skewed_out.txt")

    #r1 = pd.read_csv("smc_distribution.txt")
    #r1 = pd.read_csv("pf_distribution.txt")
    #r1_plot, = plt.plot( r1["x1"], r1["x2"],".")
    #r1["a"] = np.exp(r1["x0"])
    #r1["b"] = np.exp(r1["x1"])
    #r1["c"] = np.exp(r1["x2"])
    print (len(r1.drop_duplicates()))
    r1[r"$\lambda$"] = r1["x0"]
    #r1["b"] = r1["x1"]
    #r1["c"] = r1["x2"]
    size = (4,4)
    plt.figure(figsize=size)
    sns.distplot(r1[r"$\lambda$"])
    plt.axvline(1, color='r', linestyle='dashed', linewidth=2)
    #plt.suptitle("MCMC")
    #plt.suptitle("PMMH")
    #plt.suptitle("PMMH (skewed)")
    #plt.axvline(0.2, color='r', linestyle='dashed', linewidth=2)
    plt.draw()

    #sns.distplot(r1["a"])
    #plt.figure(figsize=size)
    #sns.distplot(r1["b"])
    #plt.axvline(0.2, color='r', linestyle='dashed', linewidth=2)
    #plt.draw()
    #plt.figure(figsize=size)
    #sns.distplot(r1["c"])
    #plt.axvline(3.0, color='r', linestyle='dashed', linewidth=2)
    #plt.draw()
    #print r1["a"].mean()
    plt.show()
    #r1_plot, = plt.scatter( r1["x1"])
    
    # adding legends
    #plt.xlabel('t')
    #plt.ylabel('error_estimate')    
    #plt.legend([r1_plot,r2_plot,r3_plot,r4_plot],\
            #["variance", "embeddeding error","sqrt(variance)","true_error"] )
    #print pd.sum(r1["x1"]*r1["x2"])/pd.sum(r1["x2"])
    #print np.dot(r1["x1"],r1["x2"])/np.sum(r1["x2"])
    #print r1["x2"].mean()
    #print r1["x3"].mean()
