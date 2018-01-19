import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import math

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

if __name__ == "__main__":
    x = []
    means = []
    std_devs = []
    means2 = []
    std_devs2 = []
    true_y = []
    df = pd.read_csv("pcn_res.txt")
    df2 = pd.read_csv("pcn_res_rts.txt")
    no_nodes=30
    for i in range(no_nodes):
        x_val=1.0/(no_nodes-1) * i
        x.append(x_val)
        means.append(df["x{}".format(i)].mean())
        std_devs.append(df["x{}".format(i)].std())

        true_y.append(1+math.sin(2*math.pi*x_val))
        if i==29:
            print means[-1]
    
    for i in range(no_nodes):
        means2.append(df2["x{}".format(i)].mean())
        std_devs2.append(df2["x{}".format(i)].std())

    fig, ax = plt.subplots()

    l1 = ax.errorbar(x,means , xerr=0.0, yerr=std_devs)
   # l2 = ax.errorbar(x,means2 , xerr=0.0, yerr=std_devs2)
    l3 = ax.errorbar(x,true_y, xerr=0.0, yerr=0.0)
    #print std_devs
    #print std_devs2
    #l3 = ax.plot(x,true_y)

    #plt.legend([l1,l2,l3],["det","stoch","true"])
    plt.legend([l1,l3],["pCN","true"])
    plt.show()
