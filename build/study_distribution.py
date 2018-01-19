import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
# read problem from command line 
if __name__ == "__main__":

    r1 = pd.read_csv("distribution.txt")
    r1_plot, = plt.plot( r1["x1"], r1["x2"],".")
    #r1_plot, = plt.scatter( r1["x1"])

    # adding legends
    #plt.xlabel('t')
    #plt.ylabel('error_estimate')    
    #plt.legend([r1_plot,r2_plot,r3_plot,r4_plot],\
            #["variance", "embeddeding error","sqrt(variance)","true_error"] )
    print r1["x1"].mean()
    plt.axvline(1, color='r', linestyle='dashed', linewidth=2)
    #print pd.sum(r1["x1"]*r1["x2"])/pd.sum(r1["x2"])
    print np.dot(r1["x1"],r1["x2"])/np.sum(r1["x2"])
    #print r1["x2"].mean()
    #print r1["x3"].mean()
    plt.show()
