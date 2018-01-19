import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
if __name__ == "__main__":
    if ( len(sys.argv) < 2 ):
        raise ValueError("You need to specify the name of a problem") 
    problem = sys.argv[1]
    forw_op_types = ["det","rts_rk"]
    
    files = os.listdir(problem+"_inverse")
    variances=[]
    max_s = 50
    for op_type in forw_op_types:
        print op_type
        a=[]
        i = 5
        while i<= max_s:
            r1 = pd.read_csv("{}_inverse/{}_dist_{}.txt".format(problem,op_type,i))
            #a.append(pow(np.log(r1["x0"].var()),2))
            a.append(r1["x0"].var()+r1["x1"].var()+r1["x2"].var())
            #a.append(r1["x0"].var())
            #a.append(r1["x0"].std())
            i+=2
            print i
        variances.append(a)
    dist_var = []
    for i in range(0,len(variances[0])):
        #dist_var.append(variances[1][i]-variances[0][i])
        dist_var.append(np.sqrt(np.abs(variances[1][i]-variances[0][i])))
        #dist_var.append(variances[1][i]-variances[0][i])
    print dist_var
    h = []
    i=5
    while i<=max_s:
        h.append(1.0/i)
        i+=2
    #r1=plt.loglog(h,variances[0],"r")[0]
    #r2=plt.loglog(h,variances[1],"b")[0]
    plt.plot(h,dist_var)
    h10 = []
    for i in h:
        h10.append(i/10.0)
    r3=plt.loglog(h,h)
    r4=plt.loglog(h,dist_var) 
    plt.xlabel("h")
    plt.ylabel("var")
    #plt.legend([r1,r2],["RK","RTS-RK"])
    plt.show()
