import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

if __name__ == "__main__":
    fig = plt.figure()
    ax = fig.gca(projection='3d')



    X = np.arange(0, 1, 1.0/30)
    Y = np.arange(0, 1, 1.0/200)
    #X, Y = np.meshgrid(X, Y)
    X = []
    Y = []
    Z = []
    r_exact = pd.read_csv("test.txt")
    no_nodes=30
    for index,row in r_exact.iterrows():
        y_tmp=[]
        x_tmp=[]
        z_tmp=[]
        for i in range(no_nodes):
            y_tmp.append(row['t'])
            x_tmp.append(1.0/(no_nodes-1)*i)
            z_tmp.append(row['v{}'.format(i+1)])
            #z_tmp.append(row['v{}'.format(i+1+no_nodes)])
        #plt.plot(x_tmp,z_tmp)
        #plt.show()
        #print row['v31']
        X.append(x_tmp)
        Y.append(y_tmp)
        Z.append(z_tmp)
        #print z_tmp
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)


    # Scalling 
    x_scale=1
    y_scale=10
    z_scale=2

    scale=np.diag([x_scale, y_scale, z_scale, 1.0])
    scale=scale*(1.0/scale.max())
    scale[3,3]=1.0

    def short_proj():
      return np.dot(Axes3D.get_proj(ax), scale)

    ax.get_proj=short_proj


    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    # Customize the z axis.
    #ax.set_zlim(-1.01, 1.01)
    #ax.zaxis.set_major_locator(LinearLocator(10))
    #ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    #ax.zaxis.set_major_locator(LinearLocator(10))
    #ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    #fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

    #print X
    #print Y

    #exact_sol = pd.read_csv("{}/exact.txt".format(problem)).iloc[-1]
    ## we also need number of $f$ 
    
    #plt.show()
