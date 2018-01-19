import matplotlib.pyplot as plt
import pylab
import matplotlib as mpl
mpl.rcParams["savefig.directory"] = ("/home/aleksa/Dropbox/master/latex/thesis/Chapter1/Figs")
circle1 = plt.Circle((-1, 0), 1, color='g',clip_on=False)
fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
pylab.xlim([-4,2])
ax.grid(True, which='both')
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
pylab.ylim([-3,3])
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()

ax.add_artist(circle1)
plt.show();
