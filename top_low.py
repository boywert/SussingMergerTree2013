from pylab import *
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

data = loadtxt('/home/boywert/sshfs/pact/home/c/cs/cs390/data_snaplist.txt')
x = data[:,3]
y = data[:,0]
fig = figure()
ax1 = fig.add_subplot(111)
ax1.scatter(x,y, s=20,lw = 0,color='k')

#ax1.set_xscale('log')
ax1.set_ylabel('Snapshot ID')
ax1.set_xlabel(r'$t/t_0$')
ax1.set_ylim([-.1,61.1])
ax1.set_xlim([-0.001,1.001])
ax2 = twiny() # ax2 is responsible for "top" axis and "right" axis
ax2.set_xticks([0.0035,0.23,0.44 , 0.62, 0.81, 1.0])

ax2.set_xticklabels(["50.0","2.0", "1.0","0.5", "0.2", "0.0"])
ax2.set_xlabel('redshift')
#ax2.axis["right"].major_ticklabels.set_visible(False)

fig.savefig("snap_info.pdf",bbox_inches='tight')
#fig.show()
