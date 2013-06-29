#!/usr/bin/python
import os 
import pylab
import sys
import numpy
import math

inputFile = 'run.list'
command = './read'

#default flag
reset_data_flag = 0
linestyles = ['','','','','','','','','','']
colors = ('b','g','r','c','m','#00FF00','k','#FFA500','#7F7F7F','#800000')
for option in sys.argv:
    if(option == '--reset-data'):
        reset_data_flag = 1
        
dat_folder = "dat/"

if(reset_data_flag == 1):
    os.system("make clean all")

if(os.path.exists('./read') < 1):
    sys.exit()

os.system("mkdir -p dat pdf100 pdf1000")

f = open(inputFile)
line = f.read().splitlines()

line.sort()

print line

nFile = len(line)

for (i,item) in enumerate(line):
    filename = item.split()
    line[i] = filename[0]
    if(reset_data_flag == 1):
        os.system(command+" "+item  )

pdf_folder = "pdf/"


os.system("mkdir -p dat "+pdf_folder+" && rm -f "+pdf_folder+"*")

if(option == '--reset-data'):
    os.system("mkdir -p dat "+pdf_folder+" && rm -f "+pdf_folder+"*")	
    os.system("mv *.dat dat/")

markersets = '+.o*psxDh^'

pylab.rc('text', usetex=True)

######################################################################################
#To convert arctan axis

def tick_function_half(X):
    V = numpy.tan(X*numpy.pi/2)
    return ["%.3f" % z for z in V]

def tick_function_full(X):
    V = numpy.tan(X*numpy.pi)
    return ["%.3f" % z for z in V]

###########################################################################
#plot mean displacement HOST
pre_ext = "_100"
host_dispmain_ext = pre_ext+"DispMainHost.dat"
print 'Start plotting',host_dispmain_ext

fig = pylab.figure()
ax1 = fig.add_subplot(111)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + host_dispmain_ext)
    ax1.plot(data[:,0],(data[:,1]+1),linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    ax1.hold(True)

ax1.hold(False)
ax1.set_xlabel(r"$\Delta_r$")
ax1.set_yscale('log')
ax1.set_xlim([0.0,1.0])
ax1.set_ylabel(r"$N+1$")

new_tick_locations = numpy.array([math.atan(0.1)*2/numpy.pi, math.atan(0.5)*2/numpy.pi, math.atan(1.)*2/numpy.pi, math.atan(2.)*2/numpy.pi,math.atan(3.)*2/numpy.pi,math.atan(10.)*2/numpy.pi])

ax1.set_xticks(new_tick_locations)
ax1.set_xticklabels(tick_function_half(new_tick_locations))

leg = ax1.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)
pylab.savefig(pdf_folder+'disp_main_host.pdf',bbox_inches='tight')
os.system('pdftops -eps '+pdf_folder+'disp_main_host.pdf')
#######################################################################################################################
# Plot length of chain
pre_ext = "_100"
linkdepth_ext = pre_ext+"LinkDepth.dat"
print 'Start plotting',linkdepth_ext
fig = pylab.figure()
ax1 = fig.add_subplot(111)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + linkdepth_ext)
    ax1.plot(data[:,0],(data[:,1]),linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
pylab.ylabel(r"$N$")
pylab.xlabel("Depth of Chains")

ax1.legend(loc='lower left', handlelength = 6,ncol=1, prop={'size':9})
pylab.savefig(pdf_folder+'linked_depth.pdf',bbox_inches='tight')

#######################################################################################################################

#######################################################################################################################
# Plot length of chains low/high mass
# Plot length of chain (low mass)

pre_ext = "_100"
linkdepth_host_ext = pre_ext+"LinkDepthHost.dat"
print 'Start plotting',linkdepth_host_ext
fig = pylab.figure()
ax1 = fig.add_subplot(311)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + linkdepth_host_ext)
    ax1.plot(data[:,0],(data[:,1])+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
#pylab.ylabel(r"$N$")
#pylab.xlabel("Depth of Chains (Low mass)")



# Plot length of chain (high mass)
pre_ext = "_500"
linkdepth_host_ext = pre_ext+"LinkDepthHost.dat"
print 'Start plotting',linkdepth_host_ext

ax2 = fig.add_subplot(312,sharex=ax1)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + linkdepth_host_ext)
    ax2.plot(data[:,0],(data[:,1])+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
pylab.ylabel(r"$N+1$")
#pylab.xlabel("Depth of Chains (High mass)")


# Plot length of chain (high mass)
pre_ext = "_1000"
linkdepth_host_ext = pre_ext+"LinkDepthHost.dat"
print 'Start plotting',linkdepth_host_ext

ax3 = fig.add_subplot(313,sharex=ax1)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + linkdepth_host_ext)
    ax3.plot(data[:,0],(data[:,1])+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
#pylab.ylabel(r"$N$")
pylab.xlabel("Depth of Chains")

ax1.legend(loc='upper right', handlelength = 4,ncol=1, fancybox=True, prop={'size':6})
pylab.savefig(pdf_folder+'linked_depth_host.pdf',bbox_inches='tight')

#######################################################################################################################
# Plot length of chains low/high mass
# Plot length of chain (low mass)
pre_ext = "_100"
linkdepth_ext = pre_ext+"LinkDepthResi.dat"
print 'Start plotting',linkdepth_ext
fig = pylab.figure()
ax1 = fig.add_subplot(311)
for (i,item) in enumerate(line):
    data = None
    data = pylab.loadtxt(dat_folder + item + linkdepth_ext)
    print dat_folder + item + linkdepth_ext
    print data
    ax1.plot(data[:,0],(data[:,1])+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
#pylab.ylabel(r"$N$")
#pylab.xlabel("Depth of Chains (Low mass)")



# Plot length of chain (high mass)
pre_ext = "_500"
linkdepth_ext = pre_ext+"LinkDepthResi.dat"
print 'Start plotting',linkdepth_ext

ax2 = fig.add_subplot(312,sharex=ax1)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + linkdepth_ext)
    ax2.plot(data[:,0],(data[:,1])+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
pylab.ylabel(r"$N+1$")
#pylab.xlabel("Depth of Chains (High mass)")


# Plot length of chain (high mass)
pre_ext = "_1000"
linkdepth_ext = pre_ext+"LinkDepthResi.dat"
print 'Start plotting',linkdepth_ext

ax3 = fig.add_subplot(313,sharex=ax1)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + linkdepth_ext)
    ax3.plot(data[:,0],(data[:,1])+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
#pylab.ylabel(r"$N$")
pylab.xlabel("Depth of Chains")

ax1.legend(loc='upper right', handlelength = 4,ncol=1, fancybox=True, prop={'size':6})
pylab.savefig(pdf_folder+'linked_depth_resi.pdf',bbox_inches='tight')
#######################################################################################################################
# Plot length of chains low/high mass
# Plot length of chain (low mass)

pre_ext = "_100"
linkdepth_ext = pre_ext+"LinkDepth.dat"
print 'Start plotting',linkdepth_ext
fig = pylab.figure()
ax1 = fig.add_subplot(311)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + linkdepth_ext)
    ax1.plot(data[:,0],(data[:,1])+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
#pylab.ylabel(r"$N$")
#pylab.xlabel("Depth of Chains (Low mass)")



# Plot length of chain (high mass)
pre_ext = "_500"
linkdepth_ext = pre_ext+"LinkDepth.dat"
print 'Start plotting',linkdepth_ext

ax2 = fig.add_subplot(312)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + linkdepth_ext)
    ax2.plot(data[:,0],(data[:,1])+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
pylab.ylabel(r"$N+1$")
#pylab.xlabel("Depth of Chains (High mass)")


# Plot length of chain (high mass)
pre_ext = "_1000"
linkdepth_ext = pre_ext+"LinkDepth.dat"
print 'Start plotting',linkdepth_ext

ax3 = fig.add_subplot(313)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + linkdepth_ext)
    ax3.plot(data[:,0],(data[:,1])+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
#pylab.ylabel(r"$N$")
pylab.xlabel(r"$l$")
ax1.set_xticklabels([])
ax2.set_xticklabels([])

ax1.set_yticks(numpy.array([1,10,100,1000,10000]))
ax2.set_yticks(numpy.array([1,10,100]))
ax3.set_yticks(numpy.array([1,10,100]))
#ax1.set_xticklabels()

pylab.subplots_adjust(hspace = 0.0001)
leg = ax1.legend(loc='upper right', handlelength = 6,ncol=1, fancybox=False, prop={'size':6})
leg.get_frame().set_linewidth(0)
pylab.savefig(pdf_folder+'linked_depth_mass.pdf',bbox_inches='tight')
os.system('pdftops -eps '+pdf_folder+'linked_depth_mass.pdf')
######################################################################################
# Plot Mass fluctuation (betacor) Highmass
pre_ext = "_1000"
h_betacorr_ext = pre_ext+"BetaCorr.dat"
print "Start plotting",h_betacorr_ext
fig = pylab.figure()
ax1 = fig.add_subplot(111)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + h_betacorr_ext)
    l = len(data)
    ax1.plot(data[:,0],data[:,1]+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
pylab.xlabel(r"$\pi^{-1}\xi_M$")
pylab.ylabel(r"$N+1$")

pylab.legend(loc='upper right', ncol=2, prop={'size':8})
pylab.savefig(pdf_folder+'betacorrhigh.pdf',bbox_inches='tight')
os.system('pdftops -eps '+pdf_folder+'betacorrhigh.pdf')
##############################################################################################

#plot tree spanning
#plot nMergers
pre_ext = "_100"
nmerger_ext = pre_ext+"Nmergers.dat"
print "Start plotting", nmerger_ext
fig = pylab.figure()
ax1 = fig.add_subplot(111)
for (i,item) in enumerate(line):
        data = pylab.loadtxt(dat_folder +item + nmerger_ext)
        ax1.plot(data[:,0],(data[:,1])+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
        pylab.hold(True)

pylab.hold(False)

#pylab.hold(False)
pylab.yscale('log')
pylab.ylabel(r"$N+1$")
pylab.xlabel(r"$N_{\rm{Mergers}}$")
leg = pylab.legend(loc='upper right', handlelength = 8,ncol=1, prop={'size':10})
leg.get_frame().set_linewidth(0)
pylab.savefig(pdf_folder+'Nmergers.pdf',bbox_inches='tight')
os.system('pdftops -eps '+pdf_folder+'betacorrhigh.pdf')


##############################################################################################
#plot splash mass  host
pre_ext = "_100"
host_lostmass_ext = pre_ext+"Lostmass_host.dat"
print 'Start plotting',host_lostmass_ext

fig = pylab.figure()
ax1 = fig.add_subplot(111)
for (i,item) in enumerate(line):
    if(os.path.exists(dat_folder +item + host_lostmass_ext) >= 1):
        data = pylab.loadtxt(dat_folder +item + host_lostmass_ext)
        ax1.plot(data[:,0]/2.,data[:,1]+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
        pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
#pylab.xlim([0.,1.])
pylab.ylabel(r"$N+1$")
pylab.xlabel(r"$\Delta N_{\rm{merge}}$")
leg = pylab.legend(loc='upper right', handlelength = 8,ncol=1, prop={'size':10})
leg.get_frame().set_linewidth(0)
pylab.savefig(pdf_folder+'Lost_N.pdf',bbox_inches='tight')
os.system('pdftops -eps '+pdf_folder+'Lost_N.pdf')
##############################################################################################
#plot splash mass  host
pre_ext = "_100"
lostmass_ext = pre_ext+"Lostmass.dat"
print 'Start plotting',lostmass_ext

fig = pylab.figure()
ax1 = fig.add_subplot(111)
for (i,item) in enumerate(line):
    if(os.path.exists(dat_folder +item + lostmass_ext) >= 1):
        data = pylab.loadtxt(dat_folder +item + lostmass_ext)
        ax1.plot(data[:,0]-1.,data[:,1]/sum(data[:,1]),linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
        pylab.hold(True)
new_tick_locations = numpy.array([math.atan(0.)*2/numpy.pi, math.atan(0.5)*2/numpy.pi, math.atan(1.)*2/numpy.pi, math.atan(2)*2/numpy.pi, math.atan(10.)*2/numpy.pi, math.atan(-0.5)*2/numpy.pi, math.atan(-1.)*2/numpy.pi, math.atan(-2.)*2/numpy.pi, math.atan(-10.)*2/numpy.pi])

ax1.set_xticks(new_tick_locations)
ax1.set_xticklabels(tick_function_half(new_tick_locations))

pylab.hold(False)
pylab.yscale('log')
pylab.ylabel(r"$P(\Delta_{\rm{merge}})$")
pylab.xlabel(r"$\Delta_{\rm{merge}}$")
pylab.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=4, prop={'size':8})
pylab.savefig(pdf_folder+'Lostmass_dlog.pdf',bbox_inches='tight')
##############################################################################################
#plot mean beta HIMASS
pre_ext = "_1000"
h_meanbeta_ext = pre_ext+"MeanBeta.dat"
print "Start plotting",h_meanbeta_ext
fig = pylab.figure()
ax1 = fig.add_subplot(111)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + h_meanbeta_ext)
    ax1.plot(data[:,0],data[:,1]+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')
#pylab.xlim([-.5,.5])
new_tick_locations = numpy.array([math.atan(0.)*2/numpy.pi, math.atan(0.5)*2/numpy.pi, math.atan(1.)*2/numpy.pi, math.atan(2)*2/numpy.pi, math.atan(10.)*2/numpy.pi, math.atan(-0.5)*2/numpy.pi, math.atan(-1.)*2/numpy.pi, math.atan(-2.)*2/numpy.pi, math.atan(-10.)*2/numpy.pi])

ax1.set_xticks(new_tick_locations)
ax1.set_xticklabels(tick_function_half(new_tick_locations))
pylab.xlabel(r"$\Delta_M$")
pylab.ylabel(r"$N+1$")
leg = pylab.legend(loc='lower center',handlelength = 8, ncol=2, prop={'size':10})
leg.get_frame().set_linewidth(0)
pylab.savefig(pdf_folder+'betahigh_hist.pdf',bbox_inches='tight')
os.system('pdftops -eps '+pdf_folder+'betahigh_hist.pdf')

############################################################################################################


##############################################################################################
#plot mean beta HIMASS
pre_ext = "_1000"
h_meanbeta_ext = pre_ext+"MeanSlope.dat"
print "Start plotting",h_meanbeta_ext
fig = pylab.figure()
ax1 = fig.add_subplot(111)
for (i,item) in enumerate(line):
    data = pylab.loadtxt(dat_folder + item + h_meanbeta_ext)
    ax1.plot(data[:,0]*2,data[:,1]+1,linestyles[i%len(linestyles)],color=colors[i%len(colors)],label = item.replace("_"," ") )
    pylab.hold(True)

pylab.hold(False)
pylab.yscale('log')

pylab.xlabel(r"$\delta M/ M$")
pylab.ylabel(r"$N+1$")
leg= pylab.legend(loc='upper left', handlelength = 6,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)
pylab.savefig(pdf_folder+'slope_hist.pdf', bbox_inches='tight')

############################################################################################################
#os.system("rm -f "+pdf_folder+'*.pdf')
