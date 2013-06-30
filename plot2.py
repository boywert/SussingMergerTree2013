import pylab 
import sys
import numpy
import math
import os

#import sets
from operator import itemgetter

sys.path.append('./lib/')
import gadgetPyIO 
import gadgetPydensity

global boxsize
global config
global cellsize
global datfolder
global runlist_file
global plot_config

runlist_file = 'run.list'
datfolder = 'dat/'
cellsize = 20.
boxsize = 62500.
f = open("config.txt")
config =f.read().splitlines()
f = None



class plotconfig(object):
    particle_color = '#505050'
    halo0_color = 'r'
    halo1_color = 'g'
    halo3_color = 'k'
    lw = 2
    halo_colours = ['r','b','g']
    colours =  ('b','g','r','c','m','#00FF00','k','#FFA500','#7F7F7F','#800000')

plot_config = plotconfig()

class output(object):
    halo = None
    trees = None
    r_misslist = None
    fluct_m_misslist = None
    zero_link = None
    lost_bin = None
    groupname = None
    mergers_bin = None
    beta_low = None
    time = None
    branch_length = None
    descendant = None
def load_group_list():
    f = open(runlist_file)
    line = f.read().splitlines()
    out = []
    for (i,item) in enumerate(line):
        filename = item.split()
        out.append(filename[0])
    out.sort()
    return out

def generate_descendant(tree):
    desc = numpy.zeros(len(tree))
    desc[:] = -1
    for i in range(len(tree)):
        for j in range(1,len(tree[i])):
            hid = int(tree[i][j])
            desc[hid] = i
    return desc

def make_plotdumbell():
    data = load('MergerTree')
    data.groupname = 'Flip-Flop'
    hid = {}
    hid[0] = 0
    plotonepanel_id_fixpos(hid, 54, data, True, 650, [46301.99609375, 26823.99609375, 8579.7548828125])
    plotonepanel_id_fixpos(hid, 55, data, True, 650, [46301.99609375, 26823.99609375, 8579.7548828125])
    plotonepanel_id_fixpos(hid, 56, data, True, 650, [46301.99609375, 26823.99609375, 8579.7548828125])

def plot_different_history():
    hid = 61000000000123
    ctid = 61000000000123
    hbtid = 61000000000887
    runlist = load_group_list()
    pylab.rc('text', usetex=True)
    fig = pylab.figure()
    ax0 = fig.add_subplot(211)
    for i in range(len(runlist)):
        print runlist[i]
        data = load(runlist[i])
        
        if(runlist[i] == 'HBT.v3'):
            haloid = UIDtoID(hbtid,data)
        elif(runlist[i] == 'CONSITENT_TREES'):
            haloid = UIDtoID(ctid,data)
        else:
            haloid = UIDtoID(hid,data)
        if(haloid != False):
            mass = []
            time = []
            cor = []
            masshalo = data.halo[haloid][2]
            snaphalo = data.halo[haloid][1]         
            mass.append(data.halo[haloid][2]/1.e12)
            time.append(data.time[data.halo[haloid][1]][4]/1e14)

            cor.append([data.time[snaphalo][4]/1.e10,masshalo/1.e12])

            while len(data.trees[haloid]) > 1:
                haloid = int(data.trees[haloid][1])
                masshalo = data.halo[haloid][2]
                snaphalo = data.halo[haloid][1]
                mass.append(masshalo/1.e12)
                time.append(data.time[snaphalo][4]/1.e10)
                cor.append([data.time[snaphalo][4]/1.e10,masshalo/1.e12])
            if(snaphalo > 0):
                mass.append(0.)
                time.append(data.time[snaphalo-1][4]/1.e10)
                cor.append([data.time[snaphalo-1][4]/1.e10,0.])

            print cor
            cor = sorted(cor, key=itemgetter(0))
            print cor
            time = [cor[ix][0] for ix in range(len(cor))]
            mass = [cor[ix][1] for ix in range(len(cor))]
            
            ax0.plot(time,mass,color=plot_config.colours[i],label = runlist[i].replace("_"," "))
            pylab.hold(True)
            mass = None
            time = None
            cor = None
    pylab.hold(False)
    leg = ax0.legend(loc='upper left', handlelength = 6,ncol=2, prop={'size':9})
    leg.get_frame().set_linewidth(0)
    #ax0.set_ylabel(r"$\mathrm{M_{vir}\, (10^{12}\, h^{-1}M_{\odot}})$")

    hid = 61000000000766
    ctid = 61000000000766
    hbtid = 61000000000888


    ax1 = fig.add_subplot(212)
    for i i