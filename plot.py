import pylab 
import sys
import numpy
import math
import os

#import sets

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
    halo4_color = 'b'
    colour = {}
    colour[0] = 'r'
    colour[1] = 'b'
    colour[2] = 'g'
    lw = 2
    
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

def make_peculiar_panels(save):
    runlist = load_group_list()
    n = 6
    zoomsize = 600.
    pos = [5600.59423828125, 28951.79296675, 25449.57421875]
    fig = {}
    fig0 = pylab.figure()
    ax0 = fig0.add_subplot(111)
    fig1 = pylab.figure()
    ID = {}
    ID[0] = 61000000000766
    ID[1] = 61000000000123
    HBTID0 = 61000000000888
    HBTID1 = 61000000000887
    haloid = {}

    colors = ('b','g','r','c','m','#00FF00','k','#FFA500','#7F7F7F','#800000')
    for i in range(len(runlist)):
        print runlist[i]
        data = load(runlist[i])
        haloid[0] = UIDtoID(ID[0],data)
        haloid[1] = UIDtoID(ID[1],data)
        if(runlist[i] == 'HBT.v3'):
            haloid[0] = UIDtoID(HBTID0,data)
            haloid[1] = UIDtoID(HBTID1,data)

        hid = haloid[1]
        mass = []
        time = []
        mass.append(data.halo[hid][2])
        time.append(data.time[data.halo[hid][1]])
    
        while len(data.trees[hid]) > 1:
            hid = int(data.trees[hid][1])
            mass.append(data.halo[hid][2])
            time.append(data.time[data.halo[hid][1]])
        ax0.plot(time,mass,color=colors[i],label = runlist[i