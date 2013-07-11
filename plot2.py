import matplotlib 
matplotlib.use('PDF')
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
    epsilon = 0.
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
        
        if(runlist[i] == 'HBT'):
            haloid = UIDtoID(hbtid,data)
        elif(runlist[i] == 'Consistent_Trees'):
            haloid = UIDtoID(ctid,data)
        else:
            haloid = UIDtoID(hid,data)
        print data.halo[haloid][0]
        if(haloid != False):
            mass = []
            time = []
            cor = []
            masshalo = data.halo[haloid][2]
            snaphalo = data.halo[haloid][1]         
            mass.append(data.halo[haloid][2]/1.e12)
            time.append(data.time[data.halo[haloid][1]][4]/1e10 +epsilon*i)

            cor.append([data.time[snaphalo][4]/1.e10,masshalo/1.e12])

            while len(data.trees[haloid]) > 1:
                haloid = int(data.trees[haloid][1])
                masshalo = data.halo[haloid][2]
                snaphalo = data.halo[haloid][1]
                mass.append(masshalo/1.e12)
                time.append(data.time[snaphalo][4]/1.e10+epsilon*i)
                cor.append([data.time[snaphalo][4]/1.e10,masshalo/1.e12])
            if(snaphalo > 0):
                mass.append(0.)
                time.append(data.time[snaphalo-1][4]/1.e10+epsilon*i)
                cor.append([data.time[snaphalo-1][4]/1.e10,0.])

            #print cor
            cor = sorted(cor, key=itemgetter(0))
            
            time = [cor[ix][0] for ix in range(len(cor))]
            mass = [cor[ix][1] for ix in range(len(cor))]
            #print time,mass
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
    for i in range(len(runlist)):
        print runlist[i]
        data = load(runlist[i])
        
        if(runlist[i] == 'HBT'):
            haloid = UIDtoID(hbtid,data)
        elif(runlist[i] == 'Consistent_Trees'):
            haloid = UIDtoID(ctid,data)
        else:
            haloid = UIDtoID(hid,data)
        print data.halo[haloid][0]
        if(haloid != False):
            mass = []
            time = []
            cor = []
            masshalo = data.halo[haloid][2]
            snaphalo = data.halo[haloid][1]         
            mass.append(data.halo[haloid][2]/1.e12)
            time.append(data.time[data.halo[haloid][1]][4]/1e10+epsilon*i)

            cor.append([data.time[snaphalo][4]/1.e10,masshalo/1.e12])

            while len(data.trees[haloid]) > 1:
                haloid = int(data.trees[haloid][1])
                masshalo = data.halo[haloid][2]
                snaphalo = data.halo[haloid][1]
                mass.append(masshalo/1.e12)
                time.append(data.time[snaphalo][4]/1.e10+epsilon*i)
                cor.append([data.time[snaphalo][4]/1.e10,masshalo/1.e12])
            if(snaphalo > 0):
                mass.append(0.)
                time.append(data.time[snaphalo-1][4]/1.e10+epsilon*i)
                cor.append([data.time[snaphalo-1][4]/1.e10,0.])

            #print cor
            cor = sorted(cor, key=itemgetter(0))
            #print cor
            time = [cor[ix][0] for ix in range(len(cor))]
            mass = [cor[ix][1] for ix in range(len(cor))]
            
            ax1.plot(time,mass,color=plot_config.colours[i],label = runlist[i].replace("_"," "))
            pylab.hold(True)
            mass = None
            time = None
            cor = None
    
    pylab.hold(False)
    pylab.subplots_adjust(hspace = 0.0001)
    ax0.set_ylim([-0.1,7])
    ax1.set_ylim([-0.1,7])
    ax0.text(0.1,3,'BLUE',color='b')
    ax1.text(0.1,3,'RED',color='r')
    ax0.set_xticklabels([])
    ax0.set_yticks(numpy.array([1,2,3,4,5,6,7]))
    ax1.set_yticks(numpy.array([0,1,2,3,4,5,6,7]))
    ax1.set_xlabel(r"$\mathrm{time\, (10^{10}\, years})$")
    yl = ax1.set_ylabel(r"$\mathrm{M_{vir}\, (10^{12}\, h^{-1}M_{\odot}})$")
    yl.set_verticalalignment('center')
    position = yl.get_position()
    #print position
    yl.set_position([0.0,1.0])
    pylab.savefig(str(hid)+'_history.pdf',bbox_inches='tight')
    os.system("pdftops -eps "+str(hid)+'_history.pdf')
    os.system("rm -f *.png *.pdf")
  
def plot_swap_mass_paper():
    hid = {}
    ctid = {}
    hbtid = {}
    hid[0] = 59000000000490
    ctid[0] = 59000000000490
    hbtid[0] = 59000000001701
    runlist = load_group_list()
    save = True
    pos = [49031.37109375, 30610.4042968, 56049.984375]
    zoomsize = 370.
    n = 3
    for i in range(len(runlist)):
        print runlist[i]
        data = load(runlist[i])
        
        if(runlist[i] == 'HBT'):
            haloid = dict(hbtid)
        elif(runlist[i] == 'Consistent_Trees'):
            haloid = dict(ctid)
        else:
            haloid = dict(hid)
        check = 1
        if(runlist[i] == 'MergerTree'):
            data.groupname = 'Particle_Identifiers'
        print haloid
        for j in haloid.keys():
            haloid[j] = UIDtoID(haloid[j],data)
            if(haloid[j] == False):
                check = 0
        if (check == 0) :
            return False
        else:
            plotserial_id_fixpos(haloid,n,data,save,zoomsize,pos)
    os.system("rm -f *.png *.pdf")

def makemultipanels_paper():
    hid = {}
    ctid = {}
    hbtid = {}
    hid[0] = 61000000000766
    hid[1] = 61000000000123
    ctid[0] = 61000000000766
    ctid[1] = 61000000000123   
    hbtid[0] = 61000000000888
    hbtid[1] = 61000000000887
    runlist = load_group_list()
    save = True
    pos = [5600.59423828125, 28951.79296675, 25449.57421875]
    zoomsize = 600.
    n=6
    for i in range(len(runlist)):
        print runlist[i]
        data = load(runlist[i])
        
        if(runlist[i] == 'HBT'):
            haloid = dict(hbtid)
        elif(runlist[i] == 'Consistent_Trees'):
            haloid = dict(ctid)
        else:
            haloid = dict(hid)
        check = 1
        print haloid
        for j in haloid.keys():
            haloid[j] = UIDtoID(haloid[j],data)
            if(haloid[j] == False):
                check = 0
        if (check == 0) :
            return False
        else:
            plotserial_id_fixpos(haloid,n,data,save,zoomsize,pos)
    os.system("rm -f *.png *.pdf")

def plot_history_compare(hid,ctid,hbtid):
    runlist = load_group_list()
    pylab.rc('text', usetex=True)
    fig = pylab.figure()
    ax0 = fig.add_subplot(111)
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
            mass = None
            time = None
            cor = None
    
    ax0.legend(loc='upper left', handlelength = 6,ncol=1, prop={'size':9})
    ax0.set_xlabel(r"$\mathrm{time\, (10^{10}\, years})$")
    ax0.set_ylabel(r"$\mathrm{M_{vir}\, (10^{12}\, h^{-1}M_{\odot}})$")
    pylab.savefig('sample_history.pdf',bbox_inches='tight')
    os.system("pdftops -eps sample_history.pdf")
    os.system("rm -f *.png *.pdf")

def plot_all_zero_link():
    runlist = load_group_list()
    for i in range(len(runlist)):
        data = load(runlist[i])
        N = len(data.zero_link)
        for j in range(N):
            plot_zero_link(j,data,True,1)

def compare_trees_table_all():
    runlist = load_group_list()
    all_data = []
    for i in range(len(runlist)):
        print 'loading data from',runlist[i]
        all_data.append(load(runlist[i]))
    table = [ [ 0 for i in range(len(runlist)) ] for j in range(len(runlist)) ]
    for i in range(len(runlist)):
        for j in range(len(runlist)):
            table[i][j] = compare_trees_table(all_data[i].trees,all_data[j].trees)
            string = ''
    string += "&"
    string += " & ".join(runlist)
    string += " \\\\ \n"
        
    for i in range(len(runlist)):
        string += runlist[i]
        string += "& "
        string_tab = []
        for j in range(len(table[i])):
            if (table[i][j] < 1.01):
                string_tab.append("%1.6f"%table[i][j] )
            else:
                string_tab.append('-')
        print string_tab
        string += " & ".join(string_tab)
        string += "\\\\ \n"
    all_data = None
    print string

def compare_progs_table_all():
    runlist = load_group_list()
    all_data = []
    for i in range(len(runlist)):
        print 'loading data from',runlist[i]
        all_data.append(load(runlist[i]))
    table = [ [ 0 for i in range(len(runlist)) ] for j in range(len(runlist)) ]
    for i in range(len(runlist)):
        for j in range(len(runlist)):
            table[i][j] = compare_progs_table(all_data[i],all_data[j])
    string = ''
    string += "&"
    string += " & ".join(runlist)
    string += " \\\\ \n"
        
    for i in range(len(runlist)):
        string += runlist[i]
        string += "& "
        string_tab = []
        for j in range(len(table[i])):
            if (table[i][j] < 1.01):
                string_tab.append("%1.6f"%table[i][j] )
            else:
                string_tab.append('-')
        print string_tab
        string += " & ".join(string_tab)
        string += "\\\\ \n"
    all_data = None
    print string

def compare_descs_table_all():
    runlist = load_group_list()
    all_data = []
    for i in range(len(runlist)):
        print 'loading data from',runlist[i]
        all_data.append(load(runlist[i]))
    table = [ [ 0 for i in range(len(runlist)) ] for j in range(len(runlist)) ]
    for i in range(len(runlist)):
        for j in range(len(runlist)):
            table[i][j] = compare_descs_table(all_data[i].descendant,all_data[j].descendant)
    string = ''
    string += "&"
    string += " & ".join(runlist)
    string += " \\\\ \n"
        
    for i in range(len(runlist)):
        string += runlist[i]
        string += "& "
        string_tab = []
        for j in range(len(table[i])):
            if (table[i][j] < 1.01):
                string_tab.append("%1.6f"%table[i][j] )
            else:
                string_tab.append('-')
        print string_tab
        string += " & ".join(string_tab)
        string += "\\\\ \n"
    all_data = None
    print string

def compare_trees_table(tree1,tree2):
    n1 = len(tree1)
    n2 = len(tree2)
    if(n1 != n2):
        print n1,'is not equal to',n2
        return 100.
    count = 0
    for i in range(n1):
        if ((len(tree1[i]) == 1) and (len(tree2[i]) == 1)):
            count += 1
        elif ((len(tree1[i]) > 1) and (len(tree2[i]) > 1)):
             if (int(tree1[i][1]) == int(tree2[i][1])):
                 count += 1
    return float(count)/float(n1)



def compare_descs_table(desc1,desc2):
    n1 = len(desc1)
    n2 = len(desc2)
    if(n1 != n2):
        print n1,'is not equal to',n2
        return 100.
    count = 0
    for i in range(n1):
        if (desc1[i] == desc2[i]):
            count += 1
    return float(count)/float(n1)




def compare_progs_table(data1,data2):
    tree1 = data1.trees
    tree2 = data2.trees
    n1 = len(tree1)
    n2 = len(tree2)
    if(n1 != n2):
        print n1,'is not equal to',n2
        return 100.
    cond = ((data1.halo['f1'] == 61) & (data1.halo['f2'] > 1.e12))
    halolist = numpy.where(cond)
    halolist = halolist[0]
    count = 0
    for i in range(len(halolist)):
        cur_id = halolist[i]
        if (len(tree1[cur_id]) == 1 and len(tree2[cur_id]) == 1):
            count += 1
        elif (len(tree1[cur_id]) > 1 and len(tree2[cur_id]) > 1):
             if int(tree1[cur_id][1]) == int(tree2[cur_id][1]):
                 count += 1
    return float(count)/len(halolist)

def load(groupname):
    data = output()
    data.groupname = groupname
    #data.halo = pylab.loadtxt(config[3]+"/"+groupname+".profile")
    dt = numpy.dtype('u8,u4,f4,f4,f4,f4,f4,f4,f4,f4')
    data.halo = numpy.fromfile(config[3]+"/"+groupname+".profile.bin", dtype=dt)
    data.halo = numpy.array(data.halo)
    
    data.time = pylab.loadtxt(config[2])
    f = open(datfolder+groupname+"_1000DispMainHost_pp.dat")
    data.r_misslist = f.read().splitlines()
    f = None
    for i in range(len(data.r_misslist)):
        temp = data.r_misslist[i].split()
        data.r_misslist[i] = temp
        temp = None

    f = open(datfolder+groupname+"_1000Beta_pp.dat")
    data.beta_low = f.read().splitlines()
    f = None
    for i in range(len(data.beta_low)):
        temp = data.beta_low[i].split()
        data.beta_low[i] = temp
        temp = None
    

    f = open(datfolder+groupname+"_1000Linkdepth_pp.dat")
    data.zero_link = f.read().splitlines()
    f = None
    for i in range(len(data.zero_link)):
        temp = data.zero_link[i].split()
        data.zero_link[i] = temp
        temp = None
    if(os.path.exists(datfolder+groupname+"_1000Lostmass_bin.dat") >= 1):
        f = open(datfolder+groupname+"_1000Lostmass_bin.dat")
        data.lost_bin = f.read().splitlines()
        f = None
        for i in range(len(data.lost_bin)):
            temp = data.lost_bin[i].split()
            data.lost_bin[i] = temp
            temp = None

    f = open(datfolder+groupname+"_1000Nmergers_bin.dat")
    data.mergers_bin = f.read().splitlines()
    f = None
    for i in range(len(data.mergers_bin)):
        temp = data.mergers_bin[i].split()
        data.mergers_bin[i] = temp
        temp = None

    f = open(config[3]+"/"+groupname+".trees")
    data.trees = f.read().splitlines()
    for i in range(len(data.trees)):
        temp = data.trees[i].split()
        data.trees[i] = temp
        if(i != int(data.trees[i][0])):
            print i," != ",data.trees[i][0]
        temp = None
    f = None
    gen_linklength(data)
    data.descendant = generate_descendant(data.trees)
    return data
def gen_linklength(data):
    data.branch_length = numpy.zeros(len(data.trees))
    cond = (data.halo['f1'] == 61)
    halo = numpy.where(cond)
    halo = halo[0]
    cond = None
    for i in range(len(halo)):
        hid = halo[i]
        n=0
        while len(data.trees[hid]) > 1:
            n += 1
            hid = int(data.trees[hid][1])
        data.branch_length[halo[i]] = n
def plot_miss_r(n,data,save):
    item = int(data.r_misslist[n][0])
    print 'data',data.r_misslist[n]
    r = float(data.r_misslist[n][1])
    pos = numpy.array([float(data.r_misslist[n][2]),float(data.r_misslist[n][3]),float(data.r_misslist[n][4])])
    plotmergers(item,data,r,pos,True,save,1.)

def plot_miss_beta(n,data,save,zoom):
    item = int(data.beta_low[n][1])
    snap1 =  data.halo[item][1]
    snap0 = data.halo[int(data.trees[item][1])][1]
    m1 = data.halo[item][2]
    m0 = data.halo[int(data.trees[item][1])][2]
    t1 = data.time[snap1][4]
    t0 = data.time[snap0][4]
    beta = (m1-m0)/(m1+m0)/(t1-t0)*(t1+t0)
    #beta = math.atan(beta) / math.asin(1.)
    print data.halo[item][2], data.halo[int(data.trees[item][1])][2], beta
    print 'data',data.beta_low[n]
    graph_mass_history(int(data.beta_low[n][0]), data,False)
    plotmergers(item,data,False,[0],False,save,zoom)

def graph_mass_history(haloid,data,save):
    hid = haloid
    mass = []
    time = []
    mass.append(data.halo[hid][2])
    time.append(data.time[data.halo[hid][1]])
    
    while len(data.trees[hid]) > 1:
        hid = int(data.trees[hid][1])
        mass.append(data.halo[hid][2])
        time.append(data.time[data.halo[hid][1]])
        #plot_general(hid,data,False,1)
    
    fig = pylab.figure()
    pylab.plot(time,mass)
    fig.show()
    
def plot_zero_link(n,data,save,zoom):
    item = int(data.zero_link[n][0])
    print 'data',data.zero_link[n]
    plotsingle(item,data,save,1.)

    

def plot_general(ihalo,data,save,zoom):
    if(ihalo > len(data.trees)):
        print ihalo,'is invalid'
        return
    if(len(data.trees[ihalo]) < 2):
        print ihalo,'does not have a progenitor'
        return
    plotmergers(ihalo,data,False,[0],False,save,zoom)



def UIDtoID(UID,data):
    list = numpy.where(data.halo['f0'] == UID)
    if(len(list[0]) > 0):
        return list[0][0]
    else:
        return False
    
def plotmergers(n,data,Radius,expected_pos,flag,save,zoom):

    item = data.trees[n]

    #convert item to integer
    for i in range(len(item)):
        item[i] = int(item[i])

    #change coordinates
    if(len(item) < 2):
        print 'No progenitor recorded'
        return False

    if (flag == False):
        if(len(item) == 2):
            r2 = numpy.array([0.,0.,0])
        else:
            h2= int(item[2])
            x2 = data.halo[h2][4]
            y2 = data.halo[h2][5]
            z2 = data.halo[h2][6]
            r2 = numpy.array([x2,y2,z2])
    else:
        r2 = expected_pos

    r2 = expected_pos
    h0= int(item[1])
    x0 = data.halo[h0][4]
    y0 = data.halo[h0][5]
    z0 = data.halo[h0][6]
    r0 = numpy.array([x0,y0,z0])


    
    x_axis = numpy.array([])
    y_axis = numpy.array([])
    z_axix = numpy.array([])
   
    h1 = int(item[0])
    x1 = data.halo[h1][4]
    y1 = data.halo[h1][5]
    z1 = data.halo[h1][6]
    r1 = numpy.array([x1,y1,z1])

    x_axis = r1 - r0
    x_axis = numpy.array([1.,0.,0.])
    y_axis = numpy.array([0.,1.,0.])
    z_axis = numpy.array([0.,0.,1.])
    for j in range(3):
        if(x_axis[j] > boxsize/2.):
            x_axis[j] -= boxsize
        if(x_axis[j] < -1.*boxsize/2.):
            x_axis[j] += boxsize
    x_axis = x_axis/(numpy.sqrt(numpy.dot(x_axis,x_axis)))
    y_axis_temp = r2 - r0
    for j in range(3):
        if(y_axis_temp[j] > boxsize/2.):
            y_axis_temp[j] -= boxsize
        if(y_axis_temp[j] < -1.*boxsize/2.):
            y_axis_temp[j] += boxsize
    y_axis_temp = y_axis_temp/(numpy.sqrt(numpy.dot(y_axis_temp,y_axis_temp)))
    z_axis = numpy.cross(x_axis,y_axis_temp)
    y_axis = numpy.cross(z_axis,x_axis)
    

    xy_plane = []
    print 'len',len(item)
    print 'r0', r0
    for (i,string) in enumerate(item):
        s = int(item[i])
        xs = data.halo[s][4]
        ys = data.halo[s][5]
        zs = data.halo[s][6]
        rs = [xs, ys, zs]
        rs = rs - r0
        print i,'dot product', numpy.dot(rs,y_axis)
        #print data.halo[s,:]
        print data.halo[s][0]
        #print 'before',rs
        for j in range(len(rs)):
            if(rs[j] > boxsize/2.):
                rs[j] -= boxsize
            if(rs[j] < -1.*boxsize/2.):
                rs[j] += boxsize
        #print 'after', rs
        l = numpy.dot(rs,z_axis)
        
        r_projected = rs - z_axis*l
        xy = [0.,0.]
        
        xy[0] = numpy.dot(r_projected,x_axis)
        xy[1] = numpy.dot(r_projected,y_axis)
        print i,xy[0],xy[1],l
        xy_plane.append(xy)


    #expected position
    rs = expected_pos - r0
    if(rs[j] > boxsize/2.):
        rs[j] -= boxsize
    if(rs[j] < -1.*boxsize/2.):
        rs[j] += boxsize

    l = numpy.dot(rs,z_axis)
        
    r_projected = rs - z_axis*l
    xy = [0.,0.]
    xy[0] = numpy.dot(r_projected,x_axis)
    xy[1] = numpy.dot(r_projected,y_axis)
    expect = [xy[0],xy[1]]
    
    ref = [xy_plane[1][0],xy_plane[1][1]]
    print "ref point", ref
    
    fig_test = pylab.figure()
    axes = fig_test.add_subplot(211)

    for (i,string) in enumerate(item):
        s = int(item[i])
        xx = xy_plane[i][0] + ref[0]
        yy = xy_plane[i][1] + ref[1]
 
        rr = data.halo[s][3]
        
        if (i==0):
            circle = pylab.Circle((xx,yy),fill=False,ls='dashed', radius=rr, alpha=.5, color='r')
            print ''
            axes.add_patch(circle)
        elif (i==1):
            circle = pylab.Circle((xx,yy),fill=False, radius=rr, alpha=.5, color=plot_config.halo0_color)
            axes.add_patch(circle)
        else:
            circle = pylab.Circle((xx,yy),fill=False, radius=rr, alpha=.5, color=plot_config.halo1_color)
            axes.add_patch(circle)

    if(flag != False):
        circle = pylab.Circle((ref[0],ref[1]),fill=False,radius=Radius,ls='dashdot', alpha=.5, color='g')
        axes.add_patch(circle)
        circle = pylab.Circle((expect[0],expect[1]),fill=False,radius=data.halo[int(item[1])][3],ls='dotted', alpha=.5, color='c')
        axes.add_patch(circle)
    #pylab.xlim([xmin,xmax])
    #pylab.ylim([ymin,ymax])
    v = pylab.axis('equal')
    #v =  pylab.axis()
    #pylab.set_aspect('equal', 'datalim')

    pylab.axis(v)
    xmin = v[0]*zoom
    xmax = v[1]*zoom
    ymin = v[2]*zoom
    ymax = v[3]*zoom

    fig = pylab.figure()
    axes = fig.add_subplot(211)

    for (i,string) in enumerate(item):
        s = int(item[i])
        xx = xy_plane[i][0] + ref[0]
        yy = xy_plane[i][1] + ref[1]
 
        rr = data.halo[s][3]
        
        if i==0:
            #circle = pylab.Circle((xx,yy),fill=False,ls='dashed', radius=rr, alpha=.5, color='r',lw=2)
            print ''
        elif i==1:
            circle = pylab.Circle((xx,yy),fill=False, radius=rr, alpha=.5, color=plot_config.halo0_color,lw=2)
            axes.add_patch(circle)
        else:
            circle = pylab.Circle((xx,yy),fill=False, radius=rr, alpha=.5, color=plot_config.halo1_color,lw=2)
            axes.add_patch(circle)

    if(flag != False):
        #circle = pylab.Circle((ref[0],ref[1]),fill=False,radius=Radius,ls='dashdot', alpha=.5, color='g',lw=2)
        #axes.add_patch(circle)
        circle = pylab.Circle((expect[0],expect[1]),fill=False,radius=data.halo[int(item[1])][3],ls='dashdot', alpha=.5, color='c',lw=2)
        axes.add_patch(circle)
    #pylab.xlim([xmin,xmax])
    #pylab.ylim([ymin,ymax])


    #xmin = -1000.
    print 'region',xmin,xmax,ymin,ymax
    snap1 = int(data.halo[int(item[1])][1])
    particle1 = getdensity2D(snap1,r0,x_axis,y_axis,z_axis,xmin,xmax,ymin,ymax)
    axes.scatter(particle1[:,0],particle1[:,1],edgecolors='none',s = 1, color=plot_config.particle_color)


    extend_radius = numpy.sqrt((xmax-xmin)**2 + (ymax-ymin)**2)

    cond = (data.halo['f1'] == snap1)
    condu = cond
    condl = cond
  # make a box to contain local halos
    for i in range(3):
        if(i==0):
            col = 'f4'
        elif(i==1):
            col = 'f5'
        else:
            col = 'f6'
        print i,col
        lower = r0[i] - extend_radius
        condu = cond
        condl = cond
        if(lower <= 0.):
            condl = condl & ((data.halo[col] >= (lower + boxsize)) | (data.halo[col] <= r0[i]))
        else:
            condl = condl & (data.halo[col] >= lower) 
            condl = condl & (data.halo[col] <= r0[i])
        upper = r0[i] + extend_radius
        if(upper >= boxsize):
            condu = condu & ((data.halo[col] <= (upper - boxsize)) | (data.halo[col] >= r0[i]))
        else:
            condu = condu & (data.halo[col] <= upper) 
            condu = condu & (data.halo[col] >= r0[i])
        cond = cond & (condl | condu)
    
    halolist = numpy.where(cond)
    print halolist
    cond = None
    condu = None
    condl = None
    #loop to plot local halo if not in the merger list
    for i in range(len(halolist[0])):
        index = halolist[0][i]
        if(index not in item):
            xs = data.halo[index][4]
            ys = data.halo[index][5]
            zs = data.halo[index][6]
            radius_s = data.halo[index][3]
            rs = numpy.array([xs,ys,zs])
            rs = rs - r0
            for j in range(3):
                if(rs[j] > boxsize/2.):
                    rs[j] -= boxsize
                if(rs[j] < -1.*boxsize/2.):
                    rs[j] += boxsize
            l = numpy.dot(rs,z_axis)
            r_projected = rs - z_axis*l
            xx = numpy.dot(r_projected,x_axis)
            yy = numpy.dot(r_projected,y_axis)
            if(xx > xmin and xx < xmax and yy > ymin and yy < ymax):
                print 'ADD',index,data.halo[index][0],xx,yy,radius_s
                circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color=plot_config.halo3_color)
                axes.add_patch(circle)
            else:
                print 'OUT OF RANGE',index,data.halo[index][0],xx,yy,radius_s

  
 


    axes2 = fig.add_subplot(212)   

    for (i,string) in enumerate(item):
        s = int(item[i])
        xx = xy_plane[i][0] + ref[0]
        yy = xy_plane[i][1] + ref[1]
 
        rr = data.halo[s][3]
        
        if i==0:
            circle = pylab.Circle((xx,yy),fill=False, radius=rr, alpha=.5, color=plot_config.halo0_color,lw=2)
            axes2.add_patch(circle)
        elif i==1:
            #circle = pylab.Circle((xx,yy),fill=False, ls='dashed',radius=rr, alpha=.5, color='m',lw=2)
            print ''
        else:
            print ''
            #circle = pylab.Circle((xx,yy),fill=False, ls='dashed',radius=rr, alpha=.5, color='g',lw=2)



    if(flag != False):
        #circle = pylab.Circle((ref[0],ref[1]),fill=False,radius=Radius,ls='dashdot', alpha=.5, color='g', lw=2)
        #axes2.add_patch(circle)
        circle = pylab.Circle((expect[0],expect[1]),fill=False,radius=data.halo[int(item[1])][3],ls='dotted', alpha=.5, color='c',lw=2)
        axes2.add_patch(circle)

    #axes.imshow(particle,interpolation='bilinear', extent=(xmin,xmax,ymin,ymax), cmap=pylab.get_cmap("binary"),origin='lower')
    snap2 = int(data.halo[int(item[0])][1])
    particle2 = getdensity2D(snap2,r0,x_axis,y_axis,z_axis,xmin,xmax,ymin,ymax)
    axes2.scatter(particle2[:,0],particle2[:,1],edgecolors='none',s = 1, color=plot_config.particle_color)
    #pylab.xlim([xmin,xmax])
    #pylab.ylim([ymin,ymax])
    #pylab.axis('equal')
    #pylab.axis(v)
    #pylab.set_aspect('equal', 'datalim')



    cond = (data.halo['f1'] == snap2)
    condu = cond
    condl = cond

  # make a box to contain local halos
    for i in range(3):
        if(i==0):
            col = 'f4'
        elif(i==1):
            col = 'f5'
        else:
            col = 'f6'
        print i,col
        lower = r0[i] - extend_radius
        condu = cond
        condl = cond
        if(lower <= 0.):
            condl = condl & ((data.halo[col] >= (lower + boxsize)) | (data.halo[col] <= r0[i]))
        else:
            condl = condl & (data.halo[col] >= lower) 
            condl = condl & (data.halo[col] <= r0[i])
        upper = r0[i] + extend_radius
        if(upper >= boxsize):
            condu = condu & ((data.halo[col] <= (upper - boxsize)) | (data.halo[col] >= r0[i]))
        else:
            condu = condu & (data.halo[col] <= upper) 
            condu = condu & (data.halo[col] >= r0[i])
        cond = cond & (condl | condu)

    halolist2 = numpy.where(cond)
    print halolist2
    cond = None
    condu = None
    condl = None
    #loop to plot local halo if not in the merger list
    print item
    for i in range(len(halolist2[0])):
        index = halolist2[0][i]
        print index
        if(index not in item):
            xs = data.halo[index][4]
            ys = data.halo[index][5]
            zs = data.halo[index][6]
            radius_s = data.halo[index][3]
            rs = numpy.array([xs,ys,zs])
            rs = rs - r0
            for j in range(3):
                if(rs[j] > boxsize/2.):
                    rs[j] -= boxsize
                if(rs[j] < -1.*boxsize/2.):
                    rs[j] += boxsize
            l = numpy.dot(rs,z_axis)
            r_projected = rs - z_axis*l
            xx = numpy.dot(r_projected,x_axis)
            yy = numpy.dot(r_projected,y_axis)
            if(xx > xmin and xx < xmax and yy > ymin and yy < ymax):
                print 'ADD',index,data.halo[index][0],xx,yy,radius_s
                circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color=plot_config.halo3_color)
                axes2.add_patch(circle)
            else:
                print 'OUT OF RANGE',index,data.halo[index][0],xx,yy,radius_s



    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    axes2.set_xlim([xmin,xmax])
    axes2.set_ylim([ymin,ymax])
    axes2.set_aspect(1.)
    axes.set_aspect(1.)
    axes2.text(xmin*(1-0.1), ymin, data.groupname+' ')
    pylab.subplots_adjust(hspace = 0.0001)

    
    axes.set_xticklabels([])
    axes.set_yticklabels([])
    axes2.set_xticklabels([])
    axes2.set_yticklabels([])

    if save is True:
        fig.savefig(data.groupname+'_'+str(data.halo[item[0]][0])+'.png',bbox_inches='tight')
    else:
        fig.show()

def plotsingle(haloid,data,save,zoom):

    item = []
    print data.halo[haloid]

    x0 = data.halo[haloid][4]
    y0 = data.halo[haloid][5]
    z0 = data.halo[haloid][6]

    r0 = numpy.array([x0,y0,z0])
    r1 = numpy.array([x0+100.,y0,z0])
    r2 = numpy.array([x0,y0+100.,z0])

    x_axis = numpy.array([])
    y_axis = numpy.array([])
    z_axix = numpy.array([])
   
    x_axis = r1 - r0
    for j in range(3):
        if(x_axis[j] > boxsize/2.):
            x_axis[j] -= boxsize
        if(x_axis[j] < -1.*boxsize/2.):
            x_axis[j] += boxsize
    x_axis = x_axis/(numpy.sqrt(numpy.dot(x_axis,x_axis)))
    y_axis_temp = r2 - r0
    for j in range(3):
        if(y_axis_temp[j] > boxsize/2.):
            y_axis_temp[j] -= boxsize
        if(y_axis_temp[j] < -1.*boxsize/2.):
            y_axis_temp[j] += boxsize
    y_axis_temp = y_axis_temp/(numpy.sqrt(numpy.dot(y_axis_temp,y_axis_temp)))
    z_axis = numpy.cross(x_axis,y_axis_temp)
    y_axis = numpy.cross(z_axis,x_axis)
    

    plotsize = data.halo[haloid][3]*3*zoom


    xmin = plotsize*-1.
    xmax = plotsize
    ymin = plotsize*-1.
    ymax = plotsize
    extend_radius = numpy.sqrt((ymax-ymin)**2 + (xmax-xmin)**2)
  
  #xmin = -1000.
    print 'region',xmin,xmax,ymin,ymax


    fig = pylab.figure()
    axes = fig.add_subplot(212)

    snap1 = data.halo[haloid][1]
    snap2 = snap1-1
    particle1 = getdensity2D(snap1,r0,x_axis,y_axis,z_axis,xmin,xmax,ymin,ymax)
    axes.scatter(particle1[:,0],particle1[:,1],edgecolors='none',s = 1, color=plot_config.particle_color)

    #extend_radius = numpy.sqrt((xmax-xmin)**2 + (ymax-ymin)**2)
    cond = (data.halo['f1'] == snap1)
    condu = cond
    condl = cond

  # make a box to contain local halos
    for i in range(3):
        if(i==0):
            col = 'f4'
        elif(i==1):
            col = 'f5'
        else:
            col = 'f6'
        print i,col
        lower = r0[i] - extend_radius
        condu = cond
        condl = cond
        if(lower <= 0.):
            condl = condl & ((data.halo[col] >= (lower + boxsize)) | (data.halo[col] <= r0[i]))
        else:
            condl = condl & (data.halo[col] >= lower) 
            condl = condl & (data.halo[col] <= r0[i])
        upper = r0[i] + extend_radius
        if(upper >= boxsize):
            condu = condu & ((data.halo[col] <= (upper - boxsize)) | (data.halo[col] >= r0[i]))
        else:
            condu = condu & (data.halo[col] <= upper) 
            condu = condu & (data.halo[col] >= r0[i])
        cond = cond & (condl | condu)
        
    
    halolist = numpy.where(cond)

    print halolist
    cond = None
    condu = None
    condl = None
    #loop to plot local halo if not in the merger list
    for i in range(len(halolist[0])):
        index = halolist[0][i]
        if(1 == 1):
            xs = data.halo[index][4]
            ys = data.halo[index][5]
            zs = data.halo[index][6]
            radius_s = data.halo[index][3]
            rs = numpy.array([xs,ys,zs])
            rs = rs - r0
            for j in range(3):
                if(rs[j] > boxsize/2.):
                    rs[j] -= boxsize
                if(rs[j] < -1.*boxsize/2.):
                    rs[j] += boxsize
            l = numpy.dot(rs,z_axis)
            r_projected = rs - z_axis*l
            xx = numpy.dot(r_projected,x_axis)
            yy = numpy.dot(r_projected,y_axis)
            if(xx > xmin and xx < xmax and yy > ymin and yy < ymax):
                print 'ADD',index,data.halo[index][0],xx,yy,radius_s
                if(index == haloid):
                    circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color=plot_config.halo0_color, lw =2)
                else:
                    circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color=plot_config.halo3_color, lw = 2)
                axes.add_patch(circle)
            else:
                print 'OUT OF RANGE',index,data.halo[index][0],xx,yy,radius_s

  
 
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    axes.set_aspect(1.)
    
   
    axes2 = fig.add_subplot(211)   

    particle2 = getdensity2D(snap2,r0,x_axis,y_axis,z_axis,xmin,xmax,ymin,ymax)
    axes2.scatter(particle2[:,0],particle2[:,1],edgecolors='none',s = 1, color=plot_config.particle_color)

    #pylab.axis('equal')
    #pylab.axis(v)
    #pylab.set_aspect('equal', 'datalim')
    snap2 = snap1-1
    cond = (data.halo['f1'] == snap2)
    condu = cond
    condl = cond

  # make a box to contain local halos
    for i in range(3):
        if(i==0):
            col = 'f4'
        elif(i==1):
            col = 'f5'
        else:
            col = 'f6'
        print i,col
        lower = r0[i] - extend_radius
        condu = cond
        condl = cond
        if(lower <= 0.):
            condl = condl & ((data.halo[col] >= (lower + boxsize)) | (data.halo[col] <= r0[i]))
        else:
            condl = condl & (data.halo[col] >= lower)  & (data.halo[col] <= r0[i])
        upper = r0[i] + extend_radius

        if(upper >= boxsize):
            condu = condu & ((data.halo[col] <= (upper - boxsize)) | (data.halo[col] >= r0[i]))
        else:
            condu = condu & (data.halo[col] <= upper) & (data.halo[col] >= r0[i])
        cond = cond & (condl | condu)
        halolist2 = numpy.where(cond)


    halolist2 = numpy.where(cond)
    print halolist2
    cond = None
    condu = None
    condl = None
    #loop to plot local halo if not in the merger list
    for i in range(len(halolist2[0])):
        index = halolist2[0][i]
        print index
        if(1==1):
            xs = data.halo[index][4]
            ys = data.halo[index][5]
            zs = data.halo[index][6]
            radius_s = data.halo[index][3]
            rs = numpy.array([xs,ys,zs])
            rs = rs - r0
            for j in range(3):
                if(rs[j] > boxsize/2.):
                    rs[j] -= boxsize
                if(rs[j] < -1.*boxsize/2.):
                    rs[j] += boxsize
            l = numpy.dot(rs,z_axis)
            r_projected = rs - z_axis*l
            xx = numpy.dot(r_projected,x_axis)
            yy = numpy.dot(r_projected,y_axis)
            if(xx > xmin and xx < xmax and yy > ymin and yy < ymax):
                print 'ADD',index,data.halo[index][0],xx,yy,radius_s
                circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color='k', lw=2)
                axes2.add_patch(circle)
            else:
                print 'OUT OF RANGE',index,data.halo[index][0],xx,yy,radius_s


    axes2.set_xlim([xmin,xmax])
    axes2.set_ylim([ymin,ymax])
    axes2.set_aspect(1.)
    pylab.subplots_adjust(hspace = 0.0001)

    
    axes.set_xticklabels([])
    axes.set_yticklabels([])
    axes2.set_xticklabels([])
    axes2.set_yticklabels([])
    axes.text(xmin*(1-0.1), ymin, data.groupname)
    
    if save is True:
        fig.savefig(data.groupname+'_'+str(data.halo[haloid][0])+'.png',bbox_inches='tight')
    else:
        fig.show()

def plotonepanel(snap,pos,data,save,zoom):

    item = []

    x0 = pos[0]
    y0 = pos[1]
    z0 = pos[2]

    r0 = numpy.array([x0,y0,z0])
    r1 = numpy.array([x0+100.,y0,z0])
    r2 = numpy.array([x0,y0+100.,z0])

    x_axis = numpy.array([])
    y_axis = numpy.array([])
    z_axix = numpy.array([])
   
    x_axis = r1 - r0
    for j in range(3):
        if(x_axis[j] > boxsize/2.):
            x_axis[j] -= boxsize
        if(x_axis[j] < -1.*boxsize/2.):
            x_axis[j] += boxsize
    x_axis = x_axis/(numpy.sqrt(numpy.dot(x_axis,x_axis)))
    y_axis_temp = r2 - r0
    for j in range(3):
        if(y_axis_temp[j] > boxsize/2.):
            y_axis_temp[j] -= boxsize
        if(y_axis_temp[j] < -1.*boxsize/2.):
            y_axis_temp[j] += boxsize
    y_axis_temp = y_axis_temp/(numpy.sqrt(numpy.dot(y_axis_temp,y_axis_temp)))
    z_axis = numpy.cross(x_axis,y_axis_temp)
    y_axis = numpy.cross(z_axis,x_axis)
    

    plotsize = 250.*zoom


    xmin = plotsize*-1.
    xmax = plotsize
    ymin = plotsize*-1.
    ymax = plotsize
    extend_radius = numpy.sqrt((ymax-ymin)**2 + (xmax-xmin)**2)
  
  #xmin = -1000.
    print 'region',xmin,xmax,ymin,ymax


    fig = pylab.figure()
    axes = fig.add_subplot(212)

    snap1 = snap
    particle1 = getdensity2D(snap1,r0,x_axis,y_axis,z_axis,xmin,xmax,ymin,ymax)
    axes.scatter(particle1[:,0],particle1[:,1],edgecolors='none',s = 1, color=plot_config.particle_color)

    #extend_radius = numpy.sqrt((xmax-xmin)**2 + (ymax-ymin)**2)
    cond = (data.halo['f1'] == snap1)
    condu = cond
    condl = cond

  # make a box to contain local halos
    for i in range(3):
        if(i==0):
            col = 'f4'
        elif(i==1):
            col = 'f5'
        else:
            col = 'f6'
        print i,col
        lower = r0[i] - extend_radius
        condu = cond
        condl = cond
        if(lower <= 0.):
            condl = condl & ((data.halo[col] >= (lower + boxsize)) | (data.halo[col] <= r0[i]))
        else:
            condl = condl & (data.halo[col] >= lower) 
            condl = condl & (data.halo[col] <= r0[i])
        upper = r0[i] + extend_radius
        if(upper >= boxsize):
            condu = condu & ((data.halo[col] <= (upper - boxsize)) | (data.halo[col] >= r0[i]))
        else:
            condu = condu & (data.halo[col] <= upper) 
            condu = condu & (data.halo[col] >= r0[i])
        cond = cond & (condl | condu)
        
    
    halolist = numpy.where(cond)

    print halolist
    cond = None
    condu = None
    condl = None
    #loop to plot local halo if not in the merger list
    for i in range(len(halolist[0])):
        index = halolist[0][i]
        if(1 == 1):
            xs = data.halo[index][4]
            ys = data.halo[index][5]
            zs = data.halo[index][6]
            radius_s = data.halo[index][3]
            rs = numpy.array([xs,ys,zs])
            rs = rs - r0
            for j in range(3):
                if(rs[j] > boxsize/2.):
                    rs[j] -= boxsize
                if(rs[j] < -1.*boxsize/2.):
                    rs[j] += boxsize
            l = numpy.dot(rs,z_axis)
            r_projected = rs - z_axis*l
            xx = numpy.dot(r_projected,x_axis)
            yy = numpy.dot(r_projected,y_axis)
            if(xx > xmin and xx < xmax and yy > ymin and yy < ymax):
                print 'ADD',index,data.halo[index][0],xx,yy,radius_s
                if(index == int(data.trees[haloid][1])):
                    circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color=plot_config.halo0_color, lw = 2)
                else:
                    circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color=plot_config.halo3_color, lw = 2)
                axes.add_patch(circle)
            else:
                print 'OUT OF RANGE',index,data.halo[index][0],xx,yy,radius_s

  
 
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    axes.set_aspect(1.)
    
    
    axes.set_xticklabels([])
    axes.set_yticklabels([])

    axes.text(xmin*(1-0.1), ymin, data.groupname)
    
    if save is True:
        fig.savefig(data.groupname+'_'+str(pos[0])+str(pos[1])+str(pos[2])+'.png',bbox_inches='tight')
    else:
        fig.show()

def plotserial_id_fixpos(haloid,n,data,save,zoomsize,pos):
    snapid = data.halo[haloid[0]][1]
    for i in range(n):
        print haloid,snapid
        plotonepanel_id_fixpos(haloid,snapid,data,save,zoomsize,pos)
        for j in haloid.keys():
            if(snapid == data.halo[haloid[j]][1]):
                if(len(data.trees[haloid[j]]) < 2):
                    haloid[j] = 0;
                else:
                    haloid[j] = int(data.trees[haloid[j]][1])

        snapid -= 1

            

def plotonepanel_id_fixpos(haloid,snapid,data,save,zoomsize,pos):

    item = []
    x0 = pos[0]
    y0 = pos[1]
    z0 = pos[2]
    #pos = [x0,y0,z0]
    r0 = numpy.array([x0,y0,z0])
    r1 = numpy.array([x0+100.,y0,z0])
    r2 = numpy.array([x0,y0+100.,z0])

    x_axis = numpy.array([])
    y_axis = numpy.array([])
    z_axix = numpy.array([])
   
    x_axis = r1 - r0
    for j in range(3):
        if(x_axis[j] > boxsize/2.):
            x_axis[j] -= boxsize
        if(x_axis[j] < -1.*boxsize/2.):
            x_axis[j] += boxsize
    x_axis = x_axis/(numpy.sqrt(numpy.dot(x_axis,x_axis)))
    y_axis_temp = r2 - r0
    for j in range(3):
        if(y_axis_temp[j] > boxsize/2.):
            y_axis_temp[j] -= boxsize
        if(y_axis_temp[j] < -1.*boxsize/2.):
            y_axis_temp[j] += boxsize
    y_axis_temp = y_axis_temp/(numpy.sqrt(numpy.dot(y_axis_temp,y_axis_temp)))
    z_axis = numpy.cross(x_axis,y_axis_temp)
    y_axis = numpy.cross(z_axis,x_axis)
    

    plotsize = zoomsize


    xmin = plotsize*-1.
    xmax = plotsize
    ymin = plotsize*-1.
    ymax = plotsize
    extend_radius = numpy.sqrt((ymax-ymin)**2 + (xmax-xmin)**2)
  
  #xmin = -1000.
    print 'region',xmin,xmax,ymin,ymax


    fig = pylab.figure()
    axes = fig.add_subplot(111)

    snap1 = snapid #data.halo[haloid][1]
    particle1 = getdensity2D(snap1,r0,x_axis,y_axis,z_axis,xmin,xmax,ymin,ymax)
    axes.scatter(particle1[:,0],particle1[:,1],edgecolors='none',s = 1, color=plot_config.particle_color)

    #extend_radius = numpy.sqrt((xmax-xmin)**2 + (ymax-ymin)**2)
    cond = (data.halo['f1'] == snap1)
    condu = cond
    condl = cond

  # make a box to contain local halos
    for i in range(3):
        if(i==0):
            col = 'f4'
        elif(i==1):
            col = 'f5'
        else:
            col = 'f6'
        print i,col
        lower = r0[i] - extend_radius
        condu = cond
        condl = cond
        if(lower <= 0.):
            condl = condl & ((data.halo[col] >= (lower + boxsize)) | (data.halo[col] <= r0[i]))
        else:
            condl = condl & (data.halo[col] >= lower) 
            condl = condl & (data.halo[col] <= r0[i])
        upper = r0[i] + extend_radius
        if(upper >= boxsize):
            condu = condu & ((data.halo[col] <= (upper - boxsize)) | (data.halo[col] >= r0[i]))
        else:
            condu = condu & (data.halo[col] <= upper) 
            condu = condu & (data.halo[col] >= r0[i])
        cond = cond & (condl | condu)
        
    
    halolist = numpy.where(cond)

    print halolist
    cond = None
    condu = None
    condl = None
    #loop to plot local halo if not in the merger list
    for i in range(len(halolist[0])):
        index = halolist[0][i]
        if(1 == 1):
            xs = data.halo[index][4]
            ys = data.halo[index][5]
            zs = data.halo[index][6]
            radius_s = data.halo[index][3]
            rs = numpy.array([xs,ys,zs])
            rs = rs - r0
            for j in range(3):
                if(rs[j] > boxsize/2.):
                    rs[j] -= boxsize
                if(rs[j] < -1.*boxsize/2.):
                    rs[j] += boxsize
            l = numpy.dot(rs,z_axis)
            r_projected = rs - z_axis*l
            xx = numpy.dot(r_projected,x_axis)
            yy = numpy.dot(r_projected,y_axis)
            if(xx > xmin and xx < xmax and yy > ymin and yy < ymax):
                print 'ADD',index,data.halo[index][0],xx,yy,radius_s
                if(index in haloid.values()):
                    key = haloid.keys()[haloid.values().index(index)]
                    circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color=plot_config.halo_colours[key], lw = 2)
                else:
                    circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color=plot_config.halo3_color, lw = 2)
                axes.add_patch(circle)
            else:
                print 'OUT OF RANGE',index,data.halo[index][0],xx,yy,radius_s

  
 
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    axes.set_aspect(1.)
    
    
    axes.set_xticklabels([])
    axes.set_yticklabels([])

    
    #axes.text(xmin*(1-0.1), ymin, data.groupname)
    axes.text(xmin*(1-0.1), ymin,'Snapshot '+str(snapid)+', '+data.groupname, fontsize=16)
    if save is True:
        fig.savefig(data.groupname+'_'+str(data.halo[haloid[0]][0])+'_'+str(snap1)+'.png',bbox_inches='tight')
        #os.system("convert -density 100 "+data.groupname+'_'+str(data.halo[haloid][0])+'_'+str(snap1)+'.pdf '+data.groupname+'_'+str(data.halo[haloid][0])+'_'+str(snap1)+'.png')
        os.system("convert -density 100 "+data.groupname+'_'+str(data.halo[haloid[0]][0])+'_'+str(snap1)+'.png EPS3:'+data.groupname+'_'+str(data.halo[haloid[0]][0])+'_'+str(snap1)+'.eps')
    else:
        fig.show()


def plotonepanel_id(haloid,data,save,zoom,shiftx,shifty):

    item = []
    x0 = data.halo[haloid][4] + shiftx
    y0 = data.halo[haloid][5] + shifty
    z0 = data.halo[haloid][6]
    pos = [x0,y0,z0]
    r0 = numpy.array([x0,y0,z0])
    r1 = numpy.array([x0+100.,y0,z0])
    r2 = numpy.array([x0,y0+100.,z0])

    x_axis = numpy.array([])
    y_axis = numpy.array([])
    z_axix = numpy.array([])
   
    x_axis = r1 - r0
    for j in range(3):
        if(x_axis[j] > boxsize/2.):
            x_axis[j] -= boxsize
        if(x_axis[j] < -1.*boxsize/2.):
            x_axis[j] += boxsize
    x_axis = x_axis/(numpy.sqrt(numpy.dot(x_axis,x_axis)))
    y_axis_temp = r2 - r0
    for j in range(3):
        if(y_axis_temp[j] > boxsize/2.):
            y_axis_temp[j] -= boxsize
        if(y_axis_temp[j] < -1.*boxsize/2.):
            y_axis_temp[j] += boxsize
    y_axis_temp = y_axis_temp/(numpy.sqrt(numpy.dot(y_axis_temp,y_axis_temp)))
    z_axis = numpy.cross(x_axis,y_axis_temp)
    y_axis = numpy.cross(z_axis,x_axis)
    

    plotsize = data.halo[haloid][3]*zoom


    xmin = plotsize*-1.
    xmax = plotsize
    ymin = plotsize*-1.
    ymax = plotsize
    extend_radius = numpy.sqrt((ymax-ymin)**2 + (xmax-xmin)**2)
  
  #xmin = -1000.
    print 'region',xmin,xmax,ymin,ymax


    fig = pylab.figure()
    axes = fig.add_subplot(111)

    snap1 = data.halo[haloid][1]
    particle1 = getdensity2D(snap1,r0,x_axis,y_axis,z_axis,xmin,xmax,ymin,ymax)
    axes.scatter(particle1[:,0],particle1[:,1],edgecolors='none',s = 1, color=plot_config.particle_color)

    #extend_radius = numpy.sqrt((xmax-xmin)**2 + (ymax-ymin)**2)
    cond = (data.halo['f1'] == snap1)
    condu = cond
    condl = cond

  # make a box to contain local halos
    for i in range(3):
        if(i==0):
            col = 'f4'
        elif(i==1):
            col = 'f5'
        else:
            col = 'f6'
        print i,col
        lower = r0[i] - extend_radius
        condu = cond
        condl = cond
        if(lower <= 0.):
            condl = condl & ((data.halo[col] >= (lower + boxsize)) | (data.halo[col] <= r0[i]))
        else:
            condl = condl & (data.halo[col] >= lower) 
            condl = condl & (data.halo[col] <= r0[i])
        upper = r0[i] + extend_radius
        if(upper >= boxsize):
            condu = condu & ((data.halo[col] <= (upper - boxsize)) | (data.halo[col] >= r0[i]))
        else:
            condu = condu & (data.halo[col] <= upper) 
            condu = condu & (data.halo[col] >= r0[i])
        cond = cond & (condl | condu)
        
    
    halolist = numpy.where(cond)

    print halolist
    cond = None
    condu = None
    condl = None
    #loop to plot local halo if not in the merger list
    for i in range(len(halolist[0])):
        index = halolist[0][i]
        if(1 == 1):
            xs = data.halo[index][4]
            ys = data.halo[index][5]
            zs = data.halo[index][6]
            radius_s = data.halo[index][3]
            rs = numpy.array([xs,ys,zs])
            rs = rs - r0
            for j in range(3):
                if(rs[j] > boxsize/2.):
                    rs[j] -= boxsize
                if(rs[j] < -1.*boxsize/2.):
                    rs[j] += boxsize
            l = numpy.dot(rs,z_axis)
            r_projected = rs - z_axis*l
            xx = numpy.dot(r_projected,x_axis)
            yy = numpy.dot(r_projected,y_axis)
            if(xx > xmin and xx < xmax and yy > ymin and yy < ymax):
                print 'ADD',index,data.halo[index][0],xx,yy,radius_s
                if(index == haloid):
                    circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color=plot_config.halo0_color, lw = 2)
                else:
                    circle = pylab.Circle((xx,yy),fill=False, radius=radius_s, alpha=.5, color=plot_config.halo3_color, lw = 2)
                axes.add_patch(circle)
            else:
                print 'OUT OF RANGE',index,data.halo[index][0],xx,yy,radius_s

  
 
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    axes.set_aspect(1.)
    
    
    axes.set_xticklabels([])
    axes.set_yticklabels([])

    #axes.text(xmin*(1-0.1), ymin, data.groupname)
    
    if save is True:
        fig.savefig(data.groupname+'_'+str(data.halo[haloid][0])+'_'+str(pos[0])+str(pos[1])+str(pos[2])+'.pdf',bbox_inches='tight')
    else:
        fig.show()




def readposition(snapid):
    filename = '%(str1)s/snapdir_%(#)03d/%(str2)s%(#)03d' %{'str1': config[4], "str2": config[1], "#" : snapid}
    print filename
    filenum = 16
    gadget = gadgetPyIO.readone(filename,filenum)
    gadget *= 1000.
    np = len(gadget)/3

    print "Total particle",np
    gadget.shape = (np,3)
    return gadget



def getdensity2D(snapid,r0,x_axis,y_axis,z_axis,min_x,max_x,min_y,max_y):
    filename = '%(str1)s/snapdir_%(#)03d/62.5_dm_%(#)03d' %{'str1': config[4], "#" : snapid}
    print filename
    filenum = 16
    x_bin = int((max_x-min_x)/cellsize)
    y_bin =  int((max_y-min_y)/cellsize)
    print "reading data with parameters",filename,filenum,r0[0],r0[1],r0[2],x_axis[0],x_axis[1],x_axis[2],y_axis[0],y_axis[1],y_axis[2],z_axis[0],z_axis[1],z_axis[2],min_x,max_x,min_y,max_y,boxsize,cellsize,x_bin,y_bin
    plane = gadgetPydensity.readone(filename,filenum,r0[0],r0[1],r0[2],x_axis[0],x_axis[1],x_axis[2],y_axis[0],y_axis[1],y_axis[2],z_axis[0],z_axis[1],z_axis[2],min_x,max_x,min_y,max_y,boxsize,cellsize,x_bin,y_bin)
    n = len(plane)/2
    plane.shape = (n,2)
    #print plane
    print "Finish reading particle data"
    return plane











def getdensity(snapid,r0,x_axis,y_axis,z_axis,min_x,max_x,min_y,max_y):
    filename = '%(str1)s/snapdir_%(#)03d/62.5_dm_%(#)03d' %{'str1': config[4],  "#" : snapid}
    print filename
    filenum = 16
    gadget = gadgetPyIO.readone(filename,filenum)
    gadget *= 1000.
    np = len(gadget)/3

    print "Total particle",np
    gadget.shape = (np,3)
    rho_av = np/(boxsize*boxsize*boxsize)
    
    gridsize = 1000.
    ngrid = int(boxsize/gridsize + 0.000001)

    #need to put particle in small grids to save time
    print 'need to put particle in small grids to save time'
    headofchain = [-1 for i in range(ngrid**3)]
    linkedlist = [-1 for i in range(np)]
    radius = numpy.sqrt((max_x-min_x)*(max_x-min_x) + (max_y-min_y)*(max_y-min_y))
    print 'run loop'
    for (i,item) in enumerate(gadget):
        
        x_block = max(min(int(item[0]/gridsize),ngrid-1),0)
        y_block = max(min(int(item[1]/gridsize),ngrid-1),0)
        z_block = max(min(int(item[2]/gridsize),ngrid-1),0)
        block = x_block + y_block*ngrid +z_block*ngrid**2
        linkedlist[i] = headofchain[block]
        headofchain[block] = i
       
    


    #make list of grid to transform particle xyz
    print 'make list of grid to transform particle xyz'
    lower_bound_x = numpy.fmod(r0[0] - radius + boxsize,boxsize)
    lower_bound_block_x = int(lower_bound_x/gridsize)
    upper_bound_x = numpy.fmod(r0[0] + radius + boxsize,boxsize)
    upper_bound_block_x = int(upper_bound_x/gridsize)
    start = lower_bound_block_x
    stop = upper_bound_block_x
    list_x = [start,stop]
    while start != stop:
        start = numpy.fmod(start+1,ngrid)
        list_x.append(start)
    list_x = sets.Set(list_x)

    lower_bound_y = numpy.fmod(r0[1] - radius + boxsize,boxsize)
    lower_bound_block_y = int(lower_bound_y/gridsize)
    upper_bound_y = numpy.fmod(r0[1] + radius + boxsize,boxsize)
    upper_bound_block_y = int(upper_bound_y/gridsize)
    start = lower_bound_block_y
    stop = upper_bound_block_y
    list_y = [start,stop]
    while start != stop:
        start = numpy.fmod(start+1,ngrid)
        list_y.append(start)
    list_y = sets.Set(list_y)

    lower_bound_z = numpy.fmod(r0[2] - radius + boxsize,boxsize)
    lower_bound_block_z = int(lower_bound_z/gridsize)
    upper_bound_z = numpy.fmod(r0[2] + radius + boxsize,boxsize)
    upper_bound_block_z = int(upper_bound_z/gridsize)
    start = lower_bound_block_z
    stop = upper_bound_block_z
    list_z = [start,stop]
    while start != stop:
        start = numpy.fmod(start+1,ngrid)
        list_z.append(start)    
    list_z = sets.Set(list_z)

    #finalise list
    print 'finalise list'
    list_all = []
    for x in list_x:
        for y in list_y:
            for z in list_z:
                block = x + y*ngrid +z*ngrid**2
                list_all.append(block)

    sigma = rho_av*2.*radius
    x_bin = int((max_x-min_x)/cellsize)
    y_bin =  int((max_y-min_y)/cellsize)
    area = [[0 for i in range(y_bin)] for j in range(x_bin)]


    #loop for all listed grids
    print 'loop for all listed grids'
    for block in list_all:
        print "checking block", block
        current = headofchain[block]
        while current != -1:
            rs = gadget[current]
            rs = rs - r0
            for j in range(len(rs)):
                if(rs[j] > boxsize/2.):
                    rs[j] -= boxsize
                if(rs[j] < -1.*boxsize/2.):
                    rs[j] += boxsize
            l = numpy.dot(rs,z_axis)
            if(numpy.fabs(l) < radius):
                r_projected = rs - z_axis*l
                x = numpy.dot(r_projected,x_axis)
                y = numpy.dot(r_projected,y_axis)
                if(x < max_x and x > min_x and y < max_y and y > min_y):
                    x_block = int((x - min_x)/cellsize)
                    y_block = int((y - min_y)/cellsize)
                    area[x_block][y_block] += 1
            current = linkedlist[current]

    gadget = None
    area /= sigma


   # for (i,item) in enumerate(xy_plane):
   #     x = xy_plane[0]
   #     y = xy_plane[1]
   #     if(x < max_x and x > min_x and y < max_y and y > min_y):
   #         x_block = int((x - min_x)/cellsize)
   #         y_block = int((y - min_y)/cellsize)
   #         area[x_block][y_block] += 1
   # area /= sigma
    return area
