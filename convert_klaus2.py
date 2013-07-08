import numpy
import pylab
mapfile = 'Mapback_order.txt'
klausfile = 'GarchingTreeMillenium_Klaus.txt'
outfile = 'GarchingTreeMillenium_Klaus'

map = {}
used = {}

f = open(mapfile)
line = f.read().splitlines()
for (i,item) in enumerate(line):
    key = item.split()
    map[str(int(key[0]))] = str(int(key[0]))
    used[str(int(key[0]))] = 0
line = None

bufout = open(outfile,'w+')
bufout.write('1'+'\n')
bufout.write('Convert from '+klausfile +'\n')
bufout.write(str(len(used))+'\n')
f = open(klausfile)
line = f.read().splitlines()
for (i,item) in enumerate(line):
    key = item.split()
    for j in range(len(key)):
        if( int(key[j]) > 1000000000000):
            key[j] = map[key[j]]
            used[key[j]] = 1
    # print key                   # 
    if(len(key) == 1):
        bufout.write(key[0]+'\n')
    elif(len(key) == 2):
        bufout.write(key[0]+'\t'+key[1]+'\n')

for item in used.keys():
    if(used[item] == 0):
        bufout.write(map[item] +'\t'+'0'+'\n')
bufout.write('END')
line = None
