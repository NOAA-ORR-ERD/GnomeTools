import numpy as np

#Headers in topo file
#TransposeArray
#BoundarySegments x
#-->Total number of bndry segs, first one listed is main outer edge
#WaterBoundaries x y
#-->Numbers in header line are the number of ow bndry pts and the total num of bndry pts
#   water boundary indexes boundary pots
#BoundaryPoints
#Topology
#DAGTree

f = open('C:\\Users\\amy.macfadyen\\Documents\\Projects\\Japan_debris\\HYCOM\\topology.dat')
lines = f.readlines()
f.close()

for ii,line in enumerate(lines):
    if line.startswith('Vertices'):
        v0 = ii + 2
    if line.startswith('BoundarySeg'):
        v1 = ii
    if line.startswith('WaterBound'):
        wb0 = ii + 1
    if line.startswith('BoundaryP'):
        wb1 = ii
        b0 = ii + 1
    if line.startswith('Topology'):
        b1 = ii
        break

vlons = []; vlats = [];
for line in lines[v0:v1]:
    vlon,vlat = (float(l) for l in line.split())
    vlons.append(vlon);vlats.append(vlat)
vlons = np.array(vlons);
vlats = np.array(vlats);

wbids = []
for line in lines[wb0:wb1]:
    wbids.append(int(line))
wbids = np.array(wbids)

bids = []
for line in lines[b0:b1]:
    bids.append(int(line))
bids = np.array(bids)

seg1_end_id = int(lines[v1+1])

blon = vlons[bids[0:seg1_end_id]]
blat = vlats[bids[0:seg1_end_id]]

wblon = blon[wbids]
wblat = blat[wbids]
