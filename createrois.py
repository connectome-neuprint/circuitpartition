"""

usage: python createrois.py connparts.json [OPT: -n NUMBER POINTS]

output:
    0.roi.json, 1.roi.json ...

ROIs can be loaded into DVID by (small change need in roi format to support this):
    curl -X POST <DVID>/api/repo/UUID/instance -d '{"roi": "name"}'
    curl -X POST <DVID>/api/node/UUID/name/roi --data-binary @0.roi.json

or can be viewed in neutu by dragging and dropping the roi files
"""

import sys
import random
from neuprint import Client
import numpy
import json

partfile=sys.argv[1]
partdata=json.load(open(partfile))

SERVER = 'emdata1.int.janelia.org:11000'
np = Client(SERVER)

#  parse partitions
partitions = {}
for entry in partdata:
    b1, b2, part = entry
    if part not in partitions:
        partitions[part] = []
    partitions[part].append((b1,b2))

MAXQUERIES = 100
points = {}

# write ROIs
for part, bodypairs in partitions.items():
    print("Retrieving partition", part)
    numqueries = MAXQUERIES
    if numqueries > len(bodypairs): # original connections were already filtered
        numqueries = len(bodypairs)
    random.shuffle(bodypairs)
    
    for iter1 in range(numqueries):
        (body1, body2) = bodypairs[iter1]
        query = "MATCH (n :`hemibrain-Neuron`)-[:Contains]->(:SynapseSet)-[:Contains]->(x :Synapse)-[:SynapsesTo]->(y :Synapse)<-[:Contains]-(:SynapseSet)<-[:Contains]-(m :`hemibrain-Neuron`) WHERE n.bodyId=%d AND m.bodyId=%d RETURN n.bodyId AS bodyId1, m.bodyId as bodyId2, x.location AS location" % (body1, body2)
        res = np.fetch_custom(query)
        for idx, row in res.iterrows():
            loc = numpy.array(row["location"]["coordinates"])/128
            loctup = (loc[0], loc[1], loc[2])
            if loctup not in points:
                points[loctup] = []
            points[loctup].append(part+1)

# find bbox
minz=miny=minx=int(999999999)
maxz=maxy=maxx=int(-999999999)

for point, partlist in points.items():
    x,y,z = point 
    if x < minx:
        minx = x
    if y < miny:
        miny = y
    if z < minz:
        minz = z
    if x > maxx:
        maxx = x
    if y > maxy:
        maxy = y
    if z > maxz:
        maxz = z

xsize = int(maxx-minx+1)
ysize = int(maxy-miny+1)
zsize = int(maxz-minz+1)

# create ROI mask and add seeds (offset by min x,y,z)
mask = numpy.zeros((xsize,ysize,zsize), numpy.int32)

for point, partlist in points.items():
    partval = max(partlist, key=partlist.count)
    coord = numpy.array(point)-numpy.array([minx,miny,minz])
    mask[int(coord[0]),int(coord[1]),int(coord[2])] = partval

# perform flood fill
import scipy.ndimage as ndimage
while 0 in mask:
    mask2 = ndimage.maximum_filter(mask, 3)
    mask2[mask != 0] = mask[mask != 0]
    mask = mask2

# make ROI from each unique seed 
# save each ROI in json with meta that explains the format
for partnum in range(len(partitions)):
    x,y,z = numpy.where(mask==(partnum+1))
    x += int(minx)
    y += int(miny)
    z += int(minz)
    
    roi = []
    for stupid in range(len(x)):
        roi.append([int(z[stupid]), int(y[stupid]), int(x[stupid])])

    roidata = {}
    roidata["type"] = "points"
    roidata["resolution"] = "128"
    roidata["order"] = "zyx"
    roidata["roi"] = roi
    fout = open(str(partnum)+".part.json", 'w')
    fout.write(json.dumps(roidata))

