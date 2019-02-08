import sys
import metis
import json
import numpy
from neuprint import Client
import random

"""
installation:

* conda create -n partition neuprint-python
* conda activate partition
* conda install metis
* pip install metis
* set NEUPRINT_APPLICATION_CREDENTIALS

usage: python createpartition.py config.yaml
"""

configfile=sys.argv[1]
data=json.load(open(configfile))

# only partition within rois (if provided)
rois = data["roiFilter"]

threshold = data["connectionthreshold"]
if threshold is None:
    threshold = 6

SERVER = 'emdata1.int.janelia.org:11000'
np = Client(SERVER)


# find nodes in the graph
roifilter = ""
roifilterset = set()
ignorelabels = set(["location", "timeStamp", "confidence", "type"])
connections = {}
points = {}

roirent = {}
if rois is not None and len(rois) > 0:
    # traverse one ROI at a time
    query = "MATCH (m :Meta:hemibrain) RETURN m.superLevelRois"
    res = np.fetch_custom(query)
    major_rois = set(res.iloc[0,0])

    for iter1 in range(len(rois)):
        print("Processing ROI:", rois[iter1])
        # warn if ROI is not a super level ROI
        if rois[iter1] not in major_rois:
            print("ROI not a major ROI, compartment stats might be wrong")
        
        # calculate pins vs compute 
        totalcompute = 0
        totalpins = 0
        query = "MATCH (n :`hemibrain-Neuron`) WHERE (n.status=\"Roughly traced\" OR n.status=\"Prelim Roughly traced\" OR n.status=\"Traced\" OR n.status=\"Leaves\") AND n." + rois[iter1] + " RETURN n.roiInfo AS roiInfo, n.bodyId AS bodyid"
    
        res = np.fetch_custom(query)
        roi2pins = set()
        for idx, row in res.iterrows():
            roidata = json.loads(row["roiInfo"])
            bodyid = row["bodyid"]
            totalcompute += roidata[rois[iter1]]["pre"]
            for troi, val in roidata.items():
                if troi != rois[iter1] and troi in major_rois:
                    totalpins += 1
                    roi2pins.add(bodyid)
                    break

        print("pins: ", totalpins, "compute:", totalcompute)
        roirent[rois[iter1]] = (totalpins, totalcompute, roi2pins)

        roifilter = "AND ("
        roifilter += ("n.`" + rois[iter1] + "`") 
        roifilter += ")"

        roifilter += " AND ("
        roifilter += ("m.`" + rois[iter1] + "`") 
        roifilter += ")"

        query = "MATCH (n :`hemibrain-Neuron`)-[:Contains]->(:SynapseSet)-[:Contains]->(x :PreSyn)-[:SynapsesTo]->(y :PostSyn)<-[:Contains]-(:SynapseSet)<-[:Contains]-(m :`hemibrain-Neuron`) WHERE (n.status=\"Roughly traced\" OR n.status=\"Prelim Roughly traced\") AND (m.status=\"Roughly traced\" OR m.status=\"Prelim Roughly traced\") "  + roifilter + " RETURN x.`" + rois[iter1] + "` AS matchroi, n.bodyId AS bodyId1, m.bodyId as bodyId2"
        res = np.fetch_custom(query)
    
        #roifilterset = set([rois[iter1]])

        for idx, row in res.iterrows():
            body1 = row["bodyId1"]
            body2 = row["bodyId2"]

            #loc = numpy.array(row["location"]["coordinates"])/128
            if body1 == body2:
                continue
            
            if (body1, body2) not in points:
                points[(body1, body2)] = []
            #points[(body1, body2)].append(loc)

            matchroi = row["matchroi"] 
            #rois2 = set(row["rois"]) - ignorelabels 
            #if len(roifilterset.intersection(rois2)) > 0:
            if matchroi:
                if (body1,body2) not in connections:
                    connections[(body1,body2)] = 1
                else:
                    connections[(body1,body2)] += 1
else:
    # basic query (no locations, otherwise no tractable)
    query = "MATCH (n :`hemibrain-Neuron`)-[z:ConnectsTo]->(m :`hemibrain-Neuron`) WHERE (n.status=\"Roughly traced\" OR n.status=\"Prelim Roughly traced\") AND (m.status=\"Roughly traced\" OR m.status=\"Prelim Roughly traced\") AND z.weight >= " + str(threshold) + " RETURN n.bodyId AS bodyId1, m.bodyId AS bodyId2, z.weight as weight"
    res = np.fetch_custom(query)

    for idx, row in res.iterrows():
        body1 = row["bodyId1"]
        body2 = row["bodyId2"]

        if body1 == body2:
            continue
        connections[(body1,body2)] = row["weight"]

# add filtered connections
connectionsfiltered = {}
for pairnode, count in connections.items():
    if count >= threshold:
        connectionsfiltered[pairnode] = count
        # make placeholder node
        body1, body2 = pairnode

        connectionsfiltered[(body1, -1)] = 0
        connectionsfiltered[(body2, -1)] = 0
        
        connectionsfiltered[(-1, body1)] = 0
        connectionsfiltered[(-1, body2)] = 0

pair2index = {}
index2pair = {}
pair2edges = {}

weightWires = data["weightWires"]
weightNodes = data["weightNodes"]
if weightWires is None:
    weightWires = False
if weightNodes is None:
    weightNodes = True

# construct input for metis
nodenum = 0
for node, count in connectionsfiltered.items():
    pair2index[node] = nodenum
    index2pair[nodenum] = node
    nodenum += 1 



for node, count in connectionsfiltered.items():
    if count != 0:
        body1, body2 = node
        midnode = pair2index[node]
        in2mid = pair2index[(body2, -1)]
        in2mid_part = pair2index[(-1, body2)]
        mid2out = pair2index[(-1, body1)]
        mid2out_part = pair2index[(body1, -1)]
        if in2mid not in pair2edges:
            pair2edges[in2mid] = []
        if mid2out not in pair2edges:
            pair2edges[mid2out] = []
        if in2mid_part not in pair2edges:
            pair2edges[in2mid_part] = []
        if mid2out_part not in pair2edges:
            pair2edges[mid2out_part] = []   
        if weightWires:
            pair2edges[midnode] = [(in2mid, count), (mid2out, count)]
            pair2edges[mid2out].append((midnode, count))
            pair2edges[in2mid].append((midnode, count))
            
            # add skeleton edge if missing
            if mid2out_part not in pair2edges[mid2out]:
                pair2edges[mid2out].append((mid2out_part, 1))
            if in2mid_part not in pair2edges[in2mid]:
                pair2edges[in2mid].append((in2mid_part, 1))

            if mid2out not in pair2edges[mid2out_part]:
                pair2edges[mid2out_part].append((mid2out, 1))
            if in2mid not in pair2edges[in2mid_part]:
                pair2edges[in2mid_part].append((in2mid, 1))
        else:
            pair2edges[midnode] = [in2mid, mid2out]
            pair2edges[mid2out].append(midnode)
            pair2edges[in2mid].append(midnode)
            
            # add skeleton edge if missing
            if mid2out_part not in pair2edges[mid2out]:
                pair2edges[mid2out].append(mid2out_part)
            if in2mid_part not in pair2edges[in2mid]:
                pair2edges[in2mid].append(in2mid_part)

            if mid2out not in pair2edges[mid2out_part]:
                pair2edges[mid2out_part].append(mid2out)
            if in2mid not in pair2edges[in2mid_part]:
                pair2edges[in2mid_part].append(in2mid)


# import graph to metis
metisadjlist = []
metisnodesz = [] # ?? node weight or node size

for iter1 in range(len(pair2index)):
    metisadjlist.append(pair2edges[iter1])
    size = int(connectionsfiltered[index2pair[iter1]])
    if size == 0:
        size = 1
    metisnodesz.append(size)
if weightNodes is None:
    metisnodesz = None

metis_data = metis.adjlist_to_metis(metisadjlist, metisnodesz)


# run partitioner 
numparts = data["numPartitions"]
if numparts is None:
    numparts = 2
partvar = data["partitionVariance"]
if partvar is None:
    partvar = 0.2
cuts, res = metis.part_graph(metis_data, numparts, ubvec=[1+partvar])

# collect results
partitions = {}
neuronpart = {}
connparts = []
computeparts = [0]*numparts
for idx, part in enumerate(res):
    body1, body2 = index2pair[idx]
    partitions[index2pair[idx]] = part
    if body1 != -1 and body2 != -1:
        computeparts[part] += connections[(body1, body2)]
        if str(body1) not in neuronpart:
            neuronpart[str(body1)] = set()
        if str(body2) not in neuronpart:
            neuronpart[str(body2)] = set()
        connparts.append([int(body1), int(body2), int(part)])
        neuronpart[str(body1)].add(part)
        neuronpart[str(body2)].add(part)

print("Partition cost (higher worse)", cuts)

print("Intrinsic neurons")
for body, partlist in neuronpart.items():
    partlist = list(partlist)
    neuronpart[body] = partlist 
    if len(partlist) == 1:
        print(body, ",", partlist[0])

print("multi-part neurons")
mpneurons = set()
for body, partlist in neuronpart.items():
    if len(partlist) > 1:
        print(body, ",", partlist)
        mpneurons.add(body)

pinparts = []
for part in range(numparts):
    pinparts.append(set())

for idx, part in enumerate(res):
    body1, body2 = index2pair[idx]
    if body1 != -1 and body2 != -1:
        if str(body1) in mpneurons:
            pinparts[part].add(body1)
        if str(body2) in mpneurons:
            pinparts[part].add(body2)

for idx, part in enumerate(res):
    body1, body2 = index2pair[idx]

for partnum in range(numparts):
    print(partnum, len(pinparts[partnum]), computeparts[partnum])

# output neuron partition mappings and connection partition mappings
fout = open("parts.json", 'w')
fout.write(json.dumps(neuronpart))

fout = open("connparts.json", 'w')
fout.write(json.dumps(connparts))


"""
# ?! !! roi extraction not tested (maybe make standalone so easy to regenerate)

# write ROIs

# sample points if non-ROI query was used and load into points
if len(roifilterset) == 0:
    numqueries = 100 * numparts
    if numqueries > len(connections): # original connections were already filtered
        numqueries = len(connections)
    keys = list(connections.keys())

    for iter1 in range(numqueries):
        (body1, body2) = random.choice(keys)
        query = "MATCH (n :`hemibrain-Neuron`)-[:Contains]->(:SynapseSet)-[:Contains]->(x :Synapse)-[:SynapsesTo]->(y :Synapse)<-[:Contains]-(:SynapseSet)<-[:Contains]-(m :`hemibrain-Neuron`) WHERE n.bodyId=%d, m.bodyId=%d RETURN n.bodyId AS bodyId1, m.bodyId as bodyId2, x.location AS location" % (body1, body2)
        res = np.fetch_custom(query)
        for idx, row in res.iterrows():
            loc = numpy.array(row["location"])/128
            if (body1, body2) not in points:
                points[(body1, body2)] = []
            points[(body1, body2)].append(loc)

# find bbox
minz=miny=minx=999999999
maxz=maxy=maxx=-999999999

for pair, pointlist in points.items()
    for point in pointlist:
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

xsize = maxx-minx+1
ysize = maxy-miny+1
zsize = maxz-minz+1

# create ROI mask and add seeds (offset by min x,y,z)
mask = numpy.zeros((xsize,ysize,zsize))

# ?! handle multiple to one??
for pair in connectionsfiltered.keys():
    b1, b2 = pair
    if b1 != -1 and b2 != -1:
        if pair in points:
            pointlist = points[pair]
            for point in pointlist:
                mask[point-numpy.array([minx,miny,minz]] = partitions[pair]

# ?! perform flood fill

# ?! make ROI from each unique seed (apply offset)
"""







