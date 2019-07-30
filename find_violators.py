import numpy
import neuprint as neu
import json

def find_violations(neuprint_addr, dataset, roilist):
    client = neu.Client('https://emdata1.int.janelia.org:11000')
    # neuron extractor query (restrict to .*raced for now)
    extractquery = f"MATCH (n :`{dataset}-Neuron`) WHERE n.status=~'.*raced' {{}} RETURN n.bodyId as bodyid"

    # get roi violation set (default all)
    checkrois = None
    if roilist is not None and len(roilist) > 0:
        checkrois = set(roilist)
        # restrict neurons to these ROIs
        roistring = "AND ("
        for roi in roilist:
            roistring += f"n.{roi} OR "
        roistring = roistring[0:-4]
        roistring += ")"
        extractquery = extractquery.format(roistring)
    else:
        roiquery = f"MATCH (m :Meta) WHERE m.dataset=\"{dataset}\" RETURN m.superLevelRois AS rois"
        roires = client.fetch_custom(roiquery)
        checkrois = set(roires["rois"].iloc[0])
        extractquery = extractquery.format("")

    # extract working set
    candbodies_res = client.fetch_custom(extractquery)
    candbodies = candbodies_res["bodyid"].tolist()
    candbodies_set = set(candbodies) 


    # query db X neurons at a time (restrict to traced neurons for the targets)
    # (connections can be one way since all we care about connections only within the body set
    # (note: m cannot be restricted to body set up front because the calls are run in minibatches)

    outputsquery=f"WITH {{}} AS TARGETS MATCH(x :`{dataset}-ConnectionSet`)-\
    [:From]->(n :`{dataset}-Neuron`) WHERE n.bodyId in TARGETS WITH collect(x) as\
    csets, n UNWIND csets as cset MATCH (cset)-[:To]->(m :`{dataset}-Neuron`) where\
    (m.status=~'.*raced') RETURN n.bodyId AS body1, cset.roiInfo AS info, m.bodyId AS body2"

    total_conns = 0
    total_splits = 0
    total_msplits = 0

    candidate_connections = [] 

    print(f"Candidate bodies: {len(candbodies)}")

    # find candidate split bodies (might as well compute
    # the number of split neurons as well)
    for iter1 in range(0, len(candbodies), 50):
    #for iter1 in range(1100, len(candbodies), 50):
        #if iter1 == 1200:
        #    break
        print(f"fetch batch: {iter1}")

        currlist = candbodies[iter1:iter1+50]
        restrictedquery = outputsquery.format(currlist) 
        
        res = client.fetch_custom(restrictedquery)

        for idx, row in res.iterrows():
            if (row["body1"] == row["body2"]) or (row["body2"] not in candbodies_set):
                continue
            roiconns = json.loads(row["info"])
            roiconns_keys = set(roiconns.keys())
            
            inter_rois = checkrois.intersection(roiconns_keys)
            if len(inter_rois) > 0:
                total_conns += 1


            if len(inter_rois) > 1:
                total_splits += 1
                if len(inter_rois) > 2:
                    total_msplits += 1
                roi2count = {}
                for roi in inter_rois:
                    roi2count[roi] = roiconns[roi]["post"]
                candidate_connections.append((row["body1"], row["body2"], roi2count))

    print(f"Total connection pairs examined: {total_conns}\nTotal split across ROIs: {total_splits}\nTotal multi-splits across ROIs: {total_msplits}")

    # distance threshold
    pixel_threshold = 400

    # global stats
    roipairs_numvio = 0
    roipairs_wgtvio = 0
    real_roi_splits = 0
    roi_scores = {}
    roi_bodies = {}


    for connpair in candidate_connections:
        (body1, body2, roi2count) = connpair
        locquery = f"MATCH (n :`{dataset}-Neuron` {{bodyId: {body1}}})<-[:From]-(x :ConnectionSet)-[:To]->(m :`{dataset}-Neuron` {{bodyId: {body2}}}) MATCH (x)-[:Contains]->(s :PostSyn) RETURN s AS synapse"
        resloc = client.fetch_custom(locquery)

        #  check each location
        roipsd = {}
        for idx, row in resloc.iterrows():
            synapse = row["synapse"]
            for (roi, count) in roi2count.items():
                if roi in synapse:
                    roipsd[roi] = synapse["location"]["coordinates"]

        # look at all pairs of ROIs
        roipsd2 = list(roipsd.items())
        for iter1 in range(0, len(roipsd2)):
            for iter2 in range(iter1+1, len(roipsd2)):
                roi1, loc1 = roipsd2[iter1]
                roi2, loc2 = roipsd2[iter2]
                loc1 = numpy.array(loc1)
                loc2 = numpy.array(loc2)

                diff = (loc1-loc2)**2
                sumsq = diff.sum()
                dist = sumsq**(1/2)

                if dist < pixel_threshold:
                    roipairs_numvio += 1
                    roiwgt = min(roi2count[roi1], roi2count[roi2])
                    roipairs_wgtvio += roiwgt
                    if roi1 > roi2:
                        roi1, roi2 = roi2, roi1
                    if (roi1, roi2) not in roi_scores:
                        roi_scores[(roi1,roi2)] = 0
                        roi_bodies[roi1 + "-" + roi2] = []
                    roi_scores[(roi1,roi2)] += roiwgt 
                    roi_bodies[roi1 + "-" + roi2].append((roiwgt,body1, body2))
                else:
                    real_roi_splits += 1

    # print stats
    print(f"Total violations: {roipairs_numvio}\nWeighted violations: {roipairs_wgtvio}\nNum real ROI splits: {real_roi_splits}\nBad ROI boundaries: {len(roi_scores)}")

    # sort roi scores and print out worst pairs as just commma separated; dump roi_bodies to file
    sorted_rois = list(zip(roi_scores.values(),roi_scores.keys()))
    sorted_rois.sort()
    sorted_rois.reverse()
    print(sorted_rois)

    fout = open("connviolations.json", 'w')
    for key, bodies in roi_bodies.items():
        bodies.sort()
        bodies.reverse()
    fout.write(json.dumps(roi_bodies))
    fout.close()

    return


#find_violations('https://emdata1.int.janelia.org:11000', "hemibrain", ["PVLP", "LO"])
find_violations('https://emdata1.int.janelia.org:11000', "hemibrain", None)
