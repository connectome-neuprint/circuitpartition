This repository contains an algorithm for partitioning the dual graph representation of the connectivity graph.
The config.json file contains an example configuration.  One can specify a list of ROIs to partition.  If no ROIs
are specified, the script will attempt to partition the entire dataset.  The user can also specify the number
of partitions and whether each wire should be weighted by the number of synapses in a given A-B connection. 

# Installation

* conda create -n partition neuprint-python
* conda activate partition
* conda install metis
* pip install metis
* pip install scipy

One must also download the neuprint private token.

# Running the command

    % export NEUPRINT_APPLICATION_CREDENTIALS="YOUR NEUPRINT TOKEN"
    % python createpartition.py config.json

This will output a list of neuron ids and their partition(s) in "parts.json".  The file "connparts.json" provids the list
of every A-B connection and the partition assignment [bodypre, bodypost, part].

# Finding ROI "violations"

find_violators.py contains a function for finding cases where a pair of neurons share connections across an ROI 
boundary.  This indicates areas of imprecise ROI boundaries.  (TODO: check for situations where a pre and post
synpse are divided by an ROI boundary, try to normalize the number of violations based on ROI size, speedup
point query)

# TODO

* Implement tool that exports a DVID ROI to enable visualization of partition
* Allow one to recursively partition and previous partitioning result
* Split connection pairs into different nodes if they are far aware
* Implement option to embed neuron's connectivity using skeleton representations


