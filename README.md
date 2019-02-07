This repository contains an algorithm for partitioning the dual graph representation of the connectivity graph.
The config.json file contains an example configuration.  One can specify a list of ROIs to partition.  If no ROIs
are specified, the script will attempt to partition the entire dataset.  The user can also specify the number
of partitions and whether each wire should be weighted by the number of synapses in a given A-B connection. 

# Installation

Currently installed in the flyem shared file system.  One must also download the neuprint private token.

# Running the command

  % export NEUPRINT_APPLICATION_CREDENTIALS="YOUR NEUPRINT TOKEN"
  % /groups/flyem/proj/create_partitions/miniconda3/envs/partition/bin/python createpartition.py config.json

This will output a list of neuron ids and their partition(s) in "parts.json".  The file "connparts.json" provids the list
of every A-B connection and the partition assignment [bodypre, bodypost, part].

# TODO

* Implement tool that exports a DVID ROI to enable visualization of partition
* Allow one to recursively partition and previous partitioning result
* Implement option to embed neuron's connectivity using skeleton representations


