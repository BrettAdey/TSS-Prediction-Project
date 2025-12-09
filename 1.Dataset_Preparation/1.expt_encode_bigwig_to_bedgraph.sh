#!/bin/bash

conda create --name bigwig ucsc-bigwigtobedgraph -c bioconda
conda activate bigwig
bigWigToBedGraph enbw_plus.bigWig enbw_plus.bedgraph
