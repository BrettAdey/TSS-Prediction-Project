#!/bin/bash

# Extract HPRT1 from experimental GROcap bedgraph
awk '$1=="chrX" && $2>=134429874 && $3<=134529874' \
  enbw_plus_perbp.bedgraph \
  > enbw_plus_HPRT1.bedgraph
