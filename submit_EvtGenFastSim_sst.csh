#!/bin/csh

#usage example: ./submit_EvtGenFastSim_sst.csh 0 40 0

rm -r LocalLibraries.package
rm LocalLibraries.zip

set prodId=`date +%F_%H-%M`

#select centrality range in %
#centrality bin edges in analysis: 0, 10, 40, 80 %
#submited for analysis: 0-40% and 40-80%
set centLow=$1
set centHigh=$2

#Two options 0 - flat pT spectrum, 1 - Levy pT spectrum (used to boost low pT D+- statistics)
set pTspectrum=$3

#number of events used for analysis is 1e5
#can set any lower number for testing
set nEvents=1e5

mkdir ./myOutput/${prodId}
mkdir ./jobs/log/${prodId}
mkdir ./jobs/err/${prodId}
mkdir ./jobs/submit/${prodId}
mkdir ./jobs/report/${prodId}
mkdir ./jobs/csh/${prodId}
mkdir ./jobs/list/${prodId}

star-submit-template -template submitToyMcEvtGen_sst.xml -entities productionId=${prodId},nEvents=${nEvents},centLow=${centLow},centHigh=${centHigh},pTspectrum=${pTspectrum}
