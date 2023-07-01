#!/bin/csh

rm -r LocalLibraries.package
rm LocalLibraries.zip

set prodId=`date +%F_%H-%M`

set centLow=$1
set centHigh=$2

set pTspectrum=$3

mkdir ./myOutput/${prodId}
mkdir ./jobs/log/${prodId}
mkdir ./jobs/err/${prodId}
mkdir ./jobs/submit/${prodId}
mkdir ./jobs/report/${prodId}
mkdir ./jobs/csh/${prodId}
mkdir ./jobs/list/${prodId}

star-submit-template -template submitToyMcEvtGen_sst.xml -entities productionId=${prodId},centLow=${centLow},centHigh=${centHigh},pTspectrum=${pTspectrum}
