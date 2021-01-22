#!/bin/csh

rm -r LocalLibraries.package
rm LocalLibraries.zip

set prodId=`date +%F_%H-%M`

mkdir ./myOutput/${prodId}
mkdir ./jobs/log/${prodId}
mkdir ./jobs/err/${prodId}
mkdir ./jobs/submit/${prodId}
mkdir ./jobs/report/${prodId}
mkdir ./jobs/csh/${prodId}
mkdir ./jobs/list/${prodId}

star-submit-template -template submitToyMcEvtGen.xml -entities productionId=${prodId}

#mv *.session.xml ./jobs/submit/${prodId}

#mv ./jobs/csh/*.csh ./jobs/csh/${prodId}
#mv ./jobs/csh/*.condor ./jobs/report/${prodId}

#mv ./jobs/report/*.log ./jobs/report/${prodId}
#mv ./jobs/report/*.report ./jobs/report/${prodId}

#mv ./jobs/list/*.list ./jobs/list/${prodId}
