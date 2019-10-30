#!/bin/bash

#If started with bash, setup environment then execute script in python
''':'
if [ ! -z "$NOTIFY_USER" ] ; then
    cat <(echo "Starting job with arguments $* on host `hostname`") $0 | mail "`whoami`"
fi
if [ -z "$CMSSW_BASE" ] ; then
    . /cvmfs/cms.cern.ch/cmsset_default.sh; eval `scramv1 runtime -sh`; . /cvmfs/cms.cern.ch/crab3/crab.sh
elif ! /usr/bin/which crab &> /dev/null; then
    . /cvmfs/cms.cern.ch/crab3/crab.sh
fi
exec python "$0" "$@"
exit 1
'''

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABAPI.RawCommand import execRaw
import sys

import trace
#trace=trace.Trace()
#trace.runfunc(execRaw,'getlog',sys.argv[1:])
#print trace.results()
execRaw('getlog',sys.argv[1:])
