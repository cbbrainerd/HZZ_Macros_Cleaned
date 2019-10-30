#!/bin/bash

#If started with bash, setup environment then execute script in python
''':'
export CRAB_BOOTSTRAP_DONE=1
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

import argparse
import os, sys

try:
    from CRABClient.UserUtilities import config, getUsernameFromSiteDB
    from CRABAPI.RawCommand import  crabCommand
except ImportError:
    try:
        os.environ['CRAB_BOOTSTRAP_DONE']
        raise
    except KeyError:
        os.execlp('/bin/bash','bash',*sys.argv)

def getConfig(analysis,sample,site,mc,data,year,listOfFiles,storageSite='T3_US_FNALLPC'):
    crabConfig=config()
    try:
        suffix=sys.argv[1]
    except IndexError:
        suffix='tmp'
    crabConfig.General.requestName='_'.join((analysis,sample,mc,data,year,'_'.join(sys.argv[1:])))
    crabConfig.JobType.scriptExe='crab_wrapper.py'
    #crabConfig.JobType.scriptExe='debug_crab.sh'
    dataset=mc if mc!='NO' else data
    crabConfig.JobType.scriptArgs=['site=%s'%site,'year=%s'%year,'mc=%s'%mc,'sample=%s'%sample,'n=1','dataset=%s'%dataset]
    crabConfig.JobType.pluginName='Analysis'
    crabConfig.Data.outputPrimaryDataset=mc if mc!='NO' else data
    crabConfig.Data.publication=False
    output_prefix='output_%s' % dataset
    crabConfig.JobType.outputFiles=[
        '%s.root' % output_prefix,
        '%s_bnn.txt' % output_prefix,
        '%s_bnn.root' % output_prefix,
        '%s_txt.txt' % output_prefix,
        '%s_txt_vbf.txt' % output_prefix
    ]
    crabConfig.JobType.psetName='PSet.py' #Don't know if needed
    crabConfig.Data.userInputFiles=listOfFiles
    #Not sure if something is needed here, or if transferred by default already
#./RunReference${1}_bkg ./BkgCards$SAMPLE$MC/bkg_input_${n}.txt 1 ./bkg_input_Fall17_AN.txt 1 ./data_input_${1}_2017_AN_miniaod.txt 1 "$SITE" "$YEAR" "$MC" "$DIRINPUT"

    crabConfig.JobType.inputFiles=[
        'crab_wrapper.py',
        'submit.sh',
        'compilereference.sh',
        'HZZ4LeptonsAnalysis_4mu.C',
        'HZZ4LeptonsAnalysis_4mu.h',
        'compilereference_4mu_signal.C',
        'compilereference_4mu_bkg.C',
        'compilereference_4mu_single.C',
        'compilereference_4mu_data.C',
        #'BkgCards%s%s/bkg_input_%s.txt' % (sample,mc,1),
        #'data_input_ ???
        'bkg_input_Fall17_AN.txt',
        'Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root',
        'ScaleFactors_mu_Moriond2018_final.root',
        'egammaEffi_txt_EGM2D_Moriond2018v1.root',
        'egammaEffi_txt_EGM2D_Moriond2018v1_gap.root',
        'egammaEffi_txt_EGM2D_runBCDEF_passingRECO_lowEt.root',
        'egammaEffi_txt_EGM2D_runBCDEF_passingRECO.root',
        'PU_Reweight_2017.root',
        'HISTOShapes2HDM_READ_ext.root',
        'HISTOShapesZpB_READ.root',
        #Make fine-grained later
        'compilereference_all_bkg.C',
        'compilereference_all_data.C',
        'HZZ4LeptonsAnalysis_4e.C',
        'HZZ4LeptonsAnalysis_4mu.C',
        'HZZ4LeptonsAnalysis_2e2mu.C',
        'RunReference2e2mu_bkg',
        'RunReference2e2mu_data',
        'HZZ4LeptonsAnalysis_all.h',
        'RunReference4e_bkg',
        'RunReference4e_data',
        'RunReference4mu_bkg',
        'RunReference4mu_data',
        'ZZMatrixElement_MEKD_src.tar.gz',
    ]
    crabConfig.Data.splitting='FileBased'
    crabConfig.Data.unitsPerJob=10
    crabConfig.Site.storageSite=storageSite
    crabConfig.General.transferOutputs=True
    crabConfig.General.transferLogs=True
    #test
    #crabConfig.Site.whitelist=[storageSite]
    return crabConfig

test_filelist=[
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_1.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_2.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_3.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_4.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_5.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_1.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_2.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_3.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_4.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_5.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_6.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_7.root'
]

class paths:
    def __init__(self):
        self.base_path='/store/user/cbrainer/'
        self.storage_site='T2_IT_Bari'

#crabConfig=getConfig('MonoHZZ','4mu','FNAL','DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8','NO','2017',test_filelist)
from get_filelist import get_filelist
from ElectronDatasets import electron_datasets
def taskname_to_lfn(taskname,primary_dataset=None):
    (number,task)=taskname.split(':',1)
    (username,taskdir)=task.split('_',1)
    if primary_dataset is None:
        (crab,analysis,primary_dataset)=taskdir.split('_')[0:3]
        assert(crab=='crab')
    return '/store/user/%s/%s/%s/%s' % (username,primary_dataset,taskdir,number)

filelist=get_filelist('DoubleEG','T2_IT_Bari',[taskname_to_lfn(task) for task in electron_datasets['DoubleEG']])

#Apparently not working with LFN now?
filelist=['root://cms-xrd-global.cern.ch/%s' % f for f in filelist]

dryrun=False
DEBUG=True
if DEBUG:
    #dryrun=True
    filelist=filelist[:10]
    import logging
    from CRABClient.ClientUtilities import setConsoleLogLevelVar
    setConsoleLogLevelVar(logging.DEBUG)

crabConfig=getConfig('MonoHZZ','4e','FNAL','NO','DoubleEG','2017',[str(f) for f in filelist])

if DEBUG:
    crabConfig.Data.unitsPerJob=1

#crabCommand('submit',config=crabConfig)

#DEBUGGING

#from submit import submit
#from raw_command_debug import crabCommandDebug

extra_args=['--wait']
if dryrun:
    extra_args.append('--dryrun')

crabCommand('submit',*extra_args,config=crabConfig)
