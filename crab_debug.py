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

import argparse

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABAPI.RawCommand import  crabCommand

def getConfig(analysis,sample,site,mc,data,year,listOfFiles,storageSite='T3_US_FNALLPC'):
    crabConfig=config()
    crabConfig.General.requestName='_'.join((analysis,sample,mc,data,year))
    crabConfig.JobType.scriptExe='crab_wrapper.py'
    crabConfig.JobType.scriptArgs=['site=%s'%site,'year=%s'%year,'mc=%s'%mc,'sample=%s'%sample,'n=1']
    crabConfig.JobType.pluginName='Analysis'
    crabConfig.Data.outputPrimaryDataset=mc if mc!='NO' else data
    crabConfig.Data.publication=False
    crabConfig.JobType.outputFiles=[
        'output.root',
        'output_bnn.txt',
        'output_bnn.root',
        'output_txt.txt',
        'output_txt_vbf.txt'
    ]
    crabConfig.JobType.psetName='PSet.py' #Don't know if needed
    crabConfig.Data.userInputFiles=listOfFiles
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
        'HISTOShapesZpB_READ.root'
    ]
    crabConfig.Data.splitting='FileBased'
    crabConfig.Data.unitsPerJob=1
    crabConfig.Site.storageSite=storageSite
    crabConfig.General.transferOutputs=True
    crabConfig.General.transferLogs=True
    return crabConfig

test_filelist=[
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_1.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_2.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_3.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_4.root',
#'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_5.root',
'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root',
'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_1.root',
'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_2.root',
'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_3.root',
'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_4.root',
'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_5.root',
'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_6.root',
'/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_7.root'
]
crabConfig=getConfig('MonoHZZ','4mu','FNAL','DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8','NO','2017',test_filelist)

crabCommand('submit',config=crabConfig)

def dryrun():
    crabCommand('submit','--dryrun',config=crabConfig)
