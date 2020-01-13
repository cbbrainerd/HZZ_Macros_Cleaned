#HZZ->4l macros

Still needs some cleaning

Setup instructions: (incomplete)
```
#!/bin/bash
#Setup CMSSW environement
version=10_2_15
. /cvmfs/cms.cern.ch/cmsset_default.sh
scram project "HZZ_macros_${version}" CMSSW "CMSSW_${version}"
cd "HZZ_macros_${version}/src"
eval `scram runtime -sh`
git cms-init
#Setup ZZMatrixElement
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
pushd ZZMatrixElement
bash setup.sh
popd
#Clone the macros
git clone https://github.com/cbbrainerd/HZZ_Macros_Cleaned macros
```

Now need to generate pileup distributions, etc.:

```
#Get golden json and pileup files, example below:
#Generate pileup histogram for 2018 data
wget https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
wget https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PileUp/pileup_latest.txt
pileupCalc.py -i Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt --inputLumiJSON pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 PU_Data_2018.root
```
