#!/bin/bash


 for file in `ls /lustre/cms/store/user/raly/MonoHiggs/MC2017_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_MonoHiggs_13TeV_merged/roottree_leptons_crab_ZZTo4L_13TeV_powheg_pythia8*`; do

#for file in `cat a.txt`; do

echo $file
cat checkentries.C | sed "s?filename?${file}?g" > tmp.C
g++ -I $ROOTSYS/include tmp.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o checkentries
./checkentries

done

