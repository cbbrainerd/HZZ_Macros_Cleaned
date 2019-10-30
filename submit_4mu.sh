#!/bin/bash
set -x -e

BASEDIR="$1"

mkdir -p "$BASEDIR/MonoHiggs_2017/MonoHiggs_SR/94X_test/jobdir"
mkdir -p "$BASEDIR/MonoHiggs_2017/MonoHiggs_SR/94X_test/histodir"

echo "Running HtoZZto4Leptons Analysis with executables RunRunReference4mu_data"
source /cvmfs/cms.cern.ch/cmsset_default.sh

export PATH
export CMSSW_BASE

melalibdir=${CMSSW_BASE}/lib/*
export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH

if [ -d "$_CONDOR_SCRATCH_DIR" ]; then
    workdir="$_CONDOR_SCRATCH_DIR";
else 
    workdir="$PWD";
fi
cd ${workdir};

savedir="$BASEDIR/MonoHiggs_2017/MonoHiggs_SR/94X_test/histodir"

echo "Working dir is $workdir"
#echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

echo "Compiling the macros"
bash compilereference.sh 4mu

./RunReference4mu_bkg ./sig_input_h150.txt 1 ./bkg_input_Fall17_AN.txt 1 ./data_input_4mu_2017_AN_miniaod.txt 1 FNAL 2017 Fall17 > ${workdir}/RunReference4mu_data_log 2>&1
cp -f ${workdir}/RunReference4mu_data_log "$BASEDIR/MonoHiggs_2017/MonoHiggs_SR/94X_test/jobdir/RunReference4mu_data_log"

mv ${workdir}/output.root    ${savedir}/.
mv ${workdir}/output_bnn.txt ${savedir}/.
mv ${workdir}/output_bnn.root ${savedir}/.
mv ${workdir}/output_txt.txt ${savedir}/.
mv ${workdir}/output_txt_vbf.txt ${savedir}/.

if [ -d "$_CONDOR_SCRATCH_DIR" ]; then
 rm -f $_CONDOR_SCRATCH_DIR/*
fi
