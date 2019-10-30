#!/bin/bash
set -e

exec 2>&1

echo "Running `printf "%q " submit.sh "$@"`"
echo "Environment is:"
env
echo "---------------"

echo "Running cmsenv"
pushd "$CMSSW_BASE/src"
eval `scramv1 runtime -sh`
popd
echo "Environment is:"
env
echo "---------------"

#if [[ "$#" -lt "6" ]]; then 
#    echo "Takes 6-7 arguments."
#    exit 1
#fi

#SITE="$1"
YEAR="$1"
IS_MC="$2"
SAMPLE="$3"
n="$4"
DATASET="$5"
FILELIST="${6:-fileList}"

BASEDIR="$CMSSW_BASE/src"

mkdir -p "$BASEDIR/MonoHiggs_2017/MonoHiggs_SR/94X_test/jobdir"
mkdir -p "$BASEDIR/MonoHiggs_2017/MonoHiggs_SR/94X_test/histodir"

echo "Running HtoZZto4Leptons Analysis with executables"
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

melalibdir="${CMSSW_BASE}/lib/${SCRAM_ARCH}"
melalibdata="${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/${SCRAM_ARCH}"
melaincdir="${CMSSW_BASE}/src/ZZMatrixElement/MELA/interface"

cmsswlibdir="$CMSSW_RELEASE_BASE/lib/${SCRAM_ARCH}"

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${melalibdir}:${melalibdata}"

#echo "Using precompiled macros"
#echo "Compiling the macros"
#bash compilereference.sh "$SAMPLE"

#Unpack madgraph source files (not included by crab by default)
echo "CMSSW_BASE: $CMSSW_BASE"
mkdir -p "$CMSSW_BASE/src/ZZMatrixElement/MEKD/"
cp ZZMatrixElement_MEKD_src.tar.gz "$CMSSW_BASE/src/ZZMatrixElement/MEKD/"
pushd "$CMSSW_BASE/src/ZZMatrixElement/MEKD/"
tar -xkf ZZMatrixElement_MEKD_src.tar.gz
rm ZZMatrixElement_MEKD_src.tar.gz
popd

#bash compilereference.sh "$SAMPLE"

if [[ "$SAMPLE" == "all" ]] ; then set 4mu 2e2mu 4e ; else set "$SAMPLE" ; fi
if [[ "$1" != "4mu" ]] && [[ "$1" != "2e2mu" ]] && [[ "$1" != "4e" ]]; then
    echo "\$SAMPLE must be one of 4mu, 2e2mu, 4e, or all, not \"${SAMPLE}\"."
    exit 1
fi
while [[ "$#" != 0 ]];
do
#./RunReference${1}_bkg ./bkg_input_${n}.txt 1 ./bkg_input_Fall17_AN.txt 1 ./data_input_${1}_2017_AN_miniaod.txt 1 "$SITE" "$YEAR" "$MC" "$DIRINPUT" "$DATASET"
#./RunReference${1}_data ./bkg_input_${n}.txt 1 ./bkg_input_Fall17_AN.txt 1 ./data_input_${1}_2017_AN_miniaod.txt 1 "$SITE" "$YEAR" "$MC" "$DATASET" "fileList"
./RunReference${1} "$IS_MC" "$DATASET" "$YEAR" "$FILELIST"
shift
done

#mv ${workdir}/output.root    ${savedir}/.
#mv ${workdir}/output_bnn.txt ${savedir}/.
#mv ${workdir}/output_bnn.root ${savedir}/.
#mv ${workdir}/output_txt.txt ${savedir}/.
#mv ${workdir}/output_txt_vbf.txt ${savedir}/.

if [ -d "$_CONDOR_SCRATCH_DIR" ]; then
 rm -f $_CONDOR_SCRATCH_DIR/*
fi
