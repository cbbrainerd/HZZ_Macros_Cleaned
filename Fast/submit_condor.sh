#!/bin/bash
set -e

#exec 2>&1

job_failed() {
#Don't want abnormal termination from condor: run it again!
    printf 'Failed with exit code %i' "$1" | tee >&2
    if [ $1 -gt 128 ]; then
        exit 11
    else
        exit $1
    fi   
}

echo "Running `printf "%q " submit.sh "$@"`"
#echo "Environment is:"
#env
#echo "---------------"

if [[ "$#" -lt "8" ]]; then 
    echo "Takes 4-6 arguments."
    exit 1
fi
#SITE="$1"
YEAR="$1"
IS_MC="$2"
FOO="$3"
n="$4"
#Replace / with . since those cause major problems
DATASET="$(tr / . <<< "$5")"
FULL_FILELIST="${6:-fileList}"
FILELIST="inputs_${8}_${4}.txt"
if [[ -n "$7" ]] ; then
    files_per_job="$7"
    echo "$files_per_job"
    jobs_processed=$(($files_per_job*$n))
    first_job=$(($jobs_processed+1))
    last_job=$(($jobs_processed+$files_per_job))
    sed -n "${first_job},${last_job}p" "$FULL_FILELIST" > "$FILELIST"
fi
CLUSTER="$8"
JOBNAME="${9-unnamed}"

#touch "fallback_${CLUSTER}.${n}.tar.gz"

#A lot of the problems seem to be with xrd, so let's just copy it and read it off disk
#Probably not a good idea to run over too many files doing it this way
#while read line
#do
#    xrdcp "$line" . || { sleep 10 ; xrdcp "$line" . ; } || { sleep 20 ; xrdcp "$line" . ; }
#done < "$FILELIST"
#
#sed 's!.*/!!' -i "$FILELIST"

#Don't need since we copy the whole ZZMatrixElement directory now
#echo "CMSSW_BASE: $CMSSW_BASE"
#mkdir -p "$CMSSW_BASE/src/ZZMatrixElement/MEKD/"
#cp ZZMatrixElement_MEKD_src.tar.gz "$CMSSW_BASE/src/ZZMatrixElement/MEKD/"
#pushd "$CMSSW_BASE/src/ZZMatrixElement/MEKD/"
#tar -xkf ZZMatrixElement_MEKD_src.tar.gz
#rm ZZMatrixElement_MEKD_src.tar.gz
#popd

. /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
export CMSV=CMSSW_10_2_15
cd /cvmfs/cms.cern.ch/${SCRAM_ARCH}/cms/cmssw/${CMSV}/src
eval `scramv1 runtime -sh`
cd -

BASEDIR="$CMSSW_BASE/src"

#Condor, pick out number of files

#mkdir -p "$BASEDIR/MonoHiggs_2017/MonoHiggs_SR/94X_test/jobdir"
#mkdir -p "$BASEDIR/MonoHiggs_2017/MonoHiggs_SR/94X_test/histodir"

#echo "Running HtoZZto4Leptons Analysis with executables"
#source /cvmfs/cms.cern.ch/cmsset_default.sh

export PATH
export CMSSW_BASE

cmsswlibdir="$CMSSW_RELEASE_BASE/lib/${SCRAM_ARCH}"

echo "Using precompiled macros"
#echo "Compiling the macros"
#bash compilereference.sh "$SAMPLE"

#fi
#bash compilereference.sh "$SAMPLE"

echo 'Running over list of files:'
cat "$FILELIST"
echo '---End Filelist---'
./ZpXFast "$IS_MC" "${DATASET}_${1}_${n}" "$YEAR" "$FILELIST" > >(tee output_"$CLUSTER"_stdout) 2> >(tee output_"$CLUSTER"_stderr >&2) || job_failed $?
#./RunReference${1}_bkg ./bkg_input_${n}.txt 1 ./bkg_input_Fall17_AN.txt 1 ./data_input_${1}_2017_AN_miniaod.txt 1 "$SITE" "$YEAR" "$MC" "$DIRINPUT" "$DATASET"
#./RunReference${1}_data ./bkg_input_${n}.txt 1 ./bkg_input_Fall17_AN.txt 1 ./data_input_${1}_2017_AN_miniaod.txt 1 "$SITE" "$YEAR" "$MC" "$DATASET" "fileList"

#Check username later
set --
xrdfs root://cmseos.fnal.gov/ mkdir -p /store/user/cbrainer/FAKE_RATE/${JOBNAME}/${CLUSTER}/
#Try to copy output to EOS. If that fails, fallback to making a tarball to transfer back to the submitting machine. If some files succeed, only the failures are added to the tarball
for output_file in output*;
do 
    xrdcp "$output_file" root://cmseos.fnal.gov//store/user/cbrainer/FAKE_RATE/${JOBNAME}/${CLUSTER}/ || set -- "$@" "$output_file"
done

#if [ "x$#" != "x0" ]; then
#    exit 1
#    rm "fallback_${CLUSTER}.${n}.tar.gz"
#    tar -czvf "fallback_${CLUSTER}.${n}.tar.gz" "$@" || exit 23
#fi

exit 0

#mv ${workdir}/output.root    ${savedir}/.
#mv ${workdir}/output_bnn.txt ${savedir}/.
#mv ${workdir}/output_bnn.root ${savedir}/.
#mv ${workdir}/output_txt.txt ${savedir}/.
#mv ${workdir}/output_txt_vbf.txt ${savedir}/.

if [ -d "$_CONDOR_SCRATCH_DIR" ]; then
 rm -f $_CONDOR_SCRATCH_DIR/*
fi
