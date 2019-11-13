#!/bin/bash

#Change to makefile at some point

melalibdir="${CMSSW_BASE}/lib/${SCRAM_ARCH}"
melalibdata="${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/${SCRAM_ARCH}"
melaincdir="${CMSSW_BASE}/src/ZZMatrixElement/MELA/interface"

cmsswlibdir="$CMSSW_RELEASE_BASE/lib/${SCRAM_ARCH}"

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${melalibdir}:${melalibdata}"
compile_x_macros_impl() {
    g++ -DSTANDALONE -I "$ROOTSYS/include" -I "$ROOFITSYS/include" -I "${melaincdir}" -I "${CMSSW_BASE}/src"  -I "${CMSSW_RELEASE_BASE}/src" compilereference_all.C HZZ4LeptonsAnalysis_${1}.C `root-config --glibs --libs --cflags` -L "$ROOFITSYS/lib"  -lRooFit -lRooFitCore -L "${CMSSW_RELEASE_BASE}/src" -L "${melalibdir}" -L "${cmsswlibdir}" -lZZMatrixElementMELA   -lZZMatrixElementMEMCalculators -lCondFormatsJetMETObjects -l JetMETCorrectionsModules -o RunReference$1 -D"product_$1" -L "${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/${SCRAM_ARCH}" ${DEBUG+-fdiagnostics-color} good_lumis.cpp HZZ4LeptonsAnalysis_common.C XYMETCorrection.C
#    g++ -DSTANDALONE -I "$ROOTSYS/include" -I "$ROOFITSYS/include" -I "${melaincdir}" -I "${CMSSW_BASE}/src"  -I "${CMSSW_RELEASE_BASE}/src" compilereference_all_data.C HZZ4LeptonsAnalysis_$1.C `root-config --glibs --libs --cflags` -L "$ROOFITSYS/lib"  -lRooFit -lRooFitCore -L "${CMSSW_RELEASE_BASE}/src" -L "${melalibdir}" -L "${cmsswlibdir}" -lZZMatrixElementMELA   -lZZMatrixElementMEMCalculators -lCondFormatsJetMETObjects -l JetMETCorrectionsModules -o RunReference$1_data -D"product_$1" -L "${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/${SCRAM_ARCH}"
    #g++ -I "$ROOTSYS/include" -I "$ROOFITSYS/include" -I "${melaincdir}" -I "${cmsswincdir}"  compilereference_$1_single.C HZZ4LeptonsAnalysis_$1.C `root-config --glibs --libs --cflags` -L "$ROOFITSYS/lib"  -lRooFit -lRooFitCore -L "${melalibdir}" -L "${cmsswlibdir}" -lZZMatrixElementMELA   -lZZMatrixElementMEMCalculators -lKaMuCaCalibration -lCondFormatsJetMETObjects -l JetMETCorrectionsModules -o RunReference$1
    #g++ -DSTANDALONE -I "$ROOTSYS/include" -I "$ROOFITSYS/include" -I "${melaincdir}" -I "${CMSSW_BASE}/src"  -I "${CMSSW_RELEASE_BASE}/src" compilereference_$1_signal.C HZZ4LeptonsAnalysis_$1.C `root-config --glibs --libs --cflags` -L "$ROOFITSYS/lib"  -lRooFit -lRooFitCore -L "${CMSSW_RELEASE_BASE}/src" -L "${melalibdir}" -L "${cmsswlibdir}" -lZZMatrixElementMELA   -lZZMatrixElementMEMCalculators -lCondFormatsJetMETObjects -l JetMETCorrectionsModules -o RunReference$1_signal
}

compile_x_macros() {
    if [ -z ${DEBUG+x} ] ; then
        compile_x_macros_impl "$@"
    else
        compile_x_macros_impl "$@" 2> >(head -n 20)
    fi
}

if [ "$1" == "4mu" ] || [ "$1" == "4e" ] || [ "$1" == "2e2mu" ]; then
    echo "Compiling $1 macros"
    compile_x_macros "$1"
elif [ "$1" == "all" ]; then    
    echo "Compiling $1 macros"
    compile_x_macros 4mu
    compile_x_macros 4e 
    compile_x_macros 2e2mu 
else 
    echo "Please provide an argument to the script: 4e, 4mu, 2e2mu or all"
    exit
fi
