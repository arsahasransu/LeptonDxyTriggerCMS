#!/bin/sh

mkdir Job_$1_$2
ls
echo "============ Created the master run folder ============"

cd Job_$1_$2
THISDIR=$(pwd)
echo $THISDIR

cd /user/asahasra/MeasureTriggerEff/CMSSW_12_6_3/src/
eval `scramv1 runtime -sh`

cd $THISDIR
echo $CMSSW_BASE
echo $1
echo $2
echo $3

cp /user/asahasra/MeasureTriggerEff/CMSSW_12_6_3/src/LeptonDxyTriggerCMS/FromNanoAOD/EgEg/robustanalyzer* ./

echo "============ Copied the cpp run files file. Listing directory contents ============"
ls
echo "======================= End listing the directory contents ========================"
echo ""

echo "Compiling the C++ files: "
g++ robustanalyzermain.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out $(( $2%20 )) 19 hists_metdata_$1_$2_ $3

echo ""
echo "Run over. Copy the output files."
ls

cp hists_metdata_$1_$2_*.root $CMSSW_BASE/src/LeptonDxyTriggerCMS/FromNanoAOD/EgEg/
