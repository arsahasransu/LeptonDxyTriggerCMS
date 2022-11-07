#!/bin/sh

mkdir Job_$1_$2
#echo "============ Created the master run folder ============"

cd Job_$1_$2
echo $TMPDIR

cd /afs/cern.ch/work/a/asahasra/private/Run3DataWithNewTriggers/CMSSW_12_4_10_patch1/src/
eval `scramv1 runtime -sh`

cd $TMPDIR/../Job_$1_$2
echo $CMSSW_BASE

cp /afs/cern.ch/work/a/asahasra/private/Run3DataWithNewTriggers/CMSSW_12_4_10_patch1/src/LeptonDxyTriggerCMS/FromNanoAOD/EgEg/robustanalyzer* ./

echo "============ Copied the cpp run files file. Listing directory contents ============"
ls
echo "======================= End listing the directory contents ========================"
echo ""

echo "Compiling the C++ files: "
g++ robustanalyzermain.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out $2 /afs/cern.ch/work/a/asahasra/private/Run3DataWithNewTriggers/CMSSW_12_4_10_patch1/src/01d4ffc5-02d3-4f41-9afb-7bb4db04a7c7.root hists_metdata_01d4ffc5_$1_ 150
#./robustanalyzer.out $2 /afs/cern.ch/work/a/asahasra/private/Run3DataWithNewTriggers/CMSSW_12_4_10_patch1/src/b4418752-9427-495f-bea3-22884abb31d0.root hists_metdata_b4418752_$1_ 150
#./robustanalyzer.out $2 /afs/cern.ch/work/a/asahasra/private/Run3DataWithNewTriggers/CMSSW_12_4_10_patch1/src/7f0e7b28-799d-487d-bdb0-348443cb4987.root hists_metdata_7f0e7b28_$1_ 150
#./robustanalyzer.out $2 /afs/cern.ch/work/a/asahasra/private/Run3DataWithNewTriggers/CMSSW_12_4_10_patch1/src/890e99ce-7d46-486e-9b4b-8866a3151496.root hists_metdata_890e99ce_$1_ 150

echo ""
echo "Run over. Copy the output files."
ls

cp hists_metdata_*_$1_$2.root $CMSSW_BASE/src/LeptonDxyTriggerCMS/FromNanoAOD/EgEg/
