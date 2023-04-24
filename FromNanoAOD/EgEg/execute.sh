echo "Compiling the C++ files: "
g++ robustanalyzermain.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out 0 7 hists_metdata_ "/pnfs/iihe/cms/store/user/asahasra/JetMET/JetMET2022C_crabRun230214Ntuplizer/230214_152148/0000/DiPhoton10_trigNtuples_10.root" &
proc0=$!
./robustanalyzer.out 1 7 hists_metdata_ "/pnfs/iihe/cms/store/user/asahasra/JetMET/JetMET2022C_crabRun230214Ntuplizer/230214_152148/0000/DiPhoton10_trigNtuples_10.root" &
proc1=$!
./robustanalyzer.out 2 7 hists_metdata_ "/pnfs/iihe/cms/store/user/asahasra/JetMET/JetMET2022C_crabRun230214Ntuplizer/230214_152148/0000/DiPhoton10_trigNtuples_10.root" & 
proc2=$!
./robustanalyzer.out 3 7 hists_metdata_ "/pnfs/iihe/cms/store/user/asahasra/JetMET/JetMET2022C_crabRun230214Ntuplizer/230214_152148/0000/DiPhoton10_trigNtuples_10.root" &
proc3=$!
./robustanalyzer.out 4 7 hists_metdata_ "/pnfs/iihe/cms/store/user/asahasra/JetMET/JetMET2022C_crabRun230214Ntuplizer/230214_152148/0000/DiPhoton10_trigNtuples_10.root" &
proc4=$!
./robustanalyzer.out 5 7 hists_metdata_ "/pnfs/iihe/cms/store/user/asahasra/JetMET/JetMET2022C_crabRun230214Ntuplizer/230214_152148/0000/DiPhoton10_trigNtuples_10.root" &
proc5=$!
./robustanalyzer.out 6 7 hists_metdata_ "/pnfs/iihe/cms/store/user/asahasra/JetMET/JetMET2022C_crabRun230214Ntuplizer/230214_152148/0000/DiPhoton10_trigNtuples_10.root" &
proc6=$!
./robustanalyzer.out 7 7 hists_metdata_ "/pnfs/iihe/cms/store/user/asahasra/JetMET/JetMET2022C_crabRun230214Ntuplizer/230214_152148/0000/DiPhoton10_trigNtuples_10.root" &
proc7=$!

while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" -o -d "/proc/${proc6}" -o -d "/proc/${proc7}" ]
do
    echo "====Still executing===="
    sleep 5
done

echo "Run over. Clean up and combine files."

rm robustanalyzer.out

hadd -f hists_metdata.root hists_metdata_?.root

rm hists_metdata_?.root
