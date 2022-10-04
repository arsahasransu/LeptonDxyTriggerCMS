echo "Compiling the C++ files: "
g++ robustanalyzermain.C data_robustanalyzer.C `root-config --cflags --glibs` -o data_robustanalyzer.out

echo "Compiling successful. Begin execution."

./data_robustanalyzer.out 0 &
proc0=$!
./data_robustanalyzer.out 1 &
proc1=$!
./data_robustanalyzer.out 2 &
proc2=$!
./data_robustanalyzer.out 3 &
proc3=$!
./data_robustanalyzer.out 4 &
proc4=$!
./data_robustanalyzer.out 5 &
proc5=$!

while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" ]
do
    echo "====Still executing===="
    sleep 5
done

echo "Run over. Clean up and combine files."

rm data_robustanalyzer.out

rm hists_M200dM20ctau3cm.root
rm hists_M200dM20ctau30cm.root
rm hists_M200dM20ctau3m.root
rm hists_M200dM20ctau1m.root
rm hists_DY.root
rm hists_Efmrl.root

hadd hists_M200dM20ctau3cm.root hists_M200dM20ctau3cm_?.root
hadd hists_M200dM20ctau30cm.root hists_M200dM20ctau30cm_?.root
hadd hists_M200dM20ctau3m.root hists_M200dM20ctau3m_?.root
hadd hists_M200dM20ctau1m.root hists_M200dM20ctau1m_?.root
hadd hists_DY.root hists_DY_?.root
hadd hists_Efmrl.root hists_Efmrl_?.root

rm hists_M200dM20ctau3cm_?.root
rm hists_M200dM20ctau30cm_?.root
rm hists_M200dM20ctau3m_?.root
rm hists_M200dM20ctau1m_?.root
rm hists_DY_?.root
rm hists_Efmrl_?.root