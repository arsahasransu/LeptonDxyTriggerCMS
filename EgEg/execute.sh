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

hadd -f hists_M200dM20ctau3cm.root hists_M200dM20ctau3cm_?.root
hadd -f hists_M200dM20ctau30cm.root hists_M200dM20ctau30cm_?.root
hadd -f hists_M200dM20ctau3m.root hists_M200dM20ctau3m_?.root
hadd -f hists_M200dM20ctau1m.root hists_M200dM20ctau1m_?.root
hadd -f hists_M200dM20ctau3cm__2.root hists_M200dM20ctau3cm__2_?.root
hadd -f hists_M200dM20ctau30cm__2.root hists_M200dM20ctau30cm__2_?.root
hadd -f hists_M200dM20ctau3m__2.root hists_M200dM20ctau3m__2_?.root
hadd -f hists_M200dM20ctau1m__2.root hists_M200dM20ctau1m__2_?.root
hadd -f hists_DY.root hists_DY_?.root
hadd -f hists_Efmrl.root hists_Efmrl?_?.root

rm hists_M200dM20ctau3cm_?.root
rm hists_M200dM20ctau30cm_?.root
rm hists_M200dM20ctau3m_?.root
rm hists_M200dM20ctau1m_?.root
rm hists_M200dM20ctau3cm__2_?.root
rm hists_M200dM20ctau30cm__2_?.root
rm hists_M200dM20ctau3m__2_?.root
rm hists_M200dM20ctau1m__2_?.root
rm hists_DY_?.root
rm hists_Efmrl?_?.root
