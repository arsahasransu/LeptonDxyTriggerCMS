echo "Compiling the C++ files: "
g++ trigFilterMain.C trigFilterAnalyzer.C `root-config --cflags --glibs` -o trigFilterAnalyzer.out

echo "Compiling successful. Begin execution."

./trigFilterAnalyzer.out 0 &
proc0=$!
./trigFilterAnalyzer.out 1 &
proc1=$!
./trigFilterAnalyzer.out 2 &
proc2=$!
./trigFilterAnalyzer.out 3 &
proc3=$!
./trigFilterAnalyzer.out 4 &
proc4=$!
./trigFilterAnalyzer.out 5 &
proc5=$!

while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" ]
do
    echo "====Still executing===="
    sleep 5
done

echo "Run over. Clean up and combine files."

rm trigFilterAnalyzer.out

hadd -f hists_M200dM20ctau3cm.root hists_M200dM20ctau3cm_?.root
hadd -f hists_M200dM20ctau30cm.root hists_M200dM20ctau30cm_?.root
hadd -f hists_M200dM20ctau1m.root hists_M200dM20ctau1m_?.root
hadd -f hists_M200dM20ctau3m.root hists_M200dM20ctau3m_?.root

rm hists_M200dM20ctau3cm_?.root
rm hists_M200dM20ctau30cm_?.root
rm hists_M200dM20ctau1m_?.root
rm hists_M200dM20ctau3m_?.root

hadd -f hists_M200dM20.root hists_M200dM20*.root
