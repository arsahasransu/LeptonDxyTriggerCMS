echo "Compiling the C++ files: "
g++ robustanalyzermain.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out 0 ./data/01d4ffc5-02d3-4f41-9afb-7bb4db04a7c7.root hists_metdata_ &
proc0=$!
./robustanalyzer.out 1 ./data/01d4ffc5-02d3-4f41-9afb-7bb4db04a7c7.root hists_metdata_ &
proc1=$!
./robustanalyzer.out 2 ./data/01d4ffc5-02d3-4f41-9afb-7bb4db04a7c7.root hists_metdata_ &
proc2=$!
./robustanalyzer.out 3 ./data/01d4ffc5-02d3-4f41-9afb-7bb4db04a7c7.root hists_metdata_ &
proc3=$!
./robustanalyzer.out 4 ./data/01d4ffc5-02d3-4f41-9afb-7bb4db04a7c7.root hists_metdata_ &
proc4=$!
./robustanalyzer.out 5 ./data/01d4ffc5-02d3-4f41-9afb-7bb4db04a7c7.root hists_metdata_ &
proc5=$!

while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" ]
do
    echo "====Still executing===="
    sleep 5
done

echo "Run over. Clean up and combine files."

rm robustanalyzer.out

hadd -f hists_metdata.root hists_metdata_?.root

rm hists_metdata_?.root
