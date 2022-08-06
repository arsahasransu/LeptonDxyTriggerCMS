echo "Compiling the C++ files: "
g++ robustanalyzermain.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out 0 &
proc0=$!
./robustanalyzer.out 1 &
proc1=$!
./robustanalyzer.out 2 &
proc2=$!
./robustanalyzer.out 3 &
proc3=$!
./robustanalyzer.out 4 &
proc4=$!
./robustanalyzer.out 5 &
proc5=$!

while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" ]
do
    echo "====Still executing===="
    sleep 5
done

echo "Run over. Clean up and combine files."

rm robustanalyzer.out

rm hists_data.root

hadd hists_data.root hists_data_?.root

rm hists_data_?.root
