echo "Compiling the C++ files: "
g++ robustanalyzermain.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out 0 ./data/890e99ce-7d46-486e-9b4b-8866a3151496.root hists_metdata_ &
proc0=$!
./robustanalyzer.out 1 ./data/890e99ce-7d46-486e-9b4b-8866a3151496.root hists_metdata_ &
proc1=$!
./robustanalyzer.out 2 ./data/890e99ce-7d46-486e-9b4b-8866a3151496.root hists_metdata_ &
proc2=$!
./robustanalyzer.out 3 ./data/890e99ce-7d46-486e-9b4b-8866a3151496.root hists_metdata_ &
proc3=$!
./robustanalyzer.out 4 ./data/890e99ce-7d46-486e-9b4b-8866a3151496.root hists_metdata_ &
proc4=$!
./robustanalyzer.out 5 ./data/890e99ce-7d46-486e-9b4b-8866a3151496.root hists_metdata_ &
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
