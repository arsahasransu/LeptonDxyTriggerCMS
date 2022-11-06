#include "robustanalyzer.hh"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {

  stringstream ss;
  ss << argv[1];
  int cnt;
  ss >> cnt;

  try {
    
    stringstream ssdata;
    ssdata<<argv[3]<<cnt<<".root";
    robustanalyzer rana_data(argv[2],ssdata.str(), true);
    rana_data.analyzersinglefile(cnt);
    
  }
  catch (char const* exc){
    cerr<<exc<<endl;
  }
  
  return -1;
}
