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

  stringstream ss2;
  ss2 << argv[2];
  int cnt2;
  ss2 >> cnt2;

  try {
    
    stringstream ssdata;
    ssdata<<argv[3]<<cnt<<".root";
    robustanalyzer rana_data(argv[4], ssdata.str(), cnt2);
    rana_data.analyzersinglefile(cnt);
    
  }
  catch (char const* exc){
    cerr<<exc<<endl;
  }
  
  return -1;
}
