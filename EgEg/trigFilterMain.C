#include "trigFilterAnalyzer.hh"
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

    stringstream ss1;
    ss1<<"hists_M200dM20ctau3cm_"<<cnt<<".root";
    trigFilterAnalyzer tfa_3cmMC("./data/DiPhoSkim220609_STHDM3cmRun3Summer21.root",ss1.str());
    tfa_3cmMC.analyzersinglefile(cnt);
    
  }
  catch (char const* exc){
    cerr<<exc<<endl;
  }
  
  return -1;
}
