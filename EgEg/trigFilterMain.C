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

    stringstream ss3cm;
    ss3cm<<"hists_M200dM20ctau3cm_"<<cnt<<".root";
    trigFilterAnalyzer tfa_3cmMC("./data/DiPhoSkim220609_STHDM3cmRun3Summer21.root",ss3cm.str());
    tfa_3cmMC.analyzersinglefile(cnt);
    
    stringstream ss30cm;
    ss30cm<<"hists_M200dM20ctau30cm_"<<cnt<<".root";
    trigFilterAnalyzer tfa_30cmMC("./data/DiPhoSkim220609_STHDM30cmRun3Summer21tuple.root",ss30cm.str());
    tfa_30cmMC.analyzersinglefile(cnt);
    
    stringstream ss1m;
    ss1m<<"hists_M200dM20ctau1m_"<<cnt<<".root";
    trigFilterAnalyzer tfa_1mMC("./data/DiPhoSkim220609_STHDM1mRun3Summer21tuple.root",ss1m.str());
    tfa_1mMC.analyzersinglefile(cnt);
    
    stringstream ss3m;
    ss3m<<"hists_M200dM20ctau3m_"<<cnt<<".root";
    trigFilterAnalyzer tfa_3mMC("./data/DiPhoSkim220609_STHDM3mRun3Summer21tuple.root",ss3m.str());
    tfa_3mMC.analyzersinglefile(cnt);
    
  }
  catch (char const* exc){
    cerr<<exc<<endl;
  }
  
  return -1;
}
