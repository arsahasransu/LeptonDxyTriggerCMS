#include "data_robustanalyzer.hh"
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
    
    stringstream ss__21;
    ss__21<<"hists_M200dM20ctau3cm__2_"<<cnt<<".root";
    data_robustanalyzer drana__2_3cmMC("./data/STHDM3cm_EgEgSkim220407__2.root",ss__21.str(), true);
    drana__2_3cmMC.analyzersinglefile(cnt);
    
    stringstream ss__24;
    ss__24<<"hists_M200dM20ctau30cm__2_"<<cnt<<".root";
    data_robustanalyzer drana__2_30cmMC("./data/STHDM30cm_EgEgSkim220407__2.root",ss__24.str(), true);
    drana__2_30cmMC.analyzersinglefile(cnt);
    
    stringstream ss__26;
    ss__26<<"hists_M200dM20ctau1m__2_"<<cnt<<".root";
    data_robustanalyzer drana__2_1mMC("./data/STHDM1m_EgEgSkim220407__2.root",ss__26.str(), true);
    drana__2_1mMC.analyzersinglefile(cnt);

    stringstream ss__25;
    ss__25<<"hists_M200dM20ctau3m__2_"<<cnt<<".root";
    data_robustanalyzer drana__2_3mMC("./data/STHDM3m_EgEgSkim220407__2.root",ss__25.str(), true);
    drana__2_3mMC.analyzersinglefile(cnt);
    
    stringstream ss1;
    ss1<<"hists_M200dM20ctau3cm_"<<cnt<<".root";
    data_robustanalyzer drana_3cmMC("./data/STHDM3cm_EgEgSkim220407.root",ss1.str(), true);
    drana_3cmMC.analyzersinglefile(cnt);
    
    stringstream ss4;
    ss4<<"hists_M200dM20ctau30cm_"<<cnt<<".root";
    data_robustanalyzer drana_30cmMC("./data/STHDM30cm_EgEgSkim220407.root",ss4.str(), true);
    drana_30cmMC.analyzersinglefile(cnt);
    
    stringstream ss6;
    ss6<<"hists_M200dM20ctau1m_"<<cnt<<".root";
    data_robustanalyzer drana_1mMC("./data/STHDM1m_EgEgSkim220407.root",ss6.str(), true);
    drana_1mMC.analyzersinglefile(cnt);

    stringstream ss5;
    ss5<<"hists_M200dM20ctau3m_"<<cnt<<".root";
    data_robustanalyzer drana_3mMC("./data/STHDM3m_EgEgSkim220407.root",ss5.str(), true);
    drana_3mMC.analyzersinglefile(cnt);

    stringstream ss3;
    ss3<<"hists_DY_"<<cnt<<".root";
    data_robustanalyzer drana_dy("./data/DYToLLM50_EgEgSkim220407.root",ss3.str(), true);
    drana_dy.analyzersinglefile(cnt);

    stringstream ss2;
    ss2<<"hists_Efmrl1_"<<cnt<<".root";
    data_robustanalyzer drana_data1("./data/Efmrl1_EgEgSkim220407.root",ss2.str(), false);
    drana_data1.analyzersinglefile(cnt);
    /*
    stringstream ss7;
    ss7<<"hists_Efmrl2_"<<cnt<<".root";
    data_robustanalyzer drana_data2("./data/Efmrl2_EgEgSkim220407.root",ss7.str(), false);
    drana_data2.analyzersinglefile(cnt);
    */
  }
  catch (char const* exc){
    cerr<<exc<<endl;
  }
  
  return -1;
}
