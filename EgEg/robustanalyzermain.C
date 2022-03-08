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
    
    stringstream ss1;
    ss1<<"hists_M200dM20ctau3cm_"<<cnt<<".root";
    data_robustanalyzer drana_3cmMC("./data/STHDM3cm_EgEgSkim22119.root",ss1.str(), true);
    drana_3cmMC.analyzersinglefile(cnt);
    
    stringstream ss4;
    ss4<<"hists_M200dM20ctau30cm_"<<cnt<<".root";
    data_robustanalyzer drana_30cmMC("./data/STHDM30cm_EgEgSkim220119.root",ss4.str(), true);
    drana_30cmMC.analyzersinglefile(cnt);
    
    stringstream ss5;
    ss5<<"hists_M200dM20ctau3m_"<<cnt<<".root";
    data_robustanalyzer drana_3mMC("./data/STHDM3m_EgEgSkim220119.root",ss5.str(), true);
    drana_3mMC.analyzersinglefile(cnt);
    
    stringstream ss6;
    ss6<<"hists_M200dM20ctau1m_"<<cnt<<".root";
    data_robustanalyzer drana_1mMC("./data/STHDM1m_EgEgSkim220119.root",ss6.str(), true);
    drana_1mMC.analyzersinglefile(cnt);
    
    stringstream ss3;
    ss3<<"hists_DY_"<<cnt<<".root";
    data_robustanalyzer drana_dy("./data/DYToLL_EgEgSkim220119.root",ss3.str(), true);
    drana_dy.analyzersinglefile(cnt);
    
    stringstream ss2;
    ss2<<"hists_Efmrl_"<<cnt<<".root";
    data_robustanalyzer drana_data1("./data/Efmrl1_EgEgSkim400.root",ss2.str(), false);
    drana_data1.analyzersinglefile(cnt);
    
    //data_robustanalyzer drana_30cmMC("./data/STHDM30cm_EgEgSkim300.root","hists_M200dM20ctau30cm.root", true);
    //drana_30cmMC.analyzersinglefile();
    
    //data_robustanalyzer drana_3cmMC("./data/STHDM3cm_EgEgSkim22119.root","hists_M200dM20ctau3cm.root", true);
    //drana_3cmMC.analyzersinglefile();
    
    //data_robustanalyzer drana_3mMC("./data/STHDM3m_EgEgSkim300.root","hists_M200dM20ctau3m.root", true);
    //drana_3mMC.analyzersinglefile();
    
    //data_robustanalyzer drana_1mMC("./data/STHDM1m_EgEgSkim300.root","hists_M200dM20ctau1m.root", true);
    //drana_1mMC.analyzersinglefile();
    
    //data_robustanalyzer drana_DYMC("./data/DYToLL_EgEgSkim300.root","hists_DY.root", true);
    //drana_DYMC.analyzersinglefile();
    
    //data_robustanalyzer drana_QCDPt20To30MC("./data/QCDPt20To30_MuEgSkim210.root","hists_QCDPt20To30.root", true);
    //drana_QCDPt20To30MC.analyzersinglefile();
    
    //data_robustanalyzer drana_QCDPt50To80MC("./data/QCDPt50To80_MuEgSkim210.root","hists_QCDPt50To80.root", true);
    //drana_QCDPt50To80MC.analyzersinglefile();
    
    //data_robustanalyzer drana_data1("./data/Efmrl1_EgEgSkim400.root","hists_Efmrl.root", false);
    //drana_data1.analyzersinglefile();
    
    //data_robustanalyzer drana_data2("./data/Efmrl2_MuEgSkim210.root","hists_Efmrl2.root", false);
    //drana_data2.analyzersinglefile();
    
    //data_robustanalyzer drana_data3("./data/Efmrl3_MuEgSkim210.root","hists_Efmrl3.root", false);
    //drana_data3.analyzersinglefile();
    
    //data_robustanalyzer drana_data4("./data/Efmrl4_MuEgSkim210.root","hists_Efmrl4.root", false);
    //drana_data4.analyzersinglefile();

  }
  catch (char const* exc){
    cerr<<exc<<endl;
  }
  
  return -1;
}
