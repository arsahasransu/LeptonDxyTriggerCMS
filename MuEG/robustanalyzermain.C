#include "data_robustanalyzer.hh"
#include <iostream>

using namespace std;

int main() {

  //data_robustanalyzer drana_30cmMC("./data/STHDM30cm_MuEgSkim200.root","hists_M200dM20ctau30cm.root");
  //drana_30cmMC.analyzersinglefile();

  data_robustanalyzer drana_3cmMC("./data/STHDM3cm_Skim210.root","hists_M200dM20ctau3cm.root", true);
  drana_3cmMC.analyzersinglefile();

  //data_robustanalyzer drana_3mMC("./data/STHDM3m_MuEgSkim200.root","hists_M200dM20ctau3m.root");
  //drana_3mMC.analyzersinglefile();

  data_robustanalyzer drana_1mMC("./data/STHDM1m_Skim210.root","hists_M200dM20ctau1m.root", true);
  drana_1mMC.analyzersinglefile();

  data_robustanalyzer drana_DYMC("./data/DYToLL_Skim210.root","hists_DY.root", true);
  drana_DYMC.analyzersinglefile();

  data_robustanalyzer drana_QCDPt20To30MC("./data/QCDPt20To30_Skim210.root","hists_QCDPt20To30.root", true);
  drana_QCDPt20To30MC.analyzersinglefile();

  data_robustanalyzer drana_QCDPt50To80MC("./data/QCDPt50To80_Skim210.root","hists_QCDPt50To80.root", true);
  drana_QCDPt50To80MC.analyzersinglefile();

  data_robustanalyzer drana_data1("./data/Efmrl1_Skim210.root","hists_Efmrl1.root", false);
  drana_data1.analyzersinglefile();

  data_robustanalyzer drana_data2("./data/Efmrl2_Skim210.root","hists_Efmrl2.root", false);
  drana_data2.analyzersinglefile();

  data_robustanalyzer drana_data3("./data/Efmrl3_Skim210.root","hists_Efmrl3.root", false);
  drana_data3.analyzersinglefile();

  data_robustanalyzer drana_data4("./data/Efmrl4_Skim210.root","hists_Efmrl4.root", false);
  drana_data4.analyzersinglefile();

  return -1;
}
