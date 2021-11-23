#include "data_robustanalyzer.hh"
#include <iostream>

using namespace std;

int main() {

  data_robustanalyzer drana_30cmMC("./data/STHDM30cm_MuMuSkim210.root","hists_M200dM20ctau30cm.root", true);
  drana_30cmMC.analyzersinglefile();

  data_robustanalyzer drana_3cmMC("./data/STHDM3cm_MuMuSkim210.root","hists_M200dM20ctau3cm.root", true);
  drana_3cmMC.analyzersinglefile();

  data_robustanalyzer drana_3mMC("./data/STHDM3m_MuMuSkim210.root","hists_M200dM20ctau3m.root", true);
  drana_3mMC.analyzersinglefile();

  data_robustanalyzer drana_1mMC("./data/STHDM1m_MuMuSkim210.root","hists_M200dM20ctau1m.root", true);
  drana_1mMC.analyzersinglefile();

  data_robustanalyzer drana_data1("./data/Efmrl1_MuMuSkim210.root","hists_Efmrl1.root", false);
  drana_data1.analyzersinglefile();

  //data_robustanalyzer drana_data2("./data/Efmrl2_MuMuSkim210.root","hists_Efmrl2.root", false);
  //drana_data2.analyzersinglefile();

  data_robustanalyzer drana_data3("./data/Efmrl3_MuMuSkim210.root","hists_Efmrl3.root", false);
  drana_data3.analyzersinglefile();

  data_robustanalyzer drana_data4("./data/Efmrl4_MuMuSkim210.root","hists_Efmrl4.root", false);
  drana_data4.analyzersinglefile();

  return -1;
}
