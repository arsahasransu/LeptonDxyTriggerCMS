#include "data_robustanalyzer.hh"
#include <iostream>

using namespace std;

int main() {

  data_robustanalyzer drana_30cmMC("./data/STHDM30cm_MuMuSkim250.root","hists_M200dM20ctau30cm", true);
  drana_30cmMC.analyzersinglefile();

  data_robustanalyzer drana_3cmMC("./data/STHDM3cm_MuMuSkim250.root","hists_M200dM20ctau3cm", true);
  drana_3cmMC.analyzersinglefile();

  data_robustanalyzer drana_3mMC("./data/STHDM3m_MuMuSkim250.root","hists_M200dM20ctau3m", true);
  drana_3mMC.analyzersinglefile();

  data_robustanalyzer drana_1mMC("./data/STHDM1m_MuMuSkim250.root","hists_M200dM20ctau1m", true);
  drana_1mMC.analyzersinglefile();

  data_robustanalyzer drana_data1("./data/Efmrl1_MuMuSkim250.root","hists_Efmrl1", false);
  drana_data1.analyzersinglefile();

  data_robustanalyzer drana_data2("./data/Efmrl2_MuMuSkim250.root","hists_Efmrl2", false);
  drana_data2.analyzersinglefile();

  data_robustanalyzer drana_data3("./data/Efmrl3_MuMuSkim250.root","hists_Efmrl3", false);
  drana_data3.analyzersinglefile();

  data_robustanalyzer drana_data4("./data/Efmrl4_MuMuSkim250.root","hists_Efmrl4", false);
  drana_data4.analyzersinglefile();

  //data_robustanalyzer drana_Mff12("./data/H2L4MuMff12_MuMuSkim250.root","hists_H2L4MuMff12", true);
  //drana_Mff12.analyzersinglefile();

  //data_robustanalyzer drana_Mff25("./data/H2L4MuMff25_MuMuSkim250.root","hists_H2L4MuMff25", true);
  //drana_Mff25.analyzersinglefile();

  //data_robustanalyzer drana_Mff50("./data/H2L4MuMff50_MuMuSkim250.root","hists_H2L4MuMff50", true);
  //drana_Mff50.analyzersinglefile();

  return -1;
}
