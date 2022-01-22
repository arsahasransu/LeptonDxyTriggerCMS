#include "DDM10_robustanalyzer.hh"
#include <iostream>

using namespace std;

int main() {

  //DDM10_robustanalyzer trial("./data/output.root","hists_trialout", true);
  //trial.analyzer();

  DDM10_robustanalyzer ddm10ana_3cmMC("./data/STHDM3cm_DDM10Skim301.root","hists_M200dM20ctau3cm", true);
  ddm10ana_3cmMC.analyzer();

  DDM10_robustanalyzer ddm10ana_30cmMC("./data/STHDM30cm_DDM10Skim301.root","hists_M200dM20ctau30cm", true);
  ddm10ana_30cmMC.analyzer();

  DDM10_robustanalyzer ddm10ana_1mMC("./data/STHDM1m_DDM10Skim301.root","hists_M200dM20ctau1m", true);
  ddm10ana_1mMC.analyzer();

  DDM10_robustanalyzer ddm10ana_3mMC("./data/STHDM3m_DDM10Skim300.root","hists_M200dM20ctau3m", true);
  ddm10ana_3mMC.analyzer();

  //DDM10_robustanalyzer ddm10ana_data1("./data/Efmrl1_DDM10Skim300.root","hists_Efmrl1", false);
  //ddm10ana_data1.analyzer();

  DDM10_robustanalyzer ddm10ana_Mff12("./data/HTo4Mu12_DDM10Skim300.root","hists_H2L4MuMff12", true);
  ddm10ana_Mff12.analyzer();

  DDM10_robustanalyzer ddm10ana_Mff25("./data/HTo4Mu25_DDM10Skim300.root","hists_H2L4MuMff25", true);
  ddm10ana_Mff25.analyzer();

  DDM10_robustanalyzer ddm10ana_Mff50("./data/HTo4Mu50_DDM10Skim300.root","hists_H2L4MuMff50", true);
  ddm10ana_Mff50.analyzer();

  return -1;
}
