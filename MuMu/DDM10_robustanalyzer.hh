#ifndef DATA_ROBUSTANALYZER_H
#define DATA_ROBUSTANALYZER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

class DDM10_robustanalyzer {

 public:
  DDM10_robustanalyzer(TString, TString, bool);
  ~DDM10_robustanalyzer();
  
  void analyzer();
  void addhistGenPart(TString);
  void fillhistpereventGenPart(TString, vector<int>);
  void sort(int*, double*, int);
  int findmin(int*, double*, int, bool, bool);
  
 private:
  
  TChain* inputChain;

  TFile* outfile;

  bool isMC;
  int genLepN, l3dim33FiltN, l2dim23FiltN, l2dim23csFiltN, l3ddm10FiltN, l2ddm10FiltN;
  int trigDoubleMu33, trigDoubleL2Mu23NV2Cha, trigDoubleL2Mu23NV2ChaCS, trigDoubleL3Mu10NoVtxDisplaced, trigDoubleL2Mu10NoVtx2ChaPromptL3Mu0Veto;
  double bsx, bsy, bsz;
  int genLepPid[100], genLepNMom[100], genLepMomPid[100];
  double genLepPt[100], genLepEta[100], genLepPhi[100], genLepVx[100], genLepVy[100], genLepVz[100];
  double l3dim33FiltPt[100], l3dim33FiltEta[100], l3dim33FiltPhi[100];
  double l2dim23FiltPt[100], l2dim23FiltEta[100], l2dim23FiltPhi[100];
  double l2dim23csFiltPt[100], l2dim23csFiltEta[100], l2dim23csFiltPhi[100];
  double l3ddm10FiltPt[100], l3ddm10FiltEta[100], l3ddm10FiltPhi[100];
  double l2ddm10FiltPt[100], l2ddm10FiltEta[100], l2ddm10FiltPhi[100];
  
  std::vector<TH1F*> all1dhists;
  TString outFileName;
};

#endif
