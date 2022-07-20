#ifndef TRIGFILTERANALYZER_H
#define TRIGFILTERANALYZER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"

using namespace std;

class trigFilterAnalyzer {

public:

  trigFilterAnalyzer(TString, TString);
  ~trigFilterAnalyzer();

  void analyzersinglefile(int);
  vector< pair<int,int> > doGenMatching(vector<int>, int, double*, double*, double*);
  void addgenhist(TString);  
  void fillgenhistinevent(TString, vector<int>, vector<bool>);
  void combineFiltmchForGen(vector<int>, vector< pair<int,int> >, vector<bool>*);
  
private:

  TChain* inputChain;

  TFile* outfile;

  int HLT_DoublePhoton33_CaloIdL, dieg33_egFiltN, HLT_DiPhoton10Time1ns, dipho10time1ns_egFiltN, HLT_DiPhoton10sminlt0p16, dipho10sminlt0p16_egFiltN;
  double dieg33_egFiltPt[100], dieg33_egFiltEta[100], dieg33_egFiltPhi[100], dipho10time1ns_egFiltPt[100], dipho10time1ns_egFiltEta[100], dipho10time1ns_egFiltPhi[100], dipho10sminlt0p16_egFiltPt[100], dipho10sminlt0p16_egFiltEta[100], dipho10sminlt0p16_egFiltPhi[100];
  int genLepN;
  double genLepPid[100], genLepPt[100], genLepEta[100], genLepPhi[100], genLepPromptEta[100], genLepPromptPhi[100], genLepVx[100], genLepVy[100], genLepVz[100], genLepDxy[100], genLepLxy[100], genLepNMom[100], genLepMomPid[100], genLepMomPt[100], genLepMomEta[100], genLepMomPhi[100], pvx[100], pvy[100], pvz[100], genLepTimeLight[100], genLepTimeAct[100], genLepTimeDiff[100];
  
  std::vector<TH1F*> all1dhists;
};

#endif
