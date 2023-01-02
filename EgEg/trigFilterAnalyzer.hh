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
  vector< pair<int,int> > doGenMatching(vector<int>,vector<int>);
  void addgenhist(TString);
  void addhist(TString);
  void addhistgenmch(TString);
  
  void fillgenhistinevent(TString, vector<int>);
  void fillhistinevent(TString, vector<int>);
  void fillhistineventgenmch(TString, vector<int>, vector<int>);

private:

  TChain* inputChain;

  TFile* outfile;

  int HLT_DoublePhoton33_CaloIdL;
  int genLepN;
  double genLepPid[100], genLepPt[100], genLepEta[100], genLepPhi[100], genLepPromptEta[100], genLepPromptPhi[100], genLepVx[100], genLepVy[100], genLepVz[100], genLepDxy[100], genLepLxy[100], genLepNMom[100], genLepMomPid[100], genLepMomPt[100], genLepMomEta[100], genLepMomPhi[100], pvx[100], pvy[100], pvz[100], genLepTimeLight[100], genLepTimeAct[100], genLepTimeDiff[100];
  
  std::vector<TH1F*> all1dhists;
};

#endif
