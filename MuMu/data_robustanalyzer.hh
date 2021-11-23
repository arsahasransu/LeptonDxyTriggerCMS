#ifndef DATA_ROBUSTANALYZER_H
#define DATA_ROBUSTANALYZER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"

using namespace std;

class data_robustanalyzer {

 public:
  data_robustanalyzer(TString, TString, bool);
  ~data_robustanalyzer();
  
  void analyzersinglefile();
  void addhist(TString);
  void addhistonce();
  void fillhistinevent(TString, vector<int>);
  void fillhistineventonce();
  void sort(int*, double*, int);
  
 private:
  
  TChain* inputChain;

  TFile* outfile;

  bool isMC;
  int genLepN, muRecoN, muFiltN, muFiltN33;
  int trigDoubleMu33;
  double bsx, bsy, bsz;
  double genLepPid[100], genLepPt[100], genLepEta[100], genLepPhi[100], genLepVx[100], genLepVy[100], genLepVz[100], genLepNMom[100], genLepMomPid[100];
  double muFiltPt[100], muFiltEta[100], muFiltPhi[100];
  double muFiltPt33[100], muFiltEta33[100], muFiltPhi33[100];
  double muRecoPt[100], muRecoEta[100], muRecoPhi[100], muRecoDxy[100], muRecoDxySig[100];

  std::vector<TH1F*> all1dhists;
};

#endif
