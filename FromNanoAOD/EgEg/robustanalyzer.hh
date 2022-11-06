#ifndef ROBUSTANALYZER_H
#define ROBUSTANALYZER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"

using namespace std;

class robustanalyzer {

 public:
  robustanalyzer(TString, TString, bool);
  ~robustanalyzer();
  
  void analyzersinglefile(int);
  void addhist(TString);
  void fillhistinevent(TString, vector<int>);
  void sort(int*, double*, int);
  
 private:
  
  TChain* inputChain;

  TFile* outfile;

  int run, lumi, bunch;
  Bool_t HLT_DiPhoton10sminlt0p12, HLT_DiPhoton10Time1p4ns, HLT_PFMET120_PFMHT120_IDTight, HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF, HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, HLT_CaloMET80_NotCleaned, HLT_PFMET200_NotCleaned, HLT_PFMET200_BeamHaloCleaned, HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;
  Float_t eventRho;
  long event;
  UInt_t eln;
  Float_t elpt[100], eleta[100], elphi[100], eldetasc[100], eldxy[100], eldxyerr[100], eldz[100], eldzerr[100], elooemoop[100], elhoe[100], elsieie[100], elconvveto[100];
  UInt_t lowpteln;
  Float_t lowptelpt[100], lowpteleta[100], lowptelphi[100];
  UInt_t phn;
  Float_t phpt[100], pheta[100], phphi[100];

  std::vector<TH1F*> all1dhists;
};

#endif
