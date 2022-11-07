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
  
  void analyzersinglefile(int, int);
  void addhist(TString);
  void fillhistinevent(TString, vector<int>);
  void addhist4LowPtElectron(TString);
  void fillhistinevent4LowPtElectron(TString, vector<int>);
  void addhist4Photon(TString);
  void fillhistinevent4Photon(TString, vector<int>);
  void sort(int*, double*, int);
  
 private:
  
  TChain* inputChain;

  TFile* outfile;

  UInt_t run, lumi, bunch;
  Bool_t HLT_DiPhoton10sminlt0p12, HLT_DiPhoton10Time1p4ns, HLT_PFMET120_PFMHT120_IDTight, HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF, HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, HLT_CaloMET80_NotCleaned, HLT_PFMET200_NotCleaned, HLT_PFMET200_BeamHaloCleaned, HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;
  Float_t eventRho;
  ULong64_t event;
  UInt_t eln;
  Float_t elpt[100], eleta[100], elphi[100], eldetasc[100], eldxy[100], eldxyerr[100], eldz[100], eldzerr[100], elooemoop[100], elhoe[100], elsieie[100], elconvveto[100];
  UInt_t lowpteln;
  Float_t lowptelpt[100], lowpteleta[100], lowptelphi[100], lowpteldxy[100], lowpteldxyerr[100], lowpteldz[100], lowpteldzerr[100], lowptelooemoop[100], lowptelhoe[100], lowptelsieie[100];
  UInt_t phn;
  Float_t phopt[100], phoeta[100], phophi[100], phosieie[100], phohoe[100];

  std::vector<TH1F*> all1dhists;
};

#endif
