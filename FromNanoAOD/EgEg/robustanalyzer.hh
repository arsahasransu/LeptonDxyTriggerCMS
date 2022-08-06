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
  long event;
  UInt_t eln;
  Float_t elpt[100], eleta[100], elphi[100];
  UInt_t lowpteln;
  Float_t lowptelpt[100], lowpteleta[100], lowptelphi[100];
  UInt_t phn;
  Float_t phpt[100], pheta[100], phphi[100];

  std::vector<TH1F*> all1dhists;
};

#endif
