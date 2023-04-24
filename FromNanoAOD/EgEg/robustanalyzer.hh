#ifndef ROBUSTANALYZER_H
#define ROBUSTANALYZER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

using namespace std;

class robustanalyzer {

 public:
  robustanalyzer(TString, TString, int);
  ~robustanalyzer();
  
  void analyzersinglefile(int);
  void addhist(TString);
  void fillhistinevent(TString, vector<int>);
  double effectivearea(double);
  void sort(int*, double*, int);
  
 private:
  
  int nC;

  TTreeReader* tree;

  TTreeReaderValue<int> *run, *lumi;
  TTreeReaderValue<double> *rho;
  TTreeReaderValue<bool> *HLT_DiPhoton10sminlt0p12, *HLT_DiPhoton10Time1p4ns, *HLT_DiPhoton10_CaloIdL, *HLTOR_METTrig;

  TTreeReaderValue<int> *eln;
  TTreeReaderValue<vector<double>> *ele, *elpt, *eleta, *elphi, *eld0, *eldz, *elseedtime, *elsmin, *elsmaj, *elsieie, *eldeta, *eldphi, *elhoe, *elchhadiso, *elneuthadiso, *elphiso, *elooemoop, *elinnerhits;
  TTreeReaderValue<vector<bool>> *elconvveto;

  TFile* outfile;

  std::vector<TH1F*> all1dhists;
};

#endif
