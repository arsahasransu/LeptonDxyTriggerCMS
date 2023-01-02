#ifndef ROBUSTANALYZER_H
#define ROBUSTANALYZER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"

using namespace std;

class robustanalyzer {
  
public:
  robustanalyzer(TString, TString, int);
  ~robustanalyzer();
  
  void analyzersinglefile(int);
  void addHLTFilterhist(TString);
  void fillHLTFilterhist(TString, vector<int>, vector<double>*, vector<double>*, vector<double>*);
  void addPhotonCollectionhist(TString);
  void fillPhotonCollectionhist(TString, vector<int>);
  void addObjectFilterAngMatchhist(TString);
  vector< pair<int,int> > fillObjectFilterAngMatchhist(TString, vector<int>, vector<double>*, vector<double>*, vector<double>*, vector<int>, vector<double>*, vector<double>*, vector<double>*);
  void sort(int*, double*, int);
  vector<int> getFiltMatchedPhoIndex(vector<pair<int,int>>, double, double);
    
private:
  
  int nC;

  TChain* inputChain;

  TFile* outfile;

  std::vector<TH1F*> all1dhists;

  bool HLT_DiPhoton10sminlt0p12, HLT_DiPhoton10Time1p4ns, HLT_DiPhoton10_CaloIdL;
  int run, lumi, dieg10sminlt0p12_usfinfilt_n, dieg10time1p4ns_usfinfilt_n, dieg10caloidl_usfinfilt_n, pho_n, pv_n;
  double bs_x, bs_y, bs_z;
  vector<double> *dieg10sminlt0p12_usfinfilt_pt, *dieg10sminlt0p12_usfinfilt_eta, *dieg10sminlt0p12_usfinfilt_phi, *dieg10time1p4ns_usfinfilt_pt, *dieg10time1p4ns_usfinfilt_eta, *dieg10time1p4ns_usfinfilt_phi, *dieg10caloidl_usfinfilt_pt, *dieg10caloidl_usfinfilt_eta, *dieg10caloidl_usfinfilt_phi;
  vector<double> *pho_pt, *pho_eta, *pho_phi, *pho_seedtime, *pho_sieie, *pho_hoe, *pho_smin, *pho_smax;
  vector<double> *pv_x, *pv_xerr, *pv_y, *pv_yerr, *pv_z, *pv_zerr;
  vector<bool> *pv_isvalid;
};

#endif
