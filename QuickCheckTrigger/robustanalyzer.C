/*
 * AUTHOR: Abanti Ranadhir Sahasransu - asahasra@cern.ch
 * The code now assumes exactly two gen electrons for a MC sample
 */

#include "robustanalyzer.hh"
#include <iostream>
#include <numeric>
#include <stdexcept>

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
robustanalyzer::robustanalyzer(TString filename, TString outfilename, int numCores) {

  nC = numCores;
  
  dieg10sminlt0p12_usfinfilt_pt = new vector<double>(0);
  dieg10sminlt0p12_usfinfilt_eta = new vector<double>(0);
  dieg10sminlt0p12_usfinfilt_phi = new vector<double>(0);
  dieg10time1p4ns_usfinfilt_pt = new vector<double>(0);
  dieg10time1p4ns_usfinfilt_eta = new vector<double>(0);
  dieg10time1p4ns_usfinfilt_phi = new vector<double>(0);
  dieg10caloidl_usfinfilt_pt = new vector<double>(0);
  dieg10caloidl_usfinfilt_eta = new vector<double>(0);
  dieg10caloidl_usfinfilt_phi = new vector<double>(0);
  
  pho_pt = new vector<double>(0);
  pho_eta = new vector<double>(0);
  pho_phi = new vector<double>(0);
  pho_seedtime = new vector<double>(0);
  pho_smin = new vector<double>(0);
  pho_smax = new vector<double>(0);
  
  pv_x = new vector<double>(0);
  pv_xerr = new vector<double>(0);
  pv_y = new vector<double>(0);
  pv_yerr = new vector<double>(0);
  pv_z = new vector<double>(0);
  pv_zerr = new vector<double>(0);
  
  inputChain = new TChain("demo/tree");
  cout<<"Initializing for file: "<<filename<<endl;
  inputChain->Add(filename);
  
  inputChain->SetBranchAddress("run", &run);
  inputChain->SetBranchAddress("lumSec", &lumi);

  inputChain->SetBranchAddress("HLT_DiPhoton10sminlt0p12", &HLT_DiPhoton10sminlt0p12);
  inputChain->SetBranchAddress("HLT_DiPhoton10Time1p4ns", &HLT_DiPhoton10Time1p4ns);
  inputChain->SetBranchAddress("HLT_DiPhoton10_CaloIdL", &HLT_DiPhoton10_CaloIdL);

  inputChain->SetBranchAddress("dieg10sminlt0p12_usfinfilt_n", &dieg10sminlt0p12_usfinfilt_n);
  inputChain->SetBranchAddress("dieg10sminlt0p12_usfinfilt_pt", &dieg10sminlt0p12_usfinfilt_pt);
  inputChain->SetBranchAddress("dieg10sminlt0p12_usfinfilt_eta", &dieg10sminlt0p12_usfinfilt_eta);
  inputChain->SetBranchAddress("dieg10sminlt0p12_usfinfilt_phi", &dieg10sminlt0p12_usfinfilt_phi);
  inputChain->SetBranchAddress("dieg10time1p4ns_usfinfilt_n", &dieg10time1p4ns_usfinfilt_n);
  inputChain->SetBranchAddress("dieg10time1p4ns_usfinfilt_pt", &dieg10time1p4ns_usfinfilt_pt);
  inputChain->SetBranchAddress("dieg10time1p4ns_usfinfilt_eta", &dieg10time1p4ns_usfinfilt_eta);
  inputChain->SetBranchAddress("dieg10time1p4ns_usfinfilt_phi", &dieg10time1p4ns_usfinfilt_phi);
  inputChain->SetBranchAddress("dieg10caloidl_usfinfilt_n", &dieg10caloidl_usfinfilt_n);
  inputChain->SetBranchAddress("dieg10caloidl_usfinfilt_pt", &dieg10caloidl_usfinfilt_pt);
  inputChain->SetBranchAddress("dieg10caloidl_usfinfilt_eta", &dieg10caloidl_usfinfilt_eta);
  inputChain->SetBranchAddress("dieg10caloidl_usfinfilt_phi", &dieg10caloidl_usfinfilt_phi);

  inputChain->SetBranchAddress("pho_n", &pho_n);
  inputChain->SetBranchAddress("pho_pt", &pho_pt);
  inputChain->SetBranchAddress("pho_eta", &pho_eta);
  inputChain->SetBranchAddress("pho_phi", &pho_phi);
  inputChain->SetBranchAddress("pho_seedtime", &pho_seedtime);
  inputChain->SetBranchAddress("pho_smin", &pho_smin);
  inputChain->SetBranchAddress("pho_smax", &pho_smax);

  inputChain->SetBranchAddress("pv_n", &pv_n);
  inputChain->SetBranchAddress("pv_x", &pv_x);
  inputChain->SetBranchAddress("pv_xerr", &pv_xerr);
  inputChain->SetBranchAddress("pv_y", &pv_y);
  inputChain->SetBranchAddress("pv_yerr", &pv_yerr);
  inputChain->SetBranchAddress("pv_z", &pv_z);
  inputChain->SetBranchAddress("pv_zerr", &pv_zerr);

  inputChain->SetBranchAddress("bs_x", &bs_x);
  inputChain->SetBranchAddress("bs_y", &bs_y);
  inputChain->SetBranchAddress("bs_z", &bs_z);

  outfile = new TFile(outfilename,"RECREATE");
}

// Fill the root file, close the root file, and handle deletions
robustanalyzer::~robustanalyzer() {
  
  inputChain->Delete();
  outfile->Write();
  outfile->Close();
}

// Analyzer for a single file
void robustanalyzer::analyzersinglefile(int splitCnt) { // Assume splitCnt to range from 0 to nCores

  int totEntries = inputChain->GetEntries();
  cout<<"Total number of entries: "<<totEntries<<endl;

  // Verfied that this logic to parallelize works
  int nCores = nC;
  // there is a lesser no.of events in the last core
  int beginevent = splitCnt*(totEntries/nCores);
  int endevent = (splitCnt+1)*(totEntries/nCores);
  if(beginevent>=totEntries) return;
  endevent = endevent<totEntries?endevent:totEntries;
  cout<<"Processing events in range: [ "<<beginevent<<" , "<<endevent<<" )"<<endl;
  int event = beginevent-1;

  addHLTFilterhist("dieg10sminlt0p12_final");
  addHLTFilterhist("dieg10time1p4ns_final");

  // Loop beginning on events
  for(unsigned int event=beginevent; event<endevent; event++) {
    
    inputChain->GetEntry(event);
    //if(event>10) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%1000==0) std::cout<<"Processed event: "<<event+1<<std::endl;
    
    if(dieg10sminlt0p12_usfinfilt_pt->size() != dieg10sminlt0p12_usfinfilt_n) throw "Error!!! Number of entries in vector \"dieg10sminlt0p12_usfinfilt_pt\" unequal entries comapred to \"dieg10sminlt0p12_usfinfilt_n\".";

    vector<int> dieg10sminlt0p12_finalidx(dieg10sminlt0p12_usfinfilt_n);
    iota(begin(dieg10sminlt0p12_finalidx), end(dieg10sminlt0p12_finalidx), 0);
    fillHLTFilterhist("dieg10sminlt0p12_final", dieg10sminlt0p12_finalidx, dieg10sminlt0p12_usfinfilt_pt, dieg10sminlt0p12_usfinfilt_eta, dieg10sminlt0p12_usfinfilt_phi);
    
    vector<int> dieg10time1p4ns_finalidx(dieg10time1p4ns_usfinfilt_n);
    iota(begin(dieg10time1p4ns_finalidx), end(dieg10time1p4ns_finalidx), 0);
    fillHLTFilterhist("dieg10time1p4ns_final", dieg10time1p4ns_finalidx, dieg10time1p4ns_usfinfilt_pt, dieg10time1p4ns_usfinfilt_eta, dieg10time1p4ns_usfinfilt_phi);

    dieg10sminlt0p12_finalidx.clear();
    dieg10time1p4ns_finalidx.clear();
  }
  
}

void robustanalyzer::fillHLTFilterhist(TString selection, vector<int> idx, vector<double> *pt, vector<double> *eta, vector<double> *phi) {

  if(idx.size()<=0) return; 
  TH1F* hmult = (TH1F*) outfile->Get(selection+"_filt_mult");
  TH1F* hpt = (TH1F*) outfile->Get(selection+"_filt_pt");
  TH1F* heta = (TH1F*) outfile->Get(selection+"_filt_eta");
  TH1F* hphi = (TH1F*) outfile->Get(selection+"_filt_phi");

  hmult->Fill(idx.size());
  for(int id : idx) {
    hpt->Fill(pt->at(id));
    heta->Fill(eta->at(id));
    hphi->Fill(phi->at(id));
  }
}

// Function to fill a set of histograms for HLT Photon Collection
void robustanalyzer::fillPhotonCollectionhist(TString selection, vector<int> idx) {

  if(idx.size()<=0) return; 
  TH1F* hmult = (TH1F*) outfile->Get(selection+"_pho_mult");
  TH1F* hpt = (TH1F*) outfile->Get(selection+"_pho_pt");
  TH1F* heta = (TH1F*) outfile->Get(selection+"_pho_eta");
  TH1F* hphi = (TH1F*) outfile->Get(selection+"_pho_phi");

}

// Function to add a set of histograms for HLT Filters
void robustanalyzer::addHLTFilterhist(TString selection) {
  all1dhists.push_back(new TH1F(selection+"_filt_mult","filter mult.",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"_filt_pt","p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_filt_eta","#eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"_filt_phi","#phi",66,-3.3,3.3));
}

// Function to add a set of histograms for HLT Photon Collection
void robustanalyzer::addPhotonCollectionhist(TString selection) {
  all1dhists.push_back(new TH1F(selection+"_pho_mult","photon mult.",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"_pho_pt","p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_pho_eta","#eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"_pho_phi","#phi",66,-3.3,3.3));
}

// Function to sort the indices based on a factor (Usually pT)
void robustanalyzer::sort(int* idx, double* factor, int n) {
  for(unsigned int i=0; i<n; i++) {
    for(unsigned int j=i+1; j<n; j++) {
      if(*(factor+*(idx+j))>*(factor+*(idx+i))) { // Sort in decreasing value of factor
	double temp = *(idx+i);
	*(idx+i) = *(idx+j);
	*(idx+j) = temp;
      }
    }
  }
}

