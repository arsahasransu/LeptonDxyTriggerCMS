/*
 * Gen selection applied for studies that need gen muons
 */

#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std;

void sort(double *pt, double *eta, double *phi, double *dxy, double *dxy_sig, int n) {

  // Sort according to pt
  for(unsigned int i=0; i<n; i++) {
    for(unsigned int j=0; j<n; j++) {
      if(*(pt+i)>*(pt+j)) {
	double temp = *(pt+i);
	*(pt+i) = *(pt+j);
	*(pt+j) = temp;
	temp = *(eta+i);
	*(eta+i) = *(eta+j);
	*(eta+j) = temp;
	temp = *(phi+i);
	*(phi+i) = *(phi+j);
	*(phi+j) = temp;
	temp = *(dxy+i);
	*(dxy+i) = *(dxy+j);
	*(dxy+j) = temp;
	temp = *(dxy_sig+i);
	*(dxy_sig+i) = *(dxy_sig+j);
	*(dxy_sig+j) = temp;
      }
    }
  }
}

int analyzer_mumu_singlefile(TString inrootfile, TString outrootfile) {

  auto chain = new TChain("events");
  chain->Add(inrootfile);
  
  int genLepN, muRecoN, muFiltN, htFiltN;
  double genLepPid[100];
  double genLepPt[100];
  double genLepEta[100];
  double genLepPhi[100];
  double genLepVx[100];
  double genLepVy[100];
  double genLepVz[100];
  double genLepNMom[100];
  double genLepMomPid[100];
  double muFiltPt[100];
  double muFiltEta[100];
  double muFiltPhi[100];
  double muRecoPt[100];
  double muRecoEta[100];
  double muRecoPhi[100];
  double muRecoDxy[100];
  double muRecoDxySig[100];
  double htFiltPt[2];
  double htFiltEta[2];
  double htFiltPhi[2];
  
  chain->SetBranchAddress("genLepn",&genLepN);
  chain->SetBranchAddress("genLepPIDarr",&genLepPid);
  chain->SetBranchAddress("genLepPtarr",&genLepPt);
  chain->SetBranchAddress("genLepEtaarr",&genLepEta);
  chain->SetBranchAddress("genLepPhiarr",&genLepPhi);
  chain->SetBranchAddress("genLepVxarr",&genLepVx);
  chain->SetBranchAddress("genLepVyarr",&genLepVy);
  chain->SetBranchAddress("genLepVzarr",&genLepVz);
  chain->SetBranchAddress("genLepNmomarr",&genLepNMom);
  chain->SetBranchAddress("genLepMomPIDarr",&genLepMomPid);
  chain->SetBranchAddress("muFiltn", &muFiltN);
  chain->SetBranchAddress("muFilt_pt", &muFiltPt);
  chain->SetBranchAddress("muFilt_eta", &muFiltEta);
  chain->SetBranchAddress("muFilt_phi", &muFiltPhi);
  chain->SetBranchAddress("mun", &muRecoN);
  chain->SetBranchAddress("mu_pt", &muRecoPt);
  chain->SetBranchAddress("mu_eta", &muRecoEta);
  chain->SetBranchAddress("mu_phi", &muRecoPhi);
  chain->SetBranchAddress("mu_dxy", &muRecoDxy);
  chain->SetBranchAddress("mu_dxy_sig", &muRecoDxySig);
  chain->SetBranchAddress("htFiltn", &htFiltN);
  chain->SetBranchAddress("htFilt_pt", &htFiltPt);
  chain->SetBranchAddress("htFilt_eta", &htFiltEta);
  chain->SetBranchAddress("htFilt_phi", &htFiltPhi);

  int totEntries = chain->GetEntries();
  std::cout<<"Events to process: "<<totEntries<<std::endl;

  auto outfile = new TFile(outrootfile,"RECREATE");

  // Genlevel objects
  auto genmumu_mult = new TH1F("genmumu_mult","#mu multiplicity",50,0,50);
  auto genmumu_pt = new TH1F("genmumu_pt","#mu p_{T} / GeV",500,0,500);
  auto genmumu_eta = new TH1F("genmumu_eta","#mu #eta",52,-2.6,2.6);
  auto genmumu_phi = new TH1F("genmumu_phi","#mu #phi",66,-3.3,3.3);
  auto genmumu_log10d0 = new TH1F("genmumu_log10d0","#mu log_{10} d_{0} / log_{10} cm",100,-5,5);

  int evtwith2genmu=0, evtwith2genmupt15=0;
  
  for(unsigned int event=0; event<totEntries; event++) {

    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;
    chain->GetEntry(event);

    // Gen-level selection
    int nMu=0, nMu15=0;
    for(int genctr=0; genctr<genLepN; genctr++) {
      if(abs(genLepMomPid[genctr])==9000007) {
	if(abs(genLepPid[genctr])==13) {
	  nMu++;
	  if(genLepPt[genctr]>15) nMu15++;
	}
      }
    }

    if(nMu<2) continue;
    if(nMu>2) cout<<"ERROR: More than 2 gen muons in an event are not possible."<<endl;
    if(nMu==2) evtwith2genmu++;
    if(nMu15==2) evtwith2genmupt15++;

    genmumu_mult->Fill(nMu);
    for(int genctr=0; genctr<genLepN; genctr++) {
      if(abs(genLepMomPid[genctr])==9000007) {
	if(abs(genLepPid[genctr])==13) {
	  genmumu_pt->Fill(genLepPt[genctr]);
	  genmumu_eta->Fill(genLepEta[genctr]);
	  genmumu_phi->Fill(genLepPhi[genctr]);
	  TLorentzVector genmu;
	  genmu.SetPtEtaPhiM(genLepPt[genctr],genLepEta[genctr],genLepPhi[genctr],0.1);
	  double d0 = genLepVx[genctr]*genmu.Py()-genLepVy[genctr]*genmu.Px();
	  d0 /= genLepPt[genctr];
	  genmumu_log10d0->Fill(TMath::Log10(TMath::Abs(d0)));
	}
      }
    }
    if(muRecoN<=0) continue;
    
  } // End of event loop

  std::cout<<"Event with 2 gen muons: "<<evtwith2genmu<<endl;
  std::cout<<"Event with 2 gen muons, pT>15: "<<evtwith2genmupt15<<endl;

  outfile->Write();
  outfile->Close();
  chain->Delete();
  return -1;
}

int main() {

  TString infile = "./data/STHDM1m_MuMu16DisplacedSkim2.root";
  TString outfile = "hists_STHDM1m_doublemu.root";
  analyzer_mumu_singlefile(infile, outfile);
  
  infile = "./data/STHDM30cm_MuMu16DisplacedSkim2.root";
  outfile = "hists_STHDM30cm_doublemu.root";
  analyzer_mumu_singlefile(infile, outfile);

  infile = "./data/STHDM3m_MuMu16DisplacedSkim2.root";
  outfile = "hists_STHDM3m_doublemu.root";
  analyzer_mumu_singlefile(infile, outfile);

  infile = "./data/STHDM3cm_MuMu16DisplacedSkim2.root";
  outfile = "hists_STHDM3cm_doublemu.root";
  analyzer_mumu_singlefile(infile, outfile);

  return -1;
}
