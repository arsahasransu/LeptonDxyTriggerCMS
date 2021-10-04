/*
 * Study the data
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

  cout<<"Processing file: "<<inrootfile<<endl;
  auto chain = new TChain("events");
  chain->Add(inrootfile);
  
  int muRecoN, muFiltN, muFiltN33, htFiltN, trigDoubleMu33;
  double muFiltPt[100];
  double muFiltEta[100];
  double muFiltPhi[100];
  double muFiltPt33[100];
  double muFiltEta33[100];
  double muFiltPhi33[100];
  double muRecoPt[100];
  double muRecoEta[100];
  double muRecoPhi[100];
  double muRecoDxy[100];
  double muRecoDxySig[100];
  double htFiltPt[2];
  double htFiltEta[2];
  double htFiltPhi[2];
  
  chain->SetBranchAddress("trig_DoubleMu33NoFiltersNoVtxDisplaced", &trigDoubleMu33);
  chain->SetBranchAddress("muFiltn", &muFiltN);
  chain->SetBranchAddress("muFilt_pt", &muFiltPt);
  chain->SetBranchAddress("muFilt_eta", &muFiltEta);
  chain->SetBranchAddress("muFilt_phi", &muFiltPhi);
  chain->SetBranchAddress("muFiltn_33", &muFiltN33);
  chain->SetBranchAddress("muFilt_pt_33", &muFiltPt33);
  chain->SetBranchAddress("muFilt_eta_33", &muFiltEta33);
  chain->SetBranchAddress("muFilt_phi_33", &muFiltPhi33);
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
  
  auto filt33mu_pt = new TH1F("filt33mu_pt","trig. #mu p_{T} / GeV",500,0,500);
  auto filt33mu_eta = new TH1F("filt33mu_eta","trig. #mu #eta",52,-2.6,2.6);
  auto filt33mu_phi = new TH1F("filt33mu_phi","trig. #mu #phi",66,-3.3,3.3);

  // parsel - pT>33, eta<2.5, d0>0.01 cm
  auto filt20parsel33_pt = new TH1F("filt20parsel33_pt","trig. #mu p_{T} / GeV",500,0,500);
  auto filt20parsel33_pt2 = new TH1F("filt20parsel33_pt2","trig. sub lead. #mu p_{T} / GeV",500,0,500);
  auto filt20parsel33_eta = new TH1F("filt20parsel33_eta","trig. #mu #eta",52,-2.6,2.6);
  auto filt20parsel33_phi = new TH1F("filt20parsel33_phi","trig. #mu #phi",66,-3.3,3.3);

  // nosel - pT>16, eta<2.5
  auto recomubasicsel_pt2 = new TH1F("recomubasicsel_pt2","trig. sub lead. #mu p_{T} / GeV",500,0,500);

  // cutdrmd0

  int cutflowgt2recomu = 0;
  int cutflowparsel = 0; // For data this is the desired rate, For signal this is fraction 1 - we wanna do better than this
  int cutflowbasicsel = 0;
  
  for(unsigned int event=0; event<totEntries; event++) {

    if(event%1000000==0) std::cout<<"Processed event: "<<event+1<<std::endl;
    //if(event>10000) break;
    chain->GetEntry(event);

    if(muRecoN<2) continue;
    cutflowgt2recomu++;
    
    // Sort the muon objects for testing the cuts
    sort(muRecoPt, muRecoEta, muRecoPhi, muRecoDxy, muRecoDxySig, muRecoN);

    // Fill the parent fiter objects
    for(unsigned int ctr=0; ctr<muFiltN33; ctr++) {
      filt33mu_pt->Fill(muFiltPt33[ctr]);
      filt33mu_eta->Fill(muFiltEta33[ctr]);
      filt33mu_phi->Fill(muFiltPhi33[ctr]);
    }
    
    // condition initialization
    // parsel
    vector<int> parselrecomuidx;
    bool parsel = true;
    // basicsel
    std::vector<unsigned int> basicselmuidx;
    bool basicselcond = true;

    for(unsigned int filtctr=0; filtctr<muFiltN; filtctr++) {

      // Find the reco object in the L3 filter
      int findpos = -1;
      for(unsigned int ctr=0; ctr<muRecoN; ctr++) {	
	if(muFiltPt[filtctr]==muRecoPt[ctr]) findpos = ctr;
      }
      if(findpos==-1) {
	cout<<"Error! Reco Object for Filter not found"<<endl;
	break;
      }

      // parselsinglemu
      bool parselsinglemu = true;
      parselsinglemu *= muRecoPt[findpos]>=33;
      parselsinglemu *= TMath::Abs(muRecoEta[findpos])<2.5;
      parselsinglemu *= TMath::Abs(muRecoDxy[findpos])>0.01;
      if(parselsinglemu) parselrecomuidx.push_back(findpos);

      // basicsel
      bool basicselsinglemucond = true;
      basicselsinglemucond *= muRecoPt[findpos]>=16;
      basicselsinglemucond *= TMath::Abs(muRecoEta[findpos])<2.5;
      if(basicselsinglemucond) basicselmuidx.push_back(findpos);

    }

    // parsel
    if(parselrecomuidx.size()!=muFiltN33) cout<<"Error: mismatching number of reco and filter muons"<<endl;
    parsel *= parselrecomuidx.size()>=1; // filter selects all the muons
    if(parsel) {
      for(unsigned int parselctr=0; parselctr<parselrecomuidx.size(); parselctr++) {
	filt20parsel33_pt->Fill(muRecoPt[parselrecomuidx[parselctr]]);
	filt20parsel33_eta->Fill(muRecoEta[parselrecomuidx[parselctr]]);
	filt20parsel33_phi->Fill(muRecoPhi[parselrecomuidx[parselctr]]);
      }
    }
    parsel *= parselrecomuidx.size()>=2; // events needs two of these electrons
    if(parsel) {
      if(trigDoubleMu33<=0) cout<<"Error: Inconsistent parent trigger match"<<endl;
      cutflowparsel++;
      filt20parsel33_pt2->Fill(muRecoPt[parselrecomuidx[1]]);
    }
    if(trigDoubleMu33>0 && !parsel) cout<<"Error: Inconsistent parent trigger match"<<endl;
      
    // basicsel
    basicselcond *= basicselmuidx.size()>=2;
    if(basicselcond) {
      cutflowbasicsel++;
      recomubasicsel_pt2->Fill(muRecoPt[basicselmuidx[1]]);
    }
    
  } // End of event loop

  cout<<cutflowgt2recomu<<"\t"<<cutflowparsel++<<"\t"<<cutflowbasicsel++<<endl;

  outfile->Write();
  outfile->Close();
  chain->Delete();
  return -1;
}

int main() {

  TString infile = "./data/Efmrl_MuMu16DisplacedSkim1.root";
  TString outfile = "hists_efmrl_MuMu16DisplacedSkim1.root";
  analyzer_mumu_singlefile(infile, outfile);
  
  infile = "./data/STHDM3cm_MuMu16DisplacedSkim2.root";
  outfile = "hists_STHDM3cm_doublemu.root";
  analyzer_mumu_singlefile(infile, outfile);
  
  infile = "./data/STHDM30cm_MuMu16DisplacedSkim2.root";
  outfile = "hists_STHDM30cm_doublemu.root";
  analyzer_mumu_singlefile(infile, outfile);
  
  infile = "./data/STHDM1m_MuMu16DisplacedSkim2.root";
  outfile = "hists_STHDM1m_doublemu.root";
  analyzer_mumu_singlefile(infile, outfile);
  
  infile = "./data/STHDM3m_MuMu16DisplacedSkim2.root";
  outfile = "hists_STHDM3m_doublemu.root";
  analyzer_mumu_singlefile(infile, outfile);
  
  return -1;
}
