#include "data_robustanalyzer.hh"
#include <iostream>
#include <numeric>

#include "TMath.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
data_robustanalyzer::data_robustanalyzer(TString filename, TString outfilename, bool issimu){

  isMC = issimu;
  
  inputChain = new TChain("events");
  cout<<"Initializing for file: "<<filename<<endl;
  inputChain->Add(filename);

  inputChain->SetBranchAddress("bsx", &bsx);
  inputChain->SetBranchAddress("bsy", &bsy);
  inputChain->SetBranchAddress("bsz", &bsz);
  inputChain->SetBranchAddress("trig_DoubleMu33NoFiltersNoVtxDisplaced", &trigDoubleMu33);
  inputChain->SetBranchAddress("muFiltn", &muFiltN);
  inputChain->SetBranchAddress("muFilt_pt", &muFiltPt);
  inputChain->SetBranchAddress("muFilt_eta", &muFiltEta);
  inputChain->SetBranchAddress("muFilt_phi", &muFiltPhi);
  inputChain->SetBranchAddress("muFiltn_33", &muFiltN33);
  inputChain->SetBranchAddress("muFilt_pt_33", &muFiltPt33);
  inputChain->SetBranchAddress("muFilt_eta_33", &muFiltEta33);
  inputChain->SetBranchAddress("muFilt_phi_33", &muFiltPhi33);
  inputChain->SetBranchAddress("mun", &muRecoN);
  inputChain->SetBranchAddress("mu_pt", &muRecoPt);
  inputChain->SetBranchAddress("mu_eta", &muRecoEta);
  inputChain->SetBranchAddress("mu_phi", &muRecoPhi);
  inputChain->SetBranchAddress("mu_dxy", &muRecoDxy);
  inputChain->SetBranchAddress("mu_dxy_sig", &muRecoDxySig);

  if(isMC) {
    inputChain->SetBranchAddress("genLepn",&genLepN);
    inputChain->SetBranchAddress("genLepPIDarr",&genLepPid);
    inputChain->SetBranchAddress("genLepPtarr",&genLepPt);
    inputChain->SetBranchAddress("genLepEtaarr",&genLepEta);
    inputChain->SetBranchAddress("genLepPhiarr",&genLepPhi);
    inputChain->SetBranchAddress("genLepVxarr",&genLepVx);
    inputChain->SetBranchAddress("genLepVyarr",&genLepVy);
    inputChain->SetBranchAddress("genLepVzarr",&genLepVz);
    inputChain->SetBranchAddress("genLepNmomarr",&genLepNMom);
    inputChain->SetBranchAddress("genLepMomPIDarr",&genLepMomPid);
  }
  
  outfile = new TFile(outfilename,"RECREATE");
}

// Fill the root file, close the root file, and handle deletions
data_robustanalyzer::~data_robustanalyzer() {

  inputChain->Delete();
  outfile->Write();
  outfile->Close();
}

// Analyzer for a single file
void data_robustanalyzer::analyzersinglefile() {

  int totEntries = inputChain->GetEntries();
  cout<<"Total number of entries: "<<totEntries<<endl;

  // Count events passing certain selections
  int nosel=0, basicsel=0, sel2=0, sel3=0, sel20=0, sel100=0;

  // Define the histograms
  addhist("nosel");
  addhist("basicsel");
  addhist("sel2"); // Selection to cross-check the trigger filter selection
  addhist("sel3"); // Selection to cross-check the parent triggers
  addhist("sel20"); // mumu DeltaR>1
  addhist("sel100"); // mumu DeltaR>1, mumuM<84 or mumuM>98

  // Add the histograms that are defined only once
  addhistonce();

  // vector of mu indices
  vector<int> noselmuidx;
  vector<int> basicselmuidx;
  vector<int> sel2muidx;
  vector<int> sel3muidx;
  vector<int> sel20muidx;
  vector<int> sel100muidx;

  // Loop beginning on events
  for(unsigned int event=0; event<totEntries; event++) {

    inputChain->GetEntry(event);
    //if(event>20) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    if(muRecoN>0) { // Do not go in here if there is not atleast one reco muon
      
      // Sort the muon objects based on their pT
      vector<int> sortedmuidx(muRecoN);
      iota(begin(sortedmuidx), end(sortedmuidx), 0);
      sort(&sortedmuidx[0], muRecoPt, muRecoN); // Verified that the algorithm works fine

      bool basicselmu = false;
      bool sel3mu = false;
      bool sel100mu = false;

      // Loop beginning on muon reco objects
      for(unsigned int muidx=0; muidx<muRecoN; muidx++) {

	unsigned int idx = sortedmuidx[muidx];
	noselmuidx.push_back(idx);

	basicselmu = true;
	basicselmu *= (TMath::Abs(muRecoEta[idx])<2.5);
	basicselmu *= (muRecoPt[idx]>=16);
	if(basicselmu) basicselmuidx.push_back(idx);

	sel3mu = true;
	sel3mu *= (muRecoPt[idx]>=33);
	sel3mu *= (TMath::Abs(muRecoEta[idx])<2.5);
	sel3mu *= (TMath::Abs(muRecoDxy[idx])>=0.01);
	if(sel3mu) sel3muidx.push_back(idx);
	
	sel100mu = true;
	sel100mu *= (muRecoPt[idx]>=20);
	sel100mu *= (TMath::Abs(muRecoEta[idx])<2.5);
	sel100mu *= (TMath::Abs(muRecoDxy[idx])>=0.01);
	sel100mu *= (TMath::Abs(muRecoDxySig[idx])>=1);
	if(sel100mu) sel100muidx.push_back(idx);
	
      } // End of loop on muon reco objects

      sel2muidx = basicselmuidx;
      sel20muidx = basicselmuidx;
      
    } // End of loop requiring one muon reco object

    if(muFiltN!=sel2muidx.size()) cout<<event<<": ***********Error! mis-match in filter and reco objects for mu16**********"<<endl;
    if(muFiltN33!=sel3muidx.size()) cout<<event<<": ***********Error! mis-match in filter and reco objects for mu33**********"<<endl;

    // Event level requirements on the selected muon objects
    bool sel20mumu = false;
    for(unsigned int idx=0; idx<sel20muidx.size(); idx++) {
      for(unsigned int idx2=idx+1; idx2<sel20muidx.size(); idx2++) {
	TLorentzVector muidx, muidx2;
	muidx.SetPtEtaPhiM(muRecoPt[sel20muidx[idx]],muRecoEta[sel20muidx[idx]],muRecoPhi[sel20muidx[idx]],0.1);
	muidx2.SetPtEtaPhiM(muRecoPt[sel20muidx[idx2]],muRecoEta[sel20muidx[idx2]],muRecoPhi[sel20muidx[idx2]],0.1);
	if(muidx.DeltaR(muidx2)>1) {
	  sel20mumu = true;
	}
      }
    }
    
    bool sel100mumu = false;
    if(sel100muidx.size()>=2) {
      TLorentzVector muidx, muidx2;
      muidx.SetPtEtaPhiM(muRecoPt[sel100muidx[0]],muRecoEta[sel100muidx[0]],muRecoPhi[sel100muidx[0]],0.1);
      muidx2.SetPtEtaPhiM(muRecoPt[sel100muidx[1]],muRecoEta[sel100muidx[1]],muRecoPhi[sel100muidx[1]],0.1);
      double invM = (muidx+muidx2).M();
      if(muidx.DeltaR(muidx2)>1 && (invM<80 || invM>100)) {
	sel100mumu = true;
      }
    }
    
    // Fill the histograms that are defined only once
    fillhistineventonce();

    // Fill histograms for the selections
    fillhistinevent("nosel", noselmuidx);
    fillhistinevent("basicsel", basicselmuidx);
    
    // Count events passing selections
    if(noselmuidx.size()>0) nosel++;
    if(basicselmuidx.size()>0) basicsel++;
    if(sel2muidx.size()>=2) {
      fillhistinevent("sel2", sel2muidx);
      sel2++;
    }
    if(sel3muidx.size()>=2) {
      fillhistinevent("sel3", sel3muidx);
      sel3++;
    }
    if(sel20mumu) {
      fillhistinevent("sel20", sel20muidx);
      sel20++;
    }
    if(sel100mumu) {
      fillhistinevent("sel100", sel100muidx);
      sel100++;
    }
        
    // Clear all the vectors
    noselmuidx.clear();
    basicselmuidx.clear();
    sel2muidx.clear();
    sel3muidx.clear();
    sel20muidx.clear();
    sel100muidx.clear();

  } // End of loop on events

  cout<<"With parent trigger: "<<sel3<<endl;
  cout<<totEntries<<"\t"<<nosel<<"\t"<<basicsel<<"\t"<<sel2<<"\t"<<sel3<<"\t"<<sel20<<"\t"<<sel100<<endl;
}

// Function to fill a set of histograms in the event
void data_robustanalyzer::fillhistinevent(TString selection, vector<int> muidx) {
  
  TH1F* mult = (TH1F*) outfile->Get(selection+"recomu_mult");

  TH1F* allpt = (TH1F*) outfile->Get(selection+"recomu_allpt");
  TH1F* alleta = (TH1F*) outfile->Get(selection+"recomu_alleta");
  TH1F* allphi = (TH1F*) outfile->Get(selection+"recomu_allphi");
  TH1F* alldxy = (TH1F*) outfile->Get(selection+"recomu_alldxy");
  TH1F* alllog10dxy = (TH1F*) outfile->Get(selection+"recomu_alllog10dxy");
  TH1F* alldxysig = (TH1F*) outfile->Get(selection+"recomu_alldxysig");
  TH1F* alllog10dxysig = (TH1F*) outfile->Get(selection+"recomu_alllog10dxysig");

  TH1F* leadpt = (TH1F*) outfile->Get(selection+"recomu_leadpt");
  TH1F* leadeta = (TH1F*) outfile->Get(selection+"recomu_leadeta");
  TH1F* leadphi = (TH1F*) outfile->Get(selection+"recomu_leadphi");
  TH1F* leaddxy = (TH1F*) outfile->Get(selection+"recomu_leaddxy");
  TH1F* leadlog10dxy = (TH1F*) outfile->Get(selection+"recomu_leadlog10dxy");
  TH1F* leaddxysig = (TH1F*) outfile->Get(selection+"recomu_leaddxysig");
  TH1F* leadlog10dxysig = (TH1F*) outfile->Get(selection+"recomu_leadlog10dxysig");

  TH1F* subleadpt = (TH1F*) outfile->Get(selection+"recomu_subleadpt");
  TH1F* subleadeta = (TH1F*) outfile->Get(selection+"recomu_subleadeta");
  TH1F* subleadphi = (TH1F*) outfile->Get(selection+"recomu_subleadphi");
  TH1F* subleaddxy = (TH1F*) outfile->Get(selection+"recomu_subleaddxy");
  TH1F* subleadlog10dxy = (TH1F*) outfile->Get(selection+"recomu_subleadlog10dxy");
  TH1F* subleaddxysig = (TH1F*) outfile->Get(selection+"recomu_subleaddxysig");
  TH1F* subleadlog10dxysig = (TH1F*) outfile->Get(selection+"recomu_subleadlog10dxysig");

  TH1F* leadsubleaddR = (TH1F*) outfile->Get(selection+"recomumu_leadsubleaddR");
  TH1F* leadsubleadM = (TH1F*) outfile->Get(selection+"recomumu_leadsubleadM");

  mult->Fill(muidx.size());
  for(unsigned int idx=0; idx<muidx.size(); idx++) {
    allpt->Fill(muRecoPt[muidx[idx]]);
    alleta->Fill(muRecoEta[muidx[idx]]);
    allphi->Fill(muRecoPhi[muidx[idx]]);
    alldxy->Fill(muRecoDxy[muidx[idx]]);
    alllog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[muidx[idx]])));
    alldxysig->Fill(muRecoDxySig[muidx[idx]]);
    alllog10dxysig->Fill(TMath::Log10(TMath::Abs(muRecoDxySig[muidx[idx]])));
  }

  if(muidx.size()>=1) {
    leadpt->Fill(muRecoPt[muidx[0]]);
    leadeta->Fill(muRecoEta[muidx[0]]);
    leadphi->Fill(muRecoPhi[muidx[0]]);
    leaddxy->Fill(muRecoDxy[muidx[0]]);
    leadlog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[muidx[0]])));
    leaddxysig->Fill(muRecoDxySig[muidx[0]]);
    leadlog10dxysig->Fill(TMath::Log10(TMath::Abs(muRecoDxySig[muidx[0]])));
  }
  
  if(muidx.size()>=2) {
    subleadpt->Fill(muRecoPt[muidx[1]]);
    subleadeta->Fill(muRecoEta[muidx[1]]);
    subleadphi->Fill(muRecoPhi[muidx[1]]);
    subleaddxy->Fill(muRecoDxy[muidx[1]]);
    subleadlog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[muidx[1]])));
    subleaddxysig->Fill(muRecoDxySig[muidx[1]]);
    subleadlog10dxysig->Fill(TMath::Log10(TMath::Abs(muRecoDxySig[muidx[1]])));
    
    TLorentzVector leadmuon, subleadmuon;
    leadmuon.SetPtEtaPhiM(muRecoPt[muidx[0]],muRecoEta[muidx[0]],muRecoPhi[muidx[0]],0.105);
    subleadmuon.SetPtEtaPhiM(muRecoPt[muidx[1]],muRecoEta[muidx[1]],muRecoPhi[muidx[1]],0.105);
    leadsubleaddR->Fill(leadmuon.DeltaR(subleadmuon));
    leadsubleadM->Fill((leadmuon+subleadmuon).M());
  }
  
}

// Function to fill a set of histograms in the event
void data_robustanalyzer::fillhistineventonce() {

  TH1F* filt16mult = (TH1F*) outfile->Get("filt16mu_mult");
  TH1F* filt16pt = (TH1F*) outfile->Get("filt16mu_allpt");
  TH1F* filt16eta = (TH1F*) outfile->Get("filt16mu_alleta");
  TH1F* filt16phi = (TH1F*) outfile->Get("filt16mu_allphi");

  TH1F* filt33mult = (TH1F*) outfile->Get("filt33mu_mult");
  TH1F* filt33pt = (TH1F*) outfile->Get("filt33mu_allpt");
  TH1F* filt33eta = (TH1F*) outfile->Get("filt33mu_alleta");
  TH1F* filt33phi = (TH1F*) outfile->Get("filt33mu_allphi");

  if(muFiltN>0) {
    filt16mult->Fill(muFiltN);
    for(unsigned int ctr=0; ctr<muFiltN; ctr++) {
      filt16pt->Fill(muFiltPt[ctr]);
      filt16eta->Fill(muFiltEta[ctr]);
      filt16phi->Fill(muFiltPhi[ctr]);
    }
  }
  if(muFiltN33>0) {
    filt33mult->Fill(muFiltN33);
    for(unsigned int ctr=0; ctr<muFiltN33; ctr++) {
      filt33pt->Fill(muFiltPt33[ctr]);
      filt33eta->Fill(muFiltEta33[ctr]);
      filt33phi->Fill(muFiltPhi33[ctr]);
    }
  }

}

// Function to add a set of histograms for a selection
void data_robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"recomu_mult","#mu multiplicity",50,-5,45));

  all1dhists.push_back(new TH1F(selection+"recomu_allpt","#mu p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recomu_alleta","#mu #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"recomu_allphi","#mu #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"recomu_alldxy","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"recomu_alllog10dxy","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"recomu_alldxysig","#mu d_{0} sig.",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"recomu_alllog10dxysig","#mu log_{10}(d_{0} sig.)",1000,-5,5));

  all1dhists.push_back(new TH1F(selection+"recomu_leadpt","#mu_{1} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recomu_leadeta","#mu_{1} #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"recomu_leadphi","#mu_{1} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"recomu_leaddxy","#mu_{1} d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"recomu_leadlog10dxy","#mu_{1} log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"recomu_leaddxysig","#mu_{1} d_{0} sig.",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"recomu_leadlog10dxysig","#mu_{1} log_{10}(d_{0} sig.)",1000,-5,5));

  all1dhists.push_back(new TH1F(selection+"recomu_subleadpt","#mu_{2} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recomu_subleadeta","#mu_{2} #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"recomu_subleadphi","#mu_{2} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"recomu_subleaddxy","#mu_{2} d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"recomu_subleadlog10dxy","#mu_{2} log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"recomu_subleaddxysig","#mu_{2} d_{0} sig.",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"recomu_subleadlog10dxysig","#mu_{2} log_{10}(d_{0} sig.)",1000,-5,5));

  all1dhists.push_back(new TH1F(selection+"recomumu_leadsubleaddR","#DeltaR(#mu_{1},#mu_{2})",1000,-1,9));
  all1dhists.push_back(new TH1F(selection+"recomumu_leadsubleadM","M(#mu_{1},#mu_{2})",1000,-50,550));
  
}

// Function to add a set of histograms once
void data_robustanalyzer::addhistonce() {

  all1dhists.push_back(new TH1F("filt16mu_mult","#mu multiplicity",50,-5,45));
  all1dhists.push_back(new TH1F("filt16mu_allpt","#mu p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F("filt16mu_alleta","#mu #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F("filt16mu_allphi","#mu #phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F("filt33mu_mult","#mu multiplicity",50,-5,45));
  all1dhists.push_back(new TH1F("filt33mu_allpt","#mu p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F("filt33mu_alleta","#mu #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F("filt33mu_allphi","#mu #phi",66,-3.3,3.3));

}

// Function to sort the indices based on a factor (Usually pT)
void data_robustanalyzer::sort(int* idx, double* factor, int n) {
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

/*
    if(muRecoN>0||muFiltN>0) {
      cout<<"Event :"<<event<<endl;
      cout<<"Reco objects"<<endl;
      for(unsigned int idx=0; idx<muRecoN; idx++) {
	cout<<muRecoPt[idx]<<"\t"<<muRecoEta[idx]<<"\t"<<muRecoPhi[idx]<<" || ";
      }
      cout<<endl;
      cout<<"Filt objects"<<endl;
      for(unsigned int idx=0; idx<muFiltN; idx++) {
	cout<<muFiltPt[idx]<<"\t"<<muFiltEta[idx]<<"\t"<<muFiltPhi[idx]<<" || ";
      }
      cout<<endl;
      cout<<"Filt 33 objects"<<endl;
      for(unsigned int idx=0; idx<muFiltN33; idx++) {
	cout<<muFiltPt33[idx]<<"\t"<<muFiltEta33[idx]<<"\t"<<muFiltPhi33[idx]<<" || ";
      }
      cout<<endl;
    }
    
*/
