#include "data_robustanalyzer.hh"
#include <iostream>
#include <numeric>

#include "TMath.h"
#include "TVector3.h"

using namespace std;

// Initialize and open the root file in the constructor
data_robustanalyzer::data_robustanalyzer(TString filename, TString outfilename, bool issimu){

  isMC = issimu;
  
  inputChain = new TChain("events");
  cout<<"Initializing for file: "<<filename<<endl;
  inputChain->Add(filename);

  inputChain->SetBranchAddress("trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL", &trigMu38Eg38);
  inputChain->SetBranchAddress("trig_HLT_Mu16NoFiltersNoVtxDisplaced_Photon20_CaloIdL", &trigMu16Eg20);
  inputChain->SetBranchAddress("muFiltn", &muFiltN);
  inputChain->SetBranchAddress("muFilt_pt", &muFiltPt);
  inputChain->SetBranchAddress("muFilt_eta", &muFiltEta);
  inputChain->SetBranchAddress("muFilt_phi", &muFiltPhi);
  inputChain->SetBranchAddress("muFiltn_38", &muFiltN_38);
  inputChain->SetBranchAddress("muFilt_pt_38", &muFiltPt_38);
  inputChain->SetBranchAddress("muFilt_eta_38", &muFiltEta_38);
  inputChain->SetBranchAddress("muFilt_phi_38", &muFiltPhi_38);
  inputChain->SetBranchAddress("l1egObjn", &l1egObjN);
  inputChain->SetBranchAddress("l1egObj_pt", &l1egObjPt);
  inputChain->SetBranchAddress("l1egObj_eta", &l1egObjEta);
  inputChain->SetBranchAddress("l1egObj_phi", &l1egObjPhi);
  inputChain->SetBranchAddress("l1Filtn", &l1FiltN);
  inputChain->SetBranchAddress("l1Filt_pt", &l1FiltPt);
  inputChain->SetBranchAddress("l1Filt_eta", &l1FiltEta);
  inputChain->SetBranchAddress("l1Filt_phi", &l1FiltPhi);
  inputChain->SetBranchAddress("mun", &muRecoN);
  inputChain->SetBranchAddress("mu_pt", &muRecoPt);
  inputChain->SetBranchAddress("mu_eta", &muRecoEta);
  inputChain->SetBranchAddress("mu_phi", &muRecoPhi);
  inputChain->SetBranchAddress("mu_dxy", &muRecoDxy);
  inputChain->SetBranchAddress("mu_dxy_sig", &muRecoDxySig);
  inputChain->SetBranchAddress("egFiltn", &egFiltN);
  inputChain->SetBranchAddress("egFilt_pt", &egFiltPt);
  inputChain->SetBranchAddress("egFilt_eta", &egFiltEta);
  inputChain->SetBranchAddress("egFilt_phi", &egFiltPhi);
  inputChain->SetBranchAddress("egFiltn_38", &egFiltN_38);
  inputChain->SetBranchAddress("egFilt_pt_38", &egFiltPt_38);
  inputChain->SetBranchAddress("egFilt_eta_38", &egFiltEta_38);
  inputChain->SetBranchAddress("egFilt_phi_38", &egFiltPhi_38);
  inputChain->SetBranchAddress("egn", &egRecoN);
  inputChain->SetBranchAddress("egptarr", &egRecoPt);
  inputChain->SetBranchAddress("egetaarr", &egRecoEta);
  inputChain->SetBranchAddress("egphiarr", &egRecoPhi);
  inputChain->SetBranchAddress("eghltEgammaClusterShapearr", &eghltEgammaClusterShape);
  inputChain->SetBranchAddress("eghltEgammaClusterShape_sigmaIEtaIEta5x5arr", &eghltEgammaClusterShape_sigmaIEtaIEta5x5);
  inputChain->SetBranchAddress("eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleanedarr", &eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned);
  inputChain->SetBranchAddress("eghltEgammaEcalPFClusterIsoarr", &eghltEgammaEcalPFClusterIso);
  inputChain->SetBranchAddress("eghltEgammaHcalPFClusterIsoarr", &eghltEgammaHcalPFClusterIso);
  inputChain->SetBranchAddress("eghltEgammaHoverEarr", &eghltEgammaHoverE);
  inputChain->SetBranchAddress("eghltEgammaPixelMatchVars_s2arr", &eghltEgammaPixelMatchVars_s2);
  inputChain->SetBranchAddress("eghltEgammaSuperClusterEnergyarr", &eghltEgammaSuperClusterEnergy);
  inputChain->SetBranchAddress("eghltEgammaEleGsfTrackIsoarr", &eghltEgammaEleGsfTrackIso);
  inputChain->SetBranchAddress("eghltEgammaGsfTrackVars_Chi2arr", &eghltEgammaGsfTrackVars_Chi2);
  inputChain->SetBranchAddress("eghltEgammaGsfTrackVars_Detaarr", &eghltEgammaGsfTrackVars_Deta);
  inputChain->SetBranchAddress("eghltEgammaGsfTrackVars_DetaSeedarr", &eghltEgammaGsfTrackVars_DetaSeed);
  inputChain->SetBranchAddress("eghltEgammaGsfTrackVars_Dphiarr", &eghltEgammaGsfTrackVars_Dphi);
  inputChain->SetBranchAddress("eghltEgammaGsfTrackVars_MissingHitsarr", &eghltEgammaGsfTrackVars_MissingHits);
  inputChain->SetBranchAddress("eghltEgammaGsfTrackVars_NLayerITarr", &eghltEgammaGsfTrackVars_NLayerIT);
  inputChain->SetBranchAddress("eghltEgammaGsfTrackVars_OneOESeedMinusOneOParr", &eghltEgammaGsfTrackVars_OneOESeedMinusOneOP);
  inputChain->SetBranchAddress("eghltEgammaGsfTrackVars_OneOESuperMinusOneOParr", &eghltEgammaGsfTrackVars_OneOESuperMinusOneOP);
  inputChain->SetBranchAddress("egEcalSeedClusterTimearr", &eghltEcalSeedClusterTime);
  //inputChain->SetBranchAddress("eghltEgammaGsfTrackVars_ValidHits", &eghltEgammaGsfTrackVars_ValidHits);
  inputChain->SetBranchAddress("egusn", &egusRecoN);
  inputChain->SetBranchAddress("egusptarr", &egusRecoPt);
  inputChain->SetBranchAddress("egusetaarr", &egusRecoEta);
  inputChain->SetBranchAddress("egusphiarr", &egusRecoPhi);
  inputChain->SetBranchAddress("egushltEgammaClusterShapearr", &egushltEgammaClusterShape);
  inputChain->SetBranchAddress("egushltEgammaClusterShape_sigmaIEtaIEta5x5arr", &egushltEgammaClusterShape_sigmaIEtaIEta5x5);
  inputChain->SetBranchAddress("egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleanedarr", &egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned);
  inputChain->SetBranchAddress("egushltEgammaEcalPFClusterIsoarr", &egushltEgammaEcalPFClusterIso);
  inputChain->SetBranchAddress("egushltEgammaHcalPFClusterIsoarr", &egushltEgammaHcalPFClusterIso);
  inputChain->SetBranchAddress("egushltEgammaHoverEarr", &egushltEgammaHoverE);
  inputChain->SetBranchAddress("egushltEgammaPixelMatchVars_s2arr", &egushltEgammaPixelMatchVars_s2);
  inputChain->SetBranchAddress("egushltEgammaSuperClusterEnergyarr", &egushltEgammaSuperClusterEnergy);
  inputChain->SetBranchAddress("egushltEgammaEleGsfTrackIsoarr", &egushltEgammaEleGsfTrackIso);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_Chi2arr", &egushltEgammaGsfTrackVars_Chi2);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_Detaarr", &egushltEgammaGsfTrackVars_Deta);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_DetaSeedarr", &egushltEgammaGsfTrackVars_DetaSeed);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_Dphiarr", &egushltEgammaGsfTrackVars_Dphi);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_MissingHitsarr", &egushltEgammaGsfTrackVars_MissingHits);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_NLayerITarr", &egushltEgammaGsfTrackVars_NLayerIT);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_OneOESeedMinusOneOParr", &egushltEgammaGsfTrackVars_OneOESeedMinusOneOP);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_OneOESuperMinusOneOParr", &egushltEgammaGsfTrackVars_OneOESuperMinusOneOP);
  inputChain->SetBranchAddress("egusEcalSeedClusterTimearr", &egushltEcalSeedClusterTime);
  //inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_ValidHits", &egushltEgammaGsfTrackVars_ValidHits);

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
  int nosel=0, basicsel=0, sel2=0, sel3=0, sel10=0;

  // Define the histograms
  addhist("nosel");
  addhist("basicsel");
  addhist("sel2"); // Selection to cross-check the trigger filter selection
  addhist("sel3"); // Selection to cross-check the parent triggers
  if(isMC) addhistcomparegenrecounseeded("sel10");
  addhistunseeded("sel10"); // Unseeded e/gamma objects, pT and eta requirements only
  
  // vector of mu indices
  vector<int> noselmuidx;
  vector<int> basicselmuidx;
  vector<int> sel2muidx;
  vector<int> sel3muidx;
  vector<int> sel10muidx;
  // vector of eg indices
  vector<int> noselegidx;
  vector<int> basicselegidx;
  vector<int> sel2egidx;
  vector<int> sel3egidx;
  vector<int> sel10egidx;
  
  // Loop beginning on events
  for(unsigned int event=0; event<totEntries; event++) {

    inputChain->GetEntry(event);
    //if(event>10000) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    if(muRecoN>0) { // Do not go in here if there is not atleast one reco muon

      // Sort the muon objects based on their pT
      vector<int> sortedmuidx(muRecoN);
      iota(begin(sortedmuidx), end(sortedmuidx), 0);
      sort(&sortedmuidx[0], muRecoPt, muRecoN); // Verified that the algorithm works fine

      bool basicselmu = false;
      bool sel3mu = false;

      // Loop beginning on muon reco objects
      for(unsigned int muidx=0; muidx<muRecoN; muidx++) {

	unsigned int idx = sortedmuidx[muidx];
	noselmuidx.push_back(idx);

	basicselmu = true;
	basicselmu *= (TMath::Abs(muRecoEta[idx])<2.5);
	basicselmu *= (muRecoPt[idx]>=16);
	if(basicselmu) basicselmuidx.push_back(idx);

	sel3mu = true;
	sel3mu *= (muRecoPt[idx]>=38);
	sel3mu *= (TMath::Abs(muRecoEta[idx])<2.5);
	sel3mu *= (TMath::Abs(muRecoDxy[idx])>=0.01);
	if(sel3mu) sel3muidx.push_back(idx);
	
      } // End of loop on muon reco objects

    } // End of loop requiring one muon reco object

    sel2muidx = basicselmuidx;
    sel10muidx = basicselmuidx;

    if(muFiltN!=sel2muidx.size()) cout<<"***********Error! mis-match in filter and reco objects for mu16**********"<<endl;
    if(muFiltN_38!=sel3muidx.size()) cout<<"***********Error! mis-match in filter and reco objects for mu38**********"<<endl;
    
    if(egRecoN>0) { // Atleast one reco eg object in the event

      // Sort the egamma objects based on their pT
      vector<int> sortedegidx(egRecoN);
      iota(begin(sortedegidx), end(sortedegidx), 0);
      sort(&sortedegidx[0], egRecoPt, egRecoN); // Verified that the algorithm works fine
      
      vector<int> sortedegusidx(egusRecoN);
      iota(begin(sortedegusidx), end(sortedegusidx), 0);
      sort(&sortedegusidx[0], egusRecoPt, egusRecoN); // Verified that the algorithm works fine

      bool basicseleg = false;
      bool sel2eg = false;
      bool sel3eg = false;
      bool sel10eg = false;
    
      // Loop beginning on egamma reco objects
      for(unsigned int egidx=0; egidx<egRecoN; egidx++) {

	unsigned int idx = sortedegidx[egidx];
	noselegidx.push_back(idx);

	basicseleg = true;
	basicseleg *= (TMath::Abs(egRecoEta[idx])<2.65);
	basicseleg *= (egRecoPt[idx]>=15);
	if(basicseleg) basicselegidx.push_back(idx);

	sel2eg = true;
	sel2eg *= (sel2muidx.size()>0);
	sel2eg *= (egRecoPt[idx]>=20);
	sel2eg *= (TMath::Abs(egRecoEta[idx])<2.65);
	bool sel2l1matchdecision = false;
	for(unsigned int l1cnt=0; l1cnt<l1FiltN; l1cnt++) {
	  bool findl1egobj = false;
	  for(unsigned int l1egcnt=0; l1egcnt<l1egObjN; l1egcnt++) {
	    if(l1egObjEta[l1egcnt]==l1FiltEta[l1cnt]) findl1egobj = true;
	  }
	  if(!findl1egobj) continue;
	  double phibinsize = 1.044;
	  double etabinsize = std::abs(egRecoEta[idx])<1.4791?0.522:1.0;
	  double etabinlow = l1FiltEta[l1cnt]-0.5*etabinsize;
	  double etabinhigh = etabinlow+etabinsize;
	  
	  double deltaphi = std::abs(egRecoPhi[idx]-l1FiltPhi[l1cnt]);
	  double pi = std::acos(-1);
	  deltaphi = (deltaphi>2*pi)?deltaphi-2*pi:deltaphi;
	  deltaphi = (deltaphi>pi)?2*pi-deltaphi:deltaphi;
	  if(egRecoEta[idx]<etabinhigh && egRecoEta[idx]>etabinlow && deltaphi<0.5*phibinsize) {
	    sel2l1matchdecision = true;
	  }
	}
	sel2eg *= sel2l1matchdecision;
	sel2eg *= abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5[idx]<0.014:eghltEgammaClusterShape_sigmaIEtaIEta5x5[idx]<0.035;
	sel2eg *= abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.15*eghltEgammaSuperClusterEnergy[idx]:eghltEgammaHoverE[idx]<0.1*eghltEgammaSuperClusterEnergy[idx];
	if(sel2eg) sel2egidx.push_back(idx);

	sel3eg = true;
	sel3eg *= (sel3muidx.size()>0);
	sel3eg *= (egRecoPt[idx]>=38);
	sel3eg *= (TMath::Abs(egRecoEta[idx])<2.65);
	bool sel3l1matchdecision = false;
	for(unsigned int l1cnt=0; l1cnt<l1FiltN; l1cnt++) {
	  bool findl1egobj = false;
	  for(unsigned int l1egcnt=0; l1egcnt<l1egObjN; l1egcnt++) {
	    if(l1egObjEta[l1egcnt]==l1FiltEta[l1cnt]) findl1egobj = true;
	  }
	  if(!findl1egobj) continue;
	  double phibinsize = 1.044;
	  double etabinsize = std::abs(egRecoEta[idx])<1.4791?0.522:1.0;
	  double etabinlow = l1FiltEta[l1cnt]-0.5*etabinsize;
	  double etabinhigh = etabinlow+etabinsize;
	  
	  double deltaphi = std::abs(egRecoPhi[idx]-l1FiltPhi[l1cnt]);
	  double pi = std::acos(-1);
	  deltaphi = (deltaphi>2*pi)?deltaphi-2*pi:deltaphi;
	  deltaphi = (deltaphi>pi)?2*pi-deltaphi:deltaphi;
	  if(egRecoEta[idx]<etabinhigh && egRecoEta[idx]>etabinlow && deltaphi<0.5*phibinsize) {
	    sel3l1matchdecision = true;
	  }
	}
	sel3eg *= sel3l1matchdecision;
	sel3eg *= abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.014:eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.035;
	sel3eg *= abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.15*eghltEgammaSuperClusterEnergy[idx]:eghltEgammaHoverE[idx]<0.1*eghltEgammaSuperClusterEnergy[idx];
	if(sel3eg) sel3egidx.push_back(idx);

      } // End of loop on egamma reco objects
            
      // Loop beginning on unseeded egamma reco objects
      for(unsigned int egidx=0; egidx<egusRecoN; egidx++) {
	
	unsigned int idx = sortedegusidx[egidx];
	
	sel10eg = true;
	sel10eg *= (sel10muidx.size()>0);
	sel10eg *= (egusRecoPt[idx]>=20);
	sel10eg *= (TMath::Abs(egusRecoEta[idx])<2.65);
	if(sel10eg) sel10egidx.push_back(idx);
	
      } // End of loop on unseeded egamma objects
      
      fillhistinevent("nosel", noselegidx);
      fillhistinevent("basicsel", basicselegidx);
      fillhistinevent("sel2", sel2egidx);
      fillhistinevent("sel3", sel3egidx);
      if(isMC && sel10egidx.size()>0) fillhistcomparegenrecounseeded("sel10", sel10egidx);
      fillhistineventunseeded("sel10", sel10egidx);
      
    } // End of condition requiring atleast one egReco object

    // If the trigger menu flow is inappropriately modelled in code, then print to undestand the cause
    if(egFiltN!=sel2egidx.size()) cout<<event<<" : ***********Error! mis-match filter and reco selection eg20**********"<<endl;
    if(egFiltN_38!=sel3egidx.size()) cout<<event<<"-----------Error! mis-match filter and reco selection eg38------------"<<endl;

    // Count events passing selections
    if(noselmuidx.size()>0 && noselegidx.size()>0) nosel++;
    if(basicselmuidx.size()>0 && basicselegidx.size()>0) basicsel++;
    if(sel2muidx.size()>0 && sel2egidx.size()>0) sel2++;
    if(sel3muidx.size()>0 && sel3egidx.size()>0) sel3++;
    if(sel10muidx.size()>0 && sel10egidx.size()>0) sel10++;
    
    // Clear all the vectors
    noselmuidx.clear();
    basicselmuidx.clear();
    sel2muidx.clear();
    sel3muidx.clear();
    sel10muidx.clear();
    noselegidx.clear();
    basicselegidx.clear();
    sel2egidx.clear();
    sel3egidx.clear();
    sel10egidx.clear();

  } // End of loop on events

  cout<<"With parent trigger: "<<sel3<<endl;
  cout<<totEntries<<"\t"<<nosel<<"\t"<<basicsel<<"\t"<<sel2<<"\t"<<sel3<<"\t"<<sel10<<endl;
}

// Function to do gen matching
int data_robustanalyzer::doGenMatchingUnseeded(vector<int> egusidx) {

  if(!isMC) {
    cout<<"Error in doing gen matching. Not MC file."<<endl;
    return -1;
  }
  int genmuidx=-1, genegidx=-1, mun=0, egn=0;
  for(unsigned int gen=0; gen<genLepN; gen++) {
    if(abs(genLepMomPid[gen])==9000007) {
      if(abs(genLepPid[gen])==11) {
	egn++;
	genegidx = gen;
      }
      if(abs(genLepPid[gen])==13) {
	mun++;
	genmuidx = gen;
      }
    }
  }

  // Require gen signature with 1mu and 1e
  int genmcheguspos = -1;
  //double mindiffpt = -1;
  if(mun==1 && egn==1 && egusidx.size()>0) {
    // Find the egusidx with the best angular match
    for(unsigned int eg=0; eg<egusidx.size(); eg++) {
      // Condition to do gen match and choose the closest in pT
      if(abs(egusRecoEta[egusidx[eg]])<1.479) {
	//if(abs(egusRecoEta[egusidx[eg]]-genLepEta[genegidx])<0.06 && (mindiffpt==-1 || abs(egusRecoPt[egusidx[eg]]-genLepPt[genegidx])<mindiffpt) ) {
	if(abs(egusRecoEta[egusidx[eg]]-genLepEta[genegidx])<0.06 ) {
	  genmcheguspos = egusidx[eg];
	  //mindiffpt = abs(egusRecoPt[egusidx[eg]]-genLepPt[genegidx]);
	  break;
	}
      }
      else {
	//if(abs(egusRecoEta[egusidx[eg]]-genLepEta[genegidx])<0.04 && (mindiffpt==-1 || abs(egusRecoPt[egusidx[eg]]-genLepPt[genegidx])<mindiffpt) ) {
	if(abs(egusRecoEta[egusidx[eg]]-genLepEta[genegidx])<0.04) {
	  genmcheguspos = egusidx[eg];
	  break;
	  //mindiffpt = abs(egusRecoPt[egusidx[eg]]-genLepPt[genegidx]);
	}
      }
    }
  }
  else {
    return -1;
  }
  
  return genmcheguspos;
}

// Function to fill a set of histograms in the event
void data_robustanalyzer::fillhistinevent(TString selection, vector<int> egidx) {

  // nothing here for now
  TH1F* egmult = (TH1F*) outfile->Get(selection+"recoeg_egmult");
  TH1F* egpt = (TH1F*) outfile->Get(selection+"recoeg_egpt");
  TH1F* egeta = (TH1F*) outfile->Get(selection+"recoeg_egeta");
  TH1F* egphi = (TH1F*) outfile->Get(selection+"recoeg_egphi");

  // Get barrel variables
  TH1F* recoeb_egclustershape = (TH1F*) outfile->Get(selection+"recoeb_egclustershape");
  TH1F* recoeb_in5x5clusshape = (TH1F*) outfile->Get(selection+"recoeb_in5x5clusshape");
  TH1F* recoeb_in5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoeb_in5x5noiseclnd");
  TH1F* recoeb_scenergy = (TH1F*) outfile->Get(selection+"recoeb_scenergy");
  TH1F* recoeb_hovere = (TH1F*) outfile->Get(selection+"recoeb_hovere");
  TH1F* recoeb_hovereoversupcluse = (TH1F*) outfile->Get(selection+"recoeb_hovereoversupcluse");
  TH1F* recoeb_ecalpfclustiso = (TH1F*) outfile->Get(selection+"recoeb_ecalpfclustiso");
  TH1F* recoeb_hcalpfclustiso = (TH1F*) outfile->Get(selection+"recoeb_hcalpfclustiso");
  TH1F* recoeb_pixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoeb_pixelmchvar_s2");
  TH1F* recoeb_trkiso = (TH1F*) outfile->Get(selection+"recoeb_trkiso");
  TH1F* recoeb_chi2 = (TH1F*) outfile->Get(selection+"recoeb_chi2");
  TH1F* recoeb_deta = (TH1F*) outfile->Get(selection+"recoeb_deta");
  TH1F* recoeb_detaseed = (TH1F*) outfile->Get(selection+"recoeb_detaseed");
  TH1F* recoeb_dphi = (TH1F*) outfile->Get(selection+"recoeb_dphi");
  TH1F* recoeb_mhits = (TH1F*) outfile->Get(selection+"recoeb_mhits");
  TH1F* recoeb_nlayerit = (TH1F*) outfile->Get(selection+"recoeb_nlayerit");
  TH1F* recoeb_ooeseedoop = (TH1F*) outfile->Get(selection+"recoeb_ooeseedoop");
  TH1F* recoeb_ooesclsoop = (TH1F*) outfile->Get(selection+"recoeb_ooesclsoop");
  TH1F* recoeb_valhits = (TH1F*) outfile->Get(selection+"recoeb_valhits");  
  TH1F* recoeb_seedclustime = (TH1F*) outfile->Get(selection+"recoeb_seedclustime");  

  // Get end-cap variables
  TH1F* recoee_egclustershape = (TH1F*) outfile->Get(selection+"recoee_egclustershape");
  TH1F* recoee_in5x5clusshape = (TH1F*) outfile->Get(selection+"recoee_in5x5clusshape");
  TH1F* recoee_in5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoee_in5x5noiseclnd");
  TH1F* recoee_scenergy = (TH1F*) outfile->Get(selection+"recoee_scenergy");
  TH1F* recoee_hovere = (TH1F*) outfile->Get(selection+"recoee_hovere");
  TH1F* recoee_hovereoversupcluse = (TH1F*) outfile->Get(selection+"recoee_hovereoversupcluse");
  TH1F* recoee_ecalpfclustiso = (TH1F*) outfile->Get(selection+"recoee_ecalpfclustiso");
  TH1F* recoee_hcalpfclustiso = (TH1F*) outfile->Get(selection+"recoee_hcalpfclustiso");
  TH1F* recoee_pixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoee_pixelmchvar_s2");
  TH1F* recoee_trkiso = (TH1F*) outfile->Get(selection+"recoee_trkiso");
  TH1F* recoee_chi2 = (TH1F*) outfile->Get(selection+"recoee_chi2");
  TH1F* recoee_deta = (TH1F*) outfile->Get(selection+"recoee_deta");
  TH1F* recoee_detaseed = (TH1F*) outfile->Get(selection+"recoee_detaseed");
  TH1F* recoee_dphi = (TH1F*) outfile->Get(selection+"recoee_dphi");
  TH1F* recoee_mhits = (TH1F*) outfile->Get(selection+"recoee_mhits");
  TH1F* recoee_nlayerit = (TH1F*) outfile->Get(selection+"recoee_nlayerit");
  TH1F* recoee_ooeseedoop = (TH1F*) outfile->Get(selection+"recoee_ooeseedoop");
  TH1F* recoee_ooesclsoop = (TH1F*) outfile->Get(selection+"recoee_ooesclsoop");
  TH1F* recoee_valhits = (TH1F*) outfile->Get(selection+"recoee_valhits");
  TH1F* recoee_seedclustime = (TH1F*) outfile->Get(selection+"recoee_seedclustime");  
  
  if(egidx.size()>0) {
    egmult->Fill(egidx.size());
    egpt->Fill(egRecoPt[egidx[0]]);
    egeta->Fill(egRecoEta[egidx[0]]);
    egphi->Fill(egRecoPhi[egidx[0]]);

    // Fill barrel variables
    if(TMath::Abs(egRecoEta[egidx[0]])<=1.479) {
      recoeb_egclustershape->Fill(eghltEgammaClusterShape[egidx[0]]);
      recoeb_in5x5clusshape->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[0]]);
      recoeb_in5x5noiseclnd->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[0]]);
      //recoeb_scenergy->Fill(eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_hovere->Fill(eghltEgammaHoverE[egidx[0]]);
      recoeb_hovereoversupcluse->Fill(eghltEgammaHoverE[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_ecalpfclustiso->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]);
      recoeb_hcalpfclustiso->Fill(eghltEgammaHcalPFClusterIso[egidx[0]]);
      if(eghltEcalSeedClusterTime[egidx[0]]!=0) recoeb_seedclustime->Fill(eghltEcalSeedClusterTime[egidx[0]]);
      if(eghltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoeb_pixelmchvar_s2->Fill(eghltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoeb_pixelmchvar_s2->Fill(-5);
      }
      if(eghltEgammaGsfTrackVars_Deta[egidx[0]]<999999) {
	recoeb_trkiso->Fill(eghltEgammaEleGsfTrackIso[egidx[0]]);
	recoeb_chi2->Fill(eghltEgammaGsfTrackVars_Chi2[egidx[0]]);
	recoeb_deta->Fill(eghltEgammaGsfTrackVars_Deta[egidx[0]]);
	recoeb_detaseed->Fill(eghltEgammaGsfTrackVars_DetaSeed[egidx[0]]);
	recoeb_dphi->Fill(eghltEgammaGsfTrackVars_Dphi[egidx[0]]);
	recoeb_mhits->Fill(eghltEgammaGsfTrackVars_MissingHits[egidx[0]]);
	recoeb_nlayerit->Fill(eghltEgammaGsfTrackVars_NLayerIT[egidx[0]]);
	recoeb_ooeseedoop->Fill(eghltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[0]]);
	recoeb_ooesclsoop->Fill(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[0]]);
	//recoeb_valhits->Fill(eghltEgammaGsfTrackVars_ValidHits[egidx[0]]);
      }
      else {
	recoeb_trkiso->Fill(-100);
	recoeb_chi2->Fill(-20);
	recoeb_deta->Fill(-5);
	recoeb_detaseed->Fill(-5);
	recoeb_dphi->Fill(-5);
	recoeb_mhits->Fill(-15);
	recoeb_nlayerit->Fill(-15);
	recoeb_ooeseedoop->Fill(-50);
	recoeb_ooesclsoop->Fill(-50);
	//recoeb_valhits->Fill(-15);
      }
    } // End of filling barrel variables
    else { // Fill end-cap variables
      recoee_egclustershape->Fill(eghltEgammaClusterShape[egidx[0]]);
      recoee_in5x5clusshape->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[0]]);
      recoee_in5x5noiseclnd->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[0]]);
      //recoee_scenergy->Fill(eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_hovere->Fill(eghltEgammaHoverE[egidx[0]]);
      recoee_hovereoversupcluse->Fill(eghltEgammaHoverE[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_ecalpfclustiso->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]);
      recoee_hcalpfclustiso->Fill(eghltEgammaHcalPFClusterIso[egidx[0]]);	
      if(eghltEcalSeedClusterTime[egidx[0]]!=0) recoee_seedclustime->Fill(eghltEcalSeedClusterTime[egidx[0]]);
      if(eghltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoee_pixelmchvar_s2->Fill(eghltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoee_pixelmchvar_s2->Fill(-5);
      }
      if(eghltEgammaGsfTrackVars_Deta[egidx[0]]<999999) {
	recoee_trkiso->Fill(eghltEgammaEleGsfTrackIso[egidx[0]]);
	recoee_chi2->Fill(eghltEgammaGsfTrackVars_Chi2[egidx[0]]);
	recoee_deta->Fill(eghltEgammaGsfTrackVars_Deta[egidx[0]]);
	recoee_detaseed->Fill(eghltEgammaGsfTrackVars_DetaSeed[egidx[0]]);
	recoee_dphi->Fill(eghltEgammaGsfTrackVars_Dphi[egidx[0]]);
	recoee_mhits->Fill(eghltEgammaGsfTrackVars_MissingHits[egidx[0]]);
	recoee_nlayerit->Fill(eghltEgammaGsfTrackVars_NLayerIT[egidx[0]]);
	recoee_ooeseedoop->Fill(eghltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[0]]);
	recoee_ooesclsoop->Fill(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[0]]);
	//recoee_valhits->Fill(eghltEgammaGsfTrackVars_ValidHits[egidx[0]]);
      }
      else {
	recoee_trkiso->Fill(-100);
	recoee_chi2->Fill(-20);
	recoee_deta->Fill(-5);
	recoee_detaseed->Fill(-5);
	recoee_dphi->Fill(-5);
	recoee_mhits->Fill(-15);
	recoee_nlayerit->Fill(-15);
	recoee_ooeseedoop->Fill(-50);
	recoee_ooesclsoop->Fill(-50);
	//recoee_valhits->Fill(-15);
      }
    } // End of filling end-cap variables

  } // End of condition requiring atleast one eg object
  
}

// Function to fill a set of histograms in the event - unseeded egamma objects
void data_robustanalyzer::fillhistineventunseeded(TString selection, vector<int> egidx) {

  // nothing here for now
  TH1F* egmult = (TH1F*) outfile->Get(selection+"recoegus_egmult");
  TH1F* egpt = (TH1F*) outfile->Get(selection+"recoegus_egpt");
  TH1F* egeta = (TH1F*) outfile->Get(selection+"recoegus_egeta");
  TH1F* egphi = (TH1F*) outfile->Get(selection+"recoegus_egphi");

  TH1F* genmchegpt = (TH1F*) outfile->Get(selection+"recoegus_genmchegpt");
  TH1F* genmchegeta = (TH1F*) outfile->Get(selection+"recoegus_genmchegeta");
  TH1F* genmchegphi = (TH1F*) outfile->Get(selection+"recoegus_genmchegphi");

  // Get barrel variables
  TH1F* recoeb_egclustershape = (TH1F*) outfile->Get(selection+"recoebus_egclustershape");
  TH1F* recoeb_in5x5clusshape = (TH1F*) outfile->Get(selection+"recoebus_in5x5clusshape");
  TH1F* recoeb_in5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoebus_in5x5noiseclnd");
  TH1F* recoeb_scenergy = (TH1F*) outfile->Get(selection+"recoebus_scenergy");
  TH1F* recoeb_hovere = (TH1F*) outfile->Get(selection+"recoebus_hovere");
  TH1F* recoeb_hovereoversupcluse = (TH1F*) outfile->Get(selection+"recoebus_hovereoversupcluse");
  TH1F* recoeb_ecalpfclustiso = (TH1F*) outfile->Get(selection+"recoebus_ecalpfclustiso");
  TH1F* recoeb_hcalpfclustiso = (TH1F*) outfile->Get(selection+"recoebus_hcalpfclustiso");
  TH1F* recoeb_pixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoebus_pixelmchvar_s2");
  TH1F* recoeb_trkiso = (TH1F*) outfile->Get(selection+"recoebus_trkiso");
  TH1F* recoeb_chi2 = (TH1F*) outfile->Get(selection+"recoebus_chi2");
  TH1F* recoeb_deta = (TH1F*) outfile->Get(selection+"recoebus_deta");
  TH1F* recoeb_detaseed = (TH1F*) outfile->Get(selection+"recoebus_detaseed");
  TH1F* recoeb_dphi = (TH1F*) outfile->Get(selection+"recoebus_dphi");
  TH1F* recoeb_mhits = (TH1F*) outfile->Get(selection+"recoebus_mhits");
  TH1F* recoeb_nlayerit = (TH1F*) outfile->Get(selection+"recoebus_nlayerit");
  TH1F* recoeb_ooeseedoop = (TH1F*) outfile->Get(selection+"recoebus_ooeseedoop");
  TH1F* recoeb_ooesclsoop = (TH1F*) outfile->Get(selection+"recoebus_ooesclsoop");
  TH1F* recoeb_valhits = (TH1F*) outfile->Get(selection+"recoebus_valhits");  
  TH1F* recoeb_seedclustime = (TH1F*) outfile->Get(selection+"recoebus_seedclustime");  

  // Get end-cap variables
  TH1F* recoee_egclustershape = (TH1F*) outfile->Get(selection+"recoeeus_egclustershape");
  TH1F* recoee_in5x5clusshape = (TH1F*) outfile->Get(selection+"recoeeus_in5x5clusshape");
  TH1F* recoee_in5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoeeus_in5x5noiseclnd");
  TH1F* recoee_scenergy = (TH1F*) outfile->Get(selection+"recoeeus_scenergy");
  TH1F* recoee_hovere = (TH1F*) outfile->Get(selection+"recoeeus_hovere");
  TH1F* recoee_hovereoversupcluse = (TH1F*) outfile->Get(selection+"recoeeus_hovereoversupcluse");
  TH1F* recoee_ecalpfclustiso = (TH1F*) outfile->Get(selection+"recoeeus_ecalpfclustiso");
  TH1F* recoee_hcalpfclustiso = (TH1F*) outfile->Get(selection+"recoeeus_hcalpfclustiso");
  TH1F* recoee_pixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoeeus_pixelmchvar_s2");
  TH1F* recoee_trkiso = (TH1F*) outfile->Get(selection+"recoeeus_trkiso");
  TH1F* recoee_chi2 = (TH1F*) outfile->Get(selection+"recoeeus_chi2");
  TH1F* recoee_deta = (TH1F*) outfile->Get(selection+"recoeeus_deta");
  TH1F* recoee_detaseed = (TH1F*) outfile->Get(selection+"recoeeus_detaseed");
  TH1F* recoee_dphi = (TH1F*) outfile->Get(selection+"recoeeus_dphi");
  TH1F* recoee_mhits = (TH1F*) outfile->Get(selection+"recoeeus_mhits");
  TH1F* recoee_nlayerit = (TH1F*) outfile->Get(selection+"recoeeus_nlayerit");
  TH1F* recoee_ooeseedoop = (TH1F*) outfile->Get(selection+"recoeeus_ooeseedoop");
  TH1F* recoee_ooesclsoop = (TH1F*) outfile->Get(selection+"recoeeus_ooesclsoop");
  TH1F* recoee_valhits = (TH1F*) outfile->Get(selection+"recoeeus_valhits");
  TH1F* recoee_seedclustime = (TH1F*) outfile->Get(selection+"recoeeus_seedclustime");  
  
  if(egidx.size()>0) {
    egmult->Fill(egidx.size());
    egpt->Fill(egRecoPt[egidx[0]]);
    egeta->Fill(egRecoEta[egidx[0]]);
    egphi->Fill(egRecoPhi[egidx[0]]);

    // Fill barrel variables
    if(TMath::Abs(egRecoEta[egidx[0]])<=1.479) {
      recoeb_egclustershape->Fill(eghltEgammaClusterShape[egidx[0]]);
      recoeb_in5x5clusshape->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[0]]);
      recoeb_in5x5noiseclnd->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[0]]);
      //recoeb_scenergy->Fill(eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_hovere->Fill(eghltEgammaHoverE[egidx[0]]);
      recoeb_hovereoversupcluse->Fill(eghltEgammaHoverE[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_ecalpfclustiso->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]);
      recoeb_hcalpfclustiso->Fill(eghltEgammaHcalPFClusterIso[egidx[0]]);
      if(eghltEcalSeedClusterTime[egidx[0]]!=0) recoeb_seedclustime->Fill(eghltEcalSeedClusterTime[egidx[0]]);
      if(eghltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoeb_pixelmchvar_s2->Fill(eghltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoeb_pixelmchvar_s2->Fill(-5);
      }
      if(eghltEgammaGsfTrackVars_Deta[egidx[0]]<999999) {
	recoeb_trkiso->Fill(eghltEgammaEleGsfTrackIso[egidx[0]]);
	recoeb_chi2->Fill(eghltEgammaGsfTrackVars_Chi2[egidx[0]]);
	recoeb_deta->Fill(eghltEgammaGsfTrackVars_Deta[egidx[0]]);
	recoeb_detaseed->Fill(eghltEgammaGsfTrackVars_DetaSeed[egidx[0]]);
	recoeb_dphi->Fill(eghltEgammaGsfTrackVars_Dphi[egidx[0]]);
	recoeb_mhits->Fill(eghltEgammaGsfTrackVars_MissingHits[egidx[0]]);
	recoeb_nlayerit->Fill(eghltEgammaGsfTrackVars_NLayerIT[egidx[0]]);
	recoeb_ooeseedoop->Fill(eghltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[0]]);
	recoeb_ooesclsoop->Fill(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[0]]);
	//recoeb_valhits->Fill(eghltEgammaGsfTrackVars_ValidHits[egidx[0]]);
      }
      else {
	recoeb_trkiso->Fill(-100);
	recoeb_chi2->Fill(-20);
	recoeb_deta->Fill(-5);
	recoeb_detaseed->Fill(-5);
	recoeb_dphi->Fill(-5);
	recoeb_mhits->Fill(-15);
	recoeb_nlayerit->Fill(-15);
	recoeb_ooeseedoop->Fill(-50);
	recoeb_ooesclsoop->Fill(-50);
	//recoeb_valhits->Fill(-15);
      }
    } // End of filling barrel variables
    else { // Fill end-cap variables
      recoee_egclustershape->Fill(eghltEgammaClusterShape[egidx[0]]);
      recoee_in5x5clusshape->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[0]]);
      recoee_in5x5noiseclnd->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[0]]);
      //recoee_scenergy->Fill(eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_hovere->Fill(eghltEgammaHoverE[egidx[0]]);
      recoee_hovereoversupcluse->Fill(eghltEgammaHoverE[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_ecalpfclustiso->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]);
      recoee_hcalpfclustiso->Fill(eghltEgammaHcalPFClusterIso[egidx[0]]);	
      if(eghltEcalSeedClusterTime[egidx[0]]!=0) recoee_seedclustime->Fill(eghltEcalSeedClusterTime[egidx[0]]);
      if(eghltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoee_pixelmchvar_s2->Fill(eghltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoee_pixelmchvar_s2->Fill(-5);
      }
      if(eghltEgammaGsfTrackVars_Deta[egidx[0]]<999999) {
	recoee_trkiso->Fill(eghltEgammaEleGsfTrackIso[egidx[0]]);
	recoee_chi2->Fill(eghltEgammaGsfTrackVars_Chi2[egidx[0]]);
	recoee_deta->Fill(eghltEgammaGsfTrackVars_Deta[egidx[0]]);
	recoee_detaseed->Fill(eghltEgammaGsfTrackVars_DetaSeed[egidx[0]]);
	recoee_dphi->Fill(eghltEgammaGsfTrackVars_Dphi[egidx[0]]);
	recoee_mhits->Fill(eghltEgammaGsfTrackVars_MissingHits[egidx[0]]);
	recoee_nlayerit->Fill(eghltEgammaGsfTrackVars_NLayerIT[egidx[0]]);
	recoee_ooeseedoop->Fill(eghltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[0]]);
	recoee_ooesclsoop->Fill(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[0]]);
	//recoee_valhits->Fill(eghltEgammaGsfTrackVars_ValidHits[egidx[0]]);
      }
      else {
	recoee_trkiso->Fill(-100);
	recoee_chi2->Fill(-20);
	recoee_deta->Fill(-5);
	recoee_detaseed->Fill(-5);
	recoee_dphi->Fill(-5);
	recoee_mhits->Fill(-15);
	recoee_nlayerit->Fill(-15);
	recoee_ooeseedoop->Fill(-50);
	recoee_ooesclsoop->Fill(-50);
	//recoee_valhits->Fill(-15);
      }
    } // End of filling end-cap variables

    if(isMC) {
      int genmchegpos = doGenMatchingUnseeded(egidx);
      if(genmchegpos!=-1) {
	genmchegpt->Fill(egRecoPt[genmchegpos]);
	genmchegeta->Fill(egRecoEta[genmchegpos]);
	genmchegphi->Fill(egRecoPhi[genmchegpos]);
      }
    }
  } // End of condition requiring atleast one eg object
  
}

// Function to fill a set of histograms for a selection - angle between gen and reco. - unseeded egamma objects
void data_robustanalyzer::fillhistcomparegenrecounseeded(TString selection, vector<int> egusidx) {

  TH1F* recoeg_genmult = (TH1F*) outfile->Get(selection+"recoegus_genmult");
  TH1F* recoeg_recomult = (TH1F*) outfile->Get(selection+"recoegus_recomult");
  TH1F* recoeb_gendEta = (TH1F*) outfile->Get(selection+"recoebus_gendEta");
  TH1F* recoeb_gendPhi = (TH1F*) outfile->Get(selection+"recoebus_gendPhi");
  TH1F* recoee_gendEta = (TH1F*) outfile->Get(selection+"recoeeus_gendEta");
  TH1F* recoee_gendPhi = (TH1F*) outfile->Get(selection+"recoeeus_gendPhi");

  if(isMC) {
    int genmuidx=-1, genegidx=-1, mun=0, egn=0;
    for(unsigned int gen=0; gen<genLepN; gen++) {
      if(abs(genLepMomPid[gen])==9000007) {
	if(abs(genLepPid[gen])==11) {
	  egn++;
	  genegidx = gen;
	}
	if(abs(genLepPid[gen])==13) {
	  mun++;
	  genmuidx = gen;
	}
      }
    }
    
    recoeg_genmult->Fill(egn);
    if(mun==1 && egn==1 && egusidx.size()>0) {
      recoeg_recomult->Fill(egusidx.size());
      for(unsigned int egctr=0; egctr<egusidx.size(); egctr++) {
	TVector3 recoeg, gene;
	recoeg.SetPtEtaPhi(egusRecoPt[egusidx[egctr]], egusRecoEta[egusidx[egctr]], egusRecoPhi[egusidx[egctr]]);
	gene.SetPtEtaPhi(genLepPt[genegidx], genLepEta[genegidx], genLepPhi[genegidx]);
	if(TMath::Abs(egusRecoEta[egusidx[egctr]])<1.479) {
	  recoeb_gendEta->Fill(recoeg.Eta()-gene.Eta());
	  recoeb_gendPhi->Fill(recoeg.DeltaPhi(gene));
	}
	else {
	  recoee_gendEta->Fill(recoeg.Eta()-gene.Eta());
	  recoee_gendPhi->Fill(recoeg.DeltaPhi(gene));
	}
      }
    }
  }
  else {
    cout<<"Error in filling angle with gen. Not MC file."<<endl;
  }
}

// Function to add a set of histograms for a selection
void data_robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"recoeg_egmult","reco N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"recoeg_egpt","reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recoeg_egeta","reco e/#gamma #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"recoeg_egphi","reco e/#gamma #phi",66,-3.3,3.3));

  // barrel variables
  all1dhists.push_back(new TH1F(selection+"recoeb_egclustershape","barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_in5x5clusshape","barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_in5x5noiseclnd","barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_hovere","barrel e/#gamma H/E",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeb_scenergy","barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeb_ecalpfclustiso","barrel e/#gamma ecal PF Iso.",10000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeb_hcalpfclustiso","barrel e/#gamma hcal PF Iso.",10000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeb_pixelmchvar_s2","barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeb_chi2","barrel e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoeb_deta","barrel e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeb_detaseed","barrel e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeb_dphi","barrel e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoeb_mhits","barrel e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_nlayerit","barrel e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_ooeseedoop","barrel e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeb_ooesclsoop","barrel e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeb_valhits","barrel e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_trkiso","barrel e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_hovereoversupcluse","barrel e/#gamma (H/E)/supClusEnergy",10000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeb_seedclustime","barrel e/#gamma_{seed} time / ns",20000,-10,10));

  // end-cap variables
  all1dhists.push_back(new TH1F(selection+"recoee_egclustershape","end-cap e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_in5x5clusshape","end-cap e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_in5x5noiseclnd","end-cap e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_hovere","end-cap e/#gamma H/E",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoee_scenergy","end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoee_ecalpfclustiso","end-cap e/#gamma ecal PF Iso.",10000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoee_hcalpfclustiso","end-cap e/#gamma hcal PF Iso.",10000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoee_pixelmchvar_s2","end-cap e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoee_chi2","end-cap e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoee_deta","end-cap e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoee_detaseed","end-cap e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoee_dphi","end-cap e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoee_mhits","end-cap e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_nlayerit","end-cap e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_ooeseedoop","end-cap e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoee_ooesclsoop","end-cap e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoee_valhits","end-cap e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_trkiso","end-cap e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_hovereoversupcluse","end-cap e/#gamma (H/E)/supClusEnergy",10000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoee_seedclustime","end-cap e/#gamma_{seed} time / ns",20000,-10,10));
}

// Function to add a set of histograms for a selection - unseeded egamma objects
void data_robustanalyzer::addhistunseeded(TString selection) {

  all1dhists.push_back(new TH1F(selection+"recoegus_egmult","reco N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"recoegus_egpt","reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recoegus_egeta","reco e/#gamma #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"recoegus_egphi","reco e/#gamma #phi",66,-3.3,3.3));

  if(isMC) {
    all1dhists.push_back(new TH1F(selection+"recoegus_genmchegpt","gen matched reco e/#gamma p_{T} / GeV",550,-50,500));
    all1dhists.push_back(new TH1F(selection+"recoegus_genmchegeta","gen matched reco e/#gamma #eta",52,-2.6,2.6));
    all1dhists.push_back(new TH1F(selection+"recoegus_genmchegphi","gen matched reco e/#gamma #phi",66,-3.3,3.3));
  }

  // barrel variables
  all1dhists.push_back(new TH1F(selection+"recoebus_egclustershape","barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_in5x5clusshape","barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_in5x5noiseclnd","barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_hovere","barrel e/#gamma H/E",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoebus_scenergy","barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoebus_ecalpfclustiso","barrel e/#gamma ecal PF Iso.",10000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoebus_hcalpfclustiso","barrel e/#gamma hcal PF Iso.",10000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoebus_pixelmchvar_s2","barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoebus_chi2","barrel e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoebus_deta","barrel e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoebus_detaseed","barrel e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoebus_dphi","barrel e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoebus_mhits","barrel e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebus_nlayerit","barrel e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebus_ooeseedoop","barrel e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoebus_ooesclsoop","barrel e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoebus_valhits","barrel e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebus_trkiso","barrel e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebus_hovereoversupcluse","barrel e/#gamma (H/E)/supClusEnergy",10000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoebus_seedclustime","barrel e/#gamma_{seed} time / ns",20000,-10,10));

  // end-cap variables
  all1dhists.push_back(new TH1F(selection+"recoeeus_egclustershape","end-cap e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_in5x5clusshape","end-cap e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_in5x5noiseclnd","end-cap e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_hovere","end-cap e/#gamma H/E",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeeus_scenergy","end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeeus_ecalpfclustiso","end-cap e/#gamma ecal PF Iso.",10000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeeus_hcalpfclustiso","end-cap e/#gamma hcal PF Iso.",10000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeeus_pixelmchvar_s2","end-cap e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeeus_chi2","end-cap e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoeeus_deta","end-cap e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeeus_detaseed","end-cap e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeeus_dphi","end-cap e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoeeus_mhits","end-cap e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeus_nlayerit","end-cap e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeus_ooeseedoop","end-cap e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeeus_ooesclsoop","end-cap e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeeus_valhits","end-cap e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeus_trkiso","end-cap e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeus_hovereoversupcluse","end-cap e/#gamma (H/E)/supClusEnergy",10000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeeus_seedclustime","end-cap e/#gamma_{seed} time / ns",20000,-10,10));
}

// Function to add a set of histograms for a selection - angle between gen and reco. - unseeded egamma objects
void data_robustanalyzer::addhistcomparegenrecounseeded(TString selection) {

  all1dhists.push_back(new TH1F(selection+"recoegus_genmult","gen e/#gamma multiplicity",10,-1,9));
  all1dhists.push_back(new TH1F(selection+"recoegus_recomult","reco. e/#gamma multiplicity",20,-1,19));
  all1dhists.push_back(new TH1F(selection+"recoebus_gendEta","#Delta#eta(reco. e/#gamma, gen e.)",20000,-1,1));
  all1dhists.push_back(new TH1F(selection+"recoebus_gendPhi","#Delta#phi(reco. e/#gamma, gen e.)",20000,-1,1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_gendEta","#Delta#eta(reco. e/#gamma, gen e.)",20000,-1,1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_gendPhi","#Delta#phi(reco. e/#gamma, gen e.)",20000,-1,1));

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
    // If the trigger menu flow is inappropriately modelled in code, then print to undestand the cause
    if(egFiltN!=sel2egidx.size()) {
      cout<<event<<" : ***********Error! mis-match filter and reco selection eg20**********"<<endl;
      cout<<"reco unselected: "<<endl;
      for(unsigned int ctr=0; ctr<egRecoN; ctr++) {
	cout<<ctr<<"\t"<<egRecoPt[ctr]<<"\t"<<egFiltEta[ctr]<<"\t"<<egFiltPhi[ctr]<<"\t"<<eghltEgammaClusterShape_sigmaIEtaIEta5x5[ctr]<<"\t"<<eghltEgammaHoverE[ctr]/eghltEgammaSuperClusterEnergy[ctr]<<" || ";
      }
      cout<<endl;
      cout<<"Filter: "<<endl;
      for(unsigned int ctr=0; ctr<egFiltN; ctr++) {
	cout<<ctr<<"\t"<<egFiltPt[ctr]<<"\t"<<egFiltEta[ctr]<<"\t"<<egFiltPhi[ctr]<<" || ";
      }
      cout<<endl;
      cout<<"L1 Eg obj: "<<endl;
      for(unsigned int l1egcnt=0; l1egcnt<l1egObjN; l1egcnt++) {
	cout<<l1egcnt<<"\t"<<l1egObjPt[l1egcnt]<<"\t"<<l1egObjEta[l1egcnt]<<"\t"<<l1egObjPhi[l1egcnt]<<" || ";
      }
      cout<<endl;
      cout<<"L1 Filter: "<<endl;
      for(unsigned int l1egcnt=0; l1egcnt<l1FiltN; l1egcnt++) {
	cout<<l1egcnt<<"\t"<<l1FiltPt[l1egcnt]<<"\t"<<l1FiltEta[l1egcnt]<<"\t"<<l1FiltPhi[l1egcnt]<<" || ";
      }
      cout<<endl;
      cout<<"reco: "<<endl;
      for(unsigned int ctr=0; ctr<sel2egidx.size(); ctr++) {
	cout<<ctr<<"\t"<<egRecoPt[sel2egidx[ctr]]<<"\t"<<egFiltEta[sel2egidx[ctr]]<<"\t"<<egFiltPhi[sel2egidx[ctr]]<<" || ";
      }
      cout<<endl;
      cout<<"-------------------------------------------------"<<endl;
    }
    if(egFiltN_38!=sel3egidx.size()) {
      cout<<event<<"-----------Error! mis-match filter and reco selection eg38------------"<<endl;
      cout<<"reco unselected: "<<endl;
      for(unsigned int ctr=0; ctr<egRecoN; ctr++) {
	cout<<ctr<<"\t"<<egRecoPt[ctr]<<"\t"<<egFiltEta[ctr]<<"\t"<<egFiltPhi[ctr]<<"\t"<<eghltEgammaClusterShape_sigmaIEtaIEta5x5[ctr]<<"\t"<<eghltEgammaHoverE[ctr]/eghltEgammaSuperClusterEnergy[ctr]<<" || ";
      }
      cout<<endl;
      cout<<"Filter 38: "<<endl;
      for(unsigned int ctr=0; ctr<egFiltN_38; ctr++) {
	cout<<ctr<<"\t"<<egFiltPt_38[ctr]<<"\t"<<egFiltEta_38[ctr]<<"\t"<<egFiltPhi_38[ctr]<<" || ";
      }
      cout<<endl;
      cout<<"L1 Eg obj: "<<endl;
      for(unsigned int l1egcnt=0; l1egcnt<l1egObjN; l1egcnt++) {
	cout<<l1egcnt<<"\t"<<l1egObjPt[l1egcnt]<<"\t"<<l1egObjEta[l1egcnt]<<"\t"<<l1egObjPhi[l1egcnt]<<" || ";
      }
      cout<<endl;
      cout<<"L1 Filter: "<<endl;
      for(unsigned int l1egcnt=0; l1egcnt<l1FiltN; l1egcnt++) {
	cout<<l1egcnt<<"\t"<<l1FiltPt[l1egcnt]<<"\t"<<l1FiltEta[l1egcnt]<<"\t"<<l1FiltPhi[l1egcnt]<<" || ";
      }
      cout<<endl;
      cout<<"reco: "<<endl;
      for(unsigned int ctr=0; ctr<sel3egidx.size(); ctr++) {
	cout<<ctr<<"\t"<<egRecoPt[sel3egidx[ctr]]<<"\t"<<egFiltEta[sel3egidx[ctr]]<<"\t"<<egFiltPhi[sel3egidx[ctr]]<<" || ";
      }
      cout<<endl;
      cout<<"-------------------------------------------------"<<endl;
    }
*/
