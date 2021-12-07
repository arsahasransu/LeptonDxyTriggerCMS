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
  inputChain->SetBranchAddress("egushltEgammaSuperClusterEnergyarr", &egushltEgammaSuperClusterEnergy);
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
  int nosel=0, basicsel=0, sel2=0, sel3=0, sel10=0, sel10us=0, sel11us=0, sel12us=0, sel13=0, sel13us=0, sel14us=0, sel20=0, sel70=0, sel80=0, sel80us=0;

  // Define the histograms
  addhist("nosel");
  addhist("basicsel");
  addhist("sel2"); // Selection to cross-check the trigger filter selection
  addhist("sel3"); // Selection to cross-check the parent triggers
  addhist("sel10"); // e/gamma objects, pT and eta requirements only
  addhist("sel13"); // sel12+mudxy>0.01
  addhist("sel80"); // new mu cuts + loose egamma cuts
  if(isMC) addhistcomparegenrecounseeded("sel10us");
  addhistunseeded("sel10us"); // Unseeded e/gamma objects, pT and eta requirements only
  addhistunseeded("sel11us"); // sel10+in5x5noiseclnd
  addhistunseeded("sel12us"); // sel11+hoe
  addhistunseeded("sel13us"); // sel12+mudxy>0.01
  addhistunseeded("sel14us"); // sel13+ecal iso.
  addhistunseeded("sel20"); // Tighter pt and eta requirements on mu and e/gamma
  addhistunseeded("sel70"); // sel20 + loose egamma cuts
  addhistunseeded("sel80us"); // new mu cuts + loose egamma cuts
  
  // vector of mu indices
  vector<int> noselmuidx;
  vector<int> basicselmuidx;
  vector<int> sel2muidx;
  vector<int> sel3muidx;
  vector<int> sel10muidx;
  vector<int> sel11muidx;
  vector<int> sel12muidx;
  vector<int> sel13muidx;
  vector<int> sel14muidx;
  vector<int> sel20muidx;
  vector<int> sel70muidx;
  vector<int> sel80muidx;
  // vector of eg indices
  vector<int> noselegidx;
  vector<int> basicselegidx;
  vector<int> sel2egidx;
  vector<int> sel3egidx;
  vector<int> sel10egidx;
  vector<int> sel10egusidx;
  vector<int> sel11egusidx;
  vector<int> sel12egusidx;
  vector<int> sel13egidx;
  vector<int> sel13egusidx;
  vector<int> sel14egusidx;
  vector<int> sel20egidx; 
  vector<int> sel70egidx;
  vector<int> sel80egidx;
  vector<int> sel80egusidx;
 
  // Loop beginning on events
  for(unsigned int event=0; event<totEntries; event++) {

    inputChain->GetEntry(event);
    //if(event>1000) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    if(muRecoN>0) { // Do not go in here if there is not atleast one reco muon

      // Sort the muon objects based on their pT
      vector<int> sortedmuidx(muRecoN);
      iota(begin(sortedmuidx), end(sortedmuidx), 0);
      sort(&sortedmuidx[0], muRecoPt, muRecoN); // Verified that the algorithm works fine

      bool basicselmu = false;
      bool sel3mu = false;
      bool sel13mu = false;
      bool sel20mu = false;
      bool sel80mu = false;

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
	
	sel13mu = true;
	sel13mu *= (muRecoPt[idx]>=16);
	sel13mu *= (TMath::Abs(muRecoEta[idx])<2.5);
	sel13mu *= (TMath::Abs(muRecoDxy[idx])>0.01);
	if(sel13mu) sel13muidx.push_back(idx);
	
	sel20mu = true;
	sel20mu *= (muRecoPt[idx]>=20);
	sel20mu *= (TMath::Abs(muRecoEta[idx])<2.5);
	sel20mu *= (TMath::Abs(muRecoDxySig[idx])>=1);
	if(sel20mu) sel20muidx.push_back(idx);
	
	sel80mu = true;
	sel80mu *= (muRecoPt[idx]>=15);
	sel80mu *= (TMath::Abs(muRecoEta[idx])<2.5);
	sel80mu *= (TMath::Abs(muRecoDxySig[idx])>=4);
	if(sel80mu) sel80muidx.push_back(idx);
	
      } // End of loop on muon reco objects

    } // End of loop requiring one muon reco object

    sel2muidx = basicselmuidx;
    sel10muidx = basicselmuidx;
    sel11muidx = basicselmuidx;
    sel12muidx = basicselmuidx;
    sel14muidx = sel13muidx;
    sel70muidx = sel20muidx;

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
      bool sel10egus = false;
      bool sel11egus = false;
      bool sel12egus = false;
      bool sel13eg = false;
      bool sel13egus = false;
      bool sel14egus = false;
      bool sel20eg = false;
      bool sel70eg = false;
      bool sel80eg = false;
      bool sel80egus = false;
    
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
	sel2eg *= isL1EgSeeded(idx);
	sel2eg *= abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5[idx]<0.014:eghltEgammaClusterShape_sigmaIEtaIEta5x5[idx]<0.035;
	sel2eg *= abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.15*eghltEgammaSuperClusterEnergy[idx]:eghltEgammaHoverE[idx]<0.1*eghltEgammaSuperClusterEnergy[idx];
	if(sel2eg) sel2egidx.push_back(idx);

	sel3eg = true;
	sel3eg *= (sel3muidx.size()>0);
	sel3eg *= (egRecoPt[idx]>=38);
	sel3eg *= (TMath::Abs(egRecoEta[idx])<2.65);
	sel3eg *= isL1EgSeeded(idx);
	sel3eg *= abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.014:eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.035;
	sel3eg *= abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.15*eghltEgammaSuperClusterEnergy[idx]:eghltEgammaHoverE[idx]<0.1*eghltEgammaSuperClusterEnergy[idx];
	if(sel3eg) sel3egidx.push_back(idx);

	sel10eg = true;
	sel10eg *= (sel10muidx.size()>0);
	sel10eg *= (egRecoPt[idx]>=15);
	sel10eg *= (TMath::Abs(egRecoEta[idx])<2.65);
	sel10eg *= isL1EgSeeded(idx);
	if(sel10eg) sel10egidx.push_back(idx);
	
	sel13eg = true;
	sel13eg *= (sel13muidx.size()>0);
	sel13eg *= (egRecoPt[idx]>=15);
	sel13eg *= (TMath::Abs(egRecoEta[idx])<2.5);
	sel13eg *= isL1EgSeeded(idx);
	sel13eg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.012:eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.03);
	sel13eg *= eghltEgammaHoverE[idx]/eghltEgammaSuperClusterEnergy[idx]<0.1;
	if(sel13eg) sel13egidx.push_back(idx);

	sel80eg = true;
	sel80eg *= (sel80muidx.size()>0);
	sel80eg *= (egRecoPt[idx]>=15);
	sel80eg *= (TMath::Abs(egRecoEta[idx])<2.5);
	sel80eg *= isL1EgSeeded(idx);
	sel80eg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.014:eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.035);
	sel80eg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.046+1.16/eghltEgammaSuperClusterEnergy[idx]:eghltEgammaHoverE[idx]<0.0324+2.52/eghltEgammaSuperClusterEnergy[idx]);
	if(sel80eg) sel80egidx.push_back(idx);
	
      } // End of loop on egamma reco objects
            
      // Loop beginning on unseeded egamma reco objects
      for(unsigned int egidx=0; egidx<egusRecoN; egidx++) {
	
	unsigned int idx = sortedegusidx[egidx];
	
	sel10egus = true;
	sel10egus *= (sel10muidx.size()>0);
	sel10egus *= (egusRecoPt[idx]>=15);
	sel10egus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	if(sel10egus) sel10egusidx.push_back(idx);
	
	sel11egus = true;
	sel11egus *= (sel11muidx.size()>0);
	sel11egus *= (egusRecoPt[idx]>=15);
	sel11egus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	sel11egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.012:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.03);
	if(sel11egus) sel11egusidx.push_back(idx);
	
	sel12egus = true;
	sel12egus *= (sel12muidx.size()>0);
	sel12egus *= (egusRecoPt[idx]>=15);
	sel12egus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	sel12egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.012:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.03);
	sel12egus *= egushltEgammaHoverE[idx]/egushltEgammaSuperClusterEnergy[idx]<0.1;
	if(sel12egus) sel12egusidx.push_back(idx);

	sel13egus = true;
	sel13egus *= (sel13muidx.size()>0);
	sel13egus *= (egusRecoPt[idx]>=15);
	sel13egus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	sel13egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.012:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.03);
	sel13egus *= egushltEgammaHoverE[idx]/egushltEgammaSuperClusterEnergy[idx]<0.1;
	if(sel13egus) sel13egusidx.push_back(idx);

	sel14egus = true;
	sel14egus *= (sel14muidx.size()>0);
	sel14egus *= (egusRecoPt[idx]>=15);
	sel14egus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	sel14egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.012:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.03);
	sel14egus *= egushltEgammaHoverE[idx]/egushltEgammaSuperClusterEnergy[idx]<0.1;
	sel14egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaEcalPFClusterIso[idx]/egushltEgammaSuperClusterEnergy[idx]<0.15:egushltEgammaEcalPFClusterIso[idx]/egushltEgammaSuperClusterEnergy[idx]<0.1);
	if(sel14egus) sel14egusidx.push_back(idx);

	sel20eg = true;
	sel20eg *= (sel20muidx.size()>0);
	sel20eg *= (egusRecoPt[idx]>=20);
	sel20eg *= (TMath::Abs(egusRecoEta[idx])<2.5);
	if(sel20eg) sel20egidx.push_back(idx);
	
	sel70eg = true;
	sel70eg *= (sel70muidx.size()>0);
	sel70eg *= (egusRecoPt[idx]>=20);
	sel70eg *= (TMath::Abs(egusRecoEta[idx])<2.5);
	sel70eg *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0106:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0387);
	sel70eg *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.046+1.16/egushltEgammaSuperClusterEnergy[idx]:egushltEgammaHoverE[idx]<0.0324+2.52/egushltEgammaSuperClusterEnergy[idx]);
	if(sel70eg) sel70egidx.push_back(idx);
	
	sel80egus = true;
	sel80egus *= (sel80muidx.size()>0);
	sel80egus *= (egusRecoPt[idx]>=15);
	sel80egus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	//sel80egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0106:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0387);
	//sel80egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.046+1.16/egushltEgammaSuperClusterEnergy[idx]:egushltEgammaHoverE[idx]<0.0324+2.52/egushltEgammaSuperClusterEnergy[idx]);
	sel80egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.014:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.035);
	sel80egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.046+1.16*egushltEgammaSuperClusterEnergy[idx]:egushltEgammaHoverE[idx]<0.0324+2.52*egushltEgammaSuperClusterEnergy[idx]);
	if(sel80egus) sel80egusidx.push_back(idx);
	
      } // End of loop on unseeded egamma objects
      
      fillhistinevent("nosel", noselegidx);
      fillhistinevent("basicsel", basicselegidx);
      fillhistinevent("sel2", sel2egidx);
      fillhistinevent("sel3", sel3egidx);
      fillhistinevent("sel10", sel10egidx);
      fillhistinevent("sel13", sel13egidx);
      fillhistinevent("sel80", sel80egidx);
      if(isMC && sel10egusidx.size()>0) fillhistcomparegenrecounseeded("sel10us", sel10egusidx);
      fillhistineventunseeded("sel10us", sel10egusidx);
      fillhistineventunseeded("sel11us", sel11egusidx);
      fillhistineventunseeded("sel12us", sel12egusidx);
      fillhistineventunseeded("sel13us", sel13egusidx);
      fillhistineventunseeded("sel14us", sel14egusidx);
      fillhistineventunseeded("sel20", sel20egidx);
      fillhistineventunseeded("sel70", sel70egidx);
      fillhistineventunseeded("sel80us", sel80egusidx);
      
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
    if(sel10muidx.size()>0 && sel10egusidx.size()>0) sel10us++;
    if(sel11muidx.size()>0 && sel11egusidx.size()>0) sel11us++;
    if(sel12muidx.size()>0 && sel12egusidx.size()>0) sel12us++;
    if(sel13muidx.size()>0 && sel13egidx.size()>0) sel13++;
    if(sel13muidx.size()>0 && sel13egusidx.size()>0) sel13us++;
    if(sel14muidx.size()>0 && sel14egusidx.size()>0) sel14us++;
    if(sel20muidx.size()>0 && sel20egidx.size()>0) sel20++;
    if(sel70muidx.size()>0 && sel70egidx.size()>0) sel70++;
    if(sel80muidx.size()>0 && sel80egidx.size()>0) sel80++;
    if(sel80muidx.size()>0 && sel80egusidx.size()>0) sel80us++;
    
    // Clear all the vectors
    noselmuidx.clear();
    basicselmuidx.clear();
    sel2muidx.clear();
    sel3muidx.clear();
    sel10muidx.clear();
    sel11muidx.clear();
    sel12muidx.clear();
    sel13muidx.clear();
    sel14muidx.clear();
    sel20muidx.clear();
    sel70muidx.clear();
    sel80muidx.clear();
    noselegidx.clear();
    basicselegidx.clear();
    sel2egidx.clear();
    sel3egidx.clear();
    sel10egidx.clear();
    sel10egusidx.clear();
    sel11egusidx.clear();
    sel12egusidx.clear();
    sel13egidx.clear();
    sel13egusidx.clear();
    sel14egusidx.clear();
    sel20egidx.clear();
    sel70egidx.clear();
    sel80egidx.clear();
    sel80egusidx.clear();

  } // End of loop on events

  cout<<"With parent trigger: "<<sel3<<endl;
  cout<<totEntries<<"\t"<<nosel<<"\t"<<basicsel<<"\t"<<sel2<<"\t"<<sel3<<"\t"<<sel10<<"\t"<<sel10us<<"\t"<<sel11us<<"\t"<<sel12us<<"\t"<<sel13<<"\t"<<sel13us<<"\t"<<sel14us<<"\t"<<sel20<<"\t"<<sel70<<"\t"<<sel80<<"\t"<<sel80us<<endl;
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
  TH1F* recoeb_hovereoverpt = (TH1F*) outfile->Get(selection+"recoeb_hovereoverpt");
  TH1F* recoeb_hovereoversupcluse = (TH1F*) outfile->Get(selection+"recoeb_hovereoversupcluse");
  TH1F* recoeb_ecalpfclustiso = (TH1F*) outfile->Get(selection+"recoeb_ecalpfclustiso");
  TH1F* recoeb_ecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeb_ecalpfclustisoovere");
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
  TH1F* recoee_hovereoverpt = (TH1F*) outfile->Get(selection+"recoee_hovereoverpt");
  TH1F* recoee_hovereoversupcluse = (TH1F*) outfile->Get(selection+"recoee_hovereoversupcluse");
  TH1F* recoee_ecalpfclustiso = (TH1F*) outfile->Get(selection+"recoee_ecalpfclustiso");
  TH1F* recoee_ecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoee_ecalpfclustisoovere");
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
      recoeb_hovereoverpt->Fill(eghltEgammaHoverE[egidx[0]]/egRecoPt[egidx[0]]);
      recoeb_hovereoversupcluse->Fill(eghltEgammaHoverE[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_ecalpfclustiso->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]);
      recoeb_ecalpfclustisoovere->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
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
      recoee_hovereoverpt->Fill(eghltEgammaHoverE[egidx[0]]/egRecoPt[egidx[0]]);
      recoee_hovereoversupcluse->Fill(eghltEgammaHoverE[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_ecalpfclustiso->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]);
      recoee_ecalpfclustisoovere->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
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

  TH1F* genmchegpt = (TH1F*) outfile->Get(selection+"recoegusgenmch_egpt");
  TH1F* genmchegeta = (TH1F*) outfile->Get(selection+"recoegusgenmch_egeta");
  TH1F* genmchegphi = (TH1F*) outfile->Get(selection+"recoegusgenmch_egphi");

  // Get barrel variables
  TH1F* recoeb_egclustershape = (TH1F*) outfile->Get(selection+"recoebus_egclustershape");
  TH1F* recoeb_in5x5clusshape = (TH1F*) outfile->Get(selection+"recoebus_in5x5clusshape");
  TH1F* recoeb_in5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoebus_in5x5noiseclnd");
  TH1F* recoeb_scenergy = (TH1F*) outfile->Get(selection+"recoebus_scenergy");
  TH1F* recoeb_hovere = (TH1F*) outfile->Get(selection+"recoebus_hovere");
  TH1F* recoeb_hovereoverpt = (TH1F*) outfile->Get(selection+"recoebus_hovereoverpt");
  TH1F* recoeb_hovereoversupcluse = (TH1F*) outfile->Get(selection+"recoebus_hovereoversupcluse");
  TH1F* recoeb_ecalpfclustiso = (TH1F*) outfile->Get(selection+"recoebus_ecalpfclustiso");
  TH1F* recoeb_ecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoebus_ecalpfclustisoovere");
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

  // Gen-matched barrel variables
  TH1F* recoebgenmch_egclustershape = (TH1F*) outfile->Get(selection+"recoebusgenmch_egclustershape");
  TH1F* recoebgenmch_in5x5clusshape = (TH1F*) outfile->Get(selection+"recoebusgenmch_in5x5clusshape");
  TH1F* recoebgenmch_in5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoebusgenmch_in5x5noiseclnd");
  TH1F* recoebgenmch_scenergy = (TH1F*) outfile->Get(selection+"recoebusgenmch_scenergy");
  TH1F* recoebgenmch_hovere = (TH1F*) outfile->Get(selection+"recoebusgenmch_hovere");
  TH1F* recoebgenmch_hovereoverpt = (TH1F*) outfile->Get(selection+"recoebusgenmch_hovereoverpt");
  TH1F* recoebgenmch_hovereoversupcluse = (TH1F*) outfile->Get(selection+"recoebusgenmch_hovereoversupcluse");
  TH1F* recoebgenmch_ecalpfclustiso = (TH1F*) outfile->Get(selection+"recoebusgenmch_ecalpfclustiso");
  TH1F* recoebgenmch_ecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoebusgenmch_ecalpfclustisoovere");
  TH1F* recoebgenmch_hcalpfclustiso = (TH1F*) outfile->Get(selection+"recoebusgenmch_hcalpfclustiso");
  TH1F* recoebgenmch_pixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoebusgenmch_pixelmchvar_s2");
  TH1F* recoebgenmch_trkiso = (TH1F*) outfile->Get(selection+"recoebusgenmch_trkiso");
  TH1F* recoebgenmch_chi2 = (TH1F*) outfile->Get(selection+"recoebusgenmch_chi2");
  TH1F* recoebgenmch_deta = (TH1F*) outfile->Get(selection+"recoebusgenmch_deta");
  TH1F* recoebgenmch_detaseed = (TH1F*) outfile->Get(selection+"recoebusgenmch_detaseed");
  TH1F* recoebgenmch_dphi = (TH1F*) outfile->Get(selection+"recoebusgenmch_dphi");
  TH1F* recoebgenmch_mhits = (TH1F*) outfile->Get(selection+"recoebusgenmch_mhits");
  TH1F* recoebgenmch_nlayerit = (TH1F*) outfile->Get(selection+"recoebusgenmch_nlayerit");
  TH1F* recoebgenmch_ooeseedoop = (TH1F*) outfile->Get(selection+"recoebusgenmch_ooeseedoop");
  TH1F* recoebgenmch_ooesclsoop = (TH1F*) outfile->Get(selection+"recoebusgenmch_ooesclsoop");
  TH1F* recoebgenmch_valhits = (TH1F*) outfile->Get(selection+"recoebusgenmch_valhits");  
  TH1F* recoebgenmch_seedclustime = (TH1F*) outfile->Get(selection+"recoebusgenmch_seedclustime");  

  // Get end-cap variables
  TH1F* recoee_egclustershape = (TH1F*) outfile->Get(selection+"recoeeus_egclustershape");
  TH1F* recoee_in5x5clusshape = (TH1F*) outfile->Get(selection+"recoeeus_in5x5clusshape");
  TH1F* recoee_in5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoeeus_in5x5noiseclnd");
  TH1F* recoee_scenergy = (TH1F*) outfile->Get(selection+"recoeeus_scenergy");
  TH1F* recoee_hovere = (TH1F*) outfile->Get(selection+"recoeeus_hovere");
  TH1F* recoee_hovereoverpt = (TH1F*) outfile->Get(selection+"recoeeus_hovereoverpt");
  TH1F* recoee_hovereoversupcluse = (TH1F*) outfile->Get(selection+"recoeeus_hovereoversupcluse");
  TH1F* recoee_ecalpfclustiso = (TH1F*) outfile->Get(selection+"recoeeus_ecalpfclustiso");
  TH1F* recoee_ecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeeus_ecalpfclustisoovere");
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

  // Gen-matched endcap variables
  TH1F* recoeegenmch_egclustershape = (TH1F*) outfile->Get(selection+"recoeeusgenmch_egclustershape");
  TH1F* recoeegenmch_in5x5clusshape = (TH1F*) outfile->Get(selection+"recoeeusgenmch_in5x5clusshape");
  TH1F* recoeegenmch_in5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoeeusgenmch_in5x5noiseclnd");
  TH1F* recoeegenmch_scenergy = (TH1F*) outfile->Get(selection+"recoeeusgenmch_scenergy");
  TH1F* recoeegenmch_hovere = (TH1F*) outfile->Get(selection+"recoeeusgenmch_hovere");
  TH1F* recoeegenmch_hovereoverpt = (TH1F*) outfile->Get(selection+"recoeeusgenmch_hovereoverpt");
  TH1F* recoeegenmch_hovereoversupcluse = (TH1F*) outfile->Get(selection+"recoeeusgenmch_hovereoversupcluse");
  TH1F* recoeegenmch_ecalpfclustiso = (TH1F*) outfile->Get(selection+"recoeeusgenmch_ecalpfclustiso");
  TH1F* recoeegenmch_ecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeeusgenmch_ecalpfclustisoovere");
  TH1F* recoeegenmch_hcalpfclustiso = (TH1F*) outfile->Get(selection+"recoeeusgenmch_hcalpfclustiso");
  TH1F* recoeegenmch_pixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoeeusgenmch_pixelmchvar_s2");
  TH1F* recoeegenmch_trkiso = (TH1F*) outfile->Get(selection+"recoeeusgenmch_trkiso");
  TH1F* recoeegenmch_chi2 = (TH1F*) outfile->Get(selection+"recoeeusgenmch_chi2");
  TH1F* recoeegenmch_deta = (TH1F*) outfile->Get(selection+"recoeeusgenmch_deta");
  TH1F* recoeegenmch_detaseed = (TH1F*) outfile->Get(selection+"recoeeusgenmch_detaseed");
  TH1F* recoeegenmch_dphi = (TH1F*) outfile->Get(selection+"recoeeusgenmch_dphi");
  TH1F* recoeegenmch_mhits = (TH1F*) outfile->Get(selection+"recoeeusgenmch_mhits");
  TH1F* recoeegenmch_nlayerit = (TH1F*) outfile->Get(selection+"recoeeusgenmch_nlayerit");
  TH1F* recoeegenmch_ooeseedoop = (TH1F*) outfile->Get(selection+"recoeeusgenmch_ooeseedoop");
  TH1F* recoeegenmch_ooesclsoop = (TH1F*) outfile->Get(selection+"recoeeusgenmch_ooesclsoop");
  TH1F* recoeegenmch_valhits = (TH1F*) outfile->Get(selection+"recoeeusgenmch_valhits");
  TH1F* recoeegenmch_seedclustime = (TH1F*) outfile->Get(selection+"recoeeusgenmch_seedclustime");  

  if(egidx.size()>0) {
    egmult->Fill(egidx.size());
    egpt->Fill(egusRecoPt[egidx[0]]);
    egeta->Fill(egusRecoEta[egidx[0]]);
    egphi->Fill(egusRecoPhi[egidx[0]]);

    // Fill barrel variables
    if(TMath::Abs(egusRecoEta[egidx[0]])<=1.479) {
      recoeb_egclustershape->Fill(egushltEgammaClusterShape[egidx[0]]);
      recoeb_in5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[0]]);
      recoeb_in5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[0]]);
      //recoeb_scenergy->Fill(egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_hovere->Fill(egushltEgammaHoverE[egidx[0]]);
      recoeb_hovereoverpt->Fill(egushltEgammaHoverE[egidx[0]]/egusRecoPt[egidx[0]]);
      recoeb_hovereoversupcluse->Fill(egushltEgammaHoverE[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_ecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]);
      recoeb_ecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_hcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[egidx[0]]);
      if(egushltEcalSeedClusterTime[egidx[0]]!=0) recoeb_seedclustime->Fill(egushltEcalSeedClusterTime[egidx[0]]);
      if(egushltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoeb_pixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoeb_pixelmchvar_s2->Fill(-5);
      }
      if(egushltEgammaGsfTrackVars_Deta[egidx[0]]<999999) {
	recoeb_trkiso->Fill(egushltEgammaEleGsfTrackIso[egidx[0]]);
	recoeb_chi2->Fill(egushltEgammaGsfTrackVars_Chi2[egidx[0]]);
	recoeb_deta->Fill(egushltEgammaGsfTrackVars_Deta[egidx[0]]);
	recoeb_detaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[egidx[0]]);
	recoeb_dphi->Fill(egushltEgammaGsfTrackVars_Dphi[egidx[0]]);
	recoeb_mhits->Fill(egushltEgammaGsfTrackVars_MissingHits[egidx[0]]);
	recoeb_nlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[egidx[0]]);
	recoeb_ooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[0]]);
	recoeb_ooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[0]]);
	//recoeb_valhits->Fill(egushltEgammaGsfTrackVars_ValidHits[egidx[0]]);
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
      recoee_egclustershape->Fill(egushltEgammaClusterShape[egidx[0]]);
      recoee_in5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[0]]);
      recoee_in5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[0]]);
      //recoee_scenergy->Fill(egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_hovere->Fill(egushltEgammaHoverE[egidx[0]]);
      recoee_hovereoverpt->Fill(egushltEgammaHoverE[egidx[0]]/egusRecoPt[egidx[0]]);
      recoee_hovereoversupcluse->Fill(egushltEgammaHoverE[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_ecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]);
      recoee_ecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_hcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[egidx[0]]);	
      if(egushltEcalSeedClusterTime[egidx[0]]!=0) recoee_seedclustime->Fill(egushltEcalSeedClusterTime[egidx[0]]);
      if(egushltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoee_pixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoee_pixelmchvar_s2->Fill(-5);
      }
      if(egushltEgammaGsfTrackVars_Deta[egidx[0]]<999999) {
	recoee_trkiso->Fill(egushltEgammaEleGsfTrackIso[egidx[0]]);
	recoee_chi2->Fill(egushltEgammaGsfTrackVars_Chi2[egidx[0]]);
	recoee_deta->Fill(egushltEgammaGsfTrackVars_Deta[egidx[0]]);
	recoee_detaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[egidx[0]]);
	recoee_dphi->Fill(egushltEgammaGsfTrackVars_Dphi[egidx[0]]);
	recoee_mhits->Fill(egushltEgammaGsfTrackVars_MissingHits[egidx[0]]);
	recoee_nlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[egidx[0]]);
	recoee_ooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[0]]);
	recoee_ooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[0]]);
	//recoee_valhits->Fill(egushltEgammaGsfTrackVars_ValidHits[egidx[0]]);
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
	// Fill gen-matched barrel variables
	if(TMath::Abs(egusRecoEta[genmchegpos])<=1.479) {
	  recoebgenmch_egclustershape->Fill(egushltEgammaClusterShape[genmchegpos]);
	  recoebgenmch_in5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[genmchegpos]);
	  recoebgenmch_in5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[genmchegpos]);
	  //recoebgenmch_scenergy->Fill(egushltEgammaSuperClusterEnergy[genmchegpos]);
	  recoebgenmch_hovere->Fill(egushltEgammaHoverE[genmchegpos]);
	  recoebgenmch_hovereoverpt->Fill(egushltEgammaHoverE[genmchegpos]/egusRecoPt[genmchegpos]);
	  recoebgenmch_hovereoversupcluse->Fill(egushltEgammaHoverE[genmchegpos]/egushltEgammaSuperClusterEnergy[genmchegpos]);
	  recoebgenmch_ecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[genmchegpos]);
	  recoebgenmch_ecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[genmchegpos]/egushltEgammaSuperClusterEnergy[genmchegpos]);
	  recoebgenmch_hcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[genmchegpos]);
	  if(egushltEcalSeedClusterTime[genmchegpos]!=0) recoebgenmch_seedclustime->Fill(egushltEcalSeedClusterTime[genmchegpos]);
	  if(egushltEgammaPixelMatchVars_s2[genmchegpos]<TMath::Power(10,36)) {
	    recoebgenmch_pixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[genmchegpos]);
	  }
	  else {
	    recoebgenmch_pixelmchvar_s2->Fill(-5);
	  }
	  if(egushltEgammaGsfTrackVars_Deta[genmchegpos]<999999) {
	    recoebgenmch_trkiso->Fill(egushltEgammaEleGsfTrackIso[genmchegpos]);
	    recoebgenmch_chi2->Fill(egushltEgammaGsfTrackVars_Chi2[genmchegpos]);
	    recoebgenmch_deta->Fill(egushltEgammaGsfTrackVars_Deta[genmchegpos]);
	    recoebgenmch_detaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[genmchegpos]);
	    recoebgenmch_dphi->Fill(egushltEgammaGsfTrackVars_Dphi[genmchegpos]);
	    recoebgenmch_mhits->Fill(egushltEgammaGsfTrackVars_MissingHits[genmchegpos]);
	    recoebgenmch_nlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[genmchegpos]);
	    recoebgenmch_ooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[genmchegpos]);
	    recoebgenmch_ooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[genmchegpos]);
	    //recoebgenmch_valhits->Fill(egushltEgammaGsfTrackVars_ValidHits[genmchegpos]);
	  }
	  else {
	    recoebgenmch_trkiso->Fill(-100);
	    recoebgenmch_chi2->Fill(-20);
	    recoebgenmch_deta->Fill(-5);
	    recoebgenmch_detaseed->Fill(-5);
	    recoebgenmch_dphi->Fill(-5);
	    recoebgenmch_mhits->Fill(-15);
	    recoebgenmch_nlayerit->Fill(-15);
	    recoebgenmch_ooeseedoop->Fill(-50);
	    recoebgenmch_ooesclsoop->Fill(-50);
	    //recoebgenmch_valhits->Fill(-15);
	  }
	} // End of filling gen-matched barrel variables
	else { // Fill gen-matched end-cap variables
	  recoeegenmch_egclustershape->Fill(egushltEgammaClusterShape[genmchegpos]);
	  recoeegenmch_in5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[genmchegpos]);
	  recoeegenmch_in5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[genmchegpos]);
	  //recoeegenmch_scenergy->Fill(egushltEgammaSuperClusterEnergy[genmchegpos]);
	  recoeegenmch_hovere->Fill(egushltEgammaHoverE[genmchegpos]);
	  recoeegenmch_hovereoverpt->Fill(egushltEgammaHoverE[genmchegpos]/egusRecoPt[genmchegpos]);
	  recoeegenmch_hovereoversupcluse->Fill(egushltEgammaHoverE[genmchegpos]/egushltEgammaSuperClusterEnergy[genmchegpos]);
	  recoeegenmch_ecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[genmchegpos]);
	  recoeegenmch_ecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[genmchegpos]/egushltEgammaSuperClusterEnergy[genmchegpos]);
	  recoeegenmch_hcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[genmchegpos]);	
	  if(egushltEcalSeedClusterTime[genmchegpos]!=0) recoeegenmch_seedclustime->Fill(egushltEcalSeedClusterTime[genmchegpos]);
	  if(egushltEgammaPixelMatchVars_s2[genmchegpos]<TMath::Power(10,36)) {
	    recoeegenmch_pixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[genmchegpos]);
	  }
	  else {
	    recoeegenmch_pixelmchvar_s2->Fill(-5);
	  }
	  if(egushltEgammaGsfTrackVars_Deta[genmchegpos]<999999) {
	    recoeegenmch_trkiso->Fill(egushltEgammaEleGsfTrackIso[genmchegpos]);
	    recoeegenmch_chi2->Fill(egushltEgammaGsfTrackVars_Chi2[genmchegpos]);
	    recoeegenmch_deta->Fill(egushltEgammaGsfTrackVars_Deta[genmchegpos]);
	    recoeegenmch_detaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[genmchegpos]);
	    recoeegenmch_dphi->Fill(egushltEgammaGsfTrackVars_Dphi[genmchegpos]);
	    recoeegenmch_mhits->Fill(egushltEgammaGsfTrackVars_MissingHits[genmchegpos]);
	    recoeegenmch_nlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[genmchegpos]);
	    recoeegenmch_ooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[genmchegpos]);
	    recoeegenmch_ooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[genmchegpos]);
	    //recoeegenmch_valhits->Fill(egushltEgammaGsfTrackVars_ValidHits[genmchegpos]);
	  }
	  else {
	    recoeegenmch_trkiso->Fill(-100);
	    recoeegenmch_chi2->Fill(-20);
	    recoeegenmch_deta->Fill(-5);
	    recoeegenmch_detaseed->Fill(-5);
	    recoeegenmch_dphi->Fill(-5);
	    recoeegenmch_mhits->Fill(-15);
	    recoeegenmch_nlayerit->Fill(-15);
	    recoeegenmch_ooeseedoop->Fill(-50);
	    recoeegenmch_ooesclsoop->Fill(-50);
	    //recoeegenmch_valhits->Fill(-15);
	  }
	} // End of filling gen-matched end-cap variables
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
  all1dhists.push_back(new TH1F(selection+"recoeb_hovere","barrel e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeb_hovereoverpt","barrel e/#gamma H/p_{T}",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoeb_scenergy","barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeb_ecalpfclustiso","barrel e/#gamma ecal PF Iso. / GeV",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_ecalpfclustisoovere","barrel e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeb_hcalpfclustiso","barrel e/#gamma hcal PF Iso.",1000,-5,95));
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
  all1dhists.push_back(new TH1F(selection+"recoee_hovere","end-cap e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoee_hovereoverpt","end-cap e/#gamma H/p_{T}",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoee_scenergy","end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoee_ecalpfclustiso","end-cap e/#gamma ecal PF Iso. / GeV",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_ecalpfclustisoovere","end-cap e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoee_hcalpfclustiso","end-cap e/#gamma hcal PF Iso.",1000,-5,95));
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
    all1dhists.push_back(new TH1F(selection+"recoegusgenmch_egpt","gen matched reco e/#gamma p_{T} / GeV",550,-50,500));
    all1dhists.push_back(new TH1F(selection+"recoegusgenmch_egeta","gen matched reco e/#gamma #eta",52,-2.6,2.6));
    all1dhists.push_back(new TH1F(selection+"recoegusgenmch_egphi","gen matched reco e/#gamma #phi",66,-3.3,3.3));
  }

  // barrel variables
  all1dhists.push_back(new TH1F(selection+"recoebus_egclustershape","barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_in5x5clusshape","barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_in5x5noiseclnd","barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_hovere","barrel e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoebus_hovereoverpt","barrel e/#gamma H/p_{T}",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoebus_scenergy","barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoebus_ecalpfclustiso","barrel e/#gamma ecal PF Iso. / GeV",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebus_ecalpfclustisoovere","barrel e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoebus_hcalpfclustiso","barrel e/#gamma hcal PF Iso.",1000,-5,95));
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

  // gen-matched barrel variables
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_egclustershape","barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_in5x5clusshape","barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_in5x5noiseclnd","barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_hovere","barrel e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_hovereoverpt","barrel e/#gamma H/p_{T}",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_scenergy","barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_ecalpfclustiso","barrel e/#gamma ecal PF Iso. / GeV",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_ecalpfclustisoovere","barrel e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_hcalpfclustiso","barrel e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_pixelmchvar_s2","barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_chi2","barrel e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_deta","barrel e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_detaseed","barrel e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_dphi","barrel e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_mhits","barrel e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_nlayerit","barrel e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_ooeseedoop","barrel e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_ooesclsoop","barrel e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_valhits","barrel e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_trkiso","barrel e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_hovereoversupcluse","barrel e/#gamma (H/E)/supClusEnergy",10000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoebusgenmch_seedclustime","barrel e/#gamma_{seed} time / ns",20000,-10,10));

  // end-cap variables
  all1dhists.push_back(new TH1F(selection+"recoeeus_egclustershape","end-cap e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_in5x5clusshape","end-cap e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_in5x5noiseclnd","end-cap e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_hovere","end-cap e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeeus_hovereoverpt","barrel e/#gamma H/p_{T}",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoeeus_scenergy","end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeeus_ecalpfclustiso","end-cap e/#gamma ecal PF Iso. / GeV",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeus_ecalpfclustisoovere","end-cap e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeeus_hcalpfclustiso","end-cap e/#gamma hcal PF Iso.",1000,-5,95));
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

  // gen-matched end-cap variables
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_egclustershape","end-cap e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_in5x5clusshape","end-cap e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_in5x5noiseclnd","end-cap e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_hovere","end-cap e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_hovereoverpt","barrel e/#gamma H/p_{T}",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_scenergy","end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_ecalpfclustiso","end-cap e/#gamma ecal PF Iso. / GeV",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_ecalpfclustisoovere","end-cap e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_hcalpfclustiso","end-cap e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_pixelmchvar_s2","end-cap e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_chi2","end-cap e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_deta","end-cap e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_detaseed","end-cap e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_dphi","end-cap e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_mhits","end-cap e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_nlayerit","end-cap e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_ooeseedoop","end-cap e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_ooesclsoop","end-cap e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_valhits","end-cap e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_trkiso","end-cap e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_hovereoversupcluse","end-cap e/#gamma (H/E)/supClusEnergy",10000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeeusgenmch_seedclustime","end-cap e/#gamma_{seed} time / ns",20000,-10,10));
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

// Function to do l1 match filtering
bool data_robustanalyzer::isL1EgSeeded(int idx) {

  bool l1matchdecision = false;
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
      l1matchdecision = true;
    }
  }

  return l1matchdecision;
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
