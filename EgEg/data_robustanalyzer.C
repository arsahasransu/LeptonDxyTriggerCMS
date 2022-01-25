#include "data_robustanalyzer.hh"
#include <iostream>
#include <numeric>

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
data_robustanalyzer::data_robustanalyzer(TString filename, TString outfilename, bool issimu){

  isMC = issimu;
  
  inputChain = new TChain("events");
  cout<<"Initializing for file: "<<filename<<endl;
  inputChain->Add(filename);

  //inputChain->SetBranchAddress("l1EgObjn", &l1egObjN);
  //inputChain->SetBranchAddress("l1EgObj_pt", &l1egObjPt);
  //inputChain->SetBranchAddress("l1EgObj_eta", &l1egObjEta);
  //inputChain->SetBranchAddress("l1EgObj_phi", &l1egObjPhi);
  inputChain->SetBranchAddress("l1Filtn", &l1FiltN);
  inputChain->SetBranchAddress("l1Filt_pt", &l1FiltPt);
  inputChain->SetBranchAddress("l1Filt_eta", &l1FiltEta);
  inputChain->SetBranchAddress("l1Filt_phi", &l1FiltPhi);
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
  inputChain->SetBranchAddress("eghltEgammaSuperClusterEnergyarr", &eghltEgammaSuperClusterEnergy);
  inputChain->SetBranchAddress("egEcalSeedClusterTimearr", &eghltEcalSeedClusterTime);
  inputChain->SetBranchAddress("eghltEgammaGsfTrackVars_ValidHitsarr", &eghltEgammaGsfTrackVars_ValidHits);
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
  inputChain->SetBranchAddress("egushltEgammaSuperClusterEnergyarr", &egushltEgammaSuperClusterEnergy);
  inputChain->SetBranchAddress("egushltEgammaPixelMatchVars_s2arr", &egushltEgammaPixelMatchVars_s2);
  inputChain->SetBranchAddress("egushltEgammaEleGsfTrackIsoarr", &egushltEgammaEleGsfTrackIso);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_Chi2arr", &egushltEgammaGsfTrackVars_Chi2);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_Detaarr", &egushltEgammaGsfTrackVars_Deta);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_DetaSeedarr", &egushltEgammaGsfTrackVars_DetaSeed);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_Dphiarr", &egushltEgammaGsfTrackVars_Dphi);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_MissingHitsarr", &egushltEgammaGsfTrackVars_MissingHits);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_NLayerITarr", &egushltEgammaGsfTrackVars_NLayerIT);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_OneOESeedMinusOneOParr", &egushltEgammaGsfTrackVars_OneOESeedMinusOneOP);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_OneOESuperMinusOneOParr", &egushltEgammaGsfTrackVars_OneOESuperMinusOneOP);
  inputChain->SetBranchAddress("egushltEgammaGsfTrackVars_ValidHitsarr", &egushltEgammaGsfTrackVars_ValidHits);
  inputChain->SetBranchAddress("egusEcalSeedClusterTimearr", &egushltEcalSeedClusterTime);
  inputChain->SetBranchAddress("egushltEgammaSuperClusterEnergyarr", &egushltEgammaSuperClusterEnergy);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaPixelMatchVars_s2arr", &eguspxlmch22hltEgammaPixelMatchVars_s2);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaEleGsfTrackIsoarr", &eguspxlmch22hltEgammaEleGsfTrackIso);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaGsfTrackVars_Chi2arr", &eguspxlmch22hltEgammaGsfTrackVars_Chi2);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaGsfTrackVars_Detaarr", &eguspxlmch22hltEgammaGsfTrackVars_Deta);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaGsfTrackVars_DetaSeedarr", &eguspxlmch22hltEgammaGsfTrackVars_DetaSeed);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaGsfTrackVars_Dphiarr", &eguspxlmch22hltEgammaGsfTrackVars_Dphi);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaGsfTrackVars_MissingHitsarr", &eguspxlmch22hltEgammaGsfTrackVars_MissingHits);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaGsfTrackVars_NLayerITarr", &eguspxlmch22hltEgammaGsfTrackVars_NLayerIT);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaGsfTrackVars_OneOESeedMinusOneOParr", &eguspxlmch22hltEgammaGsfTrackVars_OneOESeedMinusOneOP);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaGsfTrackVars_OneOESuperMinusOneOParr", &eguspxlmch22hltEgammaGsfTrackVars_OneOESuperMinusOneOP);
  inputChain->SetBranchAddress("eguspxlmch22hltEgammaGsfTrackVars_ValidHitsarr", &eguspxlmch22hltEgammaGsfTrackVars_ValidHits);

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
    inputChain->SetBranchAddress("genLepTimeLightarr",&genLepTimeLight);
    inputChain->SetBranchAddress("genLepTimeActarr",&genLepTimeAct);
    inputChain->SetBranchAddress("genLepTimeDiffarr",&genLepTimeDiff);
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
void data_robustanalyzer::analyzersinglefile(int splitCnt) { // Assume splitCnt to range from 0 to nCores

  int totEntries = inputChain->GetEntries();
  cout<<"Total number of entries: "<<totEntries<<endl;
  int nCores = 6; // Assume parallel processing over 7 cores where
  // there is a lesser no.of events in the last core
  int beginevent = splitCnt*(totEntries/nCores);
  int endevent = (splitCnt+1)*(totEntries/nCores);
  endevent = endevent<totEntries?endevent:totEntries; // Verfied that this logic to parallelize works

  // Count events passing certain selections
  int nosel=0, noselus=0, basicsel=0, basicselus=0, selelevetoid=0, selelevetoidus=0, selelevetozwindidus=0, selelevetozoppoidus=0, seleletightid=0, seleletightidus=0;

  // Define the histograms
  if(isMC) addgenhist("gennosel");
  if(isMC) addgenhist("genbarsel");
  if(isMC) addgenhist("genbarselptgt10");
  addhist("nosel");
  addhistunseeded("noselus");
  addhist("basicsel");
  addhistunseeded("basicselus");
  addhist("selelevetoid");
  addhistunseeded("selelevetoidus");
  addhistunseeded("selelevetozwindidus");
  addhistunseeded("selelevetozoppoidus");
  addhist("seleletightid");
  addhistunseeded("seleletightidus");
  
  // vector of eg indices
  vector<int> gennoselegidx;
  vector<int> genbarselegidx;
  vector<int> genbarselptgt10egidx;
  vector<int> noselegidx;
  vector<int> noselegusidx;
  vector<int> basicselegidx;
  vector<int> basicselegusidx;
  vector<int> selelevetoidegidx;
  vector<int> selelevetoidegusidx;
  vector<int> selelevetozwindidegusidx;
  vector<int> selelevetozoppoidegusidx;
  vector<int> seleletightidegidx;
  vector<int> seleletightidegusidx;
  
  // Loop beginning on events
  for(unsigned int event=beginevent; event<endevent; event++) {

    inputChain->GetEntry(event);
    //if(event>1) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    // Block for gen leptons
    if(isMC) {
      bool gennoseleg = false;
      bool genbarseleg = false;
      bool genbarselptgt10eg = false;
      for(unsigned int genCtr=0; genCtr<genLepN; genCtr++) {

	gennoseleg = true;
	gennoseleg *= abs(genLepPid[genCtr])==11;
	gennoseleg *= abs(genLepMomPid[genCtr])==9000007 || abs(genLepMomPid[genCtr])==23;
	if(gennoseleg) gennoselegidx.push_back(genCtr);

      	genbarseleg = true;
	genbarseleg *= abs(genLepPid[genCtr])==11;
	genbarseleg *= abs(genLepMomPid[genCtr])==9000007 || abs(genLepMomPid[genCtr])==23;
	genbarseleg *= abs(genLepEta[genCtr])<1.479;
	if(genbarseleg) genbarselegidx.push_back(genCtr);

	genbarselptgt10eg = true;
	genbarselptgt10eg *= abs(genLepPid[genCtr])==11;
	genbarselptgt10eg *= abs(genLepMomPid[genCtr])==9000007 || abs(genLepMomPid[genCtr])==23;
	genbarselptgt10eg *= abs(genLepEta[genCtr])<1.479;
	genbarselptgt10eg *= genLepPt[genCtr]>10;
	if(genbarselptgt10eg) genbarselptgt10egidx.push_back(genCtr);
      }

      if(gennoselegidx.size()>=2) fillgenhistinevent("gennosel",gennoselegidx);
      if(genbarselegidx.size()>=2) fillgenhistinevent("genbarsel",genbarselegidx);
      if(genbarselptgt10egidx.size()>=2) fillgenhistinevent("genbarselptgt10",genbarselptgt10egidx);
    }
    
    if(egusRecoN<0 || egRecoN<0) cout<<"Error!! Negative number of objects pas possible."<<endl;
    //if(egusRecoN<egRecoN) cout<<"Error!! number of unseeded egamma lesser than seeded egamma"<<endl;
      
    if(egRecoN>=1 && egusRecoN>=2) { // Atleast two reco eg object in the event

      // Sort the egamma objects based on their pT
      vector<int> sortedegidx(egRecoN);
      iota(begin(sortedegidx), end(sortedegidx), 0);
      sort(&sortedegidx[0], egRecoPt, egRecoN); // Verified that the algorithm works fine
      
      vector<int> sortedegusidx(egusRecoN);
      iota(begin(sortedegusidx), end(sortedegusidx), 0);
      sort(&sortedegusidx[0], egusRecoPt, egusRecoN); // Verified that the algorithm works fine

      bool basicseleg = false;
      bool basicselegus = false;
      bool selelevetoideg = false;
      bool selelevetoidegus = false;
      bool seleletightideg = false;
      bool seleletightidegus = false;
    
      // Loop beginning on egamma reco objects
      for(unsigned int egidx=0; egidx<egRecoN; egidx++) {

	unsigned int idx = sortedegidx[egidx];
	noselegidx.push_back(idx);

	basicseleg = true;
	basicseleg *= (TMath::Abs(egRecoEta[idx])<2.65);
	basicseleg *= (egRecoPt[idx]>=15);
	if(basicseleg) basicselegidx.push_back(idx);

	selelevetoideg = true;
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<2.65);
	selelevetoideg *= (egRecoPt[idx]>=15);
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0126:eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0457);
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_DetaSeed[idx])<0.00463:abs(eghltEgammaGsfTrackVars_DetaSeed[idx])<0.00814);
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_Dphi[idx])<0.148:abs(eghltEgammaGsfTrackVars_Dphi[idx])<0.19);
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.05*eghltEgammaSuperClusterEnergy[idx]+1.16:eghltEgammaHoverE[idx]<0.05*eghltEgammaSuperClusterEnergy[idx]+2.54);
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.209:abs(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.132);
	if(selelevetoideg) selelevetoidegidx.push_back(idx);

	seleletightideg = true;
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<2.65);
	seleletightideg *= (egRecoPt[idx]>=15);
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0104:eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0353);
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_DetaSeed[idx])<0.00255:abs(eghltEgammaGsfTrackVars_DetaSeed[idx])<0.00501);
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_Dphi[idx])<0.022:abs(eghltEgammaGsfTrackVars_Dphi[idx])<0.0236);
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.026*eghltEgammaSuperClusterEnergy[idx]+1.15:eghltEgammaHoverE[idx]<0.0188*eghltEgammaSuperClusterEnergy[idx]+2.06);
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.159:abs(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.0197);
	if(seleletightideg) seleletightidegidx.push_back(idx);

      } // End of loop on egamma reco objects
            
      // Loop beginning on unseeded egamma reco objects
      for(unsigned int egidx=0; egidx<egusRecoN; egidx++) {
	
	unsigned int idx = sortedegusidx[egidx];
	noselegusidx.push_back(idx);

	basicselegus = true;
	basicselegus *= (TMath::Abs(egusRecoEta[idx])<2.65);
	basicselegus *= (egusRecoPt[idx]>=15);
	if(basicselegus) basicselegusidx.push_back(idx);

      	selelevetoidegus = true;
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<2.65);
	selelevetoidegus *= (egusRecoPt[idx]>=15);
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0126:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0457);
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_DetaSeed[idx])<0.00463:abs(egushltEgammaGsfTrackVars_DetaSeed[idx])<0.00814);
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_Dphi[idx])<0.148:abs(egushltEgammaGsfTrackVars_Dphi[idx])<0.19);
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.05*egushltEgammaSuperClusterEnergy[idx]+1.16:egushltEgammaHoverE[idx]<0.05*egushltEgammaSuperClusterEnergy[idx]+2.54);
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.209:abs(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.0132);
	if(selelevetoidegus) selelevetoidegusidx.push_back(idx);

      	seleletightidegus = true;
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<2.65);
	seleletightidegus *= (egusRecoPt[idx]>=15);
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0104:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0353);
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_DetaSeed[idx])<0.00255:abs(egushltEgammaGsfTrackVars_DetaSeed[idx])<0.00501);
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_Dphi[idx])<0.022:abs(egushltEgammaGsfTrackVars_Dphi[idx])<0.0236);
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.026*egushltEgammaSuperClusterEnergy[idx]+1.15:egushltEgammaHoverE[idx]<0.0188*egushltEgammaSuperClusterEnergy[idx]+2.06);
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.159:abs(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.0197);
	if(seleletightidegus) seleletightidegusidx.push_back(idx);

      } // End of loop on unseeded egamma objects

      selelevetozwindidegusidx = selelevetoidegusidx;
      selelevetozoppoidegusidx = selelevetoidegusidx;
      
      fillhistinevent("nosel", noselegidx);
      fillhistineventunseeded("noselus", noselegusidx);
      fillhistinevent("basicsel", basicselegidx);
      fillhistineventunseeded("basicselus", basicselegusidx);
      fillhistinevent("selelevetoid", selelevetoidegidx);
      fillhistineventunseeded("selelevetoidus", selelevetoidegusidx);
      if(selelevetozwindidegusidx.size()>=2) {
	TLorentzVector elelead, elesublead;
	elelead.SetPtEtaPhiM(egusRecoPt[selelevetozwindidegusidx[0]],egusRecoEta[selelevetozwindidegusidx[0]],egusRecoPhi[selelevetozwindidegusidx[0]],0.0005);
	elesublead.SetPtEtaPhiM(egusRecoPt[selelevetozwindidegusidx[1]],egusRecoEta[selelevetozwindidegusidx[1]],egusRecoPhi[selelevetozwindidegusidx[1]],0.0005);
	if((elelead+elesublead).M()>75 && (elelead+elesublead).M()<95) {
	  fillhistineventunseeded("selelevetozwindidus", selelevetozwindidegusidx);
	}
	else {
	  fillhistineventunseeded("selelevetozoppoidus", selelevetozoppoidegusidx);
	}
      }
      else {
	fillhistineventunseeded("selelevetozoppoidus", selelevetozoppoidegusidx);
      }
      fillhistinevent("seleletightid", seleletightidegidx);
      fillhistineventunseeded("seleletightidus", seleletightidegusidx);
      
    } // End of condition requiring atleast one egReco object

    // Count events passing selections
    if(noselegidx.size()>=2) nosel++;
    if(noselegusidx.size()>=2) noselus++;
    if(basicselegidx.size()>=2) basicsel++;
    if(basicselegusidx.size()>=2) basicselus++;
    if(selelevetoidegidx.size()>=2) selelevetoid++;
    if(selelevetoidegusidx.size()>=2) selelevetoidus++;
    if(selelevetozwindidegusidx.size()>=2) {
      TLorentzVector elelead, elesublead;
      elelead.SetPtEtaPhiM(egusRecoPt[selelevetozwindidegusidx[0]],egusRecoEta[selelevetozwindidegusidx[0]],egusRecoPhi[selelevetozwindidegusidx[0]],0.0005);
      elesublead.SetPtEtaPhiM(egusRecoPt[selelevetozwindidegusidx[1]],egusRecoEta[selelevetozwindidegusidx[1]],egusRecoPhi[selelevetozwindidegusidx[1]],0.0005);
      if((elelead+elesublead).M()>75 && (elelead+elesublead).M()<95) {
        selelevetozwindidus++;
      }
      else {
        selelevetozoppoidus++;
      }
    }
    else {
      selelevetozoppoidus++;
    }
    if(seleletightidegidx.size()>=2) seleletightid++;
    if(seleletightidegusidx.size()>=2) seleletightidus++;
    
    // Clear all the vectors
    gennoselegidx.clear();
    genbarselegidx.clear();
    genbarselptgt10egidx.clear();
    noselegidx.clear();
    noselegusidx.clear();
    basicselegidx.clear();
    basicselegusidx.clear();
    selelevetoidegidx.clear();
    selelevetoidegusidx.clear();
    selelevetozwindidegusidx.clear();
    selelevetozoppoidegusidx.clear();
    seleletightidegidx.clear();
    seleletightidegusidx.clear();

  } // End of loop on events

  cout<<totEntries<<"\t"<<nosel<<"\t"<<noselus<<"\t"<<basicsel<<"\t"<<basicselus<<"\t"<<selelevetoid<<"\t"<<selelevetoidus<<"\t"<<selelevetozwindidus<<"\t"<<selelevetozoppoidus<<"\t"<<seleletightid<<"\t"<<seleletightidus<<endl;
}

// Function to do gen matching
// Currently doing genmatching for lead pT only in mueg type events
// Modify this to return two indices for egeg type events
pair<int,int> data_robustanalyzer::doGenMatchingUnseeded(vector<int> egusidx) {

  if(!isMC) {
    cout<<"Error in doingCannot do gen matching. Not MC file."<<endl;
    return make_pair(-1,-1);
  }
  int el1idx=-1, el2idx=-1, egn=0;
  for(unsigned int gen=0; gen<genLepN; gen++) {
    if(abs(genLepMomPid[gen])==9000007 || abs(genLepMomPid[gen])==23) {
      if(abs(genLepPid[gen])==11) {
	egn++;
	if(el1idx==-1) {
	  el1idx = gen;
	}
	else if(el2idx==-1) {
	  el2idx = gen;
	}
	else {
	  cout<<"Warning!!! More than two hard process gen electron in event."<<endl;
	}
      }
    }
  }
  
  // Require gen signature with 2e
  // Unseeded egamma objects are ordered with respect to pT
  int genmcheg1uspos = -1, genmcheg2uspos = -1;
  if(egn==2 && egusidx.size()>0) {
    // Find the egusidx with the best angular match
    for(unsigned int eg=0; eg<egusidx.size(); eg++) {
      if(genmcheg1uspos!=-1 && genmcheg2uspos!=-1) {
	break;
      }
      // Condition to do gen match
      if(abs(egusRecoEta[egusidx[eg]])<1.479) {
	if(abs(egusRecoEta[egusidx[eg]]-genLepEta[el1idx])<0.06 && abs(egusRecoEta[egusidx[eg]]-genLepEta[el2idx])<0.06) {
	  if(abs(egusRecoEta[egusidx[eg]]-genLepEta[el1idx])<abs(egusRecoEta[egusidx[eg]]-genLepEta[el2idx]) && genmcheg1uspos!=-1) {
	    genmcheg1uspos = egusidx[eg];
	  }
	  else {
	    genmcheg2uspos = egusidx[eg];
	  }
	}
	if(abs(egusRecoEta[egusidx[eg]]-genLepEta[el1idx])<0.06 && genmcheg1uspos!=-1) {
	  genmcheg1uspos = egusidx[eg];
	}
	if(abs(egusRecoEta[egusidx[eg]]-genLepEta[el2idx])<0.06 && genmcheg2uspos!=-1) {
	  genmcheg2uspos = egusidx[eg];
	}
      }
      else {
	if(abs(egusRecoEta[egusidx[eg]]-genLepEta[el1idx])<0.04 && abs(egusRecoEta[egusidx[eg]]-genLepEta[el2idx])<0.04) {
	  if(abs(egusRecoEta[egusidx[eg]]-genLepEta[el1idx])<abs(egusRecoEta[egusidx[eg]]-genLepEta[el2idx]) && genmcheg1uspos!=-1) {
	    genmcheg1uspos = egusidx[eg];
	  }
	  else {
	    genmcheg2uspos = egusidx[eg];
	  }
	}
	if(abs(egusRecoEta[egusidx[eg]]-genLepEta[el1idx])<0.04 && genmcheg1uspos!=-1) {
	  genmcheg1uspos = egusidx[eg];
	}
	if(abs(egusRecoEta[egusidx[eg]]-genLepEta[el2idx])<0.04 && genmcheg2uspos!=-1) {
	  genmcheg2uspos = egusidx[eg];
	}
      }
    }
  }
  else {
    return make_pair(-1,-1);
  }
  
  return make_pair(genmcheg1uspos,genmcheg2uspos);
}

// Function to fill a set of histograms for gen particles
void data_robustanalyzer::fillgenhistinevent(TString selection, vector<int> egidx) {

  TH1F* egmult = (TH1F*) outfile->Get(selection+"geneg_egmult");
  TH1F* pt = (TH1F*) outfile->Get(selection+"geneg_pt");
  TH1F* eta = (TH1F*) outfile->Get(selection+"geneg_eta");
  TH1F* phi = (TH1F*) outfile->Get(selection+"geneg_phi");
  TH1F* log10d0 = (TH1F*) outfile->Get(selection+"geneg_log10d0");
  TH1F* t1 = (TH1F*) outfile->Get(selection+"geneg_t1");
  TH1F* t0 = (TH1F*) outfile->Get(selection+"geneg_t0");
  TH1F* t1mt0 = (TH1F*) outfile->Get(selection+"geneg_t1mt0");
  TH1F* leadegpt = (TH1F*) outfile->Get(selection+"geneg_leadpt");
  TH1F* leadegeta = (TH1F*) outfile->Get(selection+"geneg_leadeta");
  TH1F* leadegphi = (TH1F*) outfile->Get(selection+"geneg_leadphi");
  TH1F* leadeglog10d0 = (TH1F*) outfile->Get(selection+"geneg_leadlog10d0");
  TH1F* subleadegpt = (TH1F*) outfile->Get(selection+"geneg_subleadpt");
  TH1F* subleadegeta = (TH1F*) outfile->Get(selection+"geneg_subleadeta");
  TH1F* subleadegphi = (TH1F*) outfile->Get(selection+"geneg_subleadphi");
  TH1F* subleadeglog10d0 = (TH1F*) outfile->Get(selection+"geneg_subleadlog10d0");

  vector<int> sortedgenegidx = egidx;
  sort(&sortedgenegidx[0], genLepPt, egidx.size()); // Verified that this works well
  
  egmult->Fill(egidx.size());
  for(unsigned int i=0; i<egidx.size(); i++) {
    TLorentzVector el;
    el.SetPtEtaPhiM(genLepPt[egidx[i]],genLepEta[egidx[i]],genLepPhi[egidx[i]],0.0005);
    double d0 = genLepVx[egidx[i]]*el.Py()-genLepVy[egidx[i]]*el.Px();
    pt->Fill(genLepPt[egidx[i]]);
    eta->Fill(genLepEta[egidx[i]]);
    phi->Fill(genLepPhi[egidx[i]]);
    log10d0->Fill(TMath::Log10(TMath::Abs(d0)));
    if(genLepTimeAct[egidx[i]]>(-1e9)) {
      t1->Fill(genLepTimeAct[egidx[i]]);
      t0->Fill(genLepTimeLight[egidx[i]]);
      t1mt0->Fill(genLepTimeDiff[egidx[i]]);
    }
  }
  TLorentzVector leadeg, subleadeg;
  leadeg.SetPtEtaPhiM(genLepPt[sortedgenegidx[0]],genLepEta[sortedgenegidx[0]],genLepPhi[sortedgenegidx[0]],0.0005);
  subleadeg.SetPtEtaPhiM(genLepPt[sortedgenegidx[1]],genLepEta[sortedgenegidx[1]],genLepPhi[sortedgenegidx[1]],0.0005);
  double leadd0 = genLepVx[sortedgenegidx[0]]*leadeg.Py()-genLepVy[sortedgenegidx[0]]*leadeg.Px();
  double subleadd0 = genLepVx[sortedgenegidx[1]]*subleadeg.Py()-genLepVy[sortedgenegidx[1]]*subleadeg.Px();
  leadd0 /= genLepPt[sortedgenegidx[0]];
  subleadd0 /= genLepPt[sortedgenegidx[1]];
  leadegpt->Fill(genLepPt[sortedgenegidx[0]]);
  leadegeta->Fill(genLepEta[sortedgenegidx[0]]);
  leadegphi->Fill(genLepPhi[sortedgenegidx[0]]);
  leadeglog10d0->Fill(TMath::Log10(TMath::Abs(leadd0)));
  subleadegpt->Fill(genLepPt[sortedgenegidx[1]]);
  subleadegeta->Fill(genLepEta[sortedgenegidx[1]]);
  subleadegphi->Fill(genLepPhi[sortedgenegidx[1]]);
  subleadeglog10d0->Fill(TMath::Log10(TMath::Abs(subleadd0)));

  sortedgenegidx.clear();
}

// Function to fill a set of histograms in the event
void data_robustanalyzer::fillhistinevent(TString selection, vector<int> egidx) {

  // nothing here for now
  TH1F* egmult = (TH1F*) outfile->Get(selection+"recoeg_egmult");
  TH1F* leadegpt = (TH1F*) outfile->Get(selection+"recoeg_leadegpt");
  TH1F* leadegeta = (TH1F*) outfile->Get(selection+"recoeg_leadegeta");
  TH1F* leadegphi = (TH1F*) outfile->Get(selection+"recoeg_leadegphi");
  TH1F* subleadegpt = (TH1F*) outfile->Get(selection+"recoeg_subleadegpt");
  TH1F* subleadegeta = (TH1F*) outfile->Get(selection+"recoeg_subleadegeta");
  TH1F* subleadegphi = (TH1F*) outfile->Get(selection+"recoeg_subleadegphi");

  // Get barrel variables - lead pT e/gamma
  TH1F* recoeb_leadegclustershape = (TH1F*) outfile->Get(selection+"recoeb_leadegclustershape");
  TH1F* recoeb_leadegin5x5clusshape = (TH1F*) outfile->Get(selection+"recoeb_leadegin5x5clusshape");
  TH1F* recoeb_leadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoeb_leadegin5x5noiseclnd");
  TH1F* recoeb_leadegscenergy = (TH1F*) outfile->Get(selection+"recoeb_leadegscenergy");
  TH1F* recoeb_leadeghovere = (TH1F*) outfile->Get(selection+"recoeb_leadeghovere");
  TH1F* recoeb_leadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"recoeb_leadeghovereoversupcluse");
  TH1F* recoeb_leadegecalpfclustiso = (TH1F*) outfile->Get(selection+"recoeb_leadegecalpfclustiso");
  TH1F* recoeb_leadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeb_leadegecalpfclustisoovere");
  TH1F* recoeb_leadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"recoeb_leadeghcalpfclustiso");
  TH1F* recoeb_leadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeb_leadeghcalpfclustisoovere");
  TH1F* recoeb_leadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoeb_leadegpixelmchvar_s2");
  TH1F* recoeb_leadegtrkiso = (TH1F*) outfile->Get(selection+"recoeb_leadegtrkiso");
  TH1F* recoeb_leadegchi2 = (TH1F*) outfile->Get(selection+"recoeb_leadegchi2");
  TH1F* recoeb_leadegdeta = (TH1F*) outfile->Get(selection+"recoeb_leadegdeta");
  TH1F* recoeb_leadegdetaseed = (TH1F*) outfile->Get(selection+"recoeb_leadegdetaseed");
  TH1F* recoeb_leadegdphi = (TH1F*) outfile->Get(selection+"recoeb_leadegdphi");
  TH1F* recoeb_leadegmhits = (TH1F*) outfile->Get(selection+"recoeb_leadegmhits");
  TH1F* recoeb_leadegnlayerit = (TH1F*) outfile->Get(selection+"recoeb_leadegnlayerit");
  TH1F* recoeb_leadegooeseedoop = (TH1F*) outfile->Get(selection+"recoeb_leadegooeseedoop");
  TH1F* recoeb_leadegooesclsoop = (TH1F*) outfile->Get(selection+"recoeb_leadegooesclsoop");
  TH1F* recoeb_leadegvalhits = (TH1F*) outfile->Get(selection+"recoeb_leadegvalhits");  
  TH1F* recoeb_leadegseedclustime = (TH1F*) outfile->Get(selection+"recoeb_leadegseedclustime");  
  
  // Get barrel variables - sublead pT e/gamma
  TH1F* recoeb_subleadegclustershape = (TH1F*) outfile->Get(selection+"recoeb_subleadegclustershape");
  TH1F* recoeb_subleadegin5x5clusshape = (TH1F*) outfile->Get(selection+"recoeb_subleadegin5x5clusshape");
  TH1F* recoeb_subleadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoeb_subleadegin5x5noiseclnd");
  TH1F* recoeb_subleadegscenergy = (TH1F*) outfile->Get(selection+"recoeb_subleadegscenergy");
  TH1F* recoeb_subleadeghovere = (TH1F*) outfile->Get(selection+"recoeb_subleadeghovere");
  TH1F* recoeb_subleadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"recoeb_subleadeghovereoversupcluse");
  TH1F* recoeb_subleadegecalpfclustiso = (TH1F*) outfile->Get(selection+"recoeb_subleadegecalpfclustiso");
  TH1F* recoeb_subleadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeb_subleadegecalpfclustisoovere");
  TH1F* recoeb_subleadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"recoeb_subleadeghcalpfclustiso");
  TH1F* recoeb_subleadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeb_subleadeghcalpfclustisoovere");
  TH1F* recoeb_subleadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoeb_subleadegpixelmchvar_s2");
  TH1F* recoeb_subleadegtrkiso = (TH1F*) outfile->Get(selection+"recoeb_subleadegtrkiso");
  TH1F* recoeb_subleadegchi2 = (TH1F*) outfile->Get(selection+"recoeb_subleadegchi2");
  TH1F* recoeb_subleadegdeta = (TH1F*) outfile->Get(selection+"recoeb_subleadegdeta");
  TH1F* recoeb_subleadegdetaseed = (TH1F*) outfile->Get(selection+"recoeb_subleadegdetaseed");
  TH1F* recoeb_subleadegdphi = (TH1F*) outfile->Get(selection+"recoeb_subleadegdphi");
  TH1F* recoeb_subleadegmhits = (TH1F*) outfile->Get(selection+"recoeb_subleadegmhits");
  TH1F* recoeb_subleadegnlayerit = (TH1F*) outfile->Get(selection+"recoeb_subleadegnlayerit");
  TH1F* recoeb_subleadegooeseedoop = (TH1F*) outfile->Get(selection+"recoeb_subleadegooeseedoop");
  TH1F* recoeb_subleadegooesclsoop = (TH1F*) outfile->Get(selection+"recoeb_subleadegooesclsoop");
  TH1F* recoeb_subleadegvalhits = (TH1F*) outfile->Get(selection+"recoeb_subleadegvalhits");  
  TH1F* recoeb_subleadegseedclustime = (TH1F*) outfile->Get(selection+"recoeb_subleadegseedclustime");  

  // invariant mass - barrel
  TH1F* recoeb_leadsubleadM = (TH1F*) outfile->Get(selection+"recoeb_leadsubleadM");

  // Get end-cap variables - lead pT e/gamma
  TH1F* recoee_leadegclustershape = (TH1F*) outfile->Get(selection+"recoee_leadegclustershape");
  TH1F* recoee_leadegin5x5clusshape = (TH1F*) outfile->Get(selection+"recoee_leadegin5x5clusshape");
  TH1F* recoee_leadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoee_leadegin5x5noiseclnd");
  TH1F* recoee_leadegscenergy = (TH1F*) outfile->Get(selection+"recoee_leadegscenergy");
  TH1F* recoee_leadeghovere = (TH1F*) outfile->Get(selection+"recoee_leadeghovere");
  TH1F* recoee_leadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"recoee_leadeghovereoversupcluse");
  TH1F* recoee_leadegecalpfclustiso = (TH1F*) outfile->Get(selection+"recoee_leadegecalpfclustiso");
  TH1F* recoee_leadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoee_leadegecalpfclustisoovere");
  TH1F* recoee_leadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"recoee_leadeghcalpfclustiso");
  TH1F* recoee_leadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoee_leadeghcalpfclustisoovere");
  TH1F* recoee_leadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoee_leadegpixelmchvar_s2");
  TH1F* recoee_leadegtrkiso = (TH1F*) outfile->Get(selection+"recoee_leadegtrkiso");
  TH1F* recoee_leadegchi2 = (TH1F*) outfile->Get(selection+"recoee_leadegchi2");
  TH1F* recoee_leadegdeta = (TH1F*) outfile->Get(selection+"recoee_leadegdeta");
  TH1F* recoee_leadegdetaseed = (TH1F*) outfile->Get(selection+"recoee_leadegdetaseed");
  TH1F* recoee_leadegdphi = (TH1F*) outfile->Get(selection+"recoee_leadegdphi");
  TH1F* recoee_leadegmhits = (TH1F*) outfile->Get(selection+"recoee_leadegmhits");
  TH1F* recoee_leadegnlayerit = (TH1F*) outfile->Get(selection+"recoee_leadegnlayerit");
  TH1F* recoee_leadegooeseedoop = (TH1F*) outfile->Get(selection+"recoee_leadegooeseedoop");
  TH1F* recoee_leadegooesclsoop = (TH1F*) outfile->Get(selection+"recoee_leadegooesclsoop");
  TH1F* recoee_leadegvalhits = (TH1F*) outfile->Get(selection+"recoee_leadegvalhits");
  TH1F* recoee_leadegseedclustime = (TH1F*) outfile->Get(selection+"recoee_leadegseedclustime");  
  
  // Get end-cap variables - sublead pT e/gamma
  TH1F* recoee_subleadegclustershape = (TH1F*) outfile->Get(selection+"recoee_subleadegclustershape");
  TH1F* recoee_subleadegin5x5clusshape = (TH1F*) outfile->Get(selection+"recoee_subleadegin5x5clusshape");
  TH1F* recoee_subleadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoee_subleadegin5x5noiseclnd");
  TH1F* recoee_subleadegscenergy = (TH1F*) outfile->Get(selection+"recoee_subleadegscenergy");
  TH1F* recoee_subleadeghovere = (TH1F*) outfile->Get(selection+"recoee_subleadeghovere");
  TH1F* recoee_subleadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"recoee_subleadeghovereoversupcluse");
  TH1F* recoee_subleadegecalpfclustiso = (TH1F*) outfile->Get(selection+"recoee_subleadegecalpfclustiso");
  TH1F* recoee_subleadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoee_subleadegecalpfclustisoovere");
  TH1F* recoee_subleadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"recoee_subleadeghcalpfclustiso");
  TH1F* recoee_subleadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoee_subleadeghcalpfclustisoovere");
  TH1F* recoee_subleadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoee_subleadegpixelmchvar_s2");
  TH1F* recoee_subleadegtrkiso = (TH1F*) outfile->Get(selection+"recoee_subleadegtrkiso");
  TH1F* recoee_subleadegchi2 = (TH1F*) outfile->Get(selection+"recoee_subleadegchi2");
  TH1F* recoee_subleadegdeta = (TH1F*) outfile->Get(selection+"recoee_subleadegdeta");
  TH1F* recoee_subleadegdetaseed = (TH1F*) outfile->Get(selection+"recoee_subleadegdetaseed");
  TH1F* recoee_subleadegdphi = (TH1F*) outfile->Get(selection+"recoee_subleadegdphi");
  TH1F* recoee_subleadegmhits = (TH1F*) outfile->Get(selection+"recoee_subleadegmhits");
  TH1F* recoee_subleadegnlayerit = (TH1F*) outfile->Get(selection+"recoee_subleadegnlayerit");
  TH1F* recoee_subleadegooeseedoop = (TH1F*) outfile->Get(selection+"recoee_subleadegooeseedoop");
  TH1F* recoee_subleadegooesclsoop = (TH1F*) outfile->Get(selection+"recoee_subleadegooesclsoop");
  TH1F* recoee_subleadegvalhits = (TH1F*) outfile->Get(selection+"recoee_subleadegvalhits");
  TH1F* recoee_subleadegseedclustime = (TH1F*) outfile->Get(selection+"recoee_subleadegseedclustime");  
  
  // invariant mass - end-cap
  TH1F* recoee_leadsubleadM = (TH1F*) outfile->Get(selection+"recoee_leadsubleadM");

  if(egidx.size()>0) {
    egmult->Fill(egidx.size());
    leadegpt->Fill(egRecoPt[egidx[0]]);
    leadegeta->Fill(egRecoEta[egidx[0]]);
    leadegphi->Fill(egRecoPhi[egidx[0]]);

    // Fill barrel variables
    if(TMath::Abs(egRecoEta[egidx[0]])<=1.479) {
      recoeb_leadegclustershape->Fill(eghltEgammaClusterShape[egidx[0]]);
      recoeb_leadegin5x5clusshape->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[0]]);
      recoeb_leadegin5x5noiseclnd->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[0]]);
      recoeb_leadegscenergy->Fill(eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_leadeghovere->Fill(eghltEgammaHoverE[egidx[0]]);
      recoeb_leadeghovereoversupcluse->Fill(eghltEgammaHoverE[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_leadegecalpfclustiso->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]);
      recoeb_leadegecalpfclustisoovere->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_leadeghcalpfclustiso->Fill(eghltEgammaHcalPFClusterIso[egidx[0]]);
      recoeb_leadeghcalpfclustisoovere->Fill(eghltEgammaHcalPFClusterIso[egidx[0]]/eghltEgammaHoverE[egidx[0]]);
      if(eghltEcalSeedClusterTime[egidx[0]]!=0) recoeb_leadegseedclustime->Fill(eghltEcalSeedClusterTime[egidx[0]]);
      if(eghltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoeb_leadegpixelmchvar_s2->Fill(eghltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoeb_leadegpixelmchvar_s2->Fill(-5);
      }
      if(eghltEgammaGsfTrackVars_Deta[egidx[0]]<999999) {
	recoeb_leadegtrkiso->Fill(eghltEgammaEleGsfTrackIso[egidx[0]]);
	recoeb_leadegchi2->Fill(eghltEgammaGsfTrackVars_Chi2[egidx[0]]);
	recoeb_leadegdeta->Fill(eghltEgammaGsfTrackVars_Deta[egidx[0]]);
	recoeb_leadegdetaseed->Fill(eghltEgammaGsfTrackVars_DetaSeed[egidx[0]]);
	recoeb_leadegdphi->Fill(eghltEgammaGsfTrackVars_Dphi[egidx[0]]);
	recoeb_leadegmhits->Fill(eghltEgammaGsfTrackVars_MissingHits[egidx[0]]);
	recoeb_leadegnlayerit->Fill(eghltEgammaGsfTrackVars_NLayerIT[egidx[0]]);
	recoeb_leadegooeseedoop->Fill(eghltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[0]]);
	recoeb_leadegooesclsoop->Fill(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[0]]);
	recoeb_leadegvalhits->Fill(eghltEgammaGsfTrackVars_ValidHits[egidx[0]]);
      }
      else {
	recoeb_leadegtrkiso->Fill(-100);
	recoeb_leadegchi2->Fill(-20);
	recoeb_leadegdeta->Fill(-5);
	recoeb_leadegdetaseed->Fill(-5);
	recoeb_leadegdphi->Fill(-5);
	recoeb_leadegmhits->Fill(-15);
	recoeb_leadegnlayerit->Fill(-15);
	recoeb_leadegooeseedoop->Fill(-50);
	recoeb_leadegooesclsoop->Fill(-50);
	recoeb_leadegvalhits->Fill(-15);
      }
    } // End of filling barrel variables
    else { // Fill end-cap variables
      recoee_leadegclustershape->Fill(eghltEgammaClusterShape[egidx[0]]);
      recoee_leadegin5x5clusshape->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[0]]);
      recoee_leadegin5x5noiseclnd->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[0]]);
      recoee_leadegscenergy->Fill(eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_leadeghovere->Fill(eghltEgammaHoverE[egidx[0]]);
      recoee_leadeghovereoversupcluse->Fill(eghltEgammaHoverE[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_leadegecalpfclustiso->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]);
      recoee_leadegecalpfclustisoovere->Fill(eghltEgammaEcalPFClusterIso[egidx[0]]/eghltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_leadeghcalpfclustiso->Fill(eghltEgammaHcalPFClusterIso[egidx[0]]);	
      recoee_leadeghcalpfclustisoovere->Fill(eghltEgammaHcalPFClusterIso[egidx[0]]/eghltEgammaHoverE[egidx[0]]);	
      if(eghltEcalSeedClusterTime[egidx[0]]!=0) recoee_leadegseedclustime->Fill(eghltEcalSeedClusterTime[egidx[0]]);
      if(eghltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoee_leadegpixelmchvar_s2->Fill(eghltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoee_leadegpixelmchvar_s2->Fill(-5);
      }
      if(eghltEgammaGsfTrackVars_Deta[egidx[0]]<999999) {
	recoee_leadegtrkiso->Fill(eghltEgammaEleGsfTrackIso[egidx[0]]);
	recoee_leadegchi2->Fill(eghltEgammaGsfTrackVars_Chi2[egidx[0]]);
	recoee_leadegdeta->Fill(eghltEgammaGsfTrackVars_Deta[egidx[0]]);
	recoee_leadegdetaseed->Fill(eghltEgammaGsfTrackVars_DetaSeed[egidx[0]]);
	recoee_leadegdphi->Fill(eghltEgammaGsfTrackVars_Dphi[egidx[0]]);
	recoee_leadegmhits->Fill(eghltEgammaGsfTrackVars_MissingHits[egidx[0]]);
	recoee_leadegnlayerit->Fill(eghltEgammaGsfTrackVars_NLayerIT[egidx[0]]);
	recoee_leadegooeseedoop->Fill(eghltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[0]]);
	recoee_leadegooesclsoop->Fill(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[0]]);
	recoee_leadegvalhits->Fill(eghltEgammaGsfTrackVars_ValidHits[egidx[0]]);
      }
      else {
	recoee_leadegtrkiso->Fill(-100);
	recoee_leadegchi2->Fill(-20);
	recoee_leadegdeta->Fill(-5);
	recoee_leadegdetaseed->Fill(-5);
	recoee_leadegdphi->Fill(-5);
	recoee_leadegmhits->Fill(-15);
	recoee_leadegnlayerit->Fill(-15);
	recoee_leadegooeseedoop->Fill(-50);
	recoee_leadegooesclsoop->Fill(-50);
	recoee_leadegvalhits->Fill(-15);
      }
    } // End of filling end-cap variables


  } // End of condition requiring atleast one eg object
  
  if(egidx.size()>=2) { // Condition requiring atleast two eg object
    TLorentzVector leadeg, subleadeg;
    leadeg.SetPtEtaPhiM(egRecoPt[egidx[0]],egRecoEta[egidx[0]],egRecoPhi[egidx[0]],0.106);
    subleadeg.SetPtEtaPhiM(egRecoPt[egidx[1]],egRecoEta[egidx[1]],egRecoPhi[egidx[1]],0.106);

    subleadegpt->Fill(egRecoPt[egidx[1]]);
    subleadegeta->Fill(egRecoEta[egidx[1]]);
    subleadegphi->Fill(egRecoPhi[egidx[1]]);

    // Fill barrel variables
    if(TMath::Abs(egRecoEta[egidx[1]])<=1.479) {
      recoeb_subleadegclustershape->Fill(eghltEgammaClusterShape[egidx[1]]);
      recoeb_subleadegin5x5clusshape->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[1]]);
      recoeb_subleadegin5x5noiseclnd->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[1]]);
      recoeb_subleadegscenergy->Fill(eghltEgammaSuperClusterEnergy[egidx[1]]);
      recoeb_subleadeghovere->Fill(eghltEgammaHoverE[egidx[1]]);
      recoeb_subleadeghovereoversupcluse->Fill(eghltEgammaHoverE[egidx[1]]/eghltEgammaSuperClusterEnergy[egidx[1]]);
      recoeb_subleadegecalpfclustiso->Fill(eghltEgammaEcalPFClusterIso[egidx[1]]);
      recoeb_subleadegecalpfclustisoovere->Fill(eghltEgammaEcalPFClusterIso[egidx[1]]/eghltEgammaSuperClusterEnergy[egidx[1]]);
      recoeb_subleadeghcalpfclustiso->Fill(eghltEgammaHcalPFClusterIso[egidx[1]]);
      recoeb_subleadeghcalpfclustisoovere->Fill(eghltEgammaHcalPFClusterIso[egidx[1]]/eghltEgammaHoverE[egidx[1]]);
      if(eghltEcalSeedClusterTime[egidx[1]]!=0) recoeb_subleadegseedclustime->Fill(eghltEcalSeedClusterTime[egidx[1]]);
      if(eghltEgammaPixelMatchVars_s2[egidx[1]]<TMath::Power(10,36)) {
	recoeb_subleadegpixelmchvar_s2->Fill(eghltEgammaPixelMatchVars_s2[egidx[1]]);
      }
      else {
	recoeb_subleadegpixelmchvar_s2->Fill(-5);
      }
      if(eghltEgammaGsfTrackVars_Deta[egidx[1]]<999999) {
	recoeb_subleadegtrkiso->Fill(eghltEgammaEleGsfTrackIso[egidx[1]]);
	recoeb_subleadegchi2->Fill(eghltEgammaGsfTrackVars_Chi2[egidx[1]]);
	recoeb_subleadegdeta->Fill(eghltEgammaGsfTrackVars_Deta[egidx[1]]);
	recoeb_subleadegdetaseed->Fill(eghltEgammaGsfTrackVars_DetaSeed[egidx[1]]);
	recoeb_subleadegdphi->Fill(eghltEgammaGsfTrackVars_Dphi[egidx[1]]);
	recoeb_subleadegmhits->Fill(eghltEgammaGsfTrackVars_MissingHits[egidx[1]]);
	recoeb_subleadegnlayerit->Fill(eghltEgammaGsfTrackVars_NLayerIT[egidx[1]]);
	recoeb_subleadegooeseedoop->Fill(eghltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[1]]);
	recoeb_subleadegooesclsoop->Fill(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[1]]);
	recoeb_subleadegvalhits->Fill(eghltEgammaGsfTrackVars_ValidHits[egidx[1]]);
      }
      else {
	recoeb_subleadegtrkiso->Fill(-100);
	recoeb_subleadegchi2->Fill(-20);
	recoeb_subleadegdeta->Fill(-5);
	recoeb_subleadegdetaseed->Fill(-5);
	recoeb_subleadegdphi->Fill(-5);
	recoeb_subleadegmhits->Fill(-15);
	recoeb_subleadegnlayerit->Fill(-15);
	recoeb_subleadegooeseedoop->Fill(-50);
	recoeb_subleadegooesclsoop->Fill(-50);
	recoeb_subleadegvalhits->Fill(-15);
      }
    } // End of filling barrel variables
    else { // Fill end-cap variables
      recoee_subleadegclustershape->Fill(eghltEgammaClusterShape[egidx[1]]);
      recoee_subleadegin5x5clusshape->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[1]]);
      recoee_subleadegin5x5noiseclnd->Fill(eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[1]]);
      recoee_subleadegscenergy->Fill(eghltEgammaSuperClusterEnergy[egidx[1]]);
      recoee_subleadeghovere->Fill(eghltEgammaHoverE[egidx[1]]);
      recoee_subleadeghovereoversupcluse->Fill(eghltEgammaHoverE[egidx[1]]/eghltEgammaSuperClusterEnergy[egidx[1]]);
      recoee_subleadegecalpfclustiso->Fill(eghltEgammaEcalPFClusterIso[egidx[1]]);
      recoee_subleadegecalpfclustisoovere->Fill(eghltEgammaEcalPFClusterIso[egidx[1]]/eghltEgammaSuperClusterEnergy[egidx[1]]);
      recoee_subleadeghcalpfclustiso->Fill(eghltEgammaHcalPFClusterIso[egidx[1]]);	
      recoee_subleadeghcalpfclustisoovere->Fill(eghltEgammaHcalPFClusterIso[egidx[1]]/eghltEgammaHoverE[egidx[1]]);	
      if(eghltEcalSeedClusterTime[egidx[1]]!=0) recoee_subleadegseedclustime->Fill(eghltEcalSeedClusterTime[egidx[1]]);
      if(eghltEgammaPixelMatchVars_s2[egidx[1]]<TMath::Power(10,36)) {
	recoee_subleadegpixelmchvar_s2->Fill(eghltEgammaPixelMatchVars_s2[egidx[1]]);
      }
      else {
	recoee_subleadegpixelmchvar_s2->Fill(-5);
      }
      if(eghltEgammaGsfTrackVars_Deta[egidx[1]]<999999) {
	recoee_subleadegtrkiso->Fill(eghltEgammaEleGsfTrackIso[egidx[1]]);
	recoee_subleadegchi2->Fill(eghltEgammaGsfTrackVars_Chi2[egidx[1]]);
	recoee_subleadegdeta->Fill(eghltEgammaGsfTrackVars_Deta[egidx[1]]);
	recoee_subleadegdetaseed->Fill(eghltEgammaGsfTrackVars_DetaSeed[egidx[1]]);
	recoee_subleadegdphi->Fill(eghltEgammaGsfTrackVars_Dphi[egidx[1]]);
	recoee_subleadegmhits->Fill(eghltEgammaGsfTrackVars_MissingHits[egidx[1]]);
	recoee_subleadegnlayerit->Fill(eghltEgammaGsfTrackVars_NLayerIT[egidx[1]]);
	recoee_subleadegooeseedoop->Fill(eghltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[1]]);
	recoee_subleadegooesclsoop->Fill(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[1]]);
	recoee_subleadegvalhits->Fill(eghltEgammaGsfTrackVars_ValidHits[egidx[1]]);
      }
      else {
	recoee_subleadegtrkiso->Fill(-100);
	recoee_subleadegchi2->Fill(-20);
	recoee_subleadegdeta->Fill(-5);
	recoee_subleadegdetaseed->Fill(-5);
	recoee_subleadegdphi->Fill(-5);
	recoee_subleadegmhits->Fill(-15);
	recoee_subleadegnlayerit->Fill(-15);
	recoee_subleadegooeseedoop->Fill(-50);
	recoee_subleadegooesclsoop->Fill(-50);
	recoee_subleadegvalhits->Fill(-15);
      }
    } // End of filling end-cap variables

    // Fill invariant mass
    if(egRecoEta[egidx[0]]<1.479 && egRecoEta[egidx[1]]<1.479) {
      recoeb_leadsubleadM->Fill((leadeg+subleadeg).M());
    }
    else {
      recoee_leadsubleadM->Fill((leadeg+subleadeg).M());
    }
  } // End of condition requiring atleast two eg object

}

// Function to fill a set of histograms in the event - unseeded egamma objects
void data_robustanalyzer::fillhistineventunseeded(TString selection, vector<int> egidx) {
  
  // nothing here for now
  TH1F* egmult = (TH1F*) outfile->Get(selection+"recoegus_egmult");
  TH1F* leadegpt = (TH1F*) outfile->Get(selection+"recoegus_leadegpt");
  TH1F* leadegeta = (TH1F*) outfile->Get(selection+"recoegus_leadegeta");
  TH1F* leadegphi = (TH1F*) outfile->Get(selection+"recoegus_leadegphi");
  TH1F* subleadegpt = (TH1F*) outfile->Get(selection+"recoegus_subleadegpt");
  TH1F* subleadegeta = (TH1F*) outfile->Get(selection+"recoegus_subleadegeta");
  TH1F* subleadegphi = (TH1F*) outfile->Get(selection+"recoegus_subleadegphi");
  
  TH1F* leadgenmchegpt = (TH1F*) outfile->Get(selection+"recoegus_leadgenmchegpt");
  TH1F* leadgenmchegeta = (TH1F*) outfile->Get(selection+"recoegus_leadgenmchegeta");
  TH1F* leadgenmchegphi = (TH1F*) outfile->Get(selection+"recoegus_leadgenmchegphi");
  TH1F* subleadgenmchegpt = (TH1F*) outfile->Get(selection+"recoegus_subleadgenmchegpt");
  TH1F* subleadgenmchegeta = (TH1F*) outfile->Get(selection+"recoegus_subleadgenmchegeta");
  TH1F* subleadgenmchegphi = (TH1F*) outfile->Get(selection+"recoegus_subleadgenmchegphi");
  
  // Get barrel variables - lead pt unseeded e/gamma
  TH1F* recoeb_leadegclustershape = (TH1F*) outfile->Get(selection+"recoebus_leadegclustershape");
  TH1F* recoeb_leadegin5x5clusshape = (TH1F*) outfile->Get(selection+"recoebus_leadegin5x5clusshape");
  TH1F* recoeb_leadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoebus_leadegin5x5noiseclnd");
  TH1F* recoeb_leadegscenergy = (TH1F*) outfile->Get(selection+"recoebus_leadegscenergy");
  TH1F* recoeb_leadeghovere = (TH1F*) outfile->Get(selection+"recoebus_leadeghovere");
  TH1F* recoeb_leadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"recoebus_leadeghovereoversupcluse");
  TH1F* recoeb_leadegecalpfclustiso = (TH1F*) outfile->Get(selection+"recoebus_leadegecalpfclustiso");
  TH1F* recoeb_leadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoebus_leadegecalpfclustisoovere");
  TH1F* recoeb_leadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"recoebus_leadeghcalpfclustiso");
  TH1F* recoeb_leadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoebus_leadeghcalpfclustisoovere");
  TH1F* recoeb_leadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoebus_leadegpixelmchvar_s2");
  TH1F* recoeb_leadegtrkiso = (TH1F*) outfile->Get(selection+"recoebus_leadegtrkiso");
  TH1F* recoeb_leadegchi2 = (TH1F*) outfile->Get(selection+"recoebus_leadegchi2");
  TH1F* recoeb_leadegdeta = (TH1F*) outfile->Get(selection+"recoebus_leadegdeta");
  TH1F* recoeb_leadegdetaseed = (TH1F*) outfile->Get(selection+"recoebus_leadegdetaseed");
  TH1F* recoeb_leadegdphi = (TH1F*) outfile->Get(selection+"recoebus_leadegdphi");
  TH1F* recoeb_leadegmhits = (TH1F*) outfile->Get(selection+"recoebus_leadegmhits");
  TH1F* recoeb_leadegnlayerit = (TH1F*) outfile->Get(selection+"recoebus_leadegnlayerit");
  TH1F* recoeb_leadegooeseedoop = (TH1F*) outfile->Get(selection+"recoebus_leadegooeseedoop");
  TH1F* recoeb_leadegooesclsoop = (TH1F*) outfile->Get(selection+"recoebus_leadegooesclsoop");
  TH1F* recoeb_leadegvalhits = (TH1F*) outfile->Get(selection+"recoebus_leadegvalhits");  
  TH1F* recoeb_leadegseedclustime = (TH1F*) outfile->Get(selection+"recoebus_leadegseedclustime");  
  
  // Get barrel variables - sublead pt unseeded e/gamma
  TH1F* recoeb_subleadegclustershape = (TH1F*) outfile->Get(selection+"recoebus_subleadegclustershape");
  TH1F* recoeb_subleadegin5x5clusshape = (TH1F*) outfile->Get(selection+"recoebus_subleadegin5x5clusshape");
  TH1F* recoeb_subleadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoebus_subleadegin5x5noiseclnd");
  TH1F* recoeb_subleadegscenergy = (TH1F*) outfile->Get(selection+"recoebus_subleadegscenergy");
  TH1F* recoeb_subleadeghovere = (TH1F*) outfile->Get(selection+"recoebus_subleadeghovere");
  TH1F* recoeb_subleadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"recoebus_subleadeghovereoversupcluse");
  TH1F* recoeb_subleadegecalpfclustiso = (TH1F*) outfile->Get(selection+"recoebus_subleadegecalpfclustiso");
  TH1F* recoeb_subleadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoebus_subleadegecalpfclustisoovere");
  TH1F* recoeb_subleadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"recoebus_subleadeghcalpfclustiso");
  TH1F* recoeb_subleadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoebus_subleadeghcalpfclustisoovere");
  TH1F* recoeb_subleadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoebus_subleadegpixelmchvar_s2");
  TH1F* recoeb_subleadegtrkiso = (TH1F*) outfile->Get(selection+"recoebus_subleadegtrkiso");
  TH1F* recoeb_subleadegchi2 = (TH1F*) outfile->Get(selection+"recoebus_subleadegchi2");
  TH1F* recoeb_subleadegdeta = (TH1F*) outfile->Get(selection+"recoebus_subleadegdeta");
  TH1F* recoeb_subleadegdetaseed = (TH1F*) outfile->Get(selection+"recoebus_subleadegdetaseed");
  TH1F* recoeb_subleadegdphi = (TH1F*) outfile->Get(selection+"recoebus_subleadegdphi");
  TH1F* recoeb_subleadegmhits = (TH1F*) outfile->Get(selection+"recoebus_subleadegmhits");
  TH1F* recoeb_subleadegnlayerit = (TH1F*) outfile->Get(selection+"recoebus_subleadegnlayerit");
  TH1F* recoeb_subleadegooeseedoop = (TH1F*) outfile->Get(selection+"recoebus_subleadegooeseedoop");
  TH1F* recoeb_subleadegooesclsoop = (TH1F*) outfile->Get(selection+"recoebus_subleadegooesclsoop");
  TH1F* recoeb_subleadegvalhits = (TH1F*) outfile->Get(selection+"recoebus_subleadegvalhits");  
  TH1F* recoeb_subleadegseedclustime = (TH1F*) outfile->Get(selection+"recoebus_subleadegseedclustime");  

  // invariant mass - barrel
  TH1F* recoeb_leadsubleadM = (TH1F*) outfile->Get(selection+"recoebus_leadsubleadM");

  // Get end-cap variables - lead pt unseeded e/gamma
  TH1F* recoee_leadegclustershape = (TH1F*) outfile->Get(selection+"recoeeus_leadegclustershape");
  TH1F* recoee_leadegin5x5clusshape = (TH1F*) outfile->Get(selection+"recoeeus_leadegin5x5clusshape");
  TH1F* recoee_leadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoeeus_leadegin5x5noiseclnd");
  TH1F* recoee_leadegscenergy = (TH1F*) outfile->Get(selection+"recoeeus_leadegscenergy");
  TH1F* recoee_leadeghovere = (TH1F*) outfile->Get(selection+"recoeeus_leadeghovere");
  TH1F* recoee_leadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"recoeeus_leadeghovereoversupcluse");
  TH1F* recoee_leadegecalpfclustiso = (TH1F*) outfile->Get(selection+"recoeeus_leadegecalpfclustiso");
  TH1F* recoee_leadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeeus_leadegecalpfclustisoovere");
  TH1F* recoee_leadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"recoeeus_leadeghcalpfclustiso");
  TH1F* recoee_leadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeeus_leadeghcalpfclustisoovere");
  TH1F* recoee_leadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoeeus_leadegpixelmchvar_s2");
  TH1F* recoee_leadegtrkiso = (TH1F*) outfile->Get(selection+"recoeeus_leadegtrkiso");
  TH1F* recoee_leadegchi2 = (TH1F*) outfile->Get(selection+"recoeeus_leadegchi2");
  TH1F* recoee_leadegdeta = (TH1F*) outfile->Get(selection+"recoeeus_leadegdeta");
  TH1F* recoee_leadegdetaseed = (TH1F*) outfile->Get(selection+"recoeeus_leadegdetaseed");
  TH1F* recoee_leadegdphi = (TH1F*) outfile->Get(selection+"recoeeus_leadegdphi");
  TH1F* recoee_leadegmhits = (TH1F*) outfile->Get(selection+"recoeeus_leadegmhits");
  TH1F* recoee_leadegnlayerit = (TH1F*) outfile->Get(selection+"recoeeus_leadegnlayerit");
  TH1F* recoee_leadegooeseedoop = (TH1F*) outfile->Get(selection+"recoeeus_leadegooeseedoop");
  TH1F* recoee_leadegooesclsoop = (TH1F*) outfile->Get(selection+"recoeeus_leadegooesclsoop");
  TH1F* recoee_leadegvalhits = (TH1F*) outfile->Get(selection+"recoeeus_leadegvalhits");
  TH1F* recoee_leadegseedclustime = (TH1F*) outfile->Get(selection+"recoeeus_leadegseedclustime");  
  
  // Get end-cap variables - sublead pt unseeded e/gamma
  TH1F* recoee_subleadegclustershape = (TH1F*) outfile->Get(selection+"recoeeus_subleadegclustershape");
  TH1F* recoee_subleadegin5x5clusshape = (TH1F*) outfile->Get(selection+"recoeeus_subleadegin5x5clusshape");
  TH1F* recoee_subleadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoeeus_subleadegin5x5noiseclnd");
  TH1F* recoee_subleadegscenergy = (TH1F*) outfile->Get(selection+"recoeeus_subleadegscenergy");
  TH1F* recoee_subleadeghovere = (TH1F*) outfile->Get(selection+"recoeeus_subleadeghovere");
  TH1F* recoee_subleadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"recoeeus_subleadeghovereoversupcluse");
  TH1F* recoee_subleadegecalpfclustiso = (TH1F*) outfile->Get(selection+"recoeeus_subleadegecalpfclustiso");
  TH1F* recoee_subleadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeeus_subleadegecalpfclustisoovere");
  TH1F* recoee_subleadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"recoeeus_subleadeghcalpfclustiso");
  TH1F* recoee_subleadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"recoeeus_subleadeghcalpfclustisoovere");
  TH1F* recoee_subleadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoeeus_subleadegpixelmchvar_s2");
  TH1F* recoee_subleadegtrkiso = (TH1F*) outfile->Get(selection+"recoeeus_subleadegtrkiso");
  TH1F* recoee_subleadegchi2 = (TH1F*) outfile->Get(selection+"recoeeus_subleadegchi2");
  TH1F* recoee_subleadegdeta = (TH1F*) outfile->Get(selection+"recoeeus_subleadegdeta");
  TH1F* recoee_subleadegdetaseed = (TH1F*) outfile->Get(selection+"recoeeus_subleadegdetaseed");
  TH1F* recoee_subleadegdphi = (TH1F*) outfile->Get(selection+"recoeeus_subleadegdphi");
  TH1F* recoee_subleadegmhits = (TH1F*) outfile->Get(selection+"recoeeus_subleadegmhits");
  TH1F* recoee_subleadegnlayerit = (TH1F*) outfile->Get(selection+"recoeeus_subleadegnlayerit");
  TH1F* recoee_subleadegooeseedoop = (TH1F*) outfile->Get(selection+"recoeeus_subleadegooeseedoop");
  TH1F* recoee_subleadegooesclsoop = (TH1F*) outfile->Get(selection+"recoeeus_subleadegooesclsoop");
  TH1F* recoee_subleadegvalhits = (TH1F*) outfile->Get(selection+"recoeeus_subleadegvalhits");
  TH1F* recoee_subleadegseedclustime = (TH1F*) outfile->Get(selection+"recoeeus_subleadegseedclustime");  
  
  // invariant mass - end-cap
  TH1F* recoee_leadsubleadM = (TH1F*) outfile->Get(selection+"recoeeus_leadsubleadM");
  
  if(egidx.size()>0) {
    egmult->Fill(egidx.size());
    leadegpt->Fill(egusRecoPt[egidx[0]]);
    leadegeta->Fill(egusRecoEta[egidx[0]]);
    leadegphi->Fill(egusRecoPhi[egidx[0]]);

    // Fill barrel variables
    if(TMath::Abs(egusRecoEta[egidx[0]])<=1.479) {
      recoeb_leadegclustershape->Fill(egushltEgammaClusterShape[egidx[0]]);
      recoeb_leadegin5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[0]]);
      recoeb_leadegin5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[0]]);
      //recoeb_leadegscenergy->Fill(egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_leadeghovere->Fill(egushltEgammaHoverE[egidx[0]]);
      recoeb_leadeghovereoversupcluse->Fill(egushltEgammaHoverE[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_leadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]);
      recoeb_leadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_leadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[egidx[0]]);
      recoeb_leadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[egidx[0]]/egushltEgammaHoverE[egidx[0]]);
      if(egushltEcalSeedClusterTime[egidx[0]]!=0) recoeb_leadegseedclustime->Fill(egushltEcalSeedClusterTime[egidx[0]]);
      if(egushltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoeb_leadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoeb_leadegpixelmchvar_s2->Fill(-5);
      }
      if(egushltEgammaGsfTrackVars_Deta[egidx[0]]<999999) {
	recoeb_leadegtrkiso->Fill(egushltEgammaEleGsfTrackIso[egidx[0]]);
	recoeb_leadegchi2->Fill(egushltEgammaGsfTrackVars_Chi2[egidx[0]]);
	recoeb_leadegdeta->Fill(egushltEgammaGsfTrackVars_Deta[egidx[0]]);
	recoeb_leadegdetaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[egidx[0]]);
	recoeb_leadegdphi->Fill(egushltEgammaGsfTrackVars_Dphi[egidx[0]]);
	recoeb_leadegmhits->Fill(egushltEgammaGsfTrackVars_MissingHits[egidx[0]]);
	recoeb_leadegnlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[egidx[0]]);
	recoeb_leadegooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[0]]);
	recoeb_leadegooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[0]]);
	//recoeb_leadegvalhits->Fill(egushltEgammaGsfTrackVars_ValidHits[egidx[0]]);
      }
      else {
	recoeb_leadegtrkiso->Fill(-100);
	recoeb_leadegchi2->Fill(-20);
	recoeb_leadegdeta->Fill(-5);
	recoeb_leadegdetaseed->Fill(-5);
	recoeb_leadegdphi->Fill(-5);
	recoeb_leadegmhits->Fill(-15);
	recoeb_leadegnlayerit->Fill(-15);
	recoeb_leadegooeseedoop->Fill(-50);
	recoeb_leadegooesclsoop->Fill(-50);
	//recoeb_leadegvalhits->Fill(-15);
      }
    } // End of filling barrel variables
    
    else { // Fill end-cap variables
      recoee_leadegclustershape->Fill(egushltEgammaClusterShape[egidx[0]]);
      recoee_leadegin5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[0]]);
      recoee_leadegin5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[0]]);
      //recoee_leadegscenergy->Fill(egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_leadeghovere->Fill(egushltEgammaHoverE[egidx[0]]);
      recoee_leadeghovereoversupcluse->Fill(egushltEgammaHoverE[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_leadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]);
      recoee_leadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_leadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[egidx[0]]);	
      recoee_leadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[egidx[0]]/egushltEgammaHoverE[egidx[0]]);	
      if(egushltEcalSeedClusterTime[egidx[0]]!=0) recoee_leadegseedclustime->Fill(egushltEcalSeedClusterTime[egidx[0]]);
      if(egushltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoee_leadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoee_leadegpixelmchvar_s2->Fill(-5);
      }
      if(egushltEgammaGsfTrackVars_Deta[egidx[0]]<999999) {
	recoee_leadegtrkiso->Fill(egushltEgammaEleGsfTrackIso[egidx[0]]);
	recoee_leadegchi2->Fill(egushltEgammaGsfTrackVars_Chi2[egidx[0]]);
	recoee_leadegdeta->Fill(egushltEgammaGsfTrackVars_Deta[egidx[0]]);
	recoee_leadegdetaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[egidx[0]]);
	recoee_leadegdphi->Fill(egushltEgammaGsfTrackVars_Dphi[egidx[0]]);
	recoee_leadegmhits->Fill(egushltEgammaGsfTrackVars_MissingHits[egidx[0]]);
	recoee_leadegnlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[egidx[0]]);
	recoee_leadegooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[0]]);
	recoee_leadegooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[0]]);
	//recoee_leadegvalhits->Fill(egushltEgammaGsfTrackVars_ValidHits[egidx[0]]);
      }
      else {
	recoee_leadegtrkiso->Fill(-100);
	recoee_leadegchi2->Fill(-20);
	recoee_leadegdeta->Fill(-5);
	recoee_leadegdetaseed->Fill(-5);
	recoee_leadegdphi->Fill(-5);
	recoee_leadegmhits->Fill(-15);
	recoee_leadegnlayerit->Fill(-15);
	recoee_leadegooeseedoop->Fill(-50);
	recoee_leadegooesclsoop->Fill(-50);
	//recoee_leadegvalhits->Fill(-15);
      }
    } // End of filling end-cap variables
    
  } // End of condition requiring atleast one eg object
  
  if(egidx.size()>=2) { // Require atleast two eg objects
    TLorentzVector leadeg, subleadeg;
    leadeg.SetPtEtaPhiM(egusRecoPt[egidx[0]],egusRecoEta[egidx[0]],egusRecoPhi[egidx[0]],0.106);
    subleadeg.SetPtEtaPhiM(egusRecoPt[egidx[1]],egusRecoEta[egidx[1]],egusRecoPhi[egidx[1]],0.106);

    subleadegpt->Fill(egusRecoPt[egidx[1]]);
    subleadegeta->Fill(egusRecoEta[egidx[1]]);
    subleadegphi->Fill(egusRecoPhi[egidx[1]]);

    // Fill barrel variables
    if(TMath::Abs(egusRecoEta[egidx[1]])<=1.479) {
      recoeb_subleadegclustershape->Fill(egushltEgammaClusterShape[egidx[1]]);
      recoeb_subleadegin5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[1]]);
      recoeb_subleadegin5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[1]]);
      recoeb_subleadegscenergy->Fill(egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoeb_subleadeghovere->Fill(egushltEgammaHoverE[egidx[1]]);
      recoeb_subleadeghovereoversupcluse->Fill(egushltEgammaHoverE[egidx[1]]/egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoeb_subleadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[egidx[1]]);
      recoeb_subleadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[egidx[1]]/egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoeb_subleadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[egidx[1]]);
      recoeb_subleadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[egidx[1]]/egushltEgammaHoverE[egidx[1]]);
      if(egushltEcalSeedClusterTime[egidx[1]]!=0) recoeb_subleadegseedclustime->Fill(egushltEcalSeedClusterTime[egidx[1]]);
      if(egushltEgammaPixelMatchVars_s2[egidx[1]]<TMath::Power(10,36)) {
	recoeb_subleadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[1]]);
      }
      else {
	recoeb_subleadegpixelmchvar_s2->Fill(-5);
      }
      if(egushltEgammaGsfTrackVars_Deta[egidx[1]]<999999) {
	recoeb_subleadegtrkiso->Fill(egushltEgammaEleGsfTrackIso[egidx[1]]);
	recoeb_subleadegchi2->Fill(egushltEgammaGsfTrackVars_Chi2[egidx[1]]);
	recoeb_subleadegdeta->Fill(egushltEgammaGsfTrackVars_Deta[egidx[1]]);
	recoeb_subleadegdetaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[egidx[1]]);
	recoeb_subleadegdphi->Fill(egushltEgammaGsfTrackVars_Dphi[egidx[1]]);
	recoeb_subleadegmhits->Fill(egushltEgammaGsfTrackVars_MissingHits[egidx[1]]);
	recoeb_subleadegnlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[egidx[1]]);
	recoeb_subleadegooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[1]]);
	recoeb_subleadegooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[1]]);
	recoeb_subleadegvalhits->Fill(egushltEgammaGsfTrackVars_ValidHits[egidx[1]]);
      }
      else {
	recoeb_subleadegtrkiso->Fill(-100);
	recoeb_subleadegchi2->Fill(-20);
	recoeb_subleadegdeta->Fill(-5);
	recoeb_subleadegdetaseed->Fill(-5);
	recoeb_subleadegdphi->Fill(-5);
	recoeb_subleadegmhits->Fill(-15);
	recoeb_subleadegnlayerit->Fill(-15);
	recoeb_subleadegooeseedoop->Fill(-50);
	recoeb_subleadegooesclsoop->Fill(-50);
	recoeb_subleadegvalhits->Fill(-15);
      }
    } // End of filling barrel variables
    else { // Fill end-cap variables
      recoee_subleadegclustershape->Fill(egushltEgammaClusterShape[egidx[1]]);
      recoee_subleadegin5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[egidx[1]]);
      recoee_subleadegin5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[egidx[1]]);
      recoee_subleadegscenergy->Fill(egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoee_subleadeghovere->Fill(egushltEgammaHoverE[egidx[1]]);
      recoee_subleadeghovereoversupcluse->Fill(egushltEgammaHoverE[egidx[1]]/egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoee_subleadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[egidx[1]]);
      recoee_subleadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[egidx[1]]/egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoee_subleadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[egidx[1]]);	
      recoee_subleadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[egidx[1]]/egushltEgammaHoverE[egidx[1]]);	
      if(egushltEcalSeedClusterTime[egidx[1]]!=0) recoee_subleadegseedclustime->Fill(egushltEcalSeedClusterTime[egidx[1]]);
      if(egushltEgammaPixelMatchVars_s2[egidx[1]]<TMath::Power(10,36)) {
	recoee_subleadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[1]]);
      }
      else {
	recoee_subleadegpixelmchvar_s2->Fill(-5);
      }
      if(egushltEgammaGsfTrackVars_Deta[egidx[1]]<999999) {
	recoee_subleadegtrkiso->Fill(egushltEgammaEleGsfTrackIso[egidx[1]]);
	recoee_subleadegchi2->Fill(egushltEgammaGsfTrackVars_Chi2[egidx[1]]);
	recoee_subleadegdeta->Fill(egushltEgammaGsfTrackVars_Deta[egidx[1]]);
	recoee_subleadegdetaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[egidx[1]]);
	recoee_subleadegdphi->Fill(egushltEgammaGsfTrackVars_Dphi[egidx[1]]);
	recoee_subleadegmhits->Fill(egushltEgammaGsfTrackVars_MissingHits[egidx[1]]);
	recoee_subleadegnlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[egidx[1]]);
	recoee_subleadegooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[egidx[1]]);
	recoee_subleadegooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[egidx[1]]);
	recoee_subleadegvalhits->Fill(egushltEgammaGsfTrackVars_ValidHits[egidx[1]]);
      }
      else {
	recoee_subleadegtrkiso->Fill(-100);
	recoee_subleadegchi2->Fill(-20);
	recoee_subleadegdeta->Fill(-5);
	recoee_subleadegdetaseed->Fill(-5);
	recoee_subleadegdphi->Fill(-5);
	recoee_subleadegmhits->Fill(-15);
	recoee_subleadegnlayerit->Fill(-15);
	recoee_subleadegooeseedoop->Fill(-50);
	recoee_subleadegooesclsoop->Fill(-50);
	recoee_subleadegvalhits->Fill(-15);
      }
    } // End of filling end-cap variables

    // Fill invariant mass
    if(egusRecoEta[egidx[0]]<1.479 && egusRecoEta[egidx[1]]<1.479) {
      recoeb_leadsubleadM->Fill((leadeg+subleadeg).M());
    }
    else {
      recoee_leadsubleadM->Fill((leadeg+subleadeg).M());
    }
  } // End of condition requiring atleast two eg object
  
    /* 
       Fix the gen matching code before uncommenting this part
    if(isMC) {
      int genmchegpos = doGenMatchingUnseeded(egidx);
      if(genmchegpos!=-1) {
	genmchegpt->Fill(egusRecoPt[genmchegpos]);
	genmchegeta->Fill(egusRecoEta[genmchegpos]);
	genmchegphi->Fill(egusRecoPhi[genmchegpos]);
      }
      }*/
}

// Function to add a set of histograms for a gen particles
void data_robustanalyzer::addgenhist(TString selection) {
  
  all1dhists.push_back(new TH1F(selection+"geneg_egmult","gen N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"geneg_pt","gen e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"geneg_eta","gen e/#gamma #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"geneg_phi","gen e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"geneg_log10d0","gen e/#gamma log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_t1","gen e/#gamma time ecal",10000,-100,900));
  all1dhists.push_back(new TH1F(selection+"geneg_t0","gen e/#gamma time ecal (equivalent prompt)",10000,-100,900));
  all1dhists.push_back(new TH1F(selection+"geneg_t1mt0","gen e/#gamma time ecal (more than prompt)",10000,-10,90));
  all1dhists.push_back(new TH1F(selection+"geneg_leadpt","gen e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"geneg_leadeta","gen e/#gamma #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"geneg_leadphi","gen e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"geneg_leadlog10d0","gen e/#gamma log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadpt","gen e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadeta","gen e/#gamma #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadphi","gen e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadlog10d0","gen e/#gamma log_{10}d_{0} / log_{10}cm",1000,-5,5));
}

// Function to add a set of histograms for a selection
void data_robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"recoeg_egmult","reco N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"recoeg_leadegpt","reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recoeg_leadegeta","reco e/#gamma #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"recoeg_leadegphi","reco e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"recoeg_subleadegpt","reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recoeg_subleadegeta","reco e/#gamma #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"recoeg_subleadegphi","reco e/#gamma #phi",66,-3.3,3.3));

  // barrel variables - lead pT e/gamma
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegclustershape","barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegin5x5clusshape","barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegin5x5noiseclnd","barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadeghovere","barrel e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegscenergy","barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegecalpfclustiso","barrel e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegecalpfclustisoovere","barrel e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadeghcalpfclustiso","barrel e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadeghcalpfclustisoovere","barrel e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegpixelmchvar_s2","barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegchi2","barrel e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegdeta","barrel e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegdetaseed","barrel e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegdphi","barrel e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegmhits","barrel e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegnlayerit","barrel e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegooeseedoop","barrel e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegooesclsoop","barrel e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegvalhits","barrel e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegtrkiso","barrel e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadeghovereoversupcluse","barrel e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegseedclustime","barrel e/#gamma_{seed} time / ns",20000,-10,10));

  // barrel variables - sub-lead pT e/gamma
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegclustershape","barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegin5x5clusshape","barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegin5x5noiseclnd","barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadeghovere","barrel e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegscenergy","barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegecalpfclustiso","barrel e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegecalpfclustisoovere","barrel e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadeghcalpfclustiso","barrel e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadeghcalpfclustisoovere","barrel e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegpixelmchvar_s2","barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegchi2","barrel e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegdeta","barrel e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegdetaseed","barrel e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegdphi","barrel e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegmhits","barrel e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegnlayerit","barrel e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegooeseedoop","barrel e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegooesclsoop","barrel e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegvalhits","barrel e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegtrkiso","barrel e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadeghovereoversupcluse","barrel e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegseedclustime","barrel e/#gamma_{seed} time / ns",20000,-10,10));

  // barrel variables - invariant mass
  all1dhists.push_back(new TH1F(selection+"recoeb_leadsubleadM","M(e/#gamma_{1},e/#gamma_{2}) / GeV",500,0,500));

  // end-cap variables - lead pT e/gamma
  all1dhists.push_back(new TH1F(selection+"recoee_leadegclustershape","end-cap e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegin5x5clusshape","end-cap e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegin5x5noiseclnd","end-cap e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_leadeghovere","end-cap e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegscenergy","end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegecalpfclustiso","end-cap e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegecalpfclustisoovere","end-cap e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoee_leadeghcalpfclustiso","end-cap e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_leadeghcalpfclustisoovere","end-cap e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegpixelmchvar_s2","end-cap e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegchi2","end-cap e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegdeta","end-cap e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegdetaseed","end-cap e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegdphi","end-cap e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegmhits","end-cap e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegnlayerit","end-cap e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegooeseedoop","end-cap e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegooesclsoop","end-cap e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegvalhits","end-cap e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegtrkiso","end-cap e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_leadeghovereoversupcluse","end-cap e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegseedclustime","end-cap e/#gamma_{seed} time / ns",20000,-10,10));

  // end-cap variables - sub-lead pT e/gamma
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegclustershape","end-cap e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegin5x5clusshape","end-cap e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegin5x5noiseclnd","end-cap e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadeghovere","end-cap e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegscenergy","end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegecalpfclustiso","end-cap e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegecalpfclustisoovere","end-cap e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadeghcalpfclustiso","end-cap e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadeghcalpfclustisoovere","end-cap e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegpixelmchvar_s2","end-cap e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegchi2","end-cap e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegdeta","end-cap e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegdetaseed","end-cap e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegdphi","end-cap e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegmhits","end-cap e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegnlayerit","end-cap e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegooeseedoop","end-cap e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegooesclsoop","end-cap e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegvalhits","end-cap e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegtrkiso","end-cap e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadeghovereoversupcluse","end-cap e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegseedclustime","end-cap e/#gamma_{seed} time / ns",20000,-10,10));

  // end-cap variables - invariant mass
  all1dhists.push_back(new TH1F(selection+"recoee_leadsubleadM","M(e/#gamma_{1},e/#gamma_{2}) / GeV",500,0,500));
}

// Function to add a set of histograms for a selection - unseeded egamma objects
void data_robustanalyzer::addhistunseeded(TString selection) {

  all1dhists.push_back(new TH1F(selection+"recoegus_egmult","reco N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"recoegus_leadegpt","reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recoegus_leadegeta","reco e/#gamma #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"recoegus_leadegphi","reco e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"recoegus_subleadegpt","reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recoegus_subleadegeta","reco e/#gamma #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"recoegus_subleadegphi","reco e/#gamma #phi",66,-3.3,3.3));

  if(isMC) {
    all1dhists.push_back(new TH1F(selection+"recoegus_leadgenmchegpt","gen matched reco e/#gamma p_{T} / GeV",550,-50,500));
    all1dhists.push_back(new TH1F(selection+"recoegus_leadgenmchegeta","gen matched reco e/#gamma #eta",52,-2.6,2.6));
    all1dhists.push_back(new TH1F(selection+"recoegus_leadgenmchegphi","gen matched reco e/#gamma #phi",66,-3.3,3.3));
    all1dhists.push_back(new TH1F(selection+"recoegus_subleadgenmchegpt","gen matched reco e/#gamma p_{T} / GeV",550,-50,500));
    all1dhists.push_back(new TH1F(selection+"recoegus_subleadgenmchegeta","gen matched reco e/#gamma #eta",52,-2.6,2.6));
    all1dhists.push_back(new TH1F(selection+"recoegus_subleadgenmchegphi","gen matched reco e/#gamma #phi",66,-3.3,3.3));
  }

  // barrel variables - lead pT unseeded e/gamma
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegclustershape","barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegin5x5clusshape","barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegin5x5noiseclnd","barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadeghovere","barrel e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegscenergy","barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegecalpfclustiso","barrel e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegecalpfclustisoovere","barrel e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadeghcalpfclustiso","barrel e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadeghcalpfclustisoovere","barrel e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegpixelmchvar_s2","barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegchi2","barrel e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegdeta","barrel e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegdetaseed","barrel e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegdphi","barrel e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegmhits","barrel e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegnlayerit","barrel e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegooeseedoop","barrel e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegooesclsoop","barrel e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegvalhits","barrel e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegtrkiso","barrel e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadeghovereoversupcluse","barrel e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegseedclustime","barrel e/#gamma_{seed} time / ns",20000,-10,10));

  // barrel variables - sub-lead pT unseeded e/gamma
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegclustershape","barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegin5x5clusshape","barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegin5x5noiseclnd","barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadeghovere","barrel e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegscenergy","barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegecalpfclustiso","barrel e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegecalpfclustisoovere","barrel e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadeghcalpfclustiso","barrel e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadeghcalpfclustisoovere","barrel e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegpixelmchvar_s2","barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegchi2","barrel e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegdeta","barrel e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegdetaseed","barrel e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegdphi","barrel e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegmhits","barrel e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegnlayerit","barrel e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegooeseedoop","barrel e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegooesclsoop","barrel e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegvalhits","barrel e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegtrkiso","barrel e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadeghovereoversupcluse","barrel e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegseedclustime","barrel e/#gamma_{seed} time / ns",20000,-10,10));

  // barrel variables - invariant mass unseeded
  all1dhists.push_back(new TH1F(selection+"recoebus_leadsubleadM","barrel M(e/#gamma_{1},e/#gamma_{2}) / GeV",500,0,500));

  // end-cap variables - lead pT unseeded e/gamma
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegclustershape","end-cap e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegin5x5clusshape","end-cap e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegin5x5noiseclnd","end-cap e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadeghovere","end-cap e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegscenergy","end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegecalpfclustiso","end-cap e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegecalpfclustisoovere","end-cap e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadeghcalpfclustiso","end-cap e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadeghcalpfclustisoovere","end-cap e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegpixelmchvar_s2","end-cap e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegchi2","end-cap e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegdeta","end-cap e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegdetaseed","end-cap e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegdphi","end-cap e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegmhits","end-cap e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegnlayerit","end-cap e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegooeseedoop","end-cap e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegooesclsoop","end-cap e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegvalhits","end-cap e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegtrkiso","end-cap e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadeghovereoversupcluse","end-cap e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegseedclustime","end-cap e/#gamma_{seed} time / ns",20000,-10,10));

  // end-cap variables - sub-lead pT unseeded e/gamma
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegclustershape","end-cap e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegin5x5clusshape","end-cap e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegin5x5noiseclnd","end-cap e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadeghovere","end-cap e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegscenergy","end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegecalpfclustiso","end-cap e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegecalpfclustisoovere","end-cap e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadeghcalpfclustiso","end-cap e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadeghcalpfclustisoovere","end-cap e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegpixelmchvar_s2","end-cap e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegchi2","end-cap e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegdeta","end-cap e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegdetaseed","end-cap e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegdphi","end-cap e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegmhits","end-cap e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegnlayerit","end-cap e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegooeseedoop","end-cap e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegooesclsoop","end-cap e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegvalhits","end-cap e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegtrkiso","end-cap e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadeghovereoversupcluse","end-cap e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegseedclustime","end-cap e/#gamma_{seed} time / ns",20000,-10,10));

  // end-cap variables - invariant mass unseeded
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadsubleadM","end-cap M(e/#gamma_{1},e/#gamma_{2}) / GeV",500,0,500));

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
