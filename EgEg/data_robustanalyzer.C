/*
 * AUTHOR: Abanti Ranadhir Sahasransu - asahasra@cern.ch
 * The code now assumes exactly two gen electrons for a MC sample
 */

#include "data_robustanalyzer.hh"
#include <iostream>
#include <numeric>
#include <stdexcept>

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
data_robustanalyzer::data_robustanalyzer(TString filename, TString outfilename, bool issimu) {

  isMC = issimu;
  
  inputChain = new TChain("events");
  cout<<"Initializing for file: "<<filename<<endl;
  inputChain->Add(filename);

  inputChain->SetBranchAddress("l1EgObjn", &l1egObjN);
  inputChain->SetBranchAddress("l1EgObj_pt", &l1egObjPt);
  inputChain->SetBranchAddress("l1EgObj_eta", &l1egObjEta);
  inputChain->SetBranchAddress("l1EgObj_phi", &l1egObjPhi);
  inputChain->SetBranchAddress("HLT_DoublePhoton33CaloIdL", &HLT_DoublePhoton33_CaloIdL);
  inputChain->SetBranchAddress("HLT_DoublePhoton70", &HLT_DoublePhoton70);
  inputChain->SetBranchAddress("l1Filtn", &l1FiltN);
  inputChain->SetBranchAddress("l1Filt_pt", &l1FiltPt);
  inputChain->SetBranchAddress("l1Filt_eta", &l1FiltEta);
  inputChain->SetBranchAddress("l1Filt_phi", &l1FiltPhi);
  inputChain->SetBranchAddress("bsx", &bsx);
  inputChain->SetBranchAddress("bsy", &bsy);
  inputChain->SetBranchAddress("bsz", &bsz);
  inputChain->SetBranchAddress("dieg70_egheusFiltn", &dieg70HeusFiltN);
  inputChain->SetBranchAddress("dieg70_egheusFilt_pt", &dieg70HeusFiltPt);
  inputChain->SetBranchAddress("dieg70_egheusFilt_eta", &dieg70HeusFiltEta);
  inputChain->SetBranchAddress("dieg70_egheusFilt_phi", &dieg70HeusFiltPhi);
  inputChain->SetBranchAddress("dieg33_egcsFiltn", &dieg33CsFiltN);
  inputChain->SetBranchAddress("dieg33_egcsFilt_pt", &dieg33CsFiltPt);
  inputChain->SetBranchAddress("dieg33_egcsFilt_eta", &dieg33CsFiltEta);
  inputChain->SetBranchAddress("dieg33_egcsFilt_phi", &dieg33CsFiltPhi);
  inputChain->SetBranchAddress("dieg33_egcsusFiltn", &dieg33CsusFiltN);
  inputChain->SetBranchAddress("dieg33_egcsusFilt_pt", &dieg33CsusFiltPt);
  inputChain->SetBranchAddress("dieg33_egcsusFilt_eta", &dieg33CsusFiltEta);
  inputChain->SetBranchAddress("dieg33_egcsusFilt_phi", &dieg33CsusFiltPhi);
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
  inputChain->SetBranchAddress("egushltEgammaClusterShape_sminarr", &egushltEgammaClusterShape_smin);
  inputChain->SetBranchAddress("egushltEgammaClusterShape_smajarr", &egushltEgammaClusterShape_smaj);
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

  if(isMC) {
    inputChain->SetBranchAddress("genLepn",&genLepN);
    inputChain->SetBranchAddress("genLepPIDarr",&genLepPid);
    inputChain->SetBranchAddress("genLepPtarr",&genLepPt);
    inputChain->SetBranchAddress("genLepEtaarr",&genLepEta);
    inputChain->SetBranchAddress("genLepPhiarr",&genLepPhi);
    inputChain->SetBranchAddress("genLepPromptEtaarr",&genLepPromptEta);
    inputChain->SetBranchAddress("genLepPromptPhiarr",&genLepPromptPhi);
    inputChain->SetBranchAddress("genLepVxarr",&genLepVx);
    inputChain->SetBranchAddress("genLepVyarr",&genLepVy);
    inputChain->SetBranchAddress("genLepVzarr",&genLepVz);
    inputChain->SetBranchAddress("genLepNmomarr",&genLepNMom);
    inputChain->SetBranchAddress("genLepMomPIDarr",&genLepMomPid);
    inputChain->SetBranchAddress("genLepMomPtarr",&genLepMomPt);
    inputChain->SetBranchAddress("genLepMomEtaarr",&genLepMomEta);
    inputChain->SetBranchAddress("genLepMomPhiarr",&genLepMomPhi);
    inputChain->SetBranchAddress("pvx", &pvx);
    inputChain->SetBranchAddress("pvy", &pvy);
    inputChain->SetBranchAddress("pvz", &pvz);
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
  int nosel=0, noselus=0, basicsel=0, basicselus=0, selelevetoid=0, selelevetoidus=0, selelevetozwindidus=0, selelevetozoppoidus=0, seleletightid=0, seleletightidus=0, dieg70id=0, dieg70idus=0, dieg33caloidl=0, cut1us=0, cut2us=0, cut3us=0, cuttimedelayonlyus=0, cuttimedelaysminus=0;
  // Trigger cross-check
  bool dipho70trig = false;
  bool dipho33caloidltrig = false;
  
  // Define the histograms
  if(isMC) addgenhist("gennosel");
  if(isMC) addgenhist("genbasicselbar");
  if(isMC) addgenhist("genbasicptgt10selbar");
  if(isMC) addgenhist("genbasicptgt10selec");
  if(isMC) addgenhist("genetabin14_16_24");
  if(isMC) addgenhist("genptgt10");
  if(isMC) addgenhist("genptgt10etalt12");
  if(isMC) addgenhist("genptgt10etabin16_24");
  if(isMC) addgenhist("genptgt10etabin14_16_24");
  if(isMC) addgenhist("genptgt10etabin14_16_24d0lt1cm");
  if(isMC) addgenhist("genbarsel");
  if(isMC) addgenhist("genbarselptgt10");
  if(isMC) addgenhist("genbasicselptgt15");
  addhist("nosel");
  addhistunseeded("noselus");
  if(isMC)addhistgenmchunseeded("gennoselAnoselus");
  if(isMC)addhistgenmchunseeded("genbasicselbarAnoselus");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selbarAnoselus");
  if(isMC)addhistgenmchunseeded("genetabin14_16_24Anoselus");
  if(isMC)addhistgenmchunseeded("genptgt10Anoselus");
  if(isMC)addhistgenmchunseeded("genptgt10etalt12Anoselus");
  if(isMC)addhistgenmchunseeded("genptgt10etabin16_24Anoselus");
  if(isMC)addhistgenmchunseeded("genptgt10etabin14_16_24Anoselus");
  if(isMC)addhistgenmchunseeded("genptgt10etabin14_16_24d0lt1cmAnoselus");
  addhist("basicsel");
  addhistunseeded("basicselus");
  if(isMC)addhistgenmchunseeded("genptgt10Abasicselus");
  if(isMC)addhistgenmchunseeded("genbasicselptgt15Abasicselus");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selbarAbasicselus");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selecAbasicselus");
  addhist("selelevetoid");
  addhistunseeded("selelevetoidus");
  addhistunseeded("selelevetozwindidus");
  addhistunseeded("selelevetozoppoidus");
  addhist("seleletightid");
  addhistunseeded("seleletightidus");
  addhistunseeded("seleletightcaloidus");
  addhist("dieg70id");
  addhistunseeded("dieg70idus");
  addhist("dieg33caloidlid");
  addhistunseeded("dieg33caloidlidus");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selbarAdieg33caloidlus");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selecAdieg33caloidlus");
  addhistunseeded("cut1us");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selbarAcut1us");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selecAcut1us");
  addhistunseeded("cut2us");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selbarAcut2us");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selecAcut2us");
  addhistunseeded("cut3us");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selbarAcut3us");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selecAcut3us");
  addhistunseeded("cuttimedelayonlyus");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selbarAcuttimedelayonlyus");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selecAcuttimedelayonlyus");
  addhistunseeded("cuttimedelaysminus");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selbarAcuttimedelaysminus");
  if(isMC)addhistgenmchunseeded("genbasicptgt10selecAcuttimedelaysminus");
  
  // Loop beginning on events
  for(unsigned int event=beginevent; event<endevent; event++) {

    // vector of eg indices
    vector<int> genelpos;
    vector<int> gennoselegidx;
    vector<int> genbasicselbaregidx;
    vector<int> genbasicptgt10selbaregidx;
    vector<int> genbasicptgt10selecegidx;
    vector<int> genetabin14_16_24egidx;
    vector<int> genptgt10egidx;
    vector<int> genptgt10etalt12egidx;
    vector<int> genptgt10etabin16_24egidx;
    vector<int> genptgt10etabin14_16_24egidx;
    vector<int> genptgt10etabin14_16_24d0lt1cmegidx;
    vector<int> genbarselegidx;
    vector<int> genbarselptgt10egidx;
    vector<int> genbasicselptgt15egidx;
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
    vector<int> seleletightcaloidegusidx;
    vector<int> dieg70idegidx;
    vector<int> dieg70idegusidx;
    vector<int> dieg33caloidlidegidx;
    vector<int> dieg33caloidlidegusidx;
    vector<int> cut1usidx;
    vector<int> cut2usidx;
    vector<int> cut3usidx;
    vector<int> cuttimedelayonlyusidx;
    vector<int> cuttimedelaysminusidx;
  
    inputChain->GetEntry(event);
    //if(event>10000) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    // Block for gen leptons
    if(isMC) {

      genelpos.push_back(-1);
      genelpos.push_back(-1);

      // Select and assign the gen electrons with pt ordering
      int numgen = 0;
      for(unsigned int genCtr=0; genCtr<genLepN; genCtr++) {

	// Calculate the impact parameter and transverse displacement
	// Use the mother of the lepton as primary vertex
	TLorentzVector part;
	part.SetPtEtaPhiM(genLepPt[genCtr],genLepEta[genCtr],genLepPhi[genCtr],0.0005);
	double d0 = (genLepVx[genCtr]-pvx[genCtr])*part.Py()-(genLepVy[genCtr]-pvy[genCtr])*part.Px();
	double lxy = (genLepVx[genCtr]-pvx[genCtr])*(genLepVx[genCtr]-pvx[genCtr]);
	lxy += (genLepVy[genCtr]-pvy[genCtr])*(genLepVy[genCtr]-pvy[genCtr]);
	lxy = TMath::Sqrt(lxy);
	d0 /= genLepPt[genCtr];
	genLepDxy[genCtr] = d0;
	genLepLxy[genCtr] = lxy;

	// Select and assign the gen electrons with pt ordering
	bool genselect = true;
	genselect *= abs(genLepPid[genCtr])==11;
	genselect *= abs(genLepMomPid[genCtr])==9000007 || abs(genLepMomPid[genCtr])==23;
	if(!genselect) continue;
	numgen++;
	if(numgen>2) throw "Error!!! Code configured for two gen electrons only";
	if(genelpos[0] == -1) genelpos[0] = genCtr;
	else {
	  if(genLepPt[genelpos[0]]<genLepPt[genCtr]) {
	    genelpos[1] = genelpos[0];
	    genelpos[0] = genCtr;
	  }
	  else {
	    genelpos[1] = genCtr;
	  }
	}
      } // End of gen loop
      if(numgen<2) continue; // Skip if event final state is not e-e
      if(numgen!=2) throw "Error!!! Code configured for exactly two gen electrons only";
      
      bool gennoseleg = false;
      bool genbasicselbareg = false;
      bool genbasicptgt10selbareg = false;
      bool genbasicptgt10seleceg = false;
      bool genetabin14_16_24eg = false;
      bool genptgt10eg = false;
      bool genptgt10etalt12eg = false;
      bool genptgt10etabin16_24eg = false;
      bool genptgt10etabin14_16_24eg = false;
      bool genptgt10etabin14_16_24d0lt1cmeg = false;
      bool genbarseleg = false;
      bool genbarselptgt10eg = false;
      bool genbasicselptgt15eg = false;

      for(int genCtr:genelpos) {

	gennoselegidx.push_back(genCtr);

        genbasicselbareg = true;
	genbasicselbareg = abs(genLepVz[genCtr])<320 || TMath::Sqrt(genLepVx[genCtr]*genLepVx[genCtr]+genLepVy[genCtr]*genLepVy[genCtr])<130;
        genbasicselbareg *= abs(genLepPromptEta[genCtr])<1.479;
	if(genbasicselbareg) genbasicselbaregidx.push_back(genCtr);
	else genbasicselbaregidx.push_back(-1);

        genbasicptgt10selbareg = true;
	genbasicptgt10selbareg = abs(genLepVz[genCtr])<320 || TMath::Sqrt(genLepVx[genCtr]*genLepVx[genCtr]+genLepVy[genCtr]*genLepVy[genCtr])<130;
        genbasicptgt10selbareg *= abs(genLepPromptEta[genCtr])<1.479;
        genbasicptgt10selbareg *= genLepPt[genCtr]>10;
	if(genbasicptgt10selbareg) genbasicptgt10selbaregidx.push_back(genCtr);
	else genbasicptgt10selbaregidx.push_back(-1);

        genbasicptgt10seleceg = true;
	genbasicptgt10seleceg = abs(genLepVz[genCtr])<320 || TMath::Sqrt(genLepVx[genCtr]*genLepVx[genCtr]+genLepVy[genCtr]*genLepVy[genCtr])<130;
        genbasicptgt10seleceg *= abs(genLepPromptEta[genCtr])>1.479;
        genbasicptgt10seleceg *= genLepPt[genCtr]>10;
	if(genbasicptgt10seleceg) genbasicptgt10selecegidx.push_back(genCtr);
	else genbasicptgt10selecegidx.push_back(-1);

	genetabin14_16_24eg = true;
	genetabin14_16_24eg *= abs(genLepEta[genCtr])<1.4 || (abs(genLepEta[genCtr])>1.6 && abs(genLepEta[genCtr])<2.4);
	if(genetabin14_16_24eg) genetabin14_16_24egidx.push_back(genCtr);
	else genetabin14_16_24egidx.push_back(-1);

	genptgt10eg = true;
	genptgt10eg *= genLepPt[genCtr]>10;
	if(genptgt10eg) genptgt10egidx.push_back(genCtr);
	else genptgt10egidx.push_back(-1);

	genptgt10etalt12eg = true;
	genptgt10etalt12eg *= genLepPt[genCtr]>10;
	genptgt10etalt12eg *= abs(genLepEta[genCtr])<1.2;
	if(genptgt10etalt12eg) genptgt10etalt12egidx.push_back(genCtr);
	else genptgt10etalt12egidx.push_back(-1);

	genptgt10etabin16_24eg = true;
	genptgt10etabin16_24eg *= genLepPt[genCtr]>10;
	genptgt10etabin16_24eg *= abs(genLepEta[genCtr])>1.6 && abs(genLepEta[genCtr])<2.4;
	if(genptgt10etabin16_24eg) genptgt10etabin16_24egidx.push_back(genCtr);
	else genptgt10etabin16_24egidx.push_back(-1);

	genptgt10etabin14_16_24eg = true;
	genptgt10etabin14_16_24eg *= genLepPt[genCtr]>10;
	genptgt10etabin14_16_24eg *= abs(genLepEta[genCtr])<1.4 || (abs(genLepEta[genCtr])>1.6 && abs(genLepEta[genCtr])<2.4);
	if(genptgt10etabin14_16_24eg) genptgt10etabin14_16_24egidx.push_back(genCtr);
	else genptgt10etabin14_16_24egidx.push_back(-1);

	genptgt10etabin14_16_24d0lt1cmeg = true;
	genptgt10etabin14_16_24d0lt1cmeg *= genLepPt[genCtr]>10;
	genptgt10etabin14_16_24d0lt1cmeg *= abs(genLepEta[genCtr])<1.4 || (abs(genLepEta[genCtr])>1.6 && abs(genLepEta[genCtr])<2.4);
	if(genptgt10etabin14_16_24d0lt1cmeg) genptgt10etabin14_16_24d0lt1cmegidx.push_back(genCtr);
	else genptgt10etabin14_16_24d0lt1cmegidx.push_back(-1);

      	genbarseleg = true;
	genbarseleg *= abs(genLepEta[genCtr])<1.479;
	if(genbarseleg) genbarselegidx.push_back(genCtr);
	else genbarselegidx.push_back(-1);

	genbarselptgt10eg = true;
	genbarselptgt10eg *= abs(genLepEta[genCtr])<1.479;
	genbarselptgt10eg *= genLepPt[genCtr]>10;
	if(genbarselptgt10eg) genbarselptgt10egidx.push_back(genCtr);
	else genbarselptgt10egidx.push_back(-1);

	genbasicselptgt15eg = true;
	genbasicselptgt15eg *= abs(genLepEta[genCtr])<2.5;
	genbasicselptgt15eg *= genLepPt[genCtr]>15;
	if(genbasicselptgt15eg) genbasicselptgt15egidx.push_back(genCtr);
        else genbasicselptgt15egidx.push_back(-1);

      }

      fillgenhistinevent("gennosel",gennoselegidx); // Verified that this is always 2 electrons
      fillgenhistinevent("genbasicselbar",genbasicselbaregidx); 
      fillgenhistinevent("genbasicptgt10selbar",genbasicptgt10selbaregidx); 
      fillgenhistinevent("genbasicptgt10selec",genbasicptgt10selecegidx); 
      fillgenhistinevent("genetabin14_16_24",genetabin14_16_24egidx);
      fillgenhistinevent("genptgt10",genptgt10egidx);
      fillgenhistinevent("genptgt10etalt12",genptgt10etalt12egidx);
      fillgenhistinevent("genptgt10etabin16_24",genptgt10etabin16_24egidx);
      fillgenhistinevent("genptgt10etabin14_16_24",genptgt10etabin14_16_24egidx);
      fillgenhistinevent("genptgt10etabin14_16_24d0lt1cm",genptgt10etabin14_16_24d0lt1cmegidx);
      fillgenhistinevent("genbarsel",genbarselegidx);
      fillgenhistinevent("genbarselptgt10",genbarselptgt10egidx);
      fillgenhistinevent("genbasicselptgt15",genbasicselptgt15egidx);
    } // End of loop on gen electrons
    
    if((HLT_DoublePhoton70==1 && dieg70HeusFiltN<2) || (HLT_DoublePhoton70==0 && dieg70HeusFiltN>=2)) throw "Error!! Inconsistent trigger result with number of objects for HLT_DoublePhoton70.";
    if((HLT_DoublePhoton33_CaloIdL==1 && dieg33CsusFiltN<2) || (HLT_DoublePhoton33_CaloIdL==0 && dieg33CsusFiltN>=2)) throw "Error!! Inconsistent trigger result with number of objects for HLT_DoublePhoton33_CaloIdL.";

    if(egusRecoN<0 || egRecoN<0) throw "Error!! Negative number of objects pas possible.";
      
    if(egusRecoN>=1) { // Atleast one reco eg us object in the event

      //if(egRecoN>egusRecoN) cout<<"Error!!! Cannot have more seeded objects than unseeded objects in an event."<<endl;
      
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
      bool seleletightcaloidegus = false;
      bool dieg70ideg = false;
      bool dieg70idegus = false;
      bool dieg33caloidlideg = false;
      bool dieg33caloidlidegus = false;
      bool cut1egus = false;
      bool cut2egus = false;
      bool cut3egus = false;
      bool cuttimedelayonlyegus = false;
      bool cuttimedelaysminegus = false;
    
      // Loop beginning on egamma reco objects
      for(unsigned int egidx=0; egidx<egRecoN; egidx++) {

	unsigned int idx = sortedegidx[egidx];
	noselegidx.push_back(idx);

	basicseleg = true;
	basicseleg *= (TMath::Abs(egRecoEta[idx])<2.5);
	basicseleg *= (egRecoPt[idx]>=15);
	if(basicseleg) basicselegidx.push_back(idx);

	selelevetoideg = true;
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<2.5);
	selelevetoideg *= (egRecoPt[idx]>=15);
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0126:eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0457);
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_DetaSeed[idx])<0.00463:abs(eghltEgammaGsfTrackVars_DetaSeed[idx])<0.00814);
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_Dphi[idx])<0.148:abs(eghltEgammaGsfTrackVars_Dphi[idx])<0.19);
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.05*eghltEgammaSuperClusterEnergy[idx]+1.16:eghltEgammaHoverE[idx]<0.05*eghltEgammaSuperClusterEnergy[idx]+2.54);
	selelevetoideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.209:abs(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.132);
	if(selelevetoideg) selelevetoidegidx.push_back(idx);

	seleletightideg = true;
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<2.5);
	seleletightideg *= (egRecoPt[idx]>=15);
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0104:eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0353);
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_DetaSeed[idx])<0.00255:abs(eghltEgammaGsfTrackVars_DetaSeed[idx])<0.00501);
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_Dphi[idx])<0.022:abs(eghltEgammaGsfTrackVars_Dphi[idx])<0.0236);
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.026*eghltEgammaSuperClusterEnergy[idx]+1.15:eghltEgammaHoverE[idx]<0.0188*eghltEgammaSuperClusterEnergy[idx]+2.06);
	seleletightideg *= (TMath::Abs(egRecoEta[idx])<1.479?abs(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.159:abs(eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.0197);
	if(seleletightideg) seleletightidegidx.push_back(idx);

	dieg70ideg = true;
	dieg70ideg *= (TMath::Abs(egRecoEta[idx])<2.65);
	dieg70ideg *= (egRecoPt[idx]>=70);
	dieg70ideg *= isL1EgSeeded(idx);
	dieg70ideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.15*eghltEgammaSuperClusterEnergy[idx]:eghltEgammaHoverE[idx]<0.1*eghltEgammaSuperClusterEnergy[idx]);
	if(dieg70ideg) dieg70idegidx.push_back(idx);

	dieg33caloidlideg = true;
	dieg33caloidlideg *= (TMath::Abs(egRecoEta[idx])<2.65);
	dieg33caloidlideg *= (egRecoPt[idx]>=33);
	dieg33caloidlideg *= isL1EgSeeded(idx);
	dieg33caloidlideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaHoverE[idx]<0.15*eghltEgammaSuperClusterEnergy[idx]:eghltEgammaHoverE[idx]<0.1*eghltEgammaSuperClusterEnergy[idx]);
	dieg33caloidlideg *= (TMath::Abs(egRecoEta[idx])<1.479?eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.014:eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.035);
	if(dieg33caloidlideg) dieg33caloidlidegidx.push_back(idx);

      } // End of loop on egamma reco objects
            
      // Loop beginning on unseeded egamma reco objects
      for(unsigned int egidx=0; egidx<egusRecoN; egidx++) {
	
	unsigned int idx = sortedegusidx[egidx];
	noselegusidx.push_back(idx);

	basicselegus = true;
	basicselegus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	basicselegus *= (egusRecoPt[idx]>=10);
	if(basicselegus) basicselegusidx.push_back(idx);

      	selelevetoidegus = true;
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	selelevetoidegus *= (egusRecoPt[idx]>=15);
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0126:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0457);
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_DetaSeed[idx])<0.00463:abs(egushltEgammaGsfTrackVars_DetaSeed[idx])<0.00814);
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_Dphi[idx])<0.148:abs(egushltEgammaGsfTrackVars_Dphi[idx])<0.19);
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.05*egushltEgammaSuperClusterEnergy[idx]+1.16:egushltEgammaHoverE[idx]<0.05*egushltEgammaSuperClusterEnergy[idx]+2.54);
	selelevetoidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.209:abs(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.0132);
	if(selelevetoidegus) selelevetoidegusidx.push_back(idx);

      	seleletightidegus = true;
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	seleletightidegus *= (egusRecoPt[idx]>=15);
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0104:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0353);
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_DetaSeed[idx])<0.00255:abs(egushltEgammaGsfTrackVars_DetaSeed[idx])<0.00501);
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_Dphi[idx])<0.022:abs(egushltEgammaGsfTrackVars_Dphi[idx])<0.0236);
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.026*egushltEgammaSuperClusterEnergy[idx]+1.15:egushltEgammaHoverE[idx]<0.0188*egushltEgammaSuperClusterEnergy[idx]+2.06);
	seleletightidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?abs(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.159:abs(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[idx])<0.0197);
	if(seleletightidegus) seleletightidegusidx.push_back(idx);

      	seleletightcaloidegus = true;
	seleletightcaloidegus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	seleletightcaloidegus *= (egusRecoPt[idx]>=15);
	seleletightcaloidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0104:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.0353);
	seleletightcaloidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.026*egushltEgammaSuperClusterEnergy[idx]+1.15:egushltEgammaHoverE[idx]<0.0188*egushltEgammaSuperClusterEnergy[idx]+2.06);
	if(seleletightcaloidegus) seleletightcaloidegusidx.push_back(idx);

      	dieg70idegus = true;
	//dieg70idegus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	dieg70idegus *= (egusRecoPt[idx]>=70);
	dieg70idegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.15*egushltEgammaSuperClusterEnergy[idx]:egushltEgammaHoverE[idx]<0.1*egushltEgammaSuperClusterEnergy[idx]);
	if(dieg70idegus) dieg70idegusidx.push_back(idx);

      	dieg33caloidlidegus = true;
	//dieg33caloidlidegus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	dieg33caloidlidegus *= (egusRecoPt[idx]>=33);
	dieg33caloidlidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.15*egushltEgammaSuperClusterEnergy[idx]:egushltEgammaHoverE[idx]<0.1*egushltEgammaSuperClusterEnergy[idx]);
	dieg33caloidlidegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.014:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.035);
	if(dieg33caloidlidegus) dieg33caloidlidegusidx.push_back(idx);

      	cut1egus = true;
	cut1egus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	cut1egus *= (egusRecoPt[idx]>=10);
	cut1egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.5*egushltEgammaSuperClusterEnergy[idx]:egushltEgammaHoverE[idx]<0.6*egushltEgammaSuperClusterEnergy[idx]);
	cut1egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.016:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.04);
	//cut1egus *= (egushltEgammaClusterShape_smin[idx]<0.4);
	if(cut1egus) cut1usidx.push_back(idx);

      	cut2egus = true;
	cut2egus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	cut2egus *= (egusRecoPt[idx]>=10);
	cut2egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.5*egushltEgammaSuperClusterEnergy[idx]:egushltEgammaHoverE[idx]<0.6*egushltEgammaSuperClusterEnergy[idx]);
	cut2egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.016:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.04);
	cut2egus *= (egushltEgammaClusterShape_smin[idx]<0.4);
	cut2egus *= (egushltEcalSeedClusterTime[idx]<2);
	if(cut2egus) cut2usidx.push_back(idx);

      	cut3egus = true;
	cut3egus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	cut3egus *= (egusRecoPt[idx]>=10);
	cut3egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.5*egushltEgammaSuperClusterEnergy[idx]:egushltEgammaHoverE[idx]<0.6*egushltEgammaSuperClusterEnergy[idx]);
	cut3egus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.016:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.04);
	cut3egus *= (egushltEgammaClusterShape_smin[idx]<0.16);
	cut3egus *= (egushltEcalSeedClusterTime[idx]<2);
	if(cut3egus) cut3usidx.push_back(idx);

      	cuttimedelayonlyegus = true;
	cuttimedelayonlyegus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	cuttimedelayonlyegus *= (egusRecoPt[idx]>=10);
	cuttimedelayonlyegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.2*egushltEgammaSuperClusterEnergy[idx]:egushltEgammaHoverE[idx]<0.2*egushltEgammaSuperClusterEnergy[idx]);
	cuttimedelayonlyegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.016:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.04);
	cuttimedelayonlyegus *= (egushltEcalSeedClusterTime[idx]>1.4);
	if(cuttimedelayonlyegus) cuttimedelayonlyusidx.push_back(idx);

      	cuttimedelaysminegus = true;
	cuttimedelaysminegus *= (TMath::Abs(egusRecoEta[idx])<2.5);
	cuttimedelaysminegus *= (egusRecoPt[idx]>=10);
	cuttimedelaysminegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaHoverE[idx]<0.2*egushltEgammaSuperClusterEnergy[idx]:egushltEgammaHoverE[idx]<0.2*egushltEgammaSuperClusterEnergy[idx]);
	cuttimedelaysminegus *= (TMath::Abs(egusRecoEta[idx])<1.479?egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.016:egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[idx]<0.04);
	cuttimedelaysminegus *= (egushltEcalSeedClusterTime[idx]<2);
	cuttimedelaysminegus *= (egushltEgammaClusterShape_smin[idx]<0.1);
	if(cuttimedelaysminegus) cuttimedelaysminusidx.push_back(idx);

      } // End of loop on unseeded egamma objects

      selelevetozwindidegusidx = selelevetoidegusidx;
      selelevetozoppoidegusidx = selelevetoidegusidx;

      fillhistinevent("nosel", noselegidx);
      fillhistineventunseeded("noselus", noselegusidx);
      fillhistinevent("basicsel", basicselegidx);
      if(basicselegusidx.size()>=2) fillhistineventunseeded("basicselus", basicselegusidx);
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
      fillhistineventunseeded("seleletightcaloidus", seleletightcaloidegusidx);
      fillhistinevent("dieg70id", dieg70idegidx);
      fillhistineventunseeded("dieg70idus", dieg70idegusidx);
      if(dieg33caloidlidegidx.size()>=1 && dieg33caloidlidegusidx.size()>=2) {
	fillhistinevent("dieg33caloidlid", dieg33caloidlidegidx);
	fillhistineventunseeded("dieg33caloidlidus", dieg33caloidlidegusidx);
      }
      if(cut1usidx.size()>=2) fillhistineventunseeded("cut1us", cut1usidx);
      if(cut2usidx.size()>=2) fillhistineventunseeded("cut2us", cut2usidx);
      if(cut3usidx.size()>=2) fillhistineventunseeded("cut3us", cut3usidx);
      if(cuttimedelayonlyusidx.size()>=2) fillhistineventunseeded("cuttimedelayonlyus", cuttimedelayonlyusidx);
      if(cuttimedelaysminusidx.size()>=2) fillhistineventunseeded("cuttimedelaysminus", cuttimedelaysminusidx);
      
    } // End of condition requiring atleast one egReco object

    // Cross-check with known triggers
    dipho70trig = false;
    if(!comparecutonobjtofilt(dieg70HeusFiltPt, dieg70HeusFiltN, egusRecoPt, dieg70idegusidx)) cout<<event<<"Mismatching objects: HLT_DoublePhoton70"<<endl;;
    if(dieg70idegidx.size()>=1 && dieg70idegusidx.size()>=2) dipho70trig = true;
    if(HLT_DoublePhoton70==true && dipho70trig==false) cout<<event<<": Type 1 - too many cuts on obj - trigger cross-check failed for: HLT_DoublePhoton70"<<endl;
    if(HLT_DoublePhoton70==false && dipho70trig==true) cout<<event<<": Type 2(accep) - not enough cuts on obj - trigger cross-check failed for: HLT_DoublePhoton70"<<endl;

    dipho33caloidltrig = false;
    if(!comparecutonobjtofilt(dieg33CsusFiltPt, dieg33CsusFiltN, egusRecoPt, dieg33caloidlidegusidx)) cout<<event<<"Mismatching objects: HLT_DoublePhoton70"<<endl;;
    if(dieg33caloidlidegidx.size()>=1 && dieg33caloidlidegusidx.size()>=2) dipho33caloidltrig = true;
    if(HLT_DoublePhoton33_CaloIdL==true && dipho33caloidltrig==false) cout<<event<<": Type 1 - too many cuts on obj - trigger cross-check failed for: HLT_DoublePhoton33"<<endl;
    if(HLT_DoublePhoton33_CaloIdL==false && dipho33caloidltrig==true) cout<<event<<": Type 2(accep) - not enough cuts on obj - trigger cross-check failed for: HLT_DoublePhoton33"<<endl;

    if(HLT_DoublePhoton70==true) dieg70idus++;
    if(HLT_DoublePhoton33_CaloIdL==true) dieg33caloidl++;
    
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
      if((elelead+elesublead).M()>70 && (elelead+elesublead).M()<100) {
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
    if(seleletightidegidx.size()>=2) seleletightid++;
    if(seleletightidegusidx.size()>=2) seleletightidus++;
  
    // Perform genmatching and fill the histograms
    if(isMC && noselegusidx.size()>=1) fillhistineventgenmchunseeded("gennoselAnoselus", gennoselegidx, noselegusidx);
    if(isMC && noselegusidx.size()>=1) fillhistineventgenmchunseeded("genbasicselbarAnoselus", genbasicselbaregidx, noselegusidx);
    if(isMC && noselegusidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selbarAnoselus", genbasicptgt10selbaregidx, noselegusidx);
    if(isMC && noselegusidx.size()>=1) fillhistineventgenmchunseeded("genetabin14_16_24Anoselus", genetabin14_16_24egidx, noselegusidx);
    if(isMC && noselegusidx.size()>=1) fillhistineventgenmchunseeded("genptgt10Anoselus", genptgt10egidx, noselegusidx);
    if(isMC && noselegusidx.size()>=1) fillhistineventgenmchunseeded("genptgt10etalt12Anoselus", genptgt10etalt12egidx, noselegusidx);
    if(isMC && noselegusidx.size()>=1) fillhistineventgenmchunseeded("genptgt10etabin16_24Anoselus", genptgt10etabin16_24egidx, noselegusidx);
    if(isMC && noselegusidx.size()>=1) fillhistineventgenmchunseeded("genptgt10etabin14_16_24Anoselus", genptgt10etabin14_16_24egidx, noselegusidx);
    if(isMC && noselegusidx.size()>=1) fillhistineventgenmchunseeded("genptgt10etabin14_16_24d0lt1cmAnoselus", genptgt10etabin14_16_24d0lt1cmegidx, noselegusidx);
    if(isMC && basicselegusidx.size()>=1) fillhistineventgenmchunseeded("genptgt10Abasicselus", genptgt10egidx, basicselegusidx);
    if(isMC && basicselegusidx.size()>=1) fillhistineventgenmchunseeded("genbasicselptgt15Abasicselus", genbasicselptgt15egidx, basicselegusidx);
    if(isMC && basicselegusidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selbarAbasicselus", genbasicptgt10selbaregidx, basicselegusidx);
    if(isMC && basicselegusidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selecAbasicselus", genbasicptgt10selecegidx, basicselegusidx);
    if(isMC && dieg33caloidlidegusidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selbarAdieg33caloidlus", genbasicptgt10selbaregidx, dieg33caloidlidegusidx);
    if(isMC && dieg33caloidlidegusidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selecAdieg33caloidlus", genbasicptgt10selecegidx, dieg33caloidlidegusidx);
    if(isMC && cut1usidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selbarAcut1us", genbasicptgt10selbaregidx, cut1usidx);
    if(isMC && cut1usidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selecAcut1us", genbasicptgt10selecegidx, cut1usidx);
    if(isMC && cut2usidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selbarAcut2us", genbasicptgt10selbaregidx, cut2usidx);
    if(isMC && cut2usidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selecAcut2us", genbasicptgt10selecegidx, cut2usidx);
    if(isMC && cut3usidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selbarAcut3us", genbasicptgt10selbaregidx, cut3usidx);
    if(isMC && cut3usidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selecAcut3us", genbasicptgt10selecegidx, cut3usidx);
    if(isMC && cuttimedelayonlyusidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selbarAcuttimedelayonlyus", genbasicptgt10selbaregidx, cuttimedelayonlyusidx);
    if(isMC && cuttimedelayonlyusidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selecAcuttimedelayonlyus", genbasicptgt10selecegidx, cuttimedelayonlyusidx);
    if(isMC && cuttimedelaysminusidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selbarAcuttimedelaysminus", genbasicptgt10selbaregidx, cuttimedelaysminusidx);
    if(isMC && cuttimedelaysminusidx.size()>=1) fillhistineventgenmchunseeded("genbasicptgt10selecAcuttimedelaysminus", genbasicptgt10selecegidx, cuttimedelaysminusidx);
    
    // Clear all the vectors
    genelpos.clear();
    gennoselegidx.clear();
    genbasicselbaregidx.clear();
    genbasicptgt10selbaregidx.clear();
    genbasicptgt10selecegidx.clear();
    genetabin14_16_24egidx.clear();
    genptgt10egidx.clear();
    genptgt10etalt12egidx.clear();
    genptgt10etabin16_24egidx.clear();
    genptgt10etabin14_16_24egidx.clear();
    genptgt10etabin14_16_24d0lt1cmegidx.clear();
    genbarselegidx.clear();
    genbarselptgt10egidx.clear();
    genbasicselptgt15egidx.clear();
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
    seleletightcaloidegusidx.clear();
    dieg70idegidx.clear();
    dieg70idegusidx.clear();
    dieg33caloidlidegidx.clear();
    dieg33caloidlidegusidx.clear();
    cut1usidx.clear();
    cut2usidx.clear();
    cut3usidx.clear();
    cuttimedelayonlyusidx.clear();
    cuttimedelaysminusidx.clear();

  } // End of loop on events

  cout<<totEntries<<"\t"<<nosel<<"\t"<<noselus<<"\t"<<basicsel<<"\t"<<basicselus<<"\t"<<selelevetoid<<"\t"<<selelevetoidus<<"\t"<<selelevetozwindidus<<"\t"<<selelevetozoppoidus<<"\t"<<dieg70idus<<"\t"<<dieg33caloidl<<endl;
}

vector< pair<int,int> > data_robustanalyzer::doGenMatchingUnseeded(vector<int> genidx, vector<int> egusidx) {

  if(genidx.size()!=2) {
    throw "Error! Analysis code only suitable for 2 gen electrons";
  }
  
  if(!isMC) {
    throw "Error! Cannot do gen matching. Not MC file.";
  }
  
  // Counter for the total no.of gen matches
  // and 1 history variable for the first gen matched position
  // Logic does not work with more than 2 gen el.
  vector< pair<int,int> > *genegusmch = new vector< pair<int,int> >;
  genegusmch->push_back(make_pair(genidx[0],-1));
  genegusmch->push_back(make_pair(genidx[1],-1));

  // Find the egusidx with the best angular match
  for(int eg : egusidx) {
    // Loop over the gen particles
    for(auto it=genegusmch->begin(); it!=genegusmch->end(); it++) {

      int genidx = (*it).first;
      int usidx = (*it).second;
      if(genidx==-1 || usidx!=-1) continue;

      // Changed to gen match with prompt equivalent of eta and phi
      double diffeta = abs(egusRecoEta[eg]-genLepPromptEta[genidx]); 
      TLorentzVector vecegus, vecgen, vecpromptgen;
      vecgen.SetPtEtaPhiM(genLepPt[genidx],genLepEta[genidx],genLepPhi[genidx],0.0005);
      vecpromptgen.SetPtEtaPhiM(vecgen.P()*TMath::Sin(2*TMath::ATan(TMath::Exp(-genLepPromptEta[genidx]))),genLepPromptEta[genidx],genLepPromptPhi[genidx],0.0005);
      vecegus.SetPtEtaPhiM(egusRecoPt[eg],egusRecoEta[eg],egusRecoPhi[eg],0.0005);
      double qdiffphi = (genLepPid[genidx]/abs(genLepPid[genidx]))*(vecpromptgen.DeltaPhi(vecegus));
      // Condition for gen matching
      if(abs(egusRecoEta[eg])<1.479) {
	if(diffeta<0.1 && qdiffphi<0.15 && qdiffphi>-0.25) {
	  (*it).first = genidx;
	  (*it).second = eg;
	}
      }
      else {
	if(diffeta<0.05 && qdiffphi<0.1 && qdiffphi>-0.15) {
	  (*it).first = genidx;
	  (*it).second = eg;
	}
      }
    } // End of gen loop
    
  } // End of unseeded egamma object loop
  
  return (*genegusmch);
}

// Function to fill a set of histograms for gen particles
void data_robustanalyzer::fillgenhistinevent(TString selection, vector<int> egidx) {

  if(egidx.size()!=2) throw "Error! Code always has 2 indices for gen electrons";

  TH1F* egmult = (TH1F*) outfile->Get(selection+"geneg_egmult");
  TH1F* egmompid = (TH1F*) outfile->Get(selection+"geneg_egmompid");
  TH1F* pt = (TH1F*) outfile->Get(selection+"geneg_pt");
  TH1F* eta = (TH1F*) outfile->Get(selection+"geneg_eta");
  TH1F* phi = (TH1F*) outfile->Get(selection+"geneg_phi");
  TH1F* prompteta = (TH1F*) outfile->Get(selection+"geneg_prompteta");
  TH1F* promptphi = (TH1F*) outfile->Get(selection+"geneg_promptphi");
  TH1F* vx= (TH1F*) outfile->Get(selection+"geneg_vx");
  TH1F* vy = (TH1F*) outfile->Get(selection+"geneg_vy");
  TH1F* vz = (TH1F*) outfile->Get(selection+"geneg_vz");
  TH1F* pvxh = (TH1F*) outfile->Get(selection+"geneg_pvx");
  TH1F* pvyh = (TH1F*) outfile->Get(selection+"geneg_pvy");
  TH1F* pvzh = (TH1F*) outfile->Get(selection+"geneg_pvz");
  TH1F* deltaetamom = (TH1F*) outfile->Get(selection+"geneg_deltaetamom");
  TH1F* deltaphimom = (TH1F*) outfile->Get(selection+"geneg_deltaphimom");
  TH1F* deltaRmom = (TH1F*) outfile->Get(selection+"geneg_deltaRmom");
  TH1F* gend0 = (TH1F*) outfile->Get(selection+"geneg_d0");
  TH1F* log10d0 = (TH1F*) outfile->Get(selection+"geneg_log10d0");
  TH1F* genlxy = (TH1F*) outfile->Get(selection+"geneg_lxy");
  TH1F* log10lxy = (TH1F*) outfile->Get(selection+"geneg_log10lxy");
  TH1F* t1 = (TH1F*) outfile->Get(selection+"geneg_t1");
  TH1F* t0 = (TH1F*) outfile->Get(selection+"geneg_t0");
  TH1F* t1mt0 = (TH1F*) outfile->Get(selection+"geneg_t1mt0");
  TH1F* leadegpt = (TH1F*) outfile->Get(selection+"geneg_leadpt");
  TH1F* leadegeta = (TH1F*) outfile->Get(selection+"geneg_leadeta");
  TH1F* leadegphi = (TH1F*) outfile->Get(selection+"geneg_leadphi");
  TH1F* leaddeltaetamom = (TH1F*) outfile->Get(selection+"geneg_leaddeltaetamom");
  TH1F* leaddeltaphimom = (TH1F*) outfile->Get(selection+"geneg_leaddeltaphimom");
  TH1F* leaddeltaRmom = (TH1F*) outfile->Get(selection+"geneg_leaddeltaRmom");
  TH1F* leadegd0 = (TH1F*) outfile->Get(selection+"geneg_leadd0");
  TH1F* leadeglog10d0 = (TH1F*) outfile->Get(selection+"geneg_leadlog10d0");
  TH1F* leadeglxy = (TH1F*) outfile->Get(selection+"geneg_leadlxy");
  TH1F* leadeglog10lxy = (TH1F*) outfile->Get(selection+"geneg_leadlog10lxy");
  TH1F* subleadegpt = (TH1F*) outfile->Get(selection+"geneg_subleadpt");
  TH1F* subleadegeta = (TH1F*) outfile->Get(selection+"geneg_subleadeta");
  TH1F* subleadegphi = (TH1F*) outfile->Get(selection+"geneg_subleadphi");
  TH1F* subleaddeltaetamom = (TH1F*) outfile->Get(selection+"geneg_subleaddeltaetamom");
  TH1F* subleaddeltaphimom = (TH1F*) outfile->Get(selection+"geneg_subleaddeltaphimom");
  TH1F* subleaddeltaRmom = (TH1F*) outfile->Get(selection+"geneg_subleaddeltaRmom");
  TH1F* subleadegd0 = (TH1F*) outfile->Get(selection+"geneg_subleadd0");
  TH1F* subleadeglog10d0 = (TH1F*) outfile->Get(selection+"geneg_subleadlog10d0");
  TH1F* subleadeglxy = (TH1F*) outfile->Get(selection+"geneg_subleadlxy");
  TH1F* subleadeglog10lxy = (TH1F*) outfile->Get(selection+"geneg_subleadlog10lxy");

  int genelmult = 0;

  if(egidx[0] != -1) {
    TVector3 el, elmom;
    el.SetPtEtaPhi(genLepPt[egidx[0]], genLepEta[egidx[0]], genLepPhi[egidx[0]]);
    elmom.SetPtEtaPhi(genLepMomPt[egidx[0]], genLepMomEta[egidx[0]], genLepMomPhi[egidx[0]]);
    egmompid->Fill(genLepMomPid[egidx[0]]);
    pt->Fill(genLepPt[egidx[0]]);
    eta->Fill(genLepEta[egidx[0]]);
    phi->Fill(genLepPhi[egidx[0]]);
    prompteta->Fill(genLepPromptEta[egidx[0]]);
    promptphi->Fill(genLepPromptPhi[egidx[0]]);
    vx->Fill(genLepVx[egidx[0]]);
    vy->Fill(genLepVy[egidx[0]]);
    vz->Fill(genLepVz[egidx[0]]);
    pvxh->Fill(pvx[egidx[0]]);
    pvyh->Fill(pvy[egidx[0]]);
    pvzh->Fill(pvz[egidx[0]]);
    deltaetamom->Fill(genLepMomEta[egidx[0]]-genLepEta[egidx[0]]);
    deltaphimom->Fill(elmom.DeltaPhi(el));
    deltaRmom->Fill(elmom.DeltaR(el));
    gend0->Fill(genLepDxy[egidx[0]]);
    log10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[egidx[0]])));
    genlxy->Fill(genLepLxy[egidx[0]]);
    log10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[egidx[0]])));
    if(genLepTimeAct[egidx[0]]>(-1e9)) {
      t1->Fill(genLepTimeAct[egidx[0]]);
      t0->Fill(genLepTimeLight[egidx[0]]);
      t1mt0->Fill(genLepTimeDiff[egidx[0]]);
    }
    leadegpt->Fill(genLepPt[egidx[0]]);
    leadegeta->Fill(genLepEta[egidx[0]]);
    leadegphi->Fill(genLepPhi[egidx[0]]);
    leaddeltaetamom->Fill(genLepMomEta[egidx[0]]-genLepEta[egidx[0]]);
    leaddeltaphimom->Fill(elmom.DeltaPhi(el));
    leaddeltaRmom->Fill(elmom.DeltaR(el));
    leadegd0->Fill(genLepDxy[egidx[0]]);
    leadeglog10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[egidx[0]])));
    leadeglxy->Fill(genLepLxy[egidx[0]]);
    leadeglog10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[egidx[0]])));
    genelmult++;
  }
  if(egidx[1] != -1) {
    TVector3 el, elmom;
    el.SetPtEtaPhi(genLepPt[egidx[1]], genLepEta[egidx[1]], genLepPhi[egidx[1]]);
    elmom.SetPtEtaPhi(genLepMomPt[egidx[1]], genLepMomEta[egidx[1]], genLepMomPhi[egidx[1]]);
    egmompid->Fill(genLepMomPid[egidx[1]]);
    pt->Fill(genLepPt[egidx[1]]);
    eta->Fill(genLepEta[egidx[1]]);
    phi->Fill(genLepPhi[egidx[1]]);
    prompteta->Fill(genLepPromptEta[egidx[1]]);
    promptphi->Fill(genLepPromptPhi[egidx[1]]);
    vx->Fill(genLepVx[egidx[1]]);
    vy->Fill(genLepVy[egidx[1]]);
    vz->Fill(genLepVz[egidx[1]]);
    deltaetamom->Fill(genLepMomEta[egidx[1]]-genLepEta[egidx[1]]);
    deltaphimom->Fill(elmom.DeltaPhi(el));
    deltaRmom->Fill(elmom.DeltaR(el));
    gend0->Fill(genLepDxy[egidx[1]]);
    log10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[egidx[1]])));
    genlxy->Fill(genLepLxy[egidx[1]]);
    log10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[egidx[1]])));
    if(genLepTimeAct[egidx[1]]>(-1e9)) {
      t1->Fill(genLepTimeAct[egidx[1]]);
      t0->Fill(genLepTimeLight[egidx[1]]);
      t1mt0->Fill(genLepTimeDiff[egidx[1]]);
    }
    subleadegpt->Fill(genLepPt[egidx[1]]);
    subleadegeta->Fill(genLepEta[egidx[1]]);
    subleadegphi->Fill(genLepPhi[egidx[1]]);
    subleaddeltaetamom->Fill(genLepMomEta[egidx[1]]-genLepEta[egidx[1]]);
    subleaddeltaphimom->Fill(elmom.DeltaPhi(el));
    subleaddeltaRmom->Fill(elmom.DeltaR(el));
    subleadegd0->Fill(genLepDxy[egidx[1]]);
    subleadeglog10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[egidx[1]])));
    subleadeglxy->Fill(genLepLxy[egidx[1]]);
    subleadeglog10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[egidx[1]])));
    genelmult++;
  }

  egmult->Fill(genelmult);
  egidx.clear();
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
	recoeb_leadegpixelmchvar_s2->Fill(-50);
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
	recoee_leadegpixelmchvar_s2->Fill(-50);
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
	recoeb_subleadegpixelmchvar_s2->Fill(-50);
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
	recoee_subleadegpixelmchvar_s2->Fill(-50);
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
  
  // Both lead and sub-lead electrons
  TH1F* recoeb_egseedclustime = (TH1F*) outfile->Get(selection+"recoebus_egseedclustime");  
  TH1F* recoeb_egpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoebus_egpixelmchvar_s2");
  TH1F* recoee_egseedclustime = (TH1F*) outfile->Get(selection+"recoeeus_egseedclustime");  
  TH1F* recoee_egpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"recoeeus_egpixelmchvar_s2");

  // Get barrel variables - lead pt unseeded e/gamma
  TH1F* recoeb_leadegclustershape = (TH1F*) outfile->Get(selection+"recoebus_leadegclustershape");
  TH1F* recoeb_leadegin5x5clusshape = (TH1F*) outfile->Get(selection+"recoebus_leadegin5x5clusshape");
  TH1F* recoeb_leadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"recoebus_leadegin5x5noiseclnd");
  TH1F* recoeb_leadegsmin = (TH1F*) outfile->Get(selection+"recoebus_leadegsmin");
  TH1F* recoeb_leadegsmaj = (TH1F*) outfile->Get(selection+"recoebus_leadegsmaj");
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
  TH1F* recoeb_subleadegsmin = (TH1F*) outfile->Get(selection+"recoebus_subleadegsmin");
  TH1F* recoeb_subleadegsmaj = (TH1F*) outfile->Get(selection+"recoebus_subleadegsmaj");
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
  TH1F* recoee_leadegsmin = (TH1F*) outfile->Get(selection+"recoeeus_leadegsmin");
  TH1F* recoee_leadegsmaj = (TH1F*) outfile->Get(selection+"recoeeus_leadegsmaj");
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
  TH1F* recoee_subleadegsmin = (TH1F*) outfile->Get(selection+"recoeeus_subleadegsmin");
  TH1F* recoee_subleadegsmaj = (TH1F*) outfile->Get(selection+"recoeeus_subleadegsmaj");
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
      recoeb_leadegsmin->Fill(egushltEgammaClusterShape_smin[egidx[0]]);
      recoeb_leadegsmaj->Fill(egushltEgammaClusterShape_smaj[egidx[0]]);
      //recoeb_leadegscenergy->Fill(egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_leadeghovere->Fill(egushltEgammaHoverE[egidx[0]]);
      recoeb_leadeghovereoversupcluse->Fill(egushltEgammaHoverE[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_leadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]);
      recoeb_leadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoeb_leadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[egidx[0]]);
      recoeb_leadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[egidx[0]]/egushltEgammaHoverE[egidx[0]]);
      if(egushltEcalSeedClusterTime[egidx[0]]!=0) recoeb_leadegseedclustime->Fill(egushltEcalSeedClusterTime[egidx[0]]);
      if(egushltEcalSeedClusterTime[egidx[0]]!=0) recoeb_egseedclustime->Fill(egushltEcalSeedClusterTime[egidx[0]]);
      if(egushltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoeb_leadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[0]]);
	recoeb_egpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoeb_leadegpixelmchvar_s2->Fill(-50);
	recoeb_egpixelmchvar_s2->Fill(-50);
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
      recoee_leadegsmin->Fill(egushltEgammaClusterShape_smin[egidx[0]]);
      recoee_leadegsmaj->Fill(egushltEgammaClusterShape_smaj[egidx[0]]);
      //recoee_leadegscenergy->Fill(egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_leadeghovere->Fill(egushltEgammaHoverE[egidx[0]]);
      recoee_leadeghovereoversupcluse->Fill(egushltEgammaHoverE[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_leadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]);
      recoee_leadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[egidx[0]]/egushltEgammaSuperClusterEnergy[egidx[0]]);
      recoee_leadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[egidx[0]]);	
      recoee_leadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[egidx[0]]/egushltEgammaHoverE[egidx[0]]);	
      if(egushltEcalSeedClusterTime[egidx[0]]!=0) recoee_leadegseedclustime->Fill(egushltEcalSeedClusterTime[egidx[0]]);
      if(egushltEcalSeedClusterTime[egidx[0]]!=0) recoee_egseedclustime->Fill(egushltEcalSeedClusterTime[egidx[0]]);
      if(egushltEgammaPixelMatchVars_s2[egidx[0]]<TMath::Power(10,36)) {
	recoee_leadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[0]]);
	recoee_egpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[0]]);
      }
      else {
	recoee_leadegpixelmchvar_s2->Fill(-50);
	recoee_egpixelmchvar_s2->Fill(-50);
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
      recoeb_subleadegsmin->Fill(egushltEgammaClusterShape_smin[egidx[1]]);
      recoeb_subleadegsmaj->Fill(egushltEgammaClusterShape_smaj[egidx[1]]);
      recoeb_subleadegscenergy->Fill(egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoeb_subleadeghovere->Fill(egushltEgammaHoverE[egidx[1]]);
      recoeb_subleadeghovereoversupcluse->Fill(egushltEgammaHoverE[egidx[1]]/egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoeb_subleadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[egidx[1]]);
      recoeb_subleadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[egidx[1]]/egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoeb_subleadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[egidx[1]]);
      recoeb_subleadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[egidx[1]]/egushltEgammaHoverE[egidx[1]]);
      if(egushltEcalSeedClusterTime[egidx[1]]!=0) recoeb_subleadegseedclustime->Fill(egushltEcalSeedClusterTime[egidx[1]]);
      if(egushltEcalSeedClusterTime[egidx[1]]!=0) recoeb_egseedclustime->Fill(egushltEcalSeedClusterTime[egidx[1]]);
      if(egushltEgammaPixelMatchVars_s2[egidx[1]]<TMath::Power(10,36)) {
	recoeb_subleadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[1]]);
	recoeb_egpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[1]]);
      }
      else {
	recoeb_subleadegpixelmchvar_s2->Fill(-50);
	recoeb_egpixelmchvar_s2->Fill(-50);
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
      recoee_subleadegsmin->Fill(egushltEgammaClusterShape_smin[egidx[1]]);
      recoee_subleadegsmaj->Fill(egushltEgammaClusterShape_smaj[egidx[1]]);
      recoee_subleadegscenergy->Fill(egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoee_subleadeghovere->Fill(egushltEgammaHoverE[egidx[1]]);
      recoee_subleadeghovereoversupcluse->Fill(egushltEgammaHoverE[egidx[1]]/egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoee_subleadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[egidx[1]]);
      recoee_subleadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[egidx[1]]/egushltEgammaSuperClusterEnergy[egidx[1]]);
      recoee_subleadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[egidx[1]]);	
      recoee_subleadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[egidx[1]]/egushltEgammaHoverE[egidx[1]]);	
      if(egushltEcalSeedClusterTime[egidx[1]]!=0) recoee_subleadegseedclustime->Fill(egushltEcalSeedClusterTime[egidx[1]]);
      if(egushltEcalSeedClusterTime[egidx[1]]!=0) recoee_egseedclustime->Fill(egushltEcalSeedClusterTime[egidx[1]]);
      if(egushltEgammaPixelMatchVars_s2[egidx[1]]<TMath::Power(10,36)) {
	recoee_subleadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[1]]);
	recoee_egpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[egidx[1]]);
      }
      else {
	recoee_subleadegpixelmchvar_s2->Fill(-50);
	recoee_egpixelmchvar_s2->Fill(-50);
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

}

// Function to fill a set of histograms in the event for gen matched objects - unseeded egamma objects
void data_robustanalyzer::fillhistineventgenmchunseeded(TString selection, vector<int> genidx, vector<int> egusidx) {

  if(genidx.size()!=2) throw "Error! Code always has 2 indices for gen electrons";
  
  // Variables before gen matching
  TH1F* geneltrigebus_dE = (TH1F*) outfile->Get(selection+"geneltrigebus_dE");
  TH1F* geneltrigebus_dPt = (TH1F*) outfile->Get(selection+"geneltrigebus_dPt");
  TH1F* geneltrigebus_dEta = (TH1F*) outfile->Get(selection+"geneltrigebus_dEta");
  TH1F* geneltrigebus_qdPhi = (TH1F*) outfile->Get(selection+"geneltrigebus_qdPhi");
  TH1F* geneltrigebus_dPromptEta = (TH1F*) outfile->Get(selection+"geneltrigebus_dPromptEta");
  TH1F* geneltrigebus_qdPromptPhi = (TH1F*) outfile->Get(selection+"geneltrigebus_qdPromptPhi");
  TH1F* geneltrigeeus_dE = (TH1F*) outfile->Get(selection+"geneltrigeeus_dE");
  TH1F* geneltrigeeus_dPt = (TH1F*) outfile->Get(selection+"geneltrigeeus_dPt");
  TH1F* geneltrigeeus_dEta = (TH1F*) outfile->Get(selection+"geneltrigeeus_dEta");
  TH1F* geneltrigeeus_qdPhi = (TH1F*) outfile->Get(selection+"geneltrigeeus_qdPhi");

  // Variables after gen matching
  TH1F* genegmult = (TH1F*) outfile->Get(selection+"recomchgenel_egmult");
  TH1F* genpt = (TH1F*) outfile->Get(selection+"recomchgenel_pt");
  TH1F* geneta = (TH1F*) outfile->Get(selection+"recomchgenel_eta");
  TH1F* genphi = (TH1F*) outfile->Get(selection+"recomchgenel_phi");
  TH1F* gend0 = (TH1F*) outfile->Get(selection+"recomchgenel_d0");
  TH1F* genlog10d0 = (TH1F*) outfile->Get(selection+"recomchgenel_log10d0");
  TH1F* genlxy = (TH1F*) outfile->Get(selection+"recomchgenel_lxy");
  TH1F* genlog10lxy = (TH1F*) outfile->Get(selection+"recomchgenel_log10lxy");
  TH1F* genleadpt = (TH1F*) outfile->Get(selection+"recomchgenel_leadpt");
  TH1F* genleadeta = (TH1F*) outfile->Get(selection+"recomchgenel_leadeta");
  TH1F* genleadphi = (TH1F*) outfile->Get(selection+"recomchgenel_leadphi");
  TH1F* genleadd0 = (TH1F*) outfile->Get(selection+"recomchgenel_leadd0");
  TH1F* genleadlog10d0 = (TH1F*) outfile->Get(selection+"recomchgenel_leadlog10d0");
  TH1F* genleadlxy = (TH1F*) outfile->Get(selection+"recomchgenel_leadlxy");
  TH1F* genleadlog10lxy = (TH1F*) outfile->Get(selection+"recomchgenel_leadlog10lxy");
  TH1F* gensubleadpt = (TH1F*) outfile->Get(selection+"recomchgenel_subleadpt");
  TH1F* gensubleadeta = (TH1F*) outfile->Get(selection+"recomchgenel_subleadeta");
  TH1F* gensubleadphi = (TH1F*) outfile->Get(selection+"recomchgenel_subleadphi");
  TH1F* gensubleadd0 = (TH1F*) outfile->Get(selection+"recomchgenel_subleadd0");
  TH1F* gensubleadlog10d0 = (TH1F*) outfile->Get(selection+"recomchgenel_subleadlog10d0");
  TH1F* gensubleadlxy = (TH1F*) outfile->Get(selection+"recomchgenel_subleadlxy");
  TH1F* gensubleadlog10lxy = (TH1F*) outfile->Get(selection+"recomchgenel_subleadlog10lxy");

  // Varibales for gen matched triiger electrons
  TH1F* leadegpt = (TH1F*) outfile->Get(selection+"genmchrecoegus_leadegpt");
  TH1F* leadegeta = (TH1F*) outfile->Get(selection+"genmchrecoegus_leadegeta");
  TH1F* leadegphi = (TH1F*) outfile->Get(selection+"genmchrecoegus_leadegphi");

  TH1F* subleadegpt = (TH1F*) outfile->Get(selection+"genmchrecoegus_subleadegpt");
  TH1F* subleadegeta = (TH1F*) outfile->Get(selection+"genmchrecoegus_subleadegeta");
  TH1F* subleadegphi = (TH1F*) outfile->Get(selection+"genmchrecoegus_subleadegphi");

  // Both lead and sub-lead electrons
  TH1F* recoeb_egseedclustime = (TH1F*) outfile->Get(selection+"genmchrecoebus_egseedclustime");  
  TH1F* recoeb_egpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"genmchrecoebus_egpixelmchvar_s2");
  TH1F* recoeb_pxlmch22_egpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"genmchrecoebus_pxlmch22_egpixelmchvar_s2");
  TH1F* recoee_egseedclustime = (TH1F*) outfile->Get(selection+"genmchrecoeeus_egseedclustime");  
  TH1F* recoee_egpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"genmchrecoeeus_egpixelmchvar_s2");
  TH1F* recoee_pxlmch22_egpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"genmchrecoeeus_pxlmch22_egpixelmchvar_s2");

  // Get barrel variables - lead pt unseeded e/gamma
  TH1F* recoeb_leadegclustershape = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegclustershape");
  TH1F* recoeb_leadegin5x5clusshape = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegin5x5clusshape");
  TH1F* recoeb_leadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegin5x5noiseclnd");
  TH1F* recoeb_leadegsmin = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegsmin");
  TH1F* recoeb_leadegsmaj = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegsmaj");
  TH1F* recoeb_leadegscenergy = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegscenergy");
  TH1F* recoeb_leadeghovere = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadeghovere");
  TH1F* recoeb_leadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadeghovereoversupcluse");
  TH1F* recoeb_leadegecalpfclustiso = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegecalpfclustiso");
  TH1F* recoeb_leadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegecalpfclustisoovere");
  TH1F* recoeb_leadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadeghcalpfclustiso");
  TH1F* recoeb_leadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadeghcalpfclustisoovere");
  TH1F* recoeb_leadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegpixelmchvar_s2");
  TH1F* recoeb_leadegtrkiso = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegtrkiso");
  TH1F* recoeb_leadegchi2 = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegchi2");
  TH1F* recoeb_leadegdeta = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegdeta");
  TH1F* recoeb_leadegdetaseed = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegdetaseed");
  TH1F* recoeb_leadegdphi = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegdphi");
  TH1F* recoeb_leadegmhits = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegmhits");
  TH1F* recoeb_leadegnlayerit = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegnlayerit");
  TH1F* recoeb_leadegooeseedoop = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegooeseedoop");
  TH1F* recoeb_leadegooesclsoop = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegooesclsoop");
  TH1F* recoeb_leadegvalhits = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegvalhits");  
  TH1F* recoeb_leadegseedclustime = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadegseedclustime");  
  
  // Get barrel variables - sublead pt unseeded e/gamma
  TH1F* recoeb_subleadegclustershape = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegclustershape");
  TH1F* recoeb_subleadegin5x5clusshape = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegin5x5clusshape");
  TH1F* recoeb_subleadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegin5x5noiseclnd");
  TH1F* recoeb_subleadegsmin = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegsmin");
  TH1F* recoeb_subleadegsmaj = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegsmaj");
  TH1F* recoeb_subleadegscenergy = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegscenergy");
  TH1F* recoeb_subleadeghovere = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadeghovere");
  TH1F* recoeb_subleadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadeghovereoversupcluse");
  TH1F* recoeb_subleadegecalpfclustiso = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegecalpfclustiso");
  TH1F* recoeb_subleadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegecalpfclustisoovere");
  TH1F* recoeb_subleadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadeghcalpfclustiso");
  TH1F* recoeb_subleadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadeghcalpfclustisoovere");
  TH1F* recoeb_subleadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegpixelmchvar_s2");
  TH1F* recoeb_subleadegtrkiso = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegtrkiso");
  TH1F* recoeb_subleadegchi2 = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegchi2");
  TH1F* recoeb_subleadegdeta = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegdeta");
  TH1F* recoeb_subleadegdetaseed = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegdetaseed");
  TH1F* recoeb_subleadegdphi = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegdphi");
  TH1F* recoeb_subleadegmhits = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegmhits");
  TH1F* recoeb_subleadegnlayerit = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegnlayerit");
  TH1F* recoeb_subleadegooeseedoop = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegooeseedoop");
  TH1F* recoeb_subleadegooesclsoop = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegooesclsoop");
  TH1F* recoeb_subleadegvalhits = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegvalhits");  
  TH1F* recoeb_subleadegseedclustime = (TH1F*) outfile->Get(selection+"genmchrecoebus_subleadegseedclustime");  

  // invariant mass - barrel
  TH1F* recoeb_leadsubleadM = (TH1F*) outfile->Get(selection+"genmchrecoebus_leadsubleadM");
  TH1F* genmchgeneltrigebus_dE = (TH1F*) outfile->Get(selection+"genmchgeneltrigebus_dE");
  TH1F* genmchgeneltrigebus_dPt = (TH1F*) outfile->Get(selection+"genmchgeneltrigebus_dPt");
  TH1F* genmchgeneltrigebus_dEta = (TH1F*) outfile->Get(selection+"genmchgeneltrigebus_dEta");
  TH1F* genmchgeneltrigebus_qdPhi = (TH1F*) outfile->Get(selection+"genmchgeneltrigebus_qdPhi");
  TH1F* genmchgeneltrigebus_dPromptE = (TH1F*) outfile->Get(selection+"genmchgeneltrigebus_dPromptE");
  TH1F* genmchgeneltrigebus_dPromptPt = (TH1F*) outfile->Get(selection+"genmchgeneltrigebus_dPromptPt");
  TH1F* genmchgeneltrigebus_dPromptEta = (TH1F*) outfile->Get(selection+"genmchgeneltrigebus_dPromptEta");
  TH1F* genmchgeneltrigebus_qdPromptPhi = (TH1F*) outfile->Get(selection+"genmchgeneltrigebus_qdPromptPhi");
  
  // Get end-cap variables - lead pt unseeded e/gamma
  TH1F* recoee_leadegclustershape = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegclustershape");
  TH1F* recoee_leadegin5x5clusshape = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegin5x5clusshape");
  TH1F* recoee_leadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegin5x5noiseclnd");
  TH1F* recoee_leadegsmin = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegsmin");
  TH1F* recoee_leadegsmaj = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegsmaj");
  TH1F* recoee_leadegscenergy = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegscenergy");
  TH1F* recoee_leadeghovere = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadeghovere");
  TH1F* recoee_leadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadeghovereoversupcluse");
  TH1F* recoee_leadegecalpfclustiso = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegecalpfclustiso");
  TH1F* recoee_leadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegecalpfclustisoovere");
  TH1F* recoee_leadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadeghcalpfclustiso");
  TH1F* recoee_leadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadeghcalpfclustisoovere");
  TH1F* recoee_leadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegpixelmchvar_s2");
  TH1F* recoee_leadegtrkiso = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegtrkiso");
  TH1F* recoee_leadegchi2 = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegchi2");
  TH1F* recoee_leadegdeta = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegdeta");
  TH1F* recoee_leadegdetaseed = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegdetaseed");
  TH1F* recoee_leadegdphi = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegdphi");
  TH1F* recoee_leadegmhits = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegmhits");
  TH1F* recoee_leadegnlayerit = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegnlayerit");
  TH1F* recoee_leadegooeseedoop = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegooeseedoop");
  TH1F* recoee_leadegooesclsoop = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegooesclsoop");
  TH1F* recoee_leadegvalhits = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegvalhits");
  TH1F* recoee_leadegseedclustime = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadegseedclustime");  
  
  // Get end-cap variables - sublead pt unseeded e/gamma
  TH1F* recoee_subleadegclustershape = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegclustershape");
  TH1F* recoee_subleadegin5x5clusshape = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegin5x5clusshape");
  TH1F* recoee_subleadegin5x5noiseclnd = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegin5x5noiseclnd");
  TH1F* recoee_subleadegsmin = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegsmin");
  TH1F* recoee_subleadegsmaj = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegsmaj");
  TH1F* recoee_subleadegscenergy = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegscenergy");
  TH1F* recoee_subleadeghovere = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadeghovere");
  TH1F* recoee_subleadeghovereoversupcluse = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadeghovereoversupcluse");
  TH1F* recoee_subleadegecalpfclustiso = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegecalpfclustiso");
  TH1F* recoee_subleadegecalpfclustisoovere = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegecalpfclustisoovere");
  TH1F* recoee_subleadeghcalpfclustiso = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadeghcalpfclustiso");
  TH1F* recoee_subleadeghcalpfclustisoovere = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadeghcalpfclustisoovere");
  TH1F* recoee_subleadegpixelmchvar_s2 = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegpixelmchvar_s2");
  TH1F* recoee_subleadegtrkiso = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegtrkiso");
  TH1F* recoee_subleadegchi2 = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegchi2");
  TH1F* recoee_subleadegdeta = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegdeta");
  TH1F* recoee_subleadegdetaseed = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegdetaseed");
  TH1F* recoee_subleadegdphi = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegdphi");
  TH1F* recoee_subleadegmhits = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegmhits");
  TH1F* recoee_subleadegnlayerit = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegnlayerit");
  TH1F* recoee_subleadegooeseedoop = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegooeseedoop");
  TH1F* recoee_subleadegooesclsoop = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegooesclsoop");
  TH1F* recoee_subleadegvalhits = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegvalhits");
  TH1F* recoee_subleadegseedclustime = (TH1F*) outfile->Get(selection+"genmchrecoeeus_subleadegseedclustime");  
  
  // invariant mass - end-cap
  TH1F* recoee_leadsubleadM = (TH1F*) outfile->Get(selection+"genmchrecoeeus_leadsubleadM");
  TH1F* genmchgeneltrigeeus_dE = (TH1F*) outfile->Get(selection+"genmchgeneltrigeeus_dE");
  TH1F* genmchgeneltrigeeus_dPt = (TH1F*) outfile->Get(selection+"genmchgeneltrigeeus_dPt");
  TH1F* genmchgeneltrigeeus_dEta = (TH1F*) outfile->Get(selection+"genmchgeneltrigeeus_dEta");
  TH1F* genmchgeneltrigeeus_qdPhi = (TH1F*) outfile->Get(selection+"genmchgeneltrigeeus_qdPhi");

  // Minimum dEta and min. qdPhi for each gen electron
  for(int gene : genidx) {
    
    if(gene==-1) continue;
    
    int genq = genLepPid[gene]/TMath::Abs(genLepPid[gene]);
    double dEtamin=1e9, qdPhimin=1e9, dEmin=1e9, dPtmin=1e9, dPromptEtamin=1e9, qdPromptPhimin=1e9;
    
    for(int egus : egusidx) {
      double dEta = genLepEta[gene]-egusRecoEta[egus];
      double dPromptEta = genLepPromptEta[gene]-egusRecoEta[egus];
      TLorentzVector vecegus, vecgen, vecpromptgen;
      vecgen.SetPtEtaPhiM(genLepPt[gene],genLepEta[gene],genLepPhi[gene],0.0005);
      vecpromptgen.SetPtEtaPhiM(vecgen.P()*TMath::Sin(2*TMath::ATan(TMath::Exp(-genLepPromptEta[gene]))),genLepPromptEta[gene],genLepPromptPhi[gene],0.0005);
      vecegus.SetPtEtaPhiM(egusRecoPt[egus],egusRecoEta[egus],egusRecoPhi[egus],0.0005);
      double qdPhi = genq*(vecgen.DeltaPhi(vecegus));
      double qdPromptPhi = genq*(vecpromptgen.DeltaPhi(vecegus));
      double dE = vecgen.E()-vecegus.E();;
      double dPt = genLepPt[gene]-egusRecoPt[egus];
      if(dEtamin==1e9 && qdPhimin==1e9 && dEmin==1e9 && dPtmin==1e9 && dPromptEtamin==1e9 && qdPromptPhimin==1e9) {
	dEtamin = dEta;
	qdPhimin = qdPhi;
	dEmin = dE;
	dPtmin = dPt;
	dPromptEtamin = dPromptEta;
	qdPromptPhimin = qdPromptPhi;
      }
      if(abs(dEta)<abs(dEtamin)) {
	dEtamin = dEta;
      }
      if(abs(qdPhi)<abs(qdPhimin)) {
	qdPhimin = qdPhi;
      }
      if(abs(dE)<abs(dEmin)) {
	dEmin = dE;
      }
      if(abs(dPt)<abs(dPtmin)) {
	dPtmin = dPt;
      }
      if(abs(dPromptEta)<abs(dPromptEtamin)) {
	dPromptEtamin = dPromptEta;
      }
      if(abs(qdPromptPhi)<abs(qdPromptPhimin)) {
	qdPromptPhimin = qdPromptPhi;
      }
    }
    
    if(abs(genLepEta[gene])<1.479) {
      geneltrigebus_dE->Fill(dEmin);
      geneltrigebus_dPt->Fill(dPtmin);
      geneltrigebus_dEta->Fill(dEtamin);
      geneltrigebus_qdPhi->Fill(qdPhimin);
      geneltrigebus_dPromptEta->Fill(dPromptEtamin);
      geneltrigebus_qdPromptPhi->Fill(qdPromptPhimin);
    }
    else {
      geneltrigeeus_dE->Fill(dEmin);
      geneltrigeeus_dPt->Fill(dPtmin);
      geneltrigeeus_dEta->Fill(dEtamin);
      geneltrigeeus_qdPhi->Fill(qdPhimin);
    }
  }

  
  // Find the gen-matched electron
  vector< pair<int,int> > genmatched = doGenMatchingUnseeded(genidx, egusidx);
  // Verified that the gen matching logic works

  int gen1 = genmatched[0].first;
  int egus1 = genmatched[0].second;
  int gen2 = genmatched[1].first;
  int egus2 = genmatched[1].second;
  
  int mchcnt = 0;

  // Fill the gen variables
  if(egus1!=-1) {
    mchcnt++;
    genpt->Fill(genLepPt[gen1]);
    geneta->Fill(genLepEta[gen1]);
    genphi->Fill(genLepPhi[gen1]);
    int genq = genLepPid[gen1]/TMath::Abs(genLepPid[gen1]);
    TLorentzVector el, promptel, vecegus;
    el.SetPtEtaPhiM(genLepPt[gen1],genLepEta[gen1],genLepPhi[gen1],0.0005);
    promptel.SetPtEtaPhiM(el.P()*TMath::Sin(2*TMath::ATan(TMath::Exp(-genLepPromptEta[gen1]))),genLepPromptEta[gen1],genLepPromptPhi[gen1],0.0005);
    vecegus.SetPtEtaPhiM(egusRecoPt[egus1],egusRecoEta[egus1],egusRecoPhi[egus1],0.0005);
    double qdPhi = genq*(el.DeltaPhi(vecegus));
    double qdPromptPhi = genq*(promptel.DeltaPhi(vecegus));
    gend0->Fill(genLepDxy[gen1]);
    genlog10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[gen1])));
    genlxy->Fill(genLepLxy[gen1]);
    genlog10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[gen1])));
    if(abs(genLepEta[gen1])<1.479) {
      genmchgeneltrigebus_dE->Fill(el.E()-vecegus.E());
      genmchgeneltrigebus_dPt->Fill(genLepPt[gen1]-egusRecoPt[egus1]);
      genmchgeneltrigebus_dEta->Fill(genLepEta[gen1]-egusRecoEta[egus1]);
      genmchgeneltrigebus_qdPhi->Fill(qdPhi);
      genmchgeneltrigebus_dPromptE->Fill(promptel.E()-vecegus.E());
      genmchgeneltrigebus_dPromptPt->Fill(promptel.Pt()-egusRecoPt[egus1]);
      genmchgeneltrigebus_dPromptEta->Fill(genLepPromptEta[gen1]-egusRecoEta[egus1]);
      genmchgeneltrigebus_qdPromptPhi->Fill(qdPromptPhi);
    }
    else {
      genmchgeneltrigeeus_dE->Fill(el.E()-vecegus.E());
      genmchgeneltrigeeus_dPt->Fill(genLepPt[gen1]-egusRecoPt[egus1]);
      genmchgeneltrigeeus_dEta->Fill(genLepEta[gen1]-egusRecoEta[egus1]);
      genmchgeneltrigeeus_qdPhi->Fill(qdPhi);
    }
    genleadpt->Fill(genLepPt[gen1]);
    genleadeta->Fill(genLepEta[gen1]);
    genleadphi->Fill(genLepPhi[gen1]);
    genleadd0->Fill(genLepDxy[gen1]);    
    genleadlog10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[gen1])));    
    genleadlxy->Fill(genLepLxy[gen1]);    
    genleadlog10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[gen1])));    
  }
  if(egus2!=-1) {
    mchcnt++;
    genpt->Fill(genLepPt[gen2]);
    geneta->Fill(genLepEta[gen2]);
    genphi->Fill(genLepPhi[gen2]);
    int genq = genLepPid[gen2]/TMath::Abs(genLepPid[gen2]);
    TLorentzVector el, promptel, vecegus;
    el.SetPtEtaPhiM(genLepPt[gen2],genLepEta[gen2],genLepPhi[gen2],0.0005);
    promptel.SetPtEtaPhiM(el.P()*TMath::Sin(2*TMath::ATan(TMath::Exp(-genLepPromptEta[gen2]))),genLepPromptEta[gen2],genLepPromptPhi[gen2],0.0005);
    vecegus.SetPtEtaPhiM(egusRecoPt[egus2],egusRecoEta[egus2],egusRecoPhi[egus2],0.0005);
    double qdPhi = genq*(el.DeltaPhi(vecegus));
    double qdPromptPhi = genq*(promptel.DeltaPhi(vecegus));
    gend0->Fill(genLepDxy[gen2]);
    genlog10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[gen2])));
    genlxy->Fill(genLepLxy[gen2]);
    genlog10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[gen2])));
    if(abs(genLepEta[gen2])<1.479) {
      genmchgeneltrigebus_dE->Fill(el.E()-vecegus.E());
      genmchgeneltrigebus_dPt->Fill(genLepPt[gen2]-egusRecoPt[egus2]);
      genmchgeneltrigebus_dEta->Fill(genLepEta[gen2]-egusRecoEta[egus2]);
      genmchgeneltrigebus_qdPhi->Fill(qdPhi);
      genmchgeneltrigebus_dPromptE->Fill(promptel.E()-vecegus.E());
      genmchgeneltrigebus_dPromptPt->Fill(promptel.Pt()-egusRecoPt[egus2]);
      genmchgeneltrigebus_dPromptEta->Fill(genLepPromptEta[gen2]-egusRecoEta[egus2]);
      genmchgeneltrigebus_qdPromptPhi->Fill(qdPromptPhi);
    }
    else {
      genmchgeneltrigeeus_dE->Fill(el.E()-vecegus.E());
      genmchgeneltrigeeus_dPt->Fill(genLepPt[gen2]-egusRecoPt[egus2]);
      genmchgeneltrigeeus_dEta->Fill(genLepEta[gen2]-egusRecoEta[egus2]);
      genmchgeneltrigeeus_qdPhi->Fill(qdPhi);
    }
    gensubleadpt->Fill(genLepPt[gen2]);
    gensubleadeta->Fill(genLepEta[gen2]);
    gensubleadphi->Fill(genLepPhi[gen2]);
    gensubleadd0->Fill(genLepDxy[gen2]);    
    gensubleadlog10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[gen2])));    
    gensubleadlxy->Fill(genLepLxy[gen2]);    
    gensubleadlog10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[gen2])));    
  }
  genegmult->Fill(mchcnt);
  
  // Fill trigger egamma objects
  if(egus1!=-1 || egus2!=-1) {
    int leadidx=-1, subleadidx=-1;
    if(egus1==-1) leadidx = egus2;
    if(egus2==-1) leadidx = egus1;
    if(egus1!=-1 && egus2!=-1) {
      leadidx = egusRecoPt[egus1]>egusRecoPt[egus2]?egus1:egus2;
      subleadidx = egusRecoPt[egus1]>egusRecoPt[egus2]?egus2:egus1;
    }
    if(leadidx==-1) {
      throw "Error!!! Logical falacy in genmatched egus object histogram filling.";
      return;
    }
    
    leadegpt->Fill(egusRecoPt[leadidx]);
    leadegeta->Fill(egusRecoEta[leadidx]);
    leadegphi->Fill(egusRecoPhi[leadidx]);
    
    // Fill barrel variables
    if(TMath::Abs(egusRecoEta[leadidx])<=1.479) {
      recoeb_leadegclustershape->Fill(egushltEgammaClusterShape[leadidx]);
      recoeb_leadegin5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[leadidx]);
      recoeb_leadegin5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[leadidx]);
      recoeb_leadegsmin->Fill(egushltEgammaClusterShape_smin[leadidx]);
      recoeb_leadegsmaj->Fill(egushltEgammaClusterShape_smaj[leadidx]);
      //recoeb_leadegscenergy->Fill(egushltEgammaSuperClusterEnergy[leadidx]);
      recoeb_leadeghovere->Fill(egushltEgammaHoverE[leadidx]);
      recoeb_leadeghovereoversupcluse->Fill(egushltEgammaHoverE[leadidx]/egushltEgammaSuperClusterEnergy[leadidx]);
      recoeb_leadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[leadidx]);
      recoeb_leadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[leadidx]/egushltEgammaSuperClusterEnergy[leadidx]);
      recoeb_leadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[leadidx]);
      recoeb_leadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[leadidx]/egushltEgammaHoverE[leadidx]);
      if(egushltEcalSeedClusterTime[leadidx]!=0) recoeb_leadegseedclustime->Fill(egushltEcalSeedClusterTime[leadidx]);
      if(egushltEcalSeedClusterTime[leadidx]!=0) recoeb_egseedclustime->Fill(egushltEcalSeedClusterTime[leadidx]);
      if(egushltEgammaPixelMatchVars_s2[leadidx]<TMath::Power(10,36)) {
	recoeb_leadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[leadidx]);
	recoeb_egpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[leadidx]);
      }
      else {
	recoeb_leadegpixelmchvar_s2->Fill(-50);
	recoeb_egpixelmchvar_s2->Fill(-50);
      }
      if(egushltEgammaGsfTrackVars_Deta[leadidx]<999999) {
	recoeb_leadegtrkiso->Fill(egushltEgammaEleGsfTrackIso[leadidx]);
	recoeb_leadegchi2->Fill(egushltEgammaGsfTrackVars_Chi2[leadidx]);
	recoeb_leadegdeta->Fill(egushltEgammaGsfTrackVars_Deta[leadidx]);
	recoeb_leadegdetaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[leadidx]);
	recoeb_leadegdphi->Fill(egushltEgammaGsfTrackVars_Dphi[leadidx]);
	recoeb_leadegmhits->Fill(egushltEgammaGsfTrackVars_MissingHits[leadidx]);
	recoeb_leadegnlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[leadidx]);
	recoeb_leadegooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[leadidx]);
	recoeb_leadegooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[leadidx]);
	//recoeb_leadegvalhits->Fill(egushltEgammaGsfTrackVars_ValidHits[leadidx]);
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
      recoee_leadegclustershape->Fill(egushltEgammaClusterShape[leadidx]);
      recoee_leadegin5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[leadidx]);
      recoee_leadegin5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[leadidx]);
      recoee_leadegsmin->Fill(egushltEgammaClusterShape_smin[leadidx]);
      recoee_leadegsmaj->Fill(egushltEgammaClusterShape_smaj[leadidx]);
      //recoee_leadegscenergy->Fill(egushltEgammaSuperClusterEnergy[leadidx]);
      recoee_leadeghovere->Fill(egushltEgammaHoverE[leadidx]);
      recoee_leadeghovereoversupcluse->Fill(egushltEgammaHoverE[leadidx]/egushltEgammaSuperClusterEnergy[leadidx]);
      recoee_leadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[leadidx]);
      recoee_leadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[leadidx]/egushltEgammaSuperClusterEnergy[leadidx]);
      recoee_leadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[leadidx]);	
      recoee_leadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[leadidx]/egushltEgammaHoverE[leadidx]);	
      if(egushltEcalSeedClusterTime[leadidx]!=0) recoee_leadegseedclustime->Fill(egushltEcalSeedClusterTime[leadidx]);
      if(egushltEcalSeedClusterTime[leadidx]!=0) recoee_egseedclustime->Fill(egushltEcalSeedClusterTime[leadidx]);
      if(egushltEgammaPixelMatchVars_s2[leadidx]<TMath::Power(10,36)) {
	recoee_leadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[leadidx]);
	recoee_egpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[leadidx]);
      }
      else {
	recoee_leadegpixelmchvar_s2->Fill(-50);
	recoee_egpixelmchvar_s2->Fill(-50);
      }
      if(egushltEgammaGsfTrackVars_Deta[leadidx]<999999) {
	recoee_leadegtrkiso->Fill(egushltEgammaEleGsfTrackIso[leadidx]);
	recoee_leadegchi2->Fill(egushltEgammaGsfTrackVars_Chi2[leadidx]);
	recoee_leadegdeta->Fill(egushltEgammaGsfTrackVars_Deta[leadidx]);
	recoee_leadegdetaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[leadidx]);
	recoee_leadegdphi->Fill(egushltEgammaGsfTrackVars_Dphi[leadidx]);
	recoee_leadegmhits->Fill(egushltEgammaGsfTrackVars_MissingHits[leadidx]);
	recoee_leadegnlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[leadidx]);
	recoee_leadegooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[leadidx]);
	recoee_leadegooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[leadidx]);
	//recoee_leadegvalhits->Fill(egushltEgammaGsfTrackVars_ValidHits[leadidx]);
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
    
    if(subleadidx!=-1) { // Require second eg object
      TLorentzVector leadeg, subleadeg;
      leadeg.SetPtEtaPhiM(egusRecoPt[leadidx],egusRecoEta[leadidx],egusRecoPhi[leadidx],0.106);
      subleadeg.SetPtEtaPhiM(egusRecoPt[subleadidx],egusRecoEta[subleadidx],egusRecoPhi[subleadidx],0.106);

      subleadegpt->Fill(egusRecoPt[subleadidx]);
      subleadegeta->Fill(egusRecoEta[subleadidx]);
      subleadegphi->Fill(egusRecoPhi[subleadidx]);
      
      // Fill barrel variables
      if(TMath::Abs(egusRecoEta[subleadidx])<=1.479) {
	recoeb_subleadegclustershape->Fill(egushltEgammaClusterShape[subleadidx]);
	recoeb_subleadegin5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[subleadidx]);
	recoeb_subleadegin5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[subleadidx]);
	recoeb_subleadegsmin->Fill(egushltEgammaClusterShape_smin[subleadidx]);
	recoeb_subleadegsmaj->Fill(egushltEgammaClusterShape_smaj[subleadidx]);
	recoeb_subleadegscenergy->Fill(egushltEgammaSuperClusterEnergy[subleadidx]);
	recoeb_subleadeghovere->Fill(egushltEgammaHoverE[subleadidx]);
	recoeb_subleadeghovereoversupcluse->Fill(egushltEgammaHoverE[subleadidx]/egushltEgammaSuperClusterEnergy[subleadidx]);
	recoeb_subleadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[subleadidx]);
	recoeb_subleadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[subleadidx]/egushltEgammaSuperClusterEnergy[subleadidx]);
	recoeb_subleadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[subleadidx]);
	recoeb_subleadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[subleadidx]/egushltEgammaHoverE[subleadidx]);
	if(egushltEcalSeedClusterTime[subleadidx]!=0) recoeb_subleadegseedclustime->Fill(egushltEcalSeedClusterTime[subleadidx]);
	if(egushltEcalSeedClusterTime[subleadidx]!=0) recoeb_egseedclustime->Fill(egushltEcalSeedClusterTime[subleadidx]);
	if(egushltEgammaPixelMatchVars_s2[subleadidx]<TMath::Power(10,36)) {
	  recoeb_subleadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[subleadidx]);
	  recoeb_egpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[subleadidx]);
	}
	else {
	  recoeb_subleadegpixelmchvar_s2->Fill(-50);
	  recoeb_egpixelmchvar_s2->Fill(-50);
	}
	if(egushltEgammaGsfTrackVars_Deta[subleadidx]<999999) {
	  recoeb_subleadegtrkiso->Fill(egushltEgammaEleGsfTrackIso[subleadidx]);
	  recoeb_subleadegchi2->Fill(egushltEgammaGsfTrackVars_Chi2[subleadidx]);
	  recoeb_subleadegdeta->Fill(egushltEgammaGsfTrackVars_Deta[subleadidx]);
	  recoeb_subleadegdetaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[subleadidx]);
	  recoeb_subleadegdphi->Fill(egushltEgammaGsfTrackVars_Dphi[subleadidx]);
	  recoeb_subleadegmhits->Fill(egushltEgammaGsfTrackVars_MissingHits[subleadidx]);
	  recoeb_subleadegnlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[subleadidx]);
	  recoeb_subleadegooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[subleadidx]);
	  recoeb_subleadegooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[subleadidx]);
	  recoeb_subleadegvalhits->Fill(egushltEgammaGsfTrackVars_ValidHits[subleadidx]);
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
	recoee_subleadegclustershape->Fill(egushltEgammaClusterShape[subleadidx]);
	recoee_subleadegin5x5clusshape->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5[subleadidx]);
	recoee_subleadegin5x5noiseclnd->Fill(egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[subleadidx]);
	recoee_subleadegsmin->Fill(egushltEgammaClusterShape_smin[subleadidx]);
	recoee_subleadegsmaj->Fill(egushltEgammaClusterShape_smaj[subleadidx]);
	recoee_subleadegscenergy->Fill(egushltEgammaSuperClusterEnergy[subleadidx]);
	recoee_subleadeghovere->Fill(egushltEgammaHoverE[subleadidx]);
	recoee_subleadeghovereoversupcluse->Fill(egushltEgammaHoverE[subleadidx]/egushltEgammaSuperClusterEnergy[subleadidx]);
	recoee_subleadegecalpfclustiso->Fill(egushltEgammaEcalPFClusterIso[subleadidx]);
	recoee_subleadegecalpfclustisoovere->Fill(egushltEgammaEcalPFClusterIso[subleadidx]/egushltEgammaSuperClusterEnergy[subleadidx]);
	recoee_subleadeghcalpfclustiso->Fill(egushltEgammaHcalPFClusterIso[subleadidx]);	
	recoee_subleadeghcalpfclustisoovere->Fill(egushltEgammaHcalPFClusterIso[subleadidx]/egushltEgammaHoverE[subleadidx]);	
	if(egushltEcalSeedClusterTime[subleadidx]!=0) recoee_subleadegseedclustime->Fill(egushltEcalSeedClusterTime[subleadidx]);
	if(egushltEcalSeedClusterTime[subleadidx]!=0) recoee_egseedclustime->Fill(egushltEcalSeedClusterTime[subleadidx]);
	if(egushltEgammaPixelMatchVars_s2[subleadidx]<TMath::Power(10,36)) {
	  recoee_subleadegpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[subleadidx]);
	  recoee_egpixelmchvar_s2->Fill(egushltEgammaPixelMatchVars_s2[subleadidx]);
	}
	else {
	  recoee_subleadegpixelmchvar_s2->Fill(-50);
	  recoee_egpixelmchvar_s2->Fill(-50);
	}
	if(egushltEgammaGsfTrackVars_Deta[subleadidx]<999999) {
	  recoee_subleadegtrkiso->Fill(egushltEgammaEleGsfTrackIso[subleadidx]);
	  recoee_subleadegchi2->Fill(egushltEgammaGsfTrackVars_Chi2[subleadidx]);
	  recoee_subleadegdeta->Fill(egushltEgammaGsfTrackVars_Deta[subleadidx]);
	  recoee_subleadegdetaseed->Fill(egushltEgammaGsfTrackVars_DetaSeed[subleadidx]);
	  recoee_subleadegdphi->Fill(egushltEgammaGsfTrackVars_Dphi[subleadidx]);
	  recoee_subleadegmhits->Fill(egushltEgammaGsfTrackVars_MissingHits[subleadidx]);
	  recoee_subleadegnlayerit->Fill(egushltEgammaGsfTrackVars_NLayerIT[subleadidx]);
	  recoee_subleadegooeseedoop->Fill(egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[subleadidx]);
	  recoee_subleadegooesclsoop->Fill(egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[subleadidx]);
	  recoee_subleadegvalhits->Fill(egushltEgammaGsfTrackVars_ValidHits[subleadidx]);
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
      if(egusRecoEta[leadidx]<1.479 && egusRecoEta[subleadidx]<1.479) {
	recoeb_leadsubleadM->Fill((leadeg+subleadeg).M());
      }
      else {
	recoee_leadsubleadM->Fill((leadeg+subleadeg).M());
	}
    } // End of sublead condition
  } // End of filling condition
  
}

// Function to add a set of histograms for a gen particles
void data_robustanalyzer::addgenhist(TString selection) {
  
  all1dhists.push_back(new TH1F(selection+"geneg_egmult","gen N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"geneg_egmompid","gen mom pdg id",100,-50,50));
  all1dhists.push_back(new TH1F(selection+"geneg_pt","gen e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"geneg_eta","gen e/#gamma #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_phi","gen e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"geneg_prompteta","gen e/#gamma corrected #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_promptphi","gen e/#gamma corrected #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"geneg_vx","gen e/#gamma v_{x} / cm",50000,-250,250));
  all1dhists.push_back(new TH1F(selection+"geneg_vy","gen e/#gamma v_{y} / cm",50000,-250,250));
  all1dhists.push_back(new TH1F(selection+"geneg_vz","gen e/#gamma v_{z} / cm",50000,-250,250));
  all1dhists.push_back(new TH1F(selection+"geneg_pvx","mother v_{x} / cm",2000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"geneg_pvy","mother v_{y} / cm",2000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"geneg_pvz","mother v_{z} / cm",50000,-250,250));
  all1dhists.push_back(new TH1F(selection+"geneg_deltaetamom","#Delta#eta(gen, mom)",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_deltaphimom","#Delta#phi(gen, mom)",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_deltaRmom","#Delta#R(gen, mom)",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"geneg_d0","gen e/#gamma d_{0} / cm",20000,-100,100));
  all1dhists.push_back(new TH1F(selection+"geneg_log10d0","gen e/#gamma log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_lxy","gen e/#gamma l_{xy} / cm",20000,-10,190));
  all1dhists.push_back(new TH1F(selection+"geneg_log10lxy","gen e/#gamma log_{10}l_{xy} / log_{10}cm",600,-1,5));
  all1dhists.push_back(new TH1F(selection+"geneg_t1","gen e/#gamma time ecal",10000,-100,900));
  all1dhists.push_back(new TH1F(selection+"geneg_t0","gen e/#gamma time ecal (equivalent prompt)",10000,-100,900));
  all1dhists.push_back(new TH1F(selection+"geneg_t1mt0","gen e/#gamma time ecal (more than prompt)",10000,-10,90));
  all1dhists.push_back(new TH1F(selection+"geneg_leadpt","gen e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"geneg_leadeta","gen e/#gamma #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_leadphi","gen e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"geneg_leaddeltaetamom","#Delta#eta(gen, mom)",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_leaddeltaphimom","#Delta#phi(gen, mom)",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_leaddeltaRmom","#Delta#R(gen, mom)",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"geneg_leadd0","gen e/#gamma d_{0} / cm",20000,-100,100));
  all1dhists.push_back(new TH1F(selection+"geneg_leadlog10d0","gen e/#gamma log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_leadlxy","gen e/#gamma l_{xy} / cm",20000,-10,190));
  all1dhists.push_back(new TH1F(selection+"geneg_leadlog10lxy","gen e/#gamma log_{10}l_{xy} / log_{10}cm",600,-1,5));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadpt","gen e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadeta","gen e/#gamma #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadphi","gen e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"geneg_subleaddeltaetamom","#Delta#eta(gen, mom)",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_subleaddeltaphimom","#Delta#phi(gen, mom)",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_subleaddeltaRmom","#Delta#R(gen, mom)",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadd0","gen e/#gamma d_{0} / cm",20000,-100,100));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadlog10d0","gen e/#gamma log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadlxy","gen e/#gamma l_{xy} / cm",20000,-10,190));
  all1dhists.push_back(new TH1F(selection+"geneg_subleadlog10lxy","gen e/#gamma log_{10}l_{xy} / log_{10}cm",600,-1,5));
}

// Function to add a set of histograms for a selection
void data_robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"recoeg_egmult","reco N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"recoeg_leadegpt","reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recoeg_leadegeta","reco e/#gamma #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"recoeg_leadegphi","reco e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"recoeg_subleadegpt","reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recoeg_subleadegeta","reco e/#gamma #eta",100,-5,5));
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
  all1dhists.push_back(new TH1F(selection+"recoegus_leadegeta","reco e/#gamma #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"recoegus_leadegphi","reco e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"recoegus_subleadegpt","reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recoegus_subleadegeta","reco e/#gamma #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"recoegus_subleadegphi","reco e/#gamma #phi",66,-3.3,3.3));

  // Both lead and sub lead electrons
  all1dhists.push_back(new TH1F(selection+"recoebus_egseedclustime","barrel e/#gamma_{seed} time / ns",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"recoeeus_egseedclustime","end-cap e/#gamma_{seed} time / ns",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"recoebus_egpixelmchvar_s2","barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeeus_egpixelmchvar_s2","end-cap e/#gamma pixelmachvar",1000,-50,950));

  // barrel variables - lead pT unseeded e/gamma
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegclustershape","barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegin5x5clusshape","barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegin5x5noiseclnd","barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegsmin","barrel e/#gamma smin",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"recoebus_leadegsmaj","barrel e/#gamma smaj",1000,0,1));
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
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegsmin","barrel e/#gamma smin",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"recoebus_subleadegsmaj","barrel e/#gamma smaj",1000,0,1));
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
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegsmin","end-cap e/#gamma smin",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_leadegsmaj","end-cap e/#gamma smaj",1000,0,1));
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
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegsmin","end-cap e/#gamma smin",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"recoeeus_subleadegsmaj","end-cap e/#gamma smaj",1000,0,1));
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

// Function to add a set of histograms for gen matched objects - unseeded egamma objects
void data_robustanalyzer::addhistgenmchunseeded(TString selection) {

  // Variables before gen match
  all1dhists.push_back(new TH1F(selection+"geneltrigebus_dE","#Delta E(gen e, trig. e/#gamma)",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"geneltrigebus_dPt","#Delta p_{T}(gen e, trig. e/#gamma)",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"geneltrigebus_dEta","#Delta#eta(gen e, trig. e/#gamma)",8000,-4,4));
  all1dhists.push_back(new TH1F(selection+"geneltrigebus_qdPhi","q#Delta#phi(gen e, trig. e/#gamma)",8000,-4,4));
  all1dhists.push_back(new TH1F(selection+"geneltrigebus_dPromptEta","#Delta#eta(gen prompt e, trig. e/#gamma)",8000,-4,4));
  all1dhists.push_back(new TH1F(selection+"geneltrigebus_qdPromptPhi","q#Delta#phi(gen prompt e, trig. e/#gamma)",8000,-4,4));
  all1dhists.push_back(new TH1F(selection+"geneltrigeeus_dE","#Delta E(gen e, trig. e/#gamma)",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"geneltrigeeus_dPt","#Delta p_{T}(gen e, trig. e/#gamma)",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"geneltrigeeus_dEta","#Delta#eta(gen e, trig. e/#gamma)",8000,-4,4));
  all1dhists.push_back(new TH1F(selection+"geneltrigeeus_qdPhi","q#Delta#phi(gen e, trig. e/#gamma)",8000,-4,4));
  
  // Variables after gen match
  all1dhists.push_back(new TH1F(selection+"recomchgenel_egmult","gen N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_pt","gen e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_eta","gen e/#gamma #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_phi","gen e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_d0","gen e/#gamma d_{0} / cm",20000,-100,100));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_log10d0","gen e/#gamma log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_lxy","gen e/#gamma l_{xy} / cm",20000,-10,190));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_log10lxy","gen e/#gamma log_{10}l_{xy} / log_{10}cm",600,-1,5));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_leadpt","gen e/#gamma_{1} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_leadeta","gen e/#gamma_{1} #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_leadphi","gen e/#gamma_{1} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_leadd0","gen e/#gamma_{1} d_{0} / cm",20000,-100,100));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_leadlog10d0","gen e/#gamma_{1} log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_leadlxy","gen e/#gamma_{1} l_{xy} / cm",20000,-10,190));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_leadlog10lxy","gen e/#gamma_{1} log_{10}l_{xy} / log_{10}cm",600,-1,5));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_subleadpt","gen e/#gamma_{2} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_subleadeta","gen e/#gamma_{2} #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_subleadphi","gen e/#gamma_{2} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_subleadd0","gen e/#gamma_{2} d_{0} / cm",20000,-100,100));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_subleadlog10d0","gen e/#gamma_{2} log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_subleadlxy","gen e/#gamma_{2} l_{xy} / cm",20000,-10,190));
  all1dhists.push_back(new TH1F(selection+"recomchgenel_subleadlog10lxy","gen e/#gamma_{2} log_{10}l_{xy} / log_{10}cm",600,-1,5));

  // Variables for gen matched trigger electron
  all1dhists.push_back(new TH1F(selection+"genmchrecoegus_leadegpt","gen matched reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"genmchrecoegus_leadegeta","gen matched reco e/#gamma #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"genmchrecoegus_leadegphi","gen matched reco e/#gamma #phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F(selection+"genmchrecoegus_subleadegpt","gen matched reco e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"genmchrecoegus_subleadegeta","gen matched reco e/#gamma #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"genmchrecoegus_subleadegphi","gen matched reco e/#gamma #phi",66,-3.3,3.3));

  // Both lead and sub lead electrons
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_egseedclustime","gen matched barrel e/#gamma_{seed} time / ns",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_egseedclustime","gen matched end-cap e/#gamma_{seed} time / ns",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_egpixelmchvar_s2","gen matched barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_egpixelmchvar_s2","gen matched end-cap e/#gamma pixelmachvar",1000,-50,950));

  // barrel variables - lead pT unseeded e/gamma
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegclustershape","gen matched barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegin5x5clusshape","gen matched barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegin5x5noiseclnd","gen matched barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegsmin","gen matched barrel e/#gamma smin",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegsmaj","gen matched barrel e/#gamma smaj",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadeghovere","gen matched barrel e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegscenergy","gen matched barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegecalpfclustiso","gen matched barrel e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegecalpfclustisoovere","gen matched barrel e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadeghcalpfclustiso","gen matched barrel e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadeghcalpfclustisoovere","gen matched barrel e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegpixelmchvar_s2","gen matched barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegchi2","gen matched barrel e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegdeta","gen matched barrel e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegdetaseed","gen matched barrel e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegdphi","gen matched barrel e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegmhits","gen matched barrel e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegnlayerit","gen matched barrel e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegooeseedoop","gen matched barrel e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegooesclsoop","gen matched barrel e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegvalhits","gen matched barrel e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegtrkiso","gen matched barrel e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadeghovereoversupcluse","gen matched barrel e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadegseedclustime","gen matched barrel e/#gamma_{seed} time / ns",20000,-10,10));

  // barrel variables - sub-lead pT unseeded e/gamma
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegclustershape","gen matched barrel e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegin5x5clusshape","gen matched barrel e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegin5x5noiseclnd","gen matched barrel e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegsmin","gen matched barrel e/#gamma smin",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegsmaj","gen matched barrel e/#gamma smaj",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadeghovere","gen matched barrel e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegscenergy","gen matched barrel e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegecalpfclustiso","gen matched barrel e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegecalpfclustisoovere","gen matched barrel e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadeghcalpfclustiso","gen matched barrel e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadeghcalpfclustisoovere","gen matched barrel e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegpixelmchvar_s2","gen matched barrel e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegchi2","gen matched barrel e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegdeta","gen matched barrel e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegdetaseed","gen matched barrel e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegdphi","gen matched barrel e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegmhits","gen matched barrel e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegnlayerit","gen matched barrel e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegooeseedoop","gen matched barrel e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegooesclsoop","gen matched barrel e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegvalhits","gen matched barrel e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegtrkiso","gen matched barrel e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadeghovereoversupcluse","gen matched barrel e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_subleadegseedclustime","gen matched barrel e/#gamma_{seed} time / ns",20000,-10,10));

  // barrel variables - invariant mass unseeded
  all1dhists.push_back(new TH1F(selection+"genmchrecoebus_leadsubleadM","gen matched barrel M(e/#gamma_{1},e/#gamma_{2}) / GeV",500,0,500));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigebus_dE","gen matched #Delta E(gen e, trig. e/#gamma)",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigebus_dPt","gen matched #Delta p_{T}(gen e, trig. e/#gamma)",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigebus_dEta","gen matched #Delta#eta(gen e, trig. e/#gamma)",8000,-4,4));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigebus_qdPhi","gen matched q#Delta#phi(gen e, trig. e/#gamma)",8000,-4,4));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigebus_dPromptE","gen matched #Delta E(gen prompt e, trig. e/#gamma)",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigebus_dPromptPt","gen matched #Delta p_{T}(gen prompt e, trig. e/#gamma)",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigebus_dPromptEta","gen matched #Delta#eta(gen prompt e, trig. e/#gamma)",8000,-4,4));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigebus_qdPromptPhi","gen matched q#Delta#phi(gen prompt e, trig. e/#gamma)",8000,-4,4));

  // end-cap variables - lead pT unseeded e/gamma
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegclustershape","gen matched end-cap e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegin5x5clusshape","gen matched end-cap e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegin5x5noiseclnd","gen matched end-cap e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegsmin","gen matched end-cap e/#gamma smin",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegsmaj","gen matched end-cap e/#gamma smaj",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadeghovere","gen matched end-cap e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegscenergy","gen matched end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegecalpfclustiso","gen matched end-cap e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegecalpfclustisoovere","gen matched end-cap e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadeghcalpfclustiso","gen matched end-cap e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadeghcalpfclustisoovere","gen matched end-cap e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegpixelmchvar_s2","gen matched end-cap e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegchi2","gen matched end-cap e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegdeta","gen matched end-cap e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegdetaseed","gen matched end-cap e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegdphi","gen matched end-cap e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegmhits","gen matched end-cap e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegnlayerit","gen matched end-cap e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegooeseedoop","gen matched end-cap e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegooesclsoop","gen matched end-cap e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegvalhits","gen matched end-cap e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegtrkiso","gen matched end-cap e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadeghovereoversupcluse","gen matched end-cap e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadegseedclustime","gen matched end-cap e/#gamma_{seed} time / ns",20000,-10,10));

  // end-cap variables - sub-lead pT unseeded e/gamma
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegclustershape","gen matched end-cap e/#gamma clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegin5x5clusshape","gen matched end-cap e/#gamma #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegin5x5noiseclnd","gen matched end-cap e/#gamma #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegsmin","gen matched end-cap e/#gamma smin",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegsmaj","gen matched end-cap e/#gamma smaj",1000,0,1));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadeghovere","gen matched end-cap e/#gamma H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegscenergy","gen matched end-cap e/#gamma SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegecalpfclustiso","gen matched end-cap e/#gamma ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegecalpfclustisoovere","gen matched end-cap e/#gamma ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadeghcalpfclustiso","gen matched end-cap e/#gamma hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadeghcalpfclustisoovere","gen matched end-cap e/#gamma hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegpixelmchvar_s2","gen matched end-cap e/#gamma pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegchi2","gen matched end-cap e/#gamma #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegdeta","gen matched end-cap e/#gamma #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegdetaseed","gen matched end-cap e/#gamma #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegdphi","gen matched end-cap e/#gamma #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegmhits","gen matched end-cap e/#gamma missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegnlayerit","gen matched end-cap e/#gamma num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegooeseedoop","gen matched end-cap e/#gamma 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegooesclsoop","gen matched end-cap e/#gamma 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegvalhits","gen matched end-cap e/#gamma valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegtrkiso","gen matched end-cap e/#gamma track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadeghovereoversupcluse","gen matched end-cap e/#gamma H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_subleadegseedclustime","gen matched end-cap e/#gamma_{seed} time / ns",20000,-10,10));

  // end-cap variables - invariant mass unseeded
  all1dhists.push_back(new TH1F(selection+"genmchrecoeeus_leadsubleadM","gen matched end-cap M(e/#gamma_{1},e/#gamma_{2}) / GeV",500,0,500));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigeeus_dE","gen matched #Delta E(gen e, trig. e/#gamma)",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigeeus_dPt","gen matched #Delta p_{T}(gen e, trig. e/#gamma)",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigeeus_dEta","gen matched #Delta#eta(gen e, trig. e/#gamma)",8000,-4,4));
  all1dhists.push_back(new TH1F(selection+"genmchgeneltrigeeus_qdPhi","gen matched q#Delta#phi(gen e, trig. e/#gamma)",8000,-4,4));
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

    TVector3 eg, l1;
    eg.SetPtEtaPhi(egRecoPt[idx], egRecoEta[idx], egRecoPhi[idx]);
    l1.SetPtEtaPhi(l1FiltPt[l1cnt], l1FiltEta[l1cnt], l1FiltPhi[l1cnt]);
    double deltaphi = std::abs(eg.DeltaPhi(l1));
    if(egRecoEta[idx]<etabinhigh && egRecoEta[idx]>etabinlow && deltaphi<0.5*phibinsize) {
      l1matchdecision = true;
    }
  }

  return l1matchdecision;
}

// Function to compare a set of cuts with filter output for comaptibility
bool data_robustanalyzer::comparecutonobjtofilt(double filt[], unsigned int filtsize, double cutobj[], vector<int> cutobjidx){
  
  //if(cutobjidx.size()!=filtsize) return false;
  //cout<<filtsize<<"\t"<<cutobjidx.size()<<endl;
  for(unsigned int filtidx=0; filtidx<filtsize; filtidx++) {
    bool foundfilt = false;
    for(int obj : cutobjidx){
      if(filt[filtidx] == cutobj[obj]) {
	foundfilt = true;
	//cout<<"Found match for: "<<filt[filtidx]<<" with obj: "<<cutobj[obj]<<endl;
      }
    }
    if(!foundfilt) return false;
  }

  return true;
}

/*
      cout<<"L1 Eg obj: ";
      for(unsigned int ctr=0; ctr<l1egObjN; ctr++) {
	cout<<l1egObjPt[ctr]<<"\t";
      }
      cout<<endl;
      cout<<"L1 Filter: ";
      for(unsigned int ctr=0; ctr<l1FiltN; ctr++) {
	cout<<l1FiltPt[ctr]<<"\t"<<l1FiltEta[ctr]<<"\t"<<l1FiltPhi[ctr]<<endl;
      }
      cout<<"Reco before cut: ";
      for(unsigned int eg=0; eg<egRecoN; eg++) {
	cout<<egRecoPt[eg]<<"\t"<<egRecoEta[eg]<<"\t"<<egRecoPhi[eg]<<"\t"<<eghltEgammaHoverE[eg]/egushltEgammaSuperClusterEnergy[eg]<<"\t"<<isL1EgSeeded(eg)<<endl;
      }
      cout<<"Reco: ";
      for(unsigned int eg : dieg70idegidx) {
	cout<<egRecoPt[eg]<<"\t";
      }
      cout<<endl;
      cout<<"UnseededFilter: ";
      for(unsigned int ctr=0; ctr<dieg70HeusFiltN; ctr++) {
	cout<<dieg33CsusFiltPt[ctr]<<"\t";
      }
      cout<<endl;
      cout<<"Unseeded Reco: ";
      for(unsigned int eg : dieg70idegusidx) {
	cout<<egusRecoPt[eg]<<"\t";
      }
      cout<<endl;
*/
