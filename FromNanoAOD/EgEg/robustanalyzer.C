/*
 * AUTHOR: Abanti Ranadhir Sahasransu - asahasra@cern.ch
 * The code was made to read NanoAOD data for preliminary trigger efficiency created by the EXO group
 * for fast processing the new Run 3 data.
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
robustanalyzer::robustanalyzer(TString filename, TString outfilename, bool issimu) {

  inputChain = new TChain("Events");
  cout<<"Initializing for file: "<<filename<<endl;
  inputChain->Add(filename);

  inputChain->SetBranchAddress("run", &run);
  inputChain->SetBranchAddress("luminosityBlock", &lumi);
  inputChain->SetBranchAddress("event", &event);
  inputChain->SetBranchAddress("Rho_fixedGridRhoFastjetAll", &eventRho);  
  inputChain->SetBranchAddress("bunchCrossing", &bunch);

  inputChain->SetBranchAddress("HLT_DiPhoton10Time1p4ns", &HLT_DiPhoton10Time1p4ns);
  inputChain->SetBranchAddress("HLT_DiPhoton10sminlt0p12", &HLT_DiPhoton10sminlt0p12);
  inputChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight);
  inputChain->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF);
  inputChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
  inputChain->SetBranchAddress("HLT_CaloMET80_NotCleaned", &HLT_CaloMET80_NotCleaned);
  inputChain->SetBranchAddress("HLT_PFMET200_NotCleaned", &HLT_PFMET200_NotCleaned);
  inputChain->SetBranchAddress("HLT_PFMET200_BeamHaloCleaned", &HLT_PFMET200_BeamHaloCleaned);
  inputChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);

  inputChain->SetBranchAddress("nElectron", &eln);
  inputChain->SetBranchAddress("Electron_deltaEtaSC", &eldetasc);
  inputChain->SetBranchAddress("Electron_dxy", &eldxy);
  inputChain->SetBranchAddress("Electron_dxyErr", &eldxyerr);
  inputChain->SetBranchAddress("Electron_dz", &eldz);
  inputChain->SetBranchAddress("Electron_dzErr", &eldzerr);
  inputChain->SetBranchAddress("Electron_eInvMinusPInv", &elooemoop);
  //inputChain->SetBranchAddress("Electron_energyErr", &elenergyerr);
  inputChain->SetBranchAddress("Electron_eta", &eleta);
  inputChain->SetBranchAddress("Electron_hoe", &elhoe);
  inputChain->SetBranchAddress("Electron_phi", &elphi);
  inputChain->SetBranchAddress("Electron_pt", &elpt);
  //inputChain->SetBranchAddress("Electron_r9", &elr9);
  //inputChain->SetBranchAddress("Electron_scEtOverPt", &elscetoverpt);
  inputChain->SetBranchAddress("Electron_sieie", &elsieie);
  inputChain->SetBranchAddress("Electron_convVeto", &elconvveto);
  //inputChain->SetBranchAddress("Electron_charge", &elcharge);
  //inputChain->SetBranchAddress("Electron_cutBased", &elcutbasedid);
  //inputChain->SetBranchAddress("Electron_pdgId", &elpdgid);
  //inputChain->SetBranchAddress("Electron_photonIdx", &elphotonIdx);

  inputChain->SetBranchAddress("nLowPtElectron", &lowpteln);
  //inputChain->SetBranchAddress("LowPtElectron_ID", &lowptelid);
  //inputChain->SetBranchAddress("LowPtElectron_deltaEtaSC", &lowpteldetasc);
  inputChain->SetBranchAddress("LowPtElectron_dxy", &lowpteldxy);
  inputChain->SetBranchAddress("LowPtElectron_dxyErr", &lowpteldxyerr);
  inputChain->SetBranchAddress("LowPtElectron_dz", &lowpteldz);
  inputChain->SetBranchAddress("LowPtElectron_dzErr", &lowpteldzerr);
  inputChain->SetBranchAddress("LowPtElectron_eInvMinusPInv", &lowptelooemoop);
  //inputChain->SetBranchAddress("LowPtElectron_energyErr", &lowptelenergyerr);
  inputChain->SetBranchAddress("LowPtElectron_eta", &lowpteleta);
  inputChain->SetBranchAddress("LowPtElectron_hoe", &lowptelhoe);
  inputChain->SetBranchAddress("LowPtElectron_phi", &lowptelphi);
  inputChain->SetBranchAddress("LowPtElectron_pt", &lowptelpt);
  //inputChain->SetBranchAddress("LowPtElectron_ptbiased", &lowptelptbias);
  //inputChain->SetBranchAddress("LowPtElectron_scEtOverPt", &lowptelscetoverpt);
  inputChain->SetBranchAddress("LowPtElectron_sieie", &lowptelsieie);
  //inputChain->SetBranchAddress("LowPtElectron_unbiased", &lowptelunbias);
  //inputChain->SetBranchAddress("LowPtElectron_charge", &lowptelcharge);
  //inputChain->SetBranchAddress("LowPtElectron_electronIdx", &lowptelid);
  //inputChain->SetBranchAddress("LowPtElectron_pdgId", &lowptelpdgid);

  inputChain->SetBranchAddress("nPhoton", &phn);
  inputChain->SetBranchAddress("Photon_pt", &phopt);
  inputChain->SetBranchAddress("Photon_eta", &phoeta);
  inputChain->SetBranchAddress("Photon_phi", &phophi);
  inputChain->SetBranchAddress("Photon_sieie", &phosieie);
  inputChain->SetBranchAddress("Photon_hoe", &phohoe);
  
  outfile = new TFile(outfilename,"RECREATE");

}

// Fill the root file, close the root file, and handle deletions
robustanalyzer::~robustanalyzer() {

  inputChain->Delete();
  outfile->Write();
  outfile->Close();
}

// Analyzer for a single file
void robustanalyzer::analyzersinglefile(int splitCnt, int numCores) { // Assume splitCnt to range from 0 to nCores

  int totEntries = inputChain->GetEntries();
  cout<<"Total number of entries: "<<totEntries<<endl;
  int nCores = numCores-1; // Assume parallel processing over numCores cores where
  // there is a lesser no.of events in the last core
  int beginevent = splitCnt*(totEntries/nCores);
  int endevent = (splitCnt+1)*(totEntries/nCores);
  endevent = endevent<totEntries?endevent:totEntries; // Verfied that this logic to parallelize works
  cout<<"Runing "<<splitCnt<<" process from "<<beginevent<<" to "<<endevent<<endl; 

  // Count events passing certain selections
  int nosel=0;

  addhist("nosel_el");
  addhist("bar_el");
  addhist("gt2_bar_el");
  addhist("id1_gt2_bar_el");
  addhist("id2_gt2_bar_el");
  addhist("mid_gt2_bar_el");
  addhist("met_mid_gt2_bar_el");
  addhist("t1p4_mid_gt2_bar_el");
  addhist("sm12_mid_gt2_bar_el");
  addhist("t1p4_met_mid_gt2_bar_el");
  addhist("sm12_met_mid_gt2_bar_el");
  addhist("ec_el");
  addhist("gt2_ec_el");
  addhist("id1_gt2_ec_el");
  addhist("id2_gt2_ec_el");
  addhist("mid_gt2_ec_el");
  addhist("met_mid_gt2_ec_el");
  addhist("t1p4_mid_gt2_ec_el");
  addhist("sm12_mid_gt2_ec_el");
  addhist("t1p4_met_mid_gt2_ec_el");
  addhist("sm12_met_mid_gt2_ec_el");

  addhist4LowPtElectron("nosel_lowptel");
  addhist4LowPtElectron("bar_lowptel");
  addhist4LowPtElectron("gt2_bar_lowptel");
  addhist4LowPtElectron("id1_gt2_bar_lowptel");
  addhist4LowPtElectron("id2_gt2_bar_lowptel");
  addhist4LowPtElectron("mid_gt2_bar_lowptel");
  addhist4LowPtElectron("met_mid_gt2_bar_lowptel");
  addhist4LowPtElectron("t1p4_mid_gt2_bar_lowptel");
  addhist4LowPtElectron("sm12_mid_gt2_bar_lowptel");
  addhist4LowPtElectron("t1p4_met_mid_gt2_bar_lowptel");
  addhist4LowPtElectron("sm12_met_mid_gt2_bar_lowptel");
  addhist4LowPtElectron("ec_lowptel");
  addhist4LowPtElectron("gt2_ec_lowptel");
  addhist4LowPtElectron("id1_gt2_ec_lowptel");
  addhist4LowPtElectron("id2_gt2_ec_lowptel");
  addhist4LowPtElectron("mid_gt2_ec_lowptel");
  addhist4LowPtElectron("met_mid_gt2_ec_lowptel");
  addhist4LowPtElectron("t1p4_mid_gt2_ec_lowptel");
  addhist4LowPtElectron("sm12_mid_gt2_ec_lowptel");
  addhist4LowPtElectron("t1p4_met_mid_gt2_ec_lowptel");
  addhist4LowPtElectron("sm12_met_mid_gt2_ec_lowptel");

  addhist4Photon("nosel_pho");
  addhist4Photon("bar_pho");
  addhist4Photon("gt2_bar_pho");
  addhist4Photon("id1_gt2_bar_pho");
  addhist4Photon("mid_gt2_bar_pho");
  addhist4Photon("met_mid_gt2_bar_pho");
  addhist4Photon("t1p4_mid_gt2_bar_pho");
  addhist4Photon("sm12_mid_gt2_bar_pho");
  addhist4Photon("t1p4_met_mid_gt2_bar_pho");
  addhist4Photon("sm12_met_mid_gt2_bar_pho");
  addhist4Photon("ec_pho");
  addhist4Photon("gt2_ec_pho");
  addhist4Photon("id1_gt2_ec_pho");
  addhist4Photon("mid_gt2_ec_pho");
  addhist4Photon("met_mid_gt2_ec_pho");
  addhist4Photon("t1p4_mid_gt2_ec_pho");
  addhist4Photon("sm12_mid_gt2_ec_pho");
  addhist4Photon("t1p4_met_mid_gt2_ec_pho");
  addhist4Photon("sm12_met_mid_gt2_ec_pho");
  
  vector<int> noselelidx;
  vector<int> barelidx;
  vector<int> id1barelidx;
  vector<int> id2barelidx;
  vector<int> midbarelidx;
  vector<int> ecelidx;
  vector<int> id1ecelidx;
  vector<int> id2ecelidx;
  vector<int> midecelidx;

  vector<int> nosellowptelidx;
  vector<int> barlowptelidx;
  vector<int> id1barlowptelidx;
  vector<int> id2barlowptelidx;
  vector<int> midbarlowptelidx;
  vector<int> eclowptelidx;
  vector<int> id1eclowptelidx;
  vector<int> id2eclowptelidx;
  vector<int> mideclowptelidx;

  vector<int> noselphoidx;
  vector<int> barphoidx;
  vector<int> id1barphoidx;
  vector<int> midbarphoidx;
  vector<int> ecphoidx;
  vector<int> id1ecphoidx;
  vector<int> midecphoidx;

  // Loop beginning on events
  for(unsigned int event=beginevent; event<endevent; event++) {

    inputChain->GetEntry(event);
    //if(event>1000) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    // Section for trigger conditions
    Bool_t mettrigs = (HLT_PFMET120_PFMHT120_IDTight ||
		       HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF ||
		       HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
		       HLT_CaloMET80_NotCleaned ||
		       HLT_PFMET200_NotCleaned ||
		       HLT_PFMET200_BeamHaloCleaned ||
		       HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);
    if((HLT_PFMET120_PFMHT120_IDTight || HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_CaloMET80_NotCleaned || HLT_PFMET200_NotCleaned || HLT_PFMET200_BeamHaloCleaned || HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight) && (!mettrigs)) throw "Logical error!! Type 1";
    if(((!HLT_PFMET120_PFMHT120_IDTight) && (!HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF) && (!HLT_PFMETNoMu120_PFMHTNoMu120_IDTight) && (!HLT_CaloMET80_NotCleaned) && (!HLT_PFMET200_NotCleaned) && (!HLT_PFMET200_BeamHaloCleaned) && (!HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight)) && (mettrigs)) throw "Logical error!! Type 2";
    Bool_t t1p4nstrig = HLT_DiPhoton10Time1p4ns;
    Bool_t sminlt0p12trig = HLT_DiPhoton10sminlt0p12;
    
    // Loop beginning on electrons
    for(unsigned int idx=0; idx<eln; idx++) {

      if(idx>0) {
	if(elpt[idx]>elpt[idx-1])                                               // Check for pt sort
	  throw "Error!!! Electrons are not pt sorted";
      }

      TLorentzVector electron;
      electron.SetPtEtaPhiM(elpt[idx], eleta[idx], elphi[idx], 0.0005);
      double energy = electron.Energy();
      
      noselelidx.push_back(idx);

      bool barelsel = true;
      barelsel *= TMath::Abs(eleta[idx])<1.479;
      if(barelsel) barelidx.push_back(idx);

      bool id1barsel = true;
      id1barsel *= TMath::Abs(eleta[idx]) < 1.479;
      id1barsel *= elsieie[idx] < 0.0103;
      if(id1barsel) id1barelidx.push_back(idx);
      
      bool id2barsel = true;
      id2barsel *= TMath::Abs(eleta[idx]) < 1.479;
      id2barsel *= elsieie[idx] < 0.0103;
      id2barsel *= elhoe[idx] < (0.0241+1.28/energy+0.042*eventRho/energy);
      if(id2barsel) id2barelidx.push_back(idx);
      
      bool midbarsel = true;
      midbarsel *= TMath::Abs(eleta[idx]) < 1.479;
      midbarsel *= elsieie[idx] < 0.0103;
      //midbarsel *= eldetasc[idx] < 0.00481;
      midbarsel *= elhoe[idx] < (0.0241+1.28/energy+0.042*eventRho/energy);
      midbarsel *= elooemoop[idx] < 0.0966;
      //midbarsel *= elconvveto[idx];
      if(midbarsel) midbarelidx.push_back(idx);
      
      bool ecelsel = true;
      ecelsel *= TMath::Abs(eleta[idx]) > 1.479;
      if(ecelsel) ecelidx.push_back(idx);

      bool id1ecsel = true;
      id1ecsel *= TMath::Abs(eleta[idx]) > 1.479;
      id1ecsel *= elsieie[idx] < 0.0278;
      if(id1ecsel) id1ecelidx.push_back(idx);
      
      bool id2ecsel = true;
      id2ecsel *= TMath::Abs(eleta[idx]) > 1.479;
      id2ecsel *= elsieie[idx] < 0.0278;
      id2ecsel *= elhoe[idx] < (0.0274+2.08/energy+0.292*eventRho/energy);
      if(id2ecsel) id2ecelidx.push_back(idx);
      
      bool midecsel = true;
      midecsel *= TMath::Abs(eleta[idx]) > 1.479;
      midecsel *= elsieie[idx] < 0.0278;
      //midecsel *= eldetasc[idx] < 0.00847;
      midecsel *= elhoe[idx] < (0.0274+2.08/energy+0.292*eventRho/energy);
      midecsel *= elooemoop[idx] < 0.0769;
      //midecsel *= elconvveto[idx];
      if(midecsel) midecelidx.push_back(idx);
      
    } // End of loop on electrons
           
    // Loop beginning on low pt electrons
    for(unsigned int idx=0; idx<lowpteln; idx++) {

      if(idx>0) {
	if(lowptelpt[idx]>lowptelpt[idx-1])                                               // Check for pt sort
	  throw "Error!!! Electrons are not pt sorted";
      }

      TLorentzVector electron;
      electron.SetPtEtaPhiM(lowptelpt[idx], lowpteleta[idx], lowptelphi[idx], 0.0005);
      double energy = electron.Energy();
      
      nosellowptelidx.push_back(idx);
      
      bool barlowptelsel = true;
      barlowptelsel *= TMath::Abs(lowpteleta[idx])<1.479;
      if(barlowptelsel) barlowptelidx.push_back(idx);

      bool id1barsel = true;
      id1barsel *= TMath::Abs(lowpteleta[idx]) < 1.479;
      id1barsel *= lowptelsieie[idx] < 0.0103;
      if(id1barsel) id1barlowptelidx.push_back(idx);
      
      bool id2barsel = true;
      id2barsel *= TMath::Abs(lowpteleta[idx]) < 1.479;
      id2barsel *= lowptelsieie[idx] < 0.0103;
      id2barsel *= lowptelhoe[idx] < (0.0241+1.28/energy+0.042*eventRho/energy);
      if(id2barsel) id2barlowptelidx.push_back(idx);
      
      bool midbarsel = true;
      midbarsel *= TMath::Abs(lowpteleta[idx]) < 1.479;
      midbarsel *= lowptelsieie[idx] < 0.0103;
      //midbarsel *= lowpteldetasc[idx] < 0.00481;
      midbarsel *= lowptelhoe[idx] < (0.0241+1.28/energy+0.042*eventRho/energy);
      midbarsel *= lowptelooemoop[idx] < 0.0966;
      //midbarsel *= elconvveto[idx];
      if(midbarsel) midbarlowptelidx.push_back(idx);
      
      bool eclowptelsel = true;
      eclowptelsel *= TMath::Abs(lowpteleta[idx]) > 1.479;
      if(eclowptelsel) eclowptelidx.push_back(idx);

      bool id1ecsel = true;
      id1ecsel *= TMath::Abs(lowpteleta[idx]) > 1.479;
      id1ecsel *= lowptelsieie[idx] < 0.0278;
      if(id1ecsel) id1eclowptelidx.push_back(idx);
      
      bool id2ecsel = true;
      id2ecsel *= TMath::Abs(lowpteleta[idx]) > 1.479;
      id2ecsel *= lowptelsieie[idx] < 0.0278;
      id2ecsel *= lowptelhoe[idx] < (0.0274+2.08/energy+0.292*eventRho/energy);
      if(id2ecsel) id2eclowptelidx.push_back(idx);
      
      bool midecsel = true;
      midecsel *= TMath::Abs(lowpteleta[idx]) > 1.479;
      midecsel *= lowptelsieie[idx] < 0.0278;
      //midecsel *= lowpteldetasc[idx] < 0.00847;
      midecsel *= lowptelhoe[idx] < (0.0274+2.08/energy+0.292*eventRho/energy);
      midecsel *= lowptelooemoop[idx] < 0.0769;
      //midecsel *= lowptelconvveto[idx];
      if(midecsel) mideclowptelidx.push_back(idx);
      
    } // End of loop on low pt electrons
       
    // Loop beginning on photons
    for(unsigned int idx=0; idx<phn; idx++) {

      if(idx>0) {
	if(phopt[idx]>phopt[idx-1])                                               // Check for pt sort
	  throw "Error!!! Electrons are not pt sorted";
      }

      TLorentzVector electron;
      electron.SetPtEtaPhiM(phopt[idx], phoeta[idx], phophi[idx], 0.0005);
      double energy = electron.Energy();
      
      noselphoidx.push_back(idx);
      
      bool barphosel = true;
      barphosel *= TMath::Abs(phoeta[idx])<1.479;
      if(barphosel) barphoidx.push_back(idx);

      bool id1barsel = true;
      id1barsel *= TMath::Abs(phoeta[idx]) < 1.479;
      id1barsel *= phosieie[idx] < 0.0103;
      if(id1barsel) id1barphoidx.push_back(idx);
      
      bool midbarsel = true;
      midbarsel *= TMath::Abs(phoeta[idx]) < 1.479;
      midbarsel *= phosieie[idx] < 0.0103;
      midbarsel *= phohoe[idx] < (0.0241+1.28/energy+0.042*eventRho/energy);
      if(midbarsel) midbarphoidx.push_back(idx);
      
      bool ecphosel = true;
      ecphosel *= TMath::Abs(phoeta[idx]) > 1.479;
      if(ecphosel) ecphoidx.push_back(idx);

      bool id1ecsel = true;
      id1ecsel *= TMath::Abs(phoeta[idx]) > 1.479;
      id1ecsel *= phosieie[idx] < 0.0278;
      if(id1ecsel) id1ecphoidx.push_back(idx);
            
      bool midecsel = true;
      midecsel *= TMath::Abs(phoeta[idx]) > 1.479;
      midecsel *= phosieie[idx] < 0.0278;
      midecsel *= phohoe[idx] < (0.0274+2.08/energy+0.292*eventRho/energy);
      if(midecsel) midecphoidx.push_back(idx);
      
    } // End of loop on photons
    
    fillhistinevent("nosel_el", noselelidx);
    fillhistinevent("bar_el", barelidx);
    if(barelidx.size()>=2) fillhistinevent("gt2_bar_el", barelidx);
    if(id1barelidx.size()>=2) fillhistinevent("id1_gt2_bar_el", id1barelidx);
    if(id2barelidx.size()>=2) fillhistinevent("id2_gt2_bar_el", id2barelidx);
    if(midbarelidx.size()>=2) fillhistinevent("mid_gt2_bar_el", midbarelidx);
    if(mettrigs && midbarelidx.size()>=2) fillhistinevent("met_mid_gt2_bar_el", midbarelidx);
    if(t1p4nstrig && midbarelidx.size()>=2) fillhistinevent("t1p4_mid_gt2_bar_el", midbarelidx);
    if(sminlt0p12trig && midbarelidx.size()>=2) fillhistinevent("sm12_mid_gt2_bar_el", midbarelidx);
    if(t1p4nstrig && mettrigs && midbarelidx.size()>=2) fillhistinevent("t1p4_met_mid_gt2_bar_el", midbarelidx);
    if(sminlt0p12trig && mettrigs && midbarelidx.size()>=2) fillhistinevent("sm12_met_mid_gt2_bar_el", midbarelidx);
    fillhistinevent("ec_el", ecelidx);
    if(ecelidx.size()>=2) fillhistinevent("gt2_ec_el", ecelidx);
    if(id1ecelidx.size()>=2) fillhistinevent("id1_gt2_ec_el", id1ecelidx);
    if(id2ecelidx.size()>=2) fillhistinevent("id2_gt2_ec_el", id2ecelidx);
    if(midecelidx.size()>=2) fillhistinevent("mid_gt2_ec_el", midecelidx);
    if(mettrigs && midecelidx.size()>=2) fillhistinevent("met_mid_gt2_ec_el", midecelidx);
    if(t1p4nstrig && midecelidx.size()>=2) fillhistinevent("t1p4_mid_gt2_ec_el", midecelidx);
    if(sminlt0p12trig && midecelidx.size()>=2) fillhistinevent("sm12_mid_gt2_ec_el", midecelidx);
    if(t1p4nstrig && mettrigs && midecelidx.size()>=2) fillhistinevent("t1p4_met_mid_gt2_ec_el", midecelidx);
    if(sminlt0p12trig && mettrigs && midecelidx.size()>=2) fillhistinevent("sm12_met_mid_gt2_ec_el", midecelidx);

    fillhistinevent4LowPtElectron("nosel_lowptel", nosellowptelidx);
    fillhistinevent4LowPtElectron("bar_lowptel", barlowptelidx);
    if(barlowptelidx.size()>=2) fillhistinevent4LowPtElectron("gt2_bar_lowptel", barlowptelidx);
    if(id1barlowptelidx.size()>=2) fillhistinevent4LowPtElectron("id1_gt2_bar_lowptel", id1barlowptelidx);
    if(id2barlowptelidx.size()>=2) fillhistinevent4LowPtElectron("id2_gt2_bar_lowptel", id2barlowptelidx);
    if(midbarlowptelidx.size()>=2) fillhistinevent4LowPtElectron("mid_gt2_bar_lowptel", midbarlowptelidx);
    if(mettrigs && midbarlowptelidx.size()>=2) fillhistinevent4LowPtElectron("met_mid_gt2_bar_lowptel", midbarlowptelidx);
    if(t1p4nstrig && midbarlowptelidx.size()>=2) fillhistinevent4LowPtElectron("t1p4_mid_gt2_bar_lowptel", midbarlowptelidx);
    if(sminlt0p12trig && midbarlowptelidx.size()>=2) fillhistinevent4LowPtElectron("sm12_mid_gt2_bar_lowptel", midbarlowptelidx);
    if(t1p4nstrig && mettrigs && midbarlowptelidx.size()>=2) fillhistinevent4LowPtElectron("t1p4_met_mid_gt2_bar_lowptel", midbarlowptelidx);
    if(sminlt0p12trig && mettrigs && midbarlowptelidx.size()>=2) fillhistinevent4LowPtElectron("sm12_met_mid_gt2_bar_lowptel", midbarlowptelidx);
    fillhistinevent4LowPtElectron("ec_lowptel", eclowptelidx);
    if(eclowptelidx.size()>=2) fillhistinevent4LowPtElectron("gt2_ec_lowptel", eclowptelidx);
    if(id1eclowptelidx.size()>=2) fillhistinevent4LowPtElectron("id1_gt2_ec_lowptel", id1eclowptelidx);
    if(id2eclowptelidx.size()>=2) fillhistinevent4LowPtElectron("id2_gt2_ec_lowptel", id2eclowptelidx);
    if(mideclowptelidx.size()>=2) fillhistinevent4LowPtElectron("mid_gt2_ec_lowptel", mideclowptelidx);
    if(mettrigs && mideclowptelidx.size()>=2) fillhistinevent4LowPtElectron("met_mid_gt2_ec_lowptel", mideclowptelidx);
    if(t1p4nstrig && mideclowptelidx.size()>=2) fillhistinevent4LowPtElectron("t1p4_mid_gt2_ec_lowptel", mideclowptelidx);
    if(sminlt0p12trig && mideclowptelidx.size()>=2) fillhistinevent4LowPtElectron("sm12_mid_gt2_ec_lowptel", mideclowptelidx);
    if(t1p4nstrig && mettrigs && mideclowptelidx.size()>=2) fillhistinevent4LowPtElectron("t1p4_met_mid_gt2_ec_lowptel", mideclowptelidx);
    if(sminlt0p12trig && mettrigs && mideclowptelidx.size()>=2) fillhistinevent4LowPtElectron("sm12_met_mid_gt2_ec_lowptel", mideclowptelidx);

    fillhistinevent4Photon("nosel_pho", noselphoidx);
    fillhistinevent4Photon("bar_pho", barphoidx);
    if(barphoidx.size()>=2) fillhistinevent4Photon("gt2_bar_pho", barphoidx);
    if(id1barphoidx.size()>=2) fillhistinevent4Photon("id1_gt2_bar_pho", id1barphoidx);
    if(midbarphoidx.size()>=2) fillhistinevent4Photon("mid_gt2_bar_pho", midbarphoidx);
    if(mettrigs && midbarphoidx.size()>=2) fillhistinevent4Photon("met_mid_gt2_bar_pho", midbarphoidx);
    if(t1p4nstrig && midbarphoidx.size()>=2) fillhistinevent4Photon("t1p4_mid_gt2_bar_pho", midbarphoidx);
    if(sminlt0p12trig && midbarphoidx.size()>=2) fillhistinevent4Photon("sm12_mid_gt2_bar_pho", midbarphoidx);
    if(t1p4nstrig && mettrigs && midbarphoidx.size()>=2) fillhistinevent4Photon("t1p4_met_mid_gt2_bar_pho", midbarphoidx);
    if(sminlt0p12trig && mettrigs && midbarphoidx.size()>=2) fillhistinevent4Photon("sm12_met_mid_gt2_bar_pho", midbarphoidx);
    fillhistinevent4Photon("ec_pho", ecphoidx);
    if(ecphoidx.size()>=2) fillhistinevent4Photon("gt2_ec_pho", ecphoidx);
    if(id1ecphoidx.size()>=2) fillhistinevent4Photon("id1_gt2_ec_pho", id1ecphoidx);
    if(midecphoidx.size()>=2) fillhistinevent4Photon("mid_gt2_ec_pho", midecphoidx);
    if(mettrigs && midecphoidx.size()>=2) fillhistinevent4Photon("met_mid_gt2_ec_pho", midecphoidx);
    if(t1p4nstrig && midecphoidx.size()>=2) fillhistinevent4Photon("t1p4_mid_gt2_ec_pho", midecphoidx);
    if(sminlt0p12trig && midecphoidx.size()>=2) fillhistinevent4Photon("sm12_mid_gt2_ec_pho", midecphoidx);
    if(t1p4nstrig && mettrigs && midecphoidx.size()>=2) fillhistinevent4Photon("t1p4_met_mid_gt2_ec_pho", midecphoidx);
    if(sminlt0p12trig && mettrigs && midecphoidx.size()>=2) fillhistinevent4Photon("sm12_met_mid_gt2_ec_pho", midecphoidx);

    // Clear all the vectors
    noselelidx.clear();
    barelidx.clear();
    id1barelidx.clear();
    id2barelidx.clear();
    midbarelidx.clear();
    ecelidx.clear();
    id1ecelidx.clear();
    id2ecelidx.clear();
    midecelidx.clear();

    nosellowptelidx.clear();
    barlowptelidx.clear();
    id1barlowptelidx.clear();
    id2barlowptelidx.clear();
    midbarlowptelidx.clear();
    eclowptelidx.clear();
    id1eclowptelidx.clear();
    id2eclowptelidx.clear();
    mideclowptelidx.clear();

    noselphoidx.clear();
    barphoidx.clear();
    id1barphoidx.clear();
    midbarphoidx.clear();
    ecphoidx.clear();
    id1ecphoidx.clear();
    midecphoidx.clear();

  } // End of loop on events

  cout<<totEntries<<"\t"<<nosel<<endl;
}

// Function to fill a set of histograms in the event
void robustanalyzer::fillhistinevent(TString selection, vector<int> elidx) {

  TH1F* evtrho = (TH1F*) outfile->Get(selection+"_event_rho");
  TH1F* evttrig = (TH1F*) outfile->Get(selection+"_event_trigdec");

  TH1F* elmult = (TH1F*) outfile->Get(selection+"_el_mult");

  TH1F* leadelpt = (TH1F*) outfile->Get(selection+"_el_lead_pt");
  TH1F* leadeleta = (TH1F*) outfile->Get(selection+"_el_lead_eta");
  TH1F* leadelphi = (TH1F*) outfile->Get(selection+"_el_lead_phi");
  TH1F* leadelsieie = (TH1F*) outfile->Get(selection+"_el_lead_sieie");
  TH1F* leadelhoe = (TH1F*) outfile->Get(selection+"_el_lead_hoe");
  TH1F* leadeldetasc = (TH1F*) outfile->Get(selection+"_el_lead_detasc");
  TH1F* leadelooemoop = (TH1F*) outfile->Get(selection+"_el_lead_ooemoop");
  TH1F* leadelconvveto = (TH1F*) outfile->Get(selection+"_el_lead_convveto");
  TH1F* leadellog10d0 = (TH1F*) outfile->Get(selection+"_el_lead_log10d0");
  TH1F* leadeld0 = (TH1F*) outfile->Get(selection+"_el_lead_d0");
  TH1F* leadeld0Err = (TH1F*) outfile->Get(selection+"_el_lead_d0Err");
  TH1F* leadellog10dz = (TH1F*) outfile->Get(selection+"_el_lead_log10dz");
  TH1F* leadeldz = (TH1F*) outfile->Get(selection+"_el_lead_dz");
  TH1F* leadeldzErr = (TH1F*) outfile->Get(selection+"_el_lead_dzErr");

  TH1F* subleadelpt = (TH1F*) outfile->Get(selection+"_el_sublead_pt");
  TH1F* subleadeleta = (TH1F*) outfile->Get(selection+"_el_sublead_eta");
  TH1F* subleadelphi = (TH1F*) outfile->Get(selection+"_el_sublead_phi");
  TH1F* subleadelsieie = (TH1F*) outfile->Get(selection+"_el_sublead_sieie");
  TH1F* subleadelhoe = (TH1F*) outfile->Get(selection+"_el_sublead_hoe");
  TH1F* subleadeldetasc = (TH1F*) outfile->Get(selection+"_el_sublead_detasc");
  TH1F* subleadelooemoop = (TH1F*) outfile->Get(selection+"_el_sublead_ooemoop");
  TH1F* subleadelconvveto = (TH1F*) outfile->Get(selection+"_el_sublead_convveto");
  TH1F* subleadellog10d0 = (TH1F*) outfile->Get(selection+"_el_sublead_log10d0");
  TH1F* subleadeld0 = (TH1F*) outfile->Get(selection+"_el_sublead_d0");
  TH1F* subleadeld0Err = (TH1F*) outfile->Get(selection+"_el_sublead_d0Err");
  TH1F* subleadellog10dz = (TH1F*) outfile->Get(selection+"_el_sublead_log10dz");
  TH1F* subleadeldz = (TH1F*) outfile->Get(selection+"_el_sublead_dz");
  TH1F* subleadeldzErr = (TH1F*) outfile->Get(selection+"_el_sublead_dzErr");

  if(elidx.size()>0) {

    // Fill the event trigger decisions
    if(HLT_DiPhoton10sminlt0p12) evttrig->Fill(2);
    if(HLT_DiPhoton10Time1p4ns) evttrig->Fill(1);
    if(HLT_PFMET120_PFMHT120_IDTight) evttrig->Fill(-1);
    if(HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF) evttrig->Fill(-2);
    if(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight) evttrig->Fill(-3);
    if(HLT_CaloMET80_NotCleaned) evttrig->Fill(-4);
    if(HLT_PFMET200_NotCleaned) evttrig->Fill(-5);
    if(HLT_PFMET200_BeamHaloCleaned) evttrig->Fill(-6);
    if(HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight) evttrig->Fill(-7);
    
    evtrho->Fill(eventRho);
    elmult->Fill(elidx.size());
    leadelpt->Fill(elpt[elidx[0]]);
    leadeleta->Fill(eleta[elidx[0]]);
    leadelphi->Fill(elphi[elidx[0]]);
    leadelsieie->Fill(elsieie[elidx[0]]);
    leadelhoe->Fill(elhoe[elidx[0]]);
    leadeldetasc->Fill(eldetasc[elidx[0]]);
    leadelooemoop->Fill(elooemoop[elidx[0]]);
    leadelconvveto->Fill(elconvveto[elidx[0]]);
    leadellog10d0->Fill(TMath::Log10(TMath::Abs(eldxy[elidx[0]])));
    leadeld0->Fill(eldxy[elidx[0]]);
    leadeld0Err->Fill(eldxyerr[elidx[0]]);
    leadellog10dz->Fill(TMath::Log10(TMath::Abs(eldz[elidx[0]])));
    leadeldz->Fill(eldz[elidx[0]]);
    leadeldzErr->Fill(eldzerr[elidx[0]]);
  } // End of condition requiring atleast one eg object
  
  if(elidx.size()>=2) { // Condition requiring atleast two eg object
    //TLorentzVector leadeg, subleadeg;
    //leadeg.SetPtEtaPhiM(egRecoPt[egidx[0]],egRecoEta[egidx[0]],egRecoPhi[egidx[0]],0.106);
    //subleadeg.SetPtEtaPhiM(egRecoPt[egidx[1]],egRecoEta[egidx[1]],egRecoPhi[egidx[1]],0.106);

    subleadelpt->Fill(elpt[elidx[1]]);
    subleadeleta->Fill(eleta[elidx[1]]);
    subleadelphi->Fill(elphi[elidx[1]]);
    subleadelsieie->Fill(elsieie[elidx[1]]);
    subleadelhoe->Fill(elhoe[elidx[1]]);
    subleadeldetasc->Fill(eldetasc[elidx[1]]);
    subleadelooemoop->Fill(elooemoop[elidx[1]]);
    subleadelconvveto->Fill(elconvveto[elidx[1]]);
    subleadellog10d0->Fill(TMath::Log10(TMath::Abs(eldxy[elidx[1]])));
    subleadeld0->Fill(eldxy[elidx[1]]);
    subleadeld0Err->Fill(eldxyerr[elidx[1]]);
    subleadellog10dz->Fill(TMath::Log10(TMath::Abs(eldz[elidx[1]])));
    subleadeldz->Fill(eldz[elidx[1]]);
    subleadeldzErr->Fill(eldzerr[elidx[1]]);
  } // End of condition requiring atleast two eg object

}

// Function to fill a set of histograms in the event
void robustanalyzer::fillhistinevent4LowPtElectron(TString selection, vector<int> elidx) {

  TH1F* evtrho = (TH1F*) outfile->Get(selection+"_event_rho");
  TH1F* evttrig = (TH1F*) outfile->Get(selection+"_event_trigdec");

  TH1F* elmult = (TH1F*) outfile->Get(selection+"_lowptel_mult");

  TH1F* leadelpt = (TH1F*) outfile->Get(selection+"_lowptel_lead_pt");
  TH1F* leadeleta = (TH1F*) outfile->Get(selection+"_lowptel_lead_eta");
  TH1F* leadelphi = (TH1F*) outfile->Get(selection+"_lowptel_lead_phi");
  TH1F* leadelsieie = (TH1F*) outfile->Get(selection+"_lowptel_lead_sieie");
  TH1F* leadelhoe = (TH1F*) outfile->Get(selection+"_lowptel_lead_hoe");
  TH1F* leadeldetasc = (TH1F*) outfile->Get(selection+"_lowptel_lead_detasc");
  TH1F* leadelooemoop = (TH1F*) outfile->Get(selection+"_lowptel_lead_ooemoop");
  TH1F* leadelconvveto = (TH1F*) outfile->Get(selection+"_lowptel_lead_convveto");
  TH1F* leadellog10d0 = (TH1F*) outfile->Get(selection+"_lowptel_lead_log10d0");
  TH1F* leadeld0 = (TH1F*) outfile->Get(selection+"_lowptel_lead_d0");
  TH1F* leadeld0Err = (TH1F*) outfile->Get(selection+"_lowptel_lead_d0Err");
  TH1F* leadellog10dz = (TH1F*) outfile->Get(selection+"_lowptel_lead_log10dz");
  TH1F* leadeldz = (TH1F*) outfile->Get(selection+"_lowptel_lead_dz");
  TH1F* leadeldzErr = (TH1F*) outfile->Get(selection+"_lowptel_lead_dzErr");

  TH1F* subleadelpt = (TH1F*) outfile->Get(selection+"_lowptel_sublead_pt");
  TH1F* subleadeleta = (TH1F*) outfile->Get(selection+"_lowptel_sublead_eta");
  TH1F* subleadelphi = (TH1F*) outfile->Get(selection+"_lowptel_sublead_phi");
  TH1F* subleadelsieie = (TH1F*) outfile->Get(selection+"_lowptel_sublead_sieie");
  TH1F* subleadelhoe = (TH1F*) outfile->Get(selection+"_lowptel_sublead_hoe");
  TH1F* subleadeldetasc = (TH1F*) outfile->Get(selection+"_lowptel_sublead_detasc");
  TH1F* subleadelooemoop = (TH1F*) outfile->Get(selection+"_lowptel_sublead_ooemoop");
  TH1F* subleadelconvveto = (TH1F*) outfile->Get(selection+"_lowptel_sublead_convveto");
  TH1F* subleadellog10d0 = (TH1F*) outfile->Get(selection+"_lowptel_sublead_log10d0");
  TH1F* subleadeld0 = (TH1F*) outfile->Get(selection+"_lowptel_sublead_d0");
  TH1F* subleadeld0Err = (TH1F*) outfile->Get(selection+"_lowptel_sublead_d0Err");
  TH1F* subleadellog10dz = (TH1F*) outfile->Get(selection+"_lowptel_sublead_log10dz");
  TH1F* subleadeldz = (TH1F*) outfile->Get(selection+"_lowptel_sublead_dz");
  TH1F* subleadeldzErr = (TH1F*) outfile->Get(selection+"_lowptel_sublead_dzErr");

  if(elidx.size()>0) {

    // Fill the event trigger decisions
    if(HLT_DiPhoton10sminlt0p12) evttrig->Fill(2);
    if(HLT_DiPhoton10Time1p4ns) evttrig->Fill(1);
    if(HLT_PFMET120_PFMHT120_IDTight) evttrig->Fill(-1);
    if(HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF) evttrig->Fill(-2);
    if(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight) evttrig->Fill(-3);
    if(HLT_CaloMET80_NotCleaned) evttrig->Fill(-4);
    if(HLT_PFMET200_NotCleaned) evttrig->Fill(-5);
    if(HLT_PFMET200_BeamHaloCleaned) evttrig->Fill(-6);
    if(HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight) evttrig->Fill(-7);
    
    evtrho->Fill(eventRho);
    elmult->Fill(elidx.size());
    leadelpt->Fill(lowptelpt[elidx[0]]);
    leadeleta->Fill(lowpteleta[elidx[0]]);
    leadelphi->Fill(lowptelphi[elidx[0]]);
    leadelsieie->Fill(lowptelsieie[elidx[0]]);
    leadelhoe->Fill(lowptelhoe[elidx[0]]);
    leadelooemoop->Fill(lowptelooemoop[elidx[0]]);
    leadellog10d0->Fill(TMath::Log10(TMath::Abs(lowpteldxy[elidx[0]])));
    leadeld0->Fill(lowpteldxy[elidx[0]]);
    leadeld0Err->Fill(lowpteldxyerr[elidx[0]]);
    leadellog10dz->Fill(TMath::Log10(TMath::Abs(lowpteldz[elidx[0]])));
    leadeldz->Fill(lowpteldz[elidx[0]]);
    leadeldzErr->Fill(lowpteldzerr[elidx[0]]);
  } // End of condition requiring atleast one eg object
  
  if(elidx.size()>=2) { // Condition requiring atleast two eg object
    //TLorentzVector leadeg, subleadeg;
    //leadeg.SetPtEtaPhiM(egRecoPt[egidx[0]],egRecoEta[egidx[0]],egRecoPhi[egidx[0]],0.106);
    //subleadeg.SetPtEtaPhiM(egRecoPt[egidx[1]],egRecoEta[egidx[1]],egRecoPhi[egidx[1]],0.106);

    subleadelpt->Fill(lowptelpt[elidx[1]]);
    subleadeleta->Fill(lowpteleta[elidx[1]]);
    subleadelphi->Fill(lowptelphi[elidx[1]]);
    subleadelsieie->Fill(lowptelsieie[elidx[1]]);
    subleadelhoe->Fill(lowptelhoe[elidx[1]]);
    subleadelooemoop->Fill(lowptelooemoop[elidx[1]]);
    subleadellog10d0->Fill(TMath::Log10(TMath::Abs(lowpteldxy[elidx[1]])));
    subleadeld0->Fill(lowpteldxy[elidx[1]]);
    subleadeld0Err->Fill(lowpteldxyerr[elidx[1]]);
    subleadellog10dz->Fill(TMath::Log10(TMath::Abs(lowpteldz[elidx[1]])));
    subleadeldz->Fill(lowpteldz[elidx[1]]);
    subleadeldzErr->Fill(lowpteldzerr[elidx[1]]);
  } // End of condition requiring atleast two eg object

}

// Function to fill a set of histograms in the event
void robustanalyzer::fillhistinevent4Photon(TString selection, vector<int> phoidx) {

  TH1F* evtrho = (TH1F*) outfile->Get(selection+"_event_rho");
  TH1F* evttrig = (TH1F*) outfile->Get(selection+"_event_trigdec");

  TH1F* phomult = (TH1F*) outfile->Get(selection+"_pho_mult");

  TH1F* leadphopt = (TH1F*) outfile->Get(selection+"_pho_lead_pt");
  TH1F* leadphoeta = (TH1F*) outfile->Get(selection+"_pho_lead_eta");
  TH1F* leadphophi = (TH1F*) outfile->Get(selection+"_pho_lead_phi");
  TH1F* leadphosieie = (TH1F*) outfile->Get(selection+"_pho_lead_sieie");
  TH1F* leadphohoe = (TH1F*) outfile->Get(selection+"_pho_lead_hoe");

  TH1F* subleadphopt = (TH1F*) outfile->Get(selection+"_pho_sublead_pt");
  TH1F* subleadphoeta = (TH1F*) outfile->Get(selection+"_pho_sublead_eta");
  TH1F* subleadphophi = (TH1F*) outfile->Get(selection+"_pho_sublead_phi");
  TH1F* subleadphosieie = (TH1F*) outfile->Get(selection+"_pho_sublead_sieie");
  TH1F* subleadphohoe = (TH1F*) outfile->Get(selection+"_pho_sublead_hoe");

  if(phoidx.size()>0) {

    // Fill the event trigger decisions
    if(HLT_DiPhoton10sminlt0p12) evttrig->Fill(2);
    if(HLT_DiPhoton10Time1p4ns) evttrig->Fill(1);
    if(HLT_PFMET120_PFMHT120_IDTight) evttrig->Fill(-1);
    if(HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF) evttrig->Fill(-2);
    if(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight) evttrig->Fill(-3);
    if(HLT_CaloMET80_NotCleaned) evttrig->Fill(-4);
    if(HLT_PFMET200_NotCleaned) evttrig->Fill(-5);
    if(HLT_PFMET200_BeamHaloCleaned) evttrig->Fill(-6);
    if(HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight) evttrig->Fill(-7);
    
    evtrho->Fill(eventRho);
    phomult->Fill(phoidx.size());
    leadphopt->Fill(phopt[phoidx[0]]);
    leadphoeta->Fill(phoeta[phoidx[0]]);
    leadphophi->Fill(phophi[phoidx[0]]);
    leadphosieie->Fill(phosieie[phoidx[0]]);
    leadphohoe->Fill(phohoe[phoidx[0]]);
  } // End of condition requiring atleast one eg object
  
  if(phoidx.size()>=2) { // Condition requiring atleast two eg object
    //TLorentzVector leadeg, subleadeg;
    //leadeg.SetPtEtaPhiM(egRecoPt[egidx[0]],egRecoEta[egidx[0]],egRecoPhi[egidx[0]],0.106);
    //subleadeg.SetPtEtaPhiM(egRecoPt[egidx[1]],egRecoEta[egidx[1]],egRecoPhi[egidx[1]],0.106);

    subleadphopt->Fill(phopt[phoidx[1]]);
    subleadphoeta->Fill(phoeta[phoidx[1]]);
    subleadphophi->Fill(phophi[phoidx[1]]);
    subleadphosieie->Fill(phosieie[phoidx[1]]);
    subleadphohoe->Fill(phohoe[phoidx[1]]);
  } // End of condition requiring atleast two eg object

}

// Function to add a set of histograms for a selection
void robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"_event_rho","rho",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"_event_trigdec","Trigger Decisions",100,-50,50));

  all1dhists.push_back(new TH1F(selection+"_el_mult","N electron",50,-5,45));

  all1dhists.push_back(new TH1F(selection+"_el_lead_pt","electron_{1} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_el_lead_eta","electron_{1} #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_lead_phi","electron_{1} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_el_lead_sieie","electron_{1} #sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"_el_lead_hoe","electron_{1} H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"_el_lead_detasc","electron_{1} #Delta#eta(SC, trk seed)",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"_el_lead_ooemoop","electron_{1} E^{-1}-p^{-1} / GeV^{-1}",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"_el_lead_convveto","electron_{1} conversion veto",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"_el_lead_log10d0","electron_{1} log_{10}d_{0} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_lead_d0","electron_{1} d_{0} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_el_lead_d0Err","electron_{1} #sigmad_{0}",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"_el_lead_log10dz","electron_{1} log_{10}d_{z} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_lead_dz","electron_{1} d_{z} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_el_lead_dzErr","electron_{1} #sigmad_{z}",5,-2,3));

  all1dhists.push_back(new TH1F(selection+"_el_sublead_pt","electron_{2} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_eta","electron_{2} #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_phi","electron_{2} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_sieie","electron_{2} #sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_hoe","electron_{2} H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_detasc","electron_{2} #Delta#eta(SC, trk seed)",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_ooemoop","electron_{2} E^{-1}-p^{-1} / GeV^{-1}",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_convveto","electron_{2} conversion veto",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_log10d0","electron_{2} log_{10}d_{0} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_d0","electron_{2} d_{0} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_d0Err","electron_{2} #sigmad_{0}",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_log10dz","electron_{2} log_{10}d_{z} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_dz","electron_{2} d_{z} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_dzErr","electron_{2} #sigmad_{z}",5,-2,3));
}

// Function to add a set of histograms for a selection
void robustanalyzer::addhist4LowPtElectron(TString selection) {

  all1dhists.push_back(new TH1F(selection+"_event_rho","rho",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"_event_trigdec","Trigger Decisions",100,-50,50));

  all1dhists.push_back(new TH1F(selection+"_lowptel_mult","N low p_{T} elec",50,-5,45));

  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_pt","low p_{T} elec_{1} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_eta","low p_{T} elec_{1} #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_phi","low p_{T} elec_{1} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_sieie","low p_{T} elec_{1} #sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_hoe","low p_{T} elec_{1} H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_detasc","low p_{T} elec_{1} #Delta#eta(SC, trk seed)",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_ooemoop","low p_{T} elec_{1} E^{-1}-p^{-1} / GeV^{-1}",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_convveto","low p_{T} elec_{1} conversion veto",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_log10d0","low p_{T} elec_{1} log_{10}d_{0} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_d0","low p_{T} elec_{1} d_{0} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_d0Err","low p_{T} elec_{1} #sigmad_{0}",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_log10dz","low p_{T} elec_{1} log_{10}d_{z} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_dz","low p_{T} elec_{1} d_{z} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_lowptel_lead_dzErr","low p_{T} elec_{1} #sigmad_{z}",5,-2,3));

  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_pt","low p_{T} elec_{2} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_eta","low p_{T} elec_{2} #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_phi","low p_{T} elec_{2} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_sieie","low p_{T} elec_{2} #sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_hoe","low p_{T} elec_{2} H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_detasc","low p_{T} elec_{2} #Delta#eta(SC, trk seed)",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_ooemoop","low p_{T} elec_{2} E^{-1}-p^{-1} / GeV^{-1}",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_convveto","low p_{T} elec_{2} conversion veto",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_log10d0","low p_{T} elec_{2} log_{10}d_{0} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_d0","low p_{T} elec_{2} d_{0} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_d0Err","low p_{T} elec_{2} #sigmad_{0}",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_log10dz","low p_{T} elec_{2} log_{10}d_{z} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_dz","low p_{T} elec_{2} d_{z} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_lowptel_sublead_dzErr","low p_{T} elec_{2} #sigmad_{z}",5,-2,3));
}

// Function to add a set of histograms for a selection
void robustanalyzer::addhist4Photon(TString selection) {

  all1dhists.push_back(new TH1F(selection+"_event_rho","rho",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"_event_trigdec","Trigger Decisions",100,-50,50));

  all1dhists.push_back(new TH1F(selection+"_pho_mult","N photon",50,-5,45));

  all1dhists.push_back(new TH1F(selection+"_pho_lead_pt","photon_{1} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_pho_lead_eta","photon_{1} #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"_pho_lead_phi","photon_{1} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_pho_lead_sieie","photon_{1} #sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"_pho_lead_hoe","photon_{1} H/E",1000,0,10));

  all1dhists.push_back(new TH1F(selection+"_pho_sublead_pt","photon_{2} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_pho_sublead_eta","photon_{2} #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"_pho_sublead_phi","photon_{2} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_pho_sublead_sieie","photon_{2} #sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"_pho_sublead_hoe","photon_{2} H/E",1000,0,10));
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
