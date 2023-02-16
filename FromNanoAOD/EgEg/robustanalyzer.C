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
robustanalyzer::robustanalyzer(TString filename, TString outfilename, int numCores) {

  nC = numCores;

  auto inputChain = new TChain("demo/tree");
  cout<<"Initializing for file: "<<filename<<endl;
  inputChain->Add(filename);

  tree = new TTreeReader(inputChain);

  run = new TTreeReaderValue<int>((*tree), "run");
  lumi = new TTreeReaderValue<int>((*tree), "lumSec");
  rho = new TTreeReaderValue<double>((*tree), "rho");
  HLT_DiPhoton10sminlt0p12 = new TTreeReaderValue<bool>((*tree), "HLT_DiPhoton10sminlt0p12");
  HLT_DiPhoton10Time1p4ns = new TTreeReaderValue<bool>((*tree), "HLT_DiPhoton10Time1p4ns");
  HLTOR_METTrig = new TTreeReaderValue<bool>((*tree), "HLTOR_METTrig");

  eln = new TTreeReaderValue<int>((*tree), "ele_n");
  ele = new TTreeReaderValue<vector<double>>((*tree), "ele_e");
  elpt = new TTreeReaderValue<vector<double>>((*tree), "ele_pt");
  eleta = new TTreeReaderValue<vector<double>>((*tree), "ele_eta");
  elphi = new TTreeReaderValue<vector<double>>((*tree), "ele_phi");
  eld0 = new TTreeReaderValue<vector<double>>((*tree), "ele_d0");
  eldz = new TTreeReaderValue<vector<double>>((*tree), "ele_dz");
  elseedtime = new TTreeReaderValue<vector<double>>((*tree), "ele_seedtime");
  elsmin = new TTreeReaderValue<vector<double>>((*tree), "ele_smin");
  elsmaj = new TTreeReaderValue<vector<double>>((*tree), "ele_smaj");
  elsieie = new TTreeReaderValue<vector<double>>((*tree), "ele_sinin_noiseclnd");
  eldeta = new TTreeReaderValue<vector<double>>((*tree), "ele_detaseed");
  eldphi = new TTreeReaderValue<vector<double>>((*tree), "ele_dphiin");
  elhoe = new TTreeReaderValue<vector<double>>((*tree), "ele_hoe");
  elchhadiso = new TTreeReaderValue<vector<double>>((*tree), "ele_chargedhadroniso");
  elneuthadiso = new TTreeReaderValue<vector<double>>((*tree), "ele_neutralhadroniso");
  elphiso = new TTreeReaderValue<vector<double>>((*tree), "ele_photoniso");
  elooemoop = new TTreeReaderValue<vector<double>>((*tree), "ele_ooemoop");
  elinnerhits = new TTreeReaderValue<vector<double>>((*tree), "ele_innerhits");
  elconvveto = new TTreeReaderValue<vector<bool>>((*tree), "ele_convveto");
  
  outfile = new TFile(outfilename,"RECREATE");

}

// Fill the root file, close the root file, and handle deletions
robustanalyzer::~robustanalyzer() {

  outfile->Write();
  outfile->Close();
}

// Analyzer for a single file
void robustanalyzer::analyzersinglefile(int splitCnt) {

  int totEntries = tree->GetEntries(true);
  cout<<"Total number of entries: "<<totEntries<<endl;
  int nCores = nC; // Assume parallel processing over numCores cores where
  // there is a lesser no.of events in the last core
  int beginevent = splitCnt*(totEntries/nCores);
  int endevent = (splitCnt+1)*(totEntries/nCores);
  endevent = endevent<totEntries?endevent:totEntries; // Verfied that this logic to parallelize works
  tree->SetEntriesRange(beginevent, endevent);
  cout<<"Runing "<<splitCnt<<" process from "<<beginevent<<" to "<<endevent<<endl;
  int event = beginevent-1;
  
  // Count events passing certain selections
  int nosel=0;
  addhist("nosel_el");
  addhist("gt2_el");
  addhist("bar_el");
  addhist("gt2_bar_el");
  addhist("ec_el");
  addhist("gt2_ec_el");
  addhist("id1_gt2_bar_el");
  addhist("id2_gt2_bar_el");
  addhist("mid_gt2_bar_el");
  addhist("met_mid_gt2_bar_el");
  addhist("t1p4_mid_gt2_bar_el");
  addhist("sm12_mid_gt2_bar_el");
  addhist("t1p4_met_mid_gt2_bar_el");
  addhist("sm12_met_mid_gt2_bar_el");
  addhist("id1_gt2_ec_el");
  addhist("id2_gt2_ec_el");
  addhist("mid_gt2_ec_el");
  addhist("met_mid_gt2_ec_el");
  addhist("t1p4_mid_gt2_ec_el");
  addhist("sm12_mid_gt2_ec_el");
  addhist("t1p4_met_mid_gt2_ec_el");
  addhist("sm12_met_mid_gt2_ec_el");
  
  vector<int> noselelidx;
  vector<int> barelidx;
  vector<int> ecelidx;
  vector<int> id1barelidx;
  vector<int> id2barelidx;
  vector<int> midbarelidx;
  vector<int> id1ecelidx;
  vector<int> id2ecelidx;
  vector<int> midecelidx;
  
  // Loop beginning on events
  while(tree->Next()) {
    
    //if(event>1000) break;
    //if(event!=283991 && event!=326114) continue;
    event++;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    // Section for trigger conditions
    Bool_t mettrigs = (*(*HLTOR_METTrig));
    Bool_t t1p4nstrig = (*(*HLT_DiPhoton10Time1p4ns));
    Bool_t sminlt0p12trig = (*(*HLT_DiPhoton10sminlt0p12));
    
    // Loop beginning on electrons
    for(unsigned int idx=0; idx<(*(*eln)); idx++) {

      if(idx>0) {
	if((*elpt)->at(idx)>(*elpt)->at(idx-1)) // Check for pt sort
	  throw "Error!!! Electrons are not pt sorted";
      }

      double energy = (*ele)->at(idx);
      
      noselelidx.push_back(idx);

      bool barelsel = true;
      barelsel *= TMath::Abs((*eleta)->at(idx))<1.479;
      if(barelsel) barelidx.push_back(idx);

      bool id1barsel = true;
      id1barsel *= TMath::Abs((*eleta)->at(idx)) < 1.479;
      id1barsel *= (*elsieie)->at(idx) < 0.0103;
      if(id1barsel) id1barelidx.push_back(idx);
      
      bool id2barsel = true;
      id2barsel *= TMath::Abs((*eleta)->at(idx)) < 1.479;
      id2barsel *= (*elsieie)->at(idx) < 0.0103;
      id2barsel *= (*elhoe)->at(idx) < (0.0241+1.28/energy+0.042*(*(*rho))/energy);
      if(id2barsel) id2barelidx.push_back(idx);
      
      bool midbarsel = true;
      midbarsel *= TMath::Abs((*eleta)->at(idx)) < 1.479;
      midbarsel *= (*elsieie)->at(idx) < 0.0103;
      //midbarsel *= eldetasc[idx] < 0.00481;
      midbarsel *= (*elhoe)->at(idx) < (0.0241+1.28/energy+0.042*(*(*rho))/energy);
      midbarsel *= (*elooemoop)->at(idx) < 0.0966;
      //midbarsel *= elconvveto[idx];
      if(midbarsel) midbarelidx.push_back(idx);
      
      bool ecelsel = true;
      ecelsel *= TMath::Abs((*eleta)->at(idx)) > 1.479;
      if(ecelsel) ecelidx.push_back(idx);
      
      bool id1ecsel = true;
      id1ecsel *= TMath::Abs((*eleta)->at(idx)) > 1.479;
      id1ecsel *= (*elsieie)->at(idx) < 0.0278;
      if(id1ecsel) id1ecelidx.push_back(idx);
      
      bool id2ecsel = true;
      id2ecsel *= TMath::Abs((*eleta)->at(idx)) > 1.479;
      id2ecsel *= (*elsieie)->at(idx) < 0.0278;
      id2ecsel *= (*elhoe)->at(idx) < (0.0274+2.08/energy+0.292*(*(*rho))/energy);
      if(id2ecsel) id2ecelidx.push_back(idx);
      
      bool midecsel = true;
      midecsel *= TMath::Abs((*eleta)->at(idx)) > 1.479;
      midecsel *= (*elsieie)->at(idx) < 0.0278;
      //midecsel *= eldetasc[idx] < 0.00847;
      midecsel *= (*elhoe)->at(idx) < (0.0274+2.08/energy+0.292*(*(*rho))/energy);
      midecsel *= (*elooemoop)->at(idx) < 0.0769;
      //midecsel *= elconvveto[idx];
      if(midecsel) midecelidx.push_back(idx);
      
    } // End of loop on electrons
               
    fillhistinevent("nosel_el", noselelidx);
    if(noselelidx.size()>=2) fillhistinevent("gt2_el", noselelidx);
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
    
    // Clear all the vectors
    noselelidx.clear();
    barelidx.clear();
    ecelidx.clear();
    id1barelidx.clear();
    id2barelidx.clear();
    midbarelidx.clear();
    id1ecelidx.clear();
    id2ecelidx.clear();
    midecelidx.clear();
    
  } // End of loop on events
  cout<<totEntries<<"\t"<<nosel<<endl;
}

// Function to add a set of histograms for a selection
void robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"_event_rho","rho",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"_event_trigdec","Trigger Decisions",20,0,20));

  all1dhists.push_back(new TH1F(selection+"_el_mult","N electron",50,-5,45));

  all1dhists.push_back(new TH1F(selection+"_el_lead_energy","electron_{1} E / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_el_lead_pt","electron_{1} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_el_lead_eta","electron_{1} #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_lead_phi","electron_{1} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_el_lead_d0","electron_{1} d_{0} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_el_lead_log10d0","electron_{1} log_{10}d_{0} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_lead_dz","electron_{1} d_{z} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_el_lead_log10dz","electron_{1} log_{10}d_{z} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_lead_seedtime","electron_{1} time_{SC} / ns",10000,-2,8));
  all1dhists.push_back(new TH1F(selection+"_el_lead_smin","electron_{1} s_{min}",20000,-2,8));
  all1dhists.push_back(new TH1F(selection+"_el_lead_smaj","electron_{1} s_{maj}",15000,-2,13));
  all1dhists.push_back(new TH1F(selection+"_el_lead_sieie","electron_{1} #sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"_el_lead_detasc","electron_{1} #Delta#eta(SC, trk seed)",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"_el_lead_dphi","electron_{1} #Delta#phi(SC, trk seed)",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"_el_lead_hoe","electron_{1} H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"_el_lead_relisowithea","electron_{1} rel.iso.",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"_el_lead_ooemoop","electron_{1} E^{-1}-p^{-1} / GeV^{-1}",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"_el_lead_innerhits","electron_{1} inner hits",50,-1,49));
  all1dhists.push_back(new TH1F(selection+"_el_lead_convveto","electron_{1} conversion veto",5,-2,3));

  all1dhists.push_back(new TH1F(selection+"_el_sublead_energy","electron_{2} E / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_pt","electron_{2} p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_eta","electron_{2} #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_phi","electron_{2} #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_d0","electron_{2} d_{0} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_log10d0","electron_{2} log_{10}d_{0} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_dz","electron_{2} d_{z} / cm",10000,-50,50));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_log10dz","electron_{2} log_{10}d_{z} / log_{10}cm",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_seedtime","electron_{2} time_{SC} / ns",10000,-2,8));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_smin","electron_{2} s_{min}",20000,-2,8));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_smaj","electron_{2} s_{maj}",15000,-2,13));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_sieie","electron_{2} #sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_detasc","electron_{2} #Delta#eta(SC, trk seed)",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_dphi","electron_{2} #Delta#phi(SC, trk seed)",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_hoe","electron_{2} H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_relisowithea","electron_{2} rel.iso.",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_ooemoop","electron_{2} E^{-1}-p^{-1} / GeV^{-1}",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_innerhits","electron_{2} inner hits",50,-1,49));
  all1dhists.push_back(new TH1F(selection+"_el_sublead_convveto","electron_{2} conversion veto",5,-2,3));
}

// Function to fill a set of histograms in the event
void robustanalyzer::fillhistinevent(TString selection, vector<int> elidx) {

  TH1F* evtrho = (TH1F*) outfile->Get(selection+"_event_rho");
  TH1F* evttrig = (TH1F*) outfile->Get(selection+"_event_trigdec");

  TH1F* elmult = (TH1F*) outfile->Get(selection+"_el_mult");

  TH1F* leadele = (TH1F*) outfile->Get(selection+"_el_lead_energy");
  TH1F* leadelpt = (TH1F*) outfile->Get(selection+"_el_lead_pt");
  TH1F* leadeleta = (TH1F*) outfile->Get(selection+"_el_lead_eta");
  TH1F* leadelphi = (TH1F*) outfile->Get(selection+"_el_lead_phi");
  TH1F* leadeld0 = (TH1F*) outfile->Get(selection+"_el_lead_d0");
  TH1F* leadellog10d0 = (TH1F*) outfile->Get(selection+"_el_lead_log10d0");
  TH1F* leadeldz = (TH1F*) outfile->Get(selection+"_el_lead_dz");
  TH1F* leadellog10dz = (TH1F*) outfile->Get(selection+"_el_lead_log10dz");
  TH1F* leadelstime = (TH1F*) outfile->Get(selection+"_el_lead_seedtime");
  TH1F* leadelsmin = (TH1F*) outfile->Get(selection+"_el_lead_smin");
  TH1F* leadelsmaj = (TH1F*) outfile->Get(selection+"_el_lead_smaj");
  TH1F* leadelsieie = (TH1F*) outfile->Get(selection+"_el_lead_sieie");
  TH1F* leadeldeta = (TH1F*) outfile->Get(selection+"_el_lead_detasc");
  TH1F* leadeldphi = (TH1F*) outfile->Get(selection+"_el_lead_dphi");
  TH1F* leadelhoe = (TH1F*) outfile->Get(selection+"_el_lead_hoe");
  TH1F* leadelrelisowea = (TH1F*) outfile->Get(selection+"_el_lead_relisowithea");
  TH1F* leadelooemooop = (TH1F*) outfile->Get(selection+"_el_lead_ooemoop");
  TH1F* leadelinnerhits = (TH1F*) outfile->Get(selection+"_el_lead_innerhits");
  TH1F* leadelconvveto = (TH1F*) outfile->Get(selection+"_el_lead_convveto");

  TH1F* subleadele = (TH1F*) outfile->Get(selection+"_el_sublead_energy");
  TH1F* subleadelpt = (TH1F*) outfile->Get(selection+"_el_sublead_pt");
  TH1F* subleadeleta = (TH1F*) outfile->Get(selection+"_el_sublead_eta");
  TH1F* subleadelphi = (TH1F*) outfile->Get(selection+"_el_sublead_phi");
  TH1F* subleadeld0 = (TH1F*) outfile->Get(selection+"_el_sublead_d0");
  TH1F* subleadellog10d0 = (TH1F*) outfile->Get(selection+"_el_sublead_log10d0");
  TH1F* subleadeldz = (TH1F*) outfile->Get(selection+"_el_sublead_dz");
  TH1F* subleadellog10dz = (TH1F*) outfile->Get(selection+"_el_sublead_log10dz");
  TH1F* subleadelstime = (TH1F*) outfile->Get(selection+"_el_sublead_seedtime");
  TH1F* subleadelsmin = (TH1F*) outfile->Get(selection+"_el_sublead_smin");
  TH1F* subleadelsmaj = (TH1F*) outfile->Get(selection+"_el_sublead_smaj");
  TH1F* subleadelsieie = (TH1F*) outfile->Get(selection+"_el_sublead_sieie");
  TH1F* subleadeldeta = (TH1F*) outfile->Get(selection+"_el_sublead_detasc");
  TH1F* subleadeldphi = (TH1F*) outfile->Get(selection+"_el_sublead_dphi");
  TH1F* subleadelhoe = (TH1F*) outfile->Get(selection+"_el_sublead_hoe");
  TH1F* subleadelrelisowea = (TH1F*) outfile->Get(selection+"_el_sublead_relisowithea");
  TH1F* subleadelooemooop = (TH1F*) outfile->Get(selection+"_el_sublead_ooemoop");
  TH1F* subleadelinnerhits = (TH1F*) outfile->Get(selection+"_el_sublead_innerhits");
  TH1F* subleadelconvveto = (TH1F*) outfile->Get(selection+"_el_sublead_convveto");

  // Fill the event trigger decisions
  if((*(*HLT_DiPhoton10sminlt0p12))) evttrig->Fill(4);
  if((*(*HLT_DiPhoton10Time1p4ns))) evttrig->Fill(2);
  if((*(*HLTOR_METTrig))) evttrig->Fill(-2);
  if(!(*(*HLTOR_METTrig)) && (*(*HLT_DiPhoton10Time1p4ns)) && (*(*HLT_DiPhoton10sminlt0p12))) evttrig->Fill(0);
  
  evtrho->Fill((*(*rho)));
    
  if(elidx.size()>0) {
    
    elmult->Fill(elidx.size());

    int leadelidx = elidx[0];
    leadele->Fill((*ele)->at(leadelidx));
    leadelpt->Fill((*elpt)->at(leadelidx));
    leadeleta->Fill((*eleta)->at(leadelidx));
    leadelphi->Fill((*elphi)->at(leadelidx));
    leadeld0->Fill((*eld0)->at(leadelidx));
    leadellog10d0->Fill(TMath::Log10(TMath::Abs((*eld0)->at(leadelidx))));
    leadeldz->Fill((*eldz)->at(leadelidx));
    leadellog10dz->Fill(TMath::Log10(TMath::Abs((*eldz)->at(leadelidx))));
    leadelstime->Fill((*elseedtime)->at(leadelidx));
    leadelsmin->Fill((*elsmin)->at(leadelidx));
    leadelsmaj->Fill((*elsmaj)->at(leadelidx));
    leadelsieie->Fill((*elsieie)->at(leadelidx));
    leadeldeta->Fill((*eldeta)->at(leadelidx));
    leadeldphi->Fill((*eldphi)->at(leadelidx));
    leadelhoe->Fill((*elhoe)->at(leadelidx));
    double ea = effectivearea((*eleta)->at(leadelidx));
    double neutiso = ((*elneuthadiso)->at(leadelidx))+((*elphiso)->at(leadelidx))-((*(*rho))*ea);
    double reliso = ((*elchhadiso)->at(leadelidx))+(neutiso>0?neutiso:0);
    leadelrelisowea->Fill(reliso);
    leadelooemooop->Fill((*elooemoop)->at(leadelidx));
    leadelinnerhits->Fill((*elinnerhits)->at(leadelidx));
    leadelconvveto->Fill((*elconvveto)->at(leadelidx));
  } // End of condition requiring atleast one eg object
  
  if(elidx.size()>=2) { // Condition requiring atleast two eg object
    //TLorentzVector leadeg, subleadeg;
    //leadeg.SetPtEtaPhiM(egRecoPt[egidx[0]],egRecoEta[egidx[0]],egRecoPhi[egidx[0]],0.106);
    //subleadeg.SetPtEtaPhiM(egRecoPt[egidx[1]],egRecoEta[egidx[1]],egRecoPhi[egidx[1]],0.106);

    int subleadelidx = elidx[1];
    subleadele->Fill((*ele)->at(subleadelidx));
    subleadelpt->Fill((*elpt)->at(subleadelidx));
    subleadeleta->Fill((*eleta)->at(subleadelidx));
    subleadelphi->Fill((*elphi)->at(subleadelidx));
    subleadeld0->Fill((*eld0)->at(subleadelidx));
    subleadellog10d0->Fill(TMath::Log10(TMath::Abs((*eld0)->at(subleadelidx))));
    subleadeldz->Fill((*eldz)->at(subleadelidx));
    subleadellog10dz->Fill(TMath::Log10(TMath::Abs((*eldz)->at(subleadelidx))));
    subleadelstime->Fill((*elseedtime)->at(subleadelidx));
    subleadelsmin->Fill((*elsmin)->at(subleadelidx));
    subleadelsmaj->Fill((*elsmaj)->at(subleadelidx));
    subleadelsieie->Fill((*elsieie)->at(subleadelidx));
    subleadeldeta->Fill((*eldeta)->at(subleadelidx));
    subleadeldphi->Fill((*eldphi)->at(subleadelidx));
    subleadelhoe->Fill((*elhoe)->at(subleadelidx));
    double ea = effectivearea((*eleta)->at(subleadelidx));
    double neutiso = ((*elneuthadiso)->at(subleadelidx))+((*elphiso)->at(subleadelidx))-((*(*rho))*ea);
    double reliso = ((*elchhadiso)->at(subleadelidx))+(neutiso>0?neutiso:0);
    subleadelrelisowea->Fill(reliso);
    subleadelooemooop->Fill((*elooemoop)->at(subleadelidx));
    subleadelinnerhits->Fill((*elinnerhits)->at(subleadelidx));
    subleadelconvveto->Fill((*elconvveto)->at(subleadelidx));

  } // End of condition requiring atleast two eg object

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

// Function to find the effective area based on slide 7 at
// https://indico.cern.ch/event/1204275/contributions/5064343/attachments/2529616/4355064/Electron_cutbasedID_preliminaryID.pdfs
double robustanalyzer::effectivearea(double eta) {

  double ea = -1.0;
  double abseta = abs(eta);

  if(abseta < 1.0) {
    ea = 0.1243;
  }
  else if(abseta >= 1.0 && abseta < 1.479) {
    ea = 0.1458;
  }
  else if(abseta >= 1.479 && abseta < 2.0) {
    ea = 0.0992;
  }
  else if(abseta >= 2.0 && abseta < 2.2) {
    ea = 0.0794;
  }
  else if(abseta >= 2.2 && abseta < 2.3) {
    ea = 0.0762;
  }
  else if(abseta >= 2.3 && abseta < 2.4) {
    ea = 0.0766;
  }
  else if(abseta >= 2.4 && abseta < 2.5) {
    ea = 0.1003;
  }
  else if(abseta >= 2.5 && abseta < 2.6) {
    ea = 0.2322;
  }
  else if(abseta >= 2.6 && abseta < 2.65) {
    ea = 0.2537;
  }
  else if(abseta >= 2.65 && abseta < 2.7) {
    ea = 0.2529;
  }
  else if(abseta >= 2.7 && abseta < 2.8) {
    ea = 0.2563;
  }
  else if(abseta >= 2.8 && abseta < 3.0) {
    ea = 0.2423;
  }
  else {
    ea = 0.0;
  }
  
  return ea;
}
