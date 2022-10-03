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

  addPhotonCollectionhist("nosel_pho_idx");
  addPhotonCollectionhist("nosel_EB_pho_idx");
  addPhotonCollectionhist("nosel_EE_pho_idx");
  addPhotonCollectionhist("sminlt0p12_EB_pho_idx");
  addPhotonCollectionhist("sminlt0p12_EE_pho_idx");
  addPhotonCollectionhist("time1p4ns_EB_pho_idx");
  addPhotonCollectionhist("time1p4ns_EE_pho_idx");

  addObjectFilterAngMatchhist("sminlt0p12_noselpho_angmch");
  addObjectFilterAngMatchhist("time1p4ns_noselpho_angmch");
  
  // Loop beginning on events
  for(unsigned int event=beginevent; event<endevent; event++) {
    
    inputChain->GetEntry(event);
    //if(event>10) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%1000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    // Cross-checks
    if(dieg10sminlt0p12_usfinfilt_pt->size() != dieg10sminlt0p12_usfinfilt_n) throw "Error!!! Number of entries in vector \"dieg10sminlt0p12_usfinfilt_pt\" unequal entries comapred to \"dieg10sminlt0p12_usfinfilt_n\".";
    if(HLT_DiPhoton10sminlt0p12 && dieg10sminlt0p12_usfinfilt_n<2) throw "Error TypeA!!! Filter dieg10sminlt0p12_usfinfilt_n<2 while HLT_DiPhoton10sminlt0p12 trigger passed.";
    if(!HLT_DiPhoton10sminlt0p12 && dieg10sminlt0p12_usfinfilt_n>=2) throw "Warning TypeB!!! Filter dieg10sminlt0p12_usfinfilt_n>=2 while HLT_DiPhoton10sminlt0p12 trigger not passed. Disable excpetion if you know what you are doing.";

    vector<int> dieg10sminlt0p12_finalidx(dieg10sminlt0p12_usfinfilt_n);
    iota(begin(dieg10sminlt0p12_finalidx), end(dieg10sminlt0p12_finalidx), 0);
    if(HLT_DiPhoton10sminlt0p12) fillHLTFilterhist("dieg10sminlt0p12_final", dieg10sminlt0p12_finalidx, dieg10sminlt0p12_usfinfilt_pt, dieg10sminlt0p12_usfinfilt_eta, dieg10sminlt0p12_usfinfilt_phi);
    
    vector<int> dieg10time1p4ns_finalidx(dieg10time1p4ns_usfinfilt_n);
    iota(begin(dieg10time1p4ns_finalidx), end(dieg10time1p4ns_finalidx), 0);
    if(HLT_DiPhoton10Time1p4ns) fillHLTFilterhist("dieg10time1p4ns_final", dieg10time1p4ns_finalidx, dieg10time1p4ns_usfinfilt_pt, dieg10time1p4ns_usfinfilt_eta, dieg10time1p4ns_usfinfilt_phi);

    vector<int> noselpho_idx(pho_pt->size());
    vector<int> noselEBpho_idx;
    vector<int> noselEEpho_idx;
    iota(begin(noselpho_idx), end(noselpho_idx), 0);
    for(int idx : noselpho_idx) {
      if(abs(pho_eta->at(idx))<1.479) {
	noselEBpho_idx.push_back(idx);
      }
      else {
	noselEEpho_idx.push_back(idx);
      }
    }
    
    fillPhotonCollectionhist("nosel_pho_idx", noselpho_idx);
    fillPhotonCollectionhist("nosel_EB_pho_idx", noselEBpho_idx);
    fillPhotonCollectionhist("nosel_EE_pho_idx", noselEEpho_idx);
    
    if(HLT_DiPhoton10sminlt0p12) {
      vector<pair<int,int>> sminlt0p12_noselpho_angmch_indices = fillObjectFilterAngMatchhist("sminlt0p12_noselpho_angmch", dieg10sminlt0p12_finalidx, dieg10sminlt0p12_usfinfilt_pt, dieg10sminlt0p12_usfinfilt_eta, dieg10sminlt0p12_usfinfilt_phi, noselpho_idx, pho_pt, pho_eta, pho_phi);
      vector<int> sminlt0p12_noselEBpho_idx = getFiltMatchedPhoIndex(sminlt0p12_noselpho_angmch_indices, 0.0, 1.479);
      vector<int> sminlt0p12_noselEEpho_idx = getFiltMatchedPhoIndex(sminlt0p12_noselpho_angmch_indices, 1.479, 4);
      if(sminlt0p12_noselEBpho_idx[0]!=-1) fillPhotonCollectionhist("sminlt0p12_EB_pho_idx", sminlt0p12_noselEBpho_idx);
      if(sminlt0p12_noselEEpho_idx[0]!=-1) fillPhotonCollectionhist("sminlt0p12_EE_pho_idx", sminlt0p12_noselEEpho_idx);
      sminlt0p12_noselEBpho_idx.clear();
      sminlt0p12_noselEEpho_idx.clear();
    }

    if(HLT_DiPhoton10Time1p4ns) {
      vector<pair<int,int>> time1p4ns_noselpho_angmch_indices = fillObjectFilterAngMatchhist("time1p4ns_noselpho_angmch", dieg10time1p4ns_finalidx, dieg10time1p4ns_usfinfilt_pt, dieg10time1p4ns_usfinfilt_eta, dieg10time1p4ns_usfinfilt_phi, noselpho_idx, pho_pt, pho_eta, pho_phi);
      vector<int> time1p4ns_noselEBpho_idx = getFiltMatchedPhoIndex(time1p4ns_noselpho_angmch_indices, 0.0, 1.479);
      vector<int> time1p4ns_noselEEpho_idx = getFiltMatchedPhoIndex(time1p4ns_noselpho_angmch_indices, 1.479, 4);
      if(time1p4ns_noselEBpho_idx[0]!=-1) fillPhotonCollectionhist("time1p4ns_EB_pho_idx", time1p4ns_noselEBpho_idx);
      if(time1p4ns_noselEEpho_idx[0]!=-1) fillPhotonCollectionhist("time1p4ns_EE_pho_idx", time1p4ns_noselEEpho_idx);
      time1p4ns_noselEBpho_idx.clear();
      time1p4ns_noselEEpho_idx.clear();
    }

    dieg10sminlt0p12_finalidx.clear();
    dieg10time1p4ns_finalidx.clear();
    noselpho_idx.clear();
    noselEBpho_idx.clear();
    noselEEpho_idx.clear();
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
  TH1F* hseedtime = (TH1F*) outfile->Get(selection+"_pho_seedtime");
  TH1F* hsmin = (TH1F*) outfile->Get(selection+"_pho_smin");
  TH1F* hsmax = (TH1F*) outfile->Get(selection+"_pho_smax");

  hmult->Fill(idx.size());
  for(int index : idx) {
    hpt->Fill(pho_pt->at(index));
    heta->Fill(pho_eta->at(index));
    hphi->Fill(pho_phi->at(index));
    hseedtime->Fill(pho_seedtime->at(index));
    hsmin->Fill(pho_smin->at(index));
    hsmax->Fill(pho_smax->at(index));
  }
}

// Function to fill a set of histograms for angular matching between filters and objects
vector< pair<int,int> > robustanalyzer::fillObjectFilterAngMatchhist(TString selection, vector<int> filtidx, vector<double> *filtpt, vector<double> *filteta, vector<double> *filtphi, vector<int> objidx, vector<double> *objpt, vector<double> *objeta, vector<double> *objphi) {

  //if(filtidx.size()==0 || objidx.size()==0) return;
  
  TH1F* premchfiltebDeta = (TH1F*) outfile->Get(selection+"_filteb_obj_premch_Deta");
  TH1F* premchfiltebDphi = (TH1F*) outfile->Get(selection+"_filteb_obj_premch_Dphi");
  TH1F* premchfiltebDr = (TH1F*) outfile->Get(selection+"_filteb_obj_premch_Dr");
  TH1F* premchfiltebDpt = (TH1F*) outfile->Get(selection+"_filteb_obj_premch_Dpt");
  TH1F* premchfilteeDeta = (TH1F*) outfile->Get(selection+"_filtee_obj_premch_Deta");
  TH1F* premchfilteeDphi = (TH1F*) outfile->Get(selection+"_filtee_obj_premch_Dphi");
  TH1F* premchfilteeDr = (TH1F*) outfile->Get(selection+"_filtee_obj_premch_Dr");
  TH1F* premchfilteeDpt = (TH1F*) outfile->Get(selection+"_filtee_obj_premch_Dpt");
  
  TH1F* aftmchfiltebDeta = (TH1F*) outfile->Get(selection+"_filteb_obj_aftmch_Deta");
  TH1F* aftmchfiltebDphi = (TH1F*) outfile->Get(selection+"_filteb_obj_aftmch_Dphi");
  TH1F* aftmchfiltebDr = (TH1F*) outfile->Get(selection+"_filteb_obj_aftmch_Dr");
  TH1F* aftmchfiltebDpt = (TH1F*) outfile->Get(selection+"_filteb_obj_aftmch_Dpt");
  TH1F* aftmchfilteeDeta = (TH1F*) outfile->Get(selection+"_filtee_obj_aftmch_Deta");
  TH1F* aftmchfilteeDphi = (TH1F*) outfile->Get(selection+"_filtee_obj_aftmch_Dphi");
  TH1F* aftmchfilteeDr = (TH1F*) outfile->Get(selection+"_filtee_obj_aftmch_Dr");
  TH1F* aftmchfilteeDpt = (TH1F*) outfile->Get(selection+"_filtee_obj_aftmch_Dpt");

  vector< pair<int,int> > matchedObjectFilterPairs;
  vector<bool> objectMatched(objidx.size(), false);
  
  // Before filter object matching
  for(int filtindex : filtidx) {
    double lowestDeta=9e9, lowestDphi=9e9, lowestDr=9e9, lowestDpt=9e9;
    TVector3 filtvec;
    filtvec.SetPtEtaPhi(filtpt->at(filtindex), filteta->at(filtindex), filtphi->at(filtindex));
    matchedObjectFilterPairs.push_back(make_pair(filtindex,-1));
    for(int objindex : objidx) {
      TVector3 objvec;
      objvec.SetPtEtaPhi(objpt->at(objindex), objeta->at(objindex), objphi->at(objindex));
      if(abs(filteta->at(filtindex)-objeta->at(objindex)) < abs(lowestDeta)) lowestDeta = filteta->at(filtindex)-objeta->at(objindex);
      if(abs(filtvec.DeltaPhi(objvec)) < abs(lowestDphi)) lowestDphi = filtvec.DeltaPhi(objvec);
      if(abs(filtvec.DeltaR(objvec)) < abs(lowestDr)) lowestDr = filtvec.DeltaR(objvec);
      if(abs(filtpt->at(filtindex)-objpt->at(objindex)) < abs(lowestDpt)) lowestDpt = filtpt->at(filtindex)-objpt->at(objindex);
    }
    if(abs(filteta->at(filtindex))<1.479) {
      if(lowestDeta!=9e9) premchfiltebDeta->Fill(lowestDeta);
      if(lowestDphi!=9e9) premchfiltebDphi->Fill(lowestDphi);
      if(lowestDr!=9e9) premchfiltebDr->Fill(lowestDr);
      if(lowestDpt!=9e9) premchfiltebDpt->Fill(lowestDpt);
    }
    else {
      if(lowestDeta!=9e9) premchfilteeDeta->Fill(lowestDeta);
      if(lowestDphi!=9e9) premchfilteeDphi->Fill(lowestDphi);
      if(lowestDr!=9e9) premchfilteeDr->Fill(lowestDr);
      if(lowestDpt!=9e9) premchfilteeDpt->Fill(lowestDpt);
    }
  }

  // Match with filter
  for(int findex=0; findex<filtidx.size(); findex++) {
    int filtindex = filtidx[findex];
    TVector3 filtvec;
    filtvec.SetPtEtaPhi(filtpt->at(filtindex), filteta->at(filtindex), filtphi->at(filtindex));
    for(int oindex=0; oindex<objidx.size(); oindex++) {
      int objindex = objidx[oindex];
      if(objectMatched[oindex]) continue;
      TVector3 objvec;
      objvec.SetPtEtaPhi(objpt->at(objindex), objeta->at(objindex), objphi->at(objindex));
      bool matchfound = false;
      // Matching condition
      if(abs(filteta->at(filtindex))<1.479) {
	if(abs(filtvec.DeltaR(objvec))<0.1) matchfound = true;
      }
      else {
	if(abs(filtvec.DeltaR(objvec))<0.075) matchfound = true;
      }
      if(matchfound) {
	objectMatched[oindex] = true;
	matchedObjectFilterPairs[findex] = make_pair(filtindex,objindex);
	break;
      }
    }
  }

  // After matching
  for(pair<int,int> matchedindex : matchedObjectFilterPairs) {
    int filtindex = matchedindex.first;
    int objindex = matchedindex.second;
    if(objindex==-1) continue;
    TVector3 filtvec;
    filtvec.SetPtEtaPhi(filtpt->at(filtindex), filteta->at(filtindex), filtphi->at(filtindex));
    TVector3 objvec;
    objvec.SetPtEtaPhi(objpt->at(objindex), objeta->at(objindex), objphi->at(objindex));
    if(abs(filteta->at(filtindex))<1.479) {
      aftmchfiltebDeta->Fill(filteta->at(filtindex)-objeta->at(objindex));
      aftmchfiltebDphi->Fill(filtvec.DeltaPhi(objvec));
      aftmchfiltebDr->Fill(filtvec.DeltaR(objvec));
      aftmchfiltebDpt->Fill(filtpt->at(filtindex)-objpt->at(objindex));
    }
    else {
      aftmchfilteeDeta->Fill(filteta->at(filtindex)-objeta->at(objindex));
      aftmchfilteeDphi->Fill(filtvec.DeltaPhi(objvec));
      aftmchfilteeDr->Fill(filtvec.DeltaR(objvec));
      aftmchfilteeDpt->Fill(filtpt->at(filtindex)-objpt->at(objindex));
    }
  }

  return matchedObjectFilterPairs;
  
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
  all1dhists.push_back(new TH1F(selection+"_pho_eta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"_pho_phi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_pho_seedtime","seed time / ns",500000,-25,25));
  all1dhists.push_back(new TH1F(selection+"_pho_smin","smin",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"_pho_smax","smax",10000,0,10));
}

// Function to add a set of histograms for angular matching between filters and objects
void robustanalyzer::addObjectFilterAngMatchhist(TString selection) {
  all1dhists.push_back(new TH1F(selection+"_filteb_obj_premch_Deta","#Delta#eta(filt, obj)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"_filteb_obj_premch_Dphi","#Delta#phi(filt, obj)",6600,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_filteb_obj_premch_Dr","#DeltaR(filt, obj)",7000,0,7));
  all1dhists.push_back(new TH1F(selection+"_filteb_obj_premch_Dpt","#Deltap_{T}(filt, obj)",200000,-1000,1000));
  all1dhists.push_back(new TH1F(selection+"_filtee_obj_premch_Deta","#Delta#eta(filt, obj)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"_filtee_obj_premch_Dphi","#Delta#phi(filt, obj)",6600,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_filtee_obj_premch_Dr","#DeltaR(filt, obj)",7000,0,7));
  all1dhists.push_back(new TH1F(selection+"_filtee_obj_premch_Dpt","#Deltap_{T}(filt, obj)",200000,-1000,1000));

  all1dhists.push_back(new TH1F(selection+"_filteb_obj_aftmch_Deta","#Delta#eta(filt, obj)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"_filteb_obj_aftmch_Dphi","#Delta#phi(filt, obj)",6600,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_filteb_obj_aftmch_Dr","#DeltaR(filt, obj)",7000,0,7));
  all1dhists.push_back(new TH1F(selection+"_filteb_obj_aftmch_Dpt","#Deltap_{T}(filt, obj)",200000,-1000,1000));
  all1dhists.push_back(new TH1F(selection+"_filtee_obj_aftmch_Deta","#Delta#eta(filt, obj)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"_filtee_obj_aftmch_Dphi","#Delta#phi(filt, obj)",6600,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"_filtee_obj_aftmch_Dr","#DeltaR(filt, obj)",7000,0,7));
  all1dhists.push_back(new TH1F(selection+"_filtee_obj_aftmch_Dpt","#Deltap_{T}(filt, obj)",200000,-1000,1000));
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

vector<int> robustanalyzer::getFiltMatchedPhoIndex(vector<pair<int,int>> matchedFilterPhotonIndices, double etamin, double etamax) {
  
  vector<int> phoIdx;
  for(pair<int,int> filtphoidx : matchedFilterPhotonIndices) {
    int phoindex = filtphoidx.second;
    if(phoindex!=-1) {
      if(abs(pho_eta->at(phoindex))>etamin && abs(pho_eta->at(phoindex))<etamax) {
	phoIdx.push_back(phoindex);
      }
    }
  }
  if(phoIdx.size()==0) phoIdx.push_back(-1);
  return phoIdx;
}
