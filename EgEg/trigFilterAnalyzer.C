/*
 * AUTHOR: Abanti Ranadhir Sahasransu - asahasra@cern.ch
 * The code now assumes exactly two gen electrons for a MC sample
 */

#include "trigFilterAnalyzer.hh"
#include <iostream>

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
trigFilterAnalyzer::trigFilterAnalyzer(TString filename, TString outfilename) {

  inputChain = new TChain("events");
  cout<<"Initializing for file: "<<filename<<endl;
  inputChain->Add(filename);

  inputChain->SetBranchAddress("HLT_DoublePhoton33CaloIdL", &HLT_DoublePhoton33_CaloIdL);
  inputChain->SetBranchAddress("dieg33_egcsusFiltn", &dieg33_egFiltN);
  inputChain->SetBranchAddress("dieg33_egcsusFilt_pt", &dieg33_egFiltPt);
  inputChain->SetBranchAddress("dieg33_egcsusFilt_eta", &dieg33_egFiltEta);
  inputChain->SetBranchAddress("dieg33_egcsusFilt_phi", &dieg33_egFiltPhi);

  inputChain->SetBranchAddress("HLT_DiPhoton10Time1ns", &HLT_DiPhoton10Time1ns);
  inputChain->SetBranchAddress("dipho10time1ns_egtimeusFiltn", &dipho10time1ns_egFiltN);
  inputChain->SetBranchAddress("dipho10time1ns_egtimeusFilt_pt", &dipho10time1ns_egFiltPt);
  inputChain->SetBranchAddress("dipho10time1ns_egtimeusFilt_eta", &dipho10time1ns_egFiltEta);
  inputChain->SetBranchAddress("dipho10time1ns_egtimeusFilt_phi", &dipho10time1ns_egFiltPhi);

  inputChain->SetBranchAddress("HLT_DiPhoton10sminlt0p16", &HLT_DiPhoton10sminlt0p16);
  inputChain->SetBranchAddress("dipho10sminlt0p16_egsminusFiltn", &dipho10sminlt0p16_egFiltN);
  inputChain->SetBranchAddress("dipho10sminlt0p16_egsminusFilt_pt", &dipho10sminlt0p16_egFiltPt);
  inputChain->SetBranchAddress("dipho10sminlt0p16_egsminusFilt_eta", &dipho10sminlt0p16_egFiltEta);
  inputChain->SetBranchAddress("dipho10sminlt0p16_egsminusFilt_phi", &dipho10sminlt0p16_egFiltPhi);

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
  
  outfile = new TFile(outfilename,"RECREATE");
}

// Fill the root file, close the root file, and handle deletions
trigFilterAnalyzer::~trigFilterAnalyzer() {

  inputChain->Delete();
  outfile->Write();
  outfile->Close();
}

// Analyzer for a single file
void trigFilterAnalyzer::analyzersinglefile(int splitCnt) { // Assume splitCnt to range from 0 to nCores

  int totEntries = inputChain->GetEntries();
  cout<<"Total number of entries: "<<totEntries<<endl;
  int nCores = 6; // Assume parallel processing over 7 cores where
  // there is a lesser no.of events in the last core
  int beginevent = splitCnt*(totEntries/nCores);
  int endevent = (splitCnt+1)*(totEntries/nCores);
  endevent = endevent<totEntries?endevent:totEntries; // Verfied that this logic to parallelize works

  // Define the histograms
  addgenhist("gennosel");
  addgenhist("gennosel_dieg33mch");
  addgenhist("gennosel_time1nsAsminlt0p16mch");
  addgenhist("genptgt10barsel");
  addgenhist("genptgt10barsel_dieg33mch");
  addgenhist("genptgt10barsel_time1nsAsminlt0p16mch");

  // Loop beginning on events
  for(unsigned int event=beginevent; event<endevent; event++) {

    inputChain->GetEntry(event);
    //if(event>10000) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    // vector of indices
    vector<int> genelpos;
    genelpos.push_back(-1);
    genelpos.push_back(-1);
    vector<int> gennoselegidx;
    vector<bool> gennoselegidxdec;
    vector<int> genptgt10barselegidx;
    vector<bool> genptgt10barselegidxdec;
    
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

    bool genptgt10barseleg = false;
    
    for(int genCtr:genelpos) {
      
      //gennoselegidx.push_back(genCtr);
      //gennoselegidxdec.push_back(true);

      genptgt10barseleg = true;
      genptgt10barseleg *= abs(genLepPid[genCtr])==11;
      genptgt10barseleg *= abs(genLepVz[genCtr])<320;
      genptgt10barseleg *= TMath::Sqrt(genLepVx[genCtr]*genLepVx[genCtr]+genLepVy[genCtr]*genLepVy[genCtr])<130;
      genptgt10barseleg *= abs(genLepPromptEta[genCtr])<1.479;
      genptgt10barseleg *= genLepPt[genCtr]>10;;
      if(genptgt10barseleg) {
	genptgt10barselegidx.push_back(genCtr);
	genptgt10barselegidxdec.push_back(true);
      }
	
    }
    //fillgenhistinevent("gennosel",gennoselegidx, gennoselegidxdec); // Verified that this is always 2 electrons
    fillgenhistinevent("genptgt10barsel",genptgt10barselegidx, genptgt10barselegidxdec);

    // Cross check for the trigger filters with trigger result
    if(HLT_DoublePhoton33_CaloIdL>0 && dieg33_egFiltN<2) throw "Type 1 Error!!! HLT_DoublePhoton33_CaloIdL true but not enough filter egamma";
    if(HLT_DoublePhoton33_CaloIdL==0 && dieg33_egFiltN>=2) throw "Type 2 Error!!! HLT_DoublePhoton33_CaloIdL false with required filter egamma objects";
    if(HLT_DiPhoton10Time1ns>0 && dipho10time1ns_egFiltN<2) throw "Type 1 Error!!! HLT_DiPhoton10Time1ns true but not enough filter egamma";
    if(HLT_DiPhoton10Time1ns==0 && dipho10time1ns_egFiltN>=2) throw "Type 2 Error!!! HLT_DiPhoton10Time1ns false with required filter egamma objects";
    if(HLT_DiPhoton10sminlt0p16>0 && dipho10sminlt0p16_egFiltN<2) throw "Type 1 Error!!! HLT_DiPhoton10sminlt0p16 true but not enough filter egamma";
    if(HLT_DiPhoton10sminlt0p16==0 && dipho10sminlt0p16_egFiltN>=2) throw "Type 2 Error!!! HLT_DiPhoton10sminlt0p16 false with required filter egamma objects";

    // Vector for gen matching eg filters
    vector< pair<int, int> > dieg33mchgennoselidx;
    vector< pair<int, int> > time1nsmchgennoselidx;
    vector< pair<int, int> > sminlt0p16mchgennoselidx;
    vector<bool> gennosel_dieg33mchegidxdec = gennoselegidxdec;
    vector<bool> gennosel_time1nsAsminlt0p16mchegidxdec = gennoselegidxdec;
    vector< pair<int, int> > dieg33mchgenptgt10barselidx;
    vector< pair<int, int> > time1nsmchgenptgt10barselidx;
    vector< pair<int, int> > sminlt0p16mchgenptgt10barselidx;
    vector<bool> genptgt10barsel_dieg33mchegidxdec = genptgt10barselegidxdec;
    vector<bool> genptgt10barsel_time1nsAsminlt0p16mchegidxdec = genptgt10barselegidxdec;

    // Find the filters matching to gen
    if(HLT_DoublePhoton33_CaloIdL>0) {
      //dieg33mchgennoselidx = doGenMatching(gennoselegidx, dieg33_egFiltN, &dieg33_egFiltPt[0], &dieg33_egFiltEta[0], &dieg33_egFiltPhi[0]);
      dieg33mchgenptgt10barselidx = doGenMatching(genptgt10barselegidx, dieg33_egFiltN, &dieg33_egFiltPt[0], &dieg33_egFiltEta[0], &dieg33_egFiltPhi[0]);
    }
    if(HLT_DiPhoton10Time1ns>0) {
      //time1nsmchgennoselidx = doGenMatching(gennoselegidx, dipho10time1ns_egFiltN, &dipho10time1ns_egFiltPt[0], &dipho10time1ns_egFiltEta[0], &dipho10time1ns_egFiltPhi[0]);
      time1nsmchgenptgt10barselidx = doGenMatching(genptgt10barselegidx, dipho10time1ns_egFiltN, &dipho10time1ns_egFiltPt[0], &dipho10time1ns_egFiltEta[0], &dipho10time1ns_egFiltPhi[0]);
    }
    if(HLT_DiPhoton10sminlt0p16>0) {
      //sminlt0p16mchgennoselidx = doGenMatching(gennoselegidx, dipho10sminlt0p16_egFiltN, &dipho10sminlt0p16_egFiltPt[0], &dipho10sminlt0p16_egFiltEta[0], &dipho10sminlt0p16_egFiltPhi[0]);
      sminlt0p16mchgenptgt10barselidx = doGenMatching(genptgt10barselegidx, dipho10sminlt0p16_egFiltN, &dipho10sminlt0p16_egFiltPt[0], &dipho10sminlt0p16_egFiltEta[0], &dipho10sminlt0p16_egFiltPhi[0]);
    }

    // Histograms after gen matching
    //combineFiltmchForGen(gennoselegidx, dieg33mchgennoselidx, &gennosel_dieg33mchegidxdec);
    //fillgenhistinevent("gennosel_dieg33mch",gennoselegidx, gennosel_dieg33mchegidxdec);

    //combineFiltmchForGen(gennoselegidx, time1nsmchgennoselidx, &gennosel_time1nsAsminlt0p16mchegidxdec);
    //combineFiltmchForGen(gennoselegidx, sminlt0p16mchgennoselidx, &gennosel_time1nsAsminlt0p16mchegidxdec);
    //fillgenhistinevent("gennosel_time1nsAsminlt0p16mch",gennoselegidx, gennosel_time1nsAsminlt0p16mchegidxdec);
    
    combineFiltmchForGen(genptgt10barselegidx, dieg33mchgenptgt10barselidx, &genptgt10barsel_dieg33mchegidxdec);
    fillgenhistinevent("genptgt10barsel_dieg33mch",genptgt10barselegidx, genptgt10barsel_dieg33mchegidxdec);

    combineFiltmchForGen(genptgt10barselegidx, time1nsmchgenptgt10barselidx, &genptgt10barsel_time1nsAsminlt0p16mchegidxdec);
    combineFiltmchForGen(genptgt10barselegidx, sminlt0p16mchgenptgt10barselidx, &genptgt10barsel_time1nsAsminlt0p16mchegidxdec);
    fillgenhistinevent("genptgt10barsel_time1nsAsminlt0p16mch",genptgt10barselegidx, genptgt10barsel_time1nsAsminlt0p16mchegidxdec);
    
    // Clear all the vectors
    genelpos.clear();
    gennoselegidx.clear();
    gennoselegidxdec.clear();
    dieg33mchgennoselidx.clear();
    time1nsmchgennoselidx.clear();
    sminlt0p16mchgennoselidx.clear();
    gennosel_dieg33mchegidxdec.clear();
    gennosel_time1nsAsminlt0p16mchegidxdec.clear();
    genptgt10barselegidx.clear();
    genptgt10barselegidxdec.clear();
    dieg33mchgenptgt10barselidx.clear();
    time1nsmchgenptgt10barselidx.clear();
    sminlt0p16mchgenptgt10barselidx.clear();
    genptgt10barsel_dieg33mchegidxdec.clear();
    genptgt10barsel_time1nsAsminlt0p16mchegidxdec.clear();
    
  } // End of event loop

  cout<<totEntries<<endl;

}

vector< pair<int,int> > trigFilterAnalyzer::doGenMatching(vector<int> genidx, int egnum, double* egpt, double* egeta, double* egphi) {
  
  vector< pair<int, int> > filtmchgen;

  for(unsigned int egctr=0; egctr<egnum; egctr++) {

    for(int gen : genidx) {

      // Gen match with prompt equivalent of gen eta and gen phi
      bool foundfiltmchgen = false;
      double diffeta = abs((*(egeta+egctr))-genLepPromptEta[gen]); 
      TLorentzVector veceg, vecgen, vecpromptgen;
      vecgen.SetPtEtaPhiM(genLepPt[gen],genLepEta[gen],genLepPhi[gen],0.0005);
      vecpromptgen.SetPtEtaPhiM(vecgen.P()*TMath::Sin(2*TMath::ATan(TMath::Exp(-genLepPromptEta[gen]))),genLepPromptEta[gen],genLepPromptPhi[gen],0.0005);
      veceg.SetPtEtaPhiM((*(egpt+egctr)), (*(egeta+egctr)), (*(egphi+egctr)), 0.0005);
      double qdiffphi = (genLepPid[gen]/abs(genLepPid[gen]))*(vecpromptgen.DeltaPhi(veceg));
      // Condition for gen matching
      if(abs((*(egeta+egctr)))<1.479) {
	if(diffeta<0.1 && qdiffphi<0.15 && qdiffphi>-0.25) {
	  foundfiltmchgen = true;
	}
      }
      else {
	if(diffeta<0.05 && qdiffphi<0.1 && qdiffphi>-0.15) {
	  foundfiltmchgen = true;
	}
      } // End of gen matching

      // Check if gen already has a match
      bool genhasmch = false;
      for(unsigned int mchctr=0; mchctr<filtmchgen.size(); mchctr++) {
	if(filtmchgen[mchctr].first == gen) {
	  genhasmch = true;
	  break;
	}
      }
      if(genhasmch) break;

      if(foundfiltmchgen && !genhasmch) {
	filtmchgen.push_back(make_pair(gen, egctr));
	break;
      }
      
    } // End of loop on gen els
    
  } // End of loop on eg filters
    
  return (filtmchgen);
}

// Function to fill a set of histograms for gen particles
void trigFilterAnalyzer::fillgenhistinevent(TString selection, vector<int> egidx, vector<bool> egdecidx) {

  if(egidx.size() != egdecidx.size()) throw "Incompatible sizes for gen electron vector and gen decision vector";
  
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

  int leadegpos=-1, subleadegpos=-1, dummy=-1;
  for(unsigned int egidxpos=0; egidxpos<egidx.size(); egidxpos++) {

    // Find the lead and sublead electron by pt
    if(leadegpos==-1) {
      leadegpos = egidxpos;
    }
    else if(genLepPt[egidx[leadegpos]]<genLepPt[egidx[egidxpos]]) {
      subleadegpos = leadegpos;
      leadegpos = egidxpos;
    }
    else if(subleadegpos==-1) {
      subleadegpos = egidxpos;
    }
    else if(genLepPt[egidx[subleadegpos]]<genLepPt[egidx[egidxpos]]) {
      subleadegpos = egidxpos;
    }
    else {
      dummy = egidx[egidxpos];
    }
    
    if(egdecidx[egidxpos]) {
      TVector3 el, elmom;
      el.SetPtEtaPhi(genLepPt[egidx[egidxpos]], genLepEta[egidx[egidxpos]], genLepPhi[egidx[egidxpos]]);
      elmom.SetPtEtaPhi(genLepMomPt[egidx[egidxpos]], genLepMomEta[egidx[egidxpos]], genLepMomPhi[egidx[egidxpos]]);
      egmompid->Fill(genLepMomPid[egidx[egidxpos]]);
      pt->Fill(genLepPt[egidx[egidxpos]]);
      eta->Fill(genLepEta[egidx[egidxpos]]);
      phi->Fill(genLepPhi[egidx[egidxpos]]);
      prompteta->Fill(genLepPromptEta[egidx[egidxpos]]);
      promptphi->Fill(genLepPromptPhi[egidx[egidxpos]]);
      vx->Fill(genLepVx[egidx[egidxpos]]);
      vy->Fill(genLepVy[egidx[egidxpos]]);
      vz->Fill(genLepVz[egidx[egidxpos]]);
      pvxh->Fill(pvx[egidx[egidxpos]]);
      pvyh->Fill(pvy[egidx[egidxpos]]);
      pvzh->Fill(pvz[egidx[egidxpos]]);
      deltaetamom->Fill(genLepMomEta[egidx[egidxpos]]-genLepEta[egidx[egidxpos]]);
      deltaphimom->Fill(elmom.DeltaPhi(el));
      deltaRmom->Fill(elmom.DeltaR(el));
      gend0->Fill(genLepDxy[egidx[egidxpos]]);
      log10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[egidx[egidxpos]])));
      genlxy->Fill(genLepLxy[egidx[egidxpos]]);
      log10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[egidx[egidxpos]])));
      if(genLepTimeAct[egidx[egidxpos]]>(-1e9)) {
	t1->Fill(genLepTimeAct[egidx[egidxpos]]);
	t0->Fill(genLepTimeLight[egidx[egidxpos]]);
	t1mt0->Fill(genLepTimeDiff[egidx[egidxpos]]);
      }
      genelmult++;
    }
  }

  if(leadegpos!=-1 && egdecidx[leadegpos]) {
    TVector3 el, elmom;
    el.SetPtEtaPhi(genLepPt[egidx[leadegpos]], genLepEta[egidx[leadegpos]], genLepPhi[egidx[leadegpos]]);
    elmom.SetPtEtaPhi(genLepMomPt[egidx[leadegpos]], genLepMomEta[egidx[leadegpos]], genLepMomPhi[egidx[leadegpos]]);
    leadegpt->Fill(genLepPt[egidx[leadegpos]]);
    leadegeta->Fill(genLepEta[egidx[leadegpos]]);
    leadegphi->Fill(genLepPhi[egidx[leadegpos]]);
    leaddeltaetamom->Fill(genLepMomEta[egidx[leadegpos]]-genLepEta[egidx[leadegpos]]);
    leaddeltaphimom->Fill(elmom.DeltaPhi(el));
    leaddeltaRmom->Fill(elmom.DeltaR(el));
    leadegd0->Fill(genLepDxy[egidx[leadegpos]]);
    leadeglog10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[egidx[leadegpos]])));
    leadeglxy->Fill(genLepLxy[egidx[leadegpos]]);
    leadeglog10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[egidx[leadegpos]])));
  }

  if(subleadegpos!=-1 && egdecidx[subleadegpos]) {
    TVector3 el, elmom;
    el.SetPtEtaPhi(genLepPt[egidx[subleadegpos]], genLepEta[egidx[subleadegpos]], genLepPhi[egidx[subleadegpos]]);
    elmom.SetPtEtaPhi(genLepMomPt[egidx[subleadegpos]], genLepMomEta[egidx[subleadegpos]], genLepMomPhi[egidx[subleadegpos]]);
    leadegpt->Fill(genLepPt[egidx[subleadegpos]]);
    leadegeta->Fill(genLepEta[egidx[subleadegpos]]);
    leadegphi->Fill(genLepPhi[egidx[subleadegpos]]);
    leaddeltaetamom->Fill(genLepMomEta[egidx[subleadegpos]]-genLepEta[egidx[subleadegpos]]);
    leaddeltaphimom->Fill(elmom.DeltaPhi(el));
    leaddeltaRmom->Fill(elmom.DeltaR(el));
    leadegd0->Fill(genLepDxy[egidx[subleadegpos]]);
    leadeglog10d0->Fill(TMath::Log10(TMath::Abs(genLepDxy[egidx[subleadegpos]])));
    leadeglxy->Fill(genLepLxy[egidx[subleadegpos]]);
    leadeglog10lxy->Fill(TMath::Log10(TMath::Abs(genLepLxy[egidx[subleadegpos]])));
  }

  egmult->Fill(genelmult);
  egidx.clear();
}

// Function to add a set of histograms for a gen particles
void trigFilterAnalyzer::addgenhist(TString selection) {
  
  all1dhists.push_back(new TH1F(selection+"geneg_egmult","gen N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"geneg_egmompid","gen mom pdg id",100,-50,50));
  all1dhists.push_back(new TH1F(selection+"geneg_pt","gen e/#gamma p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"geneg_eta","gen e/#gamma #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_phi","gen e/#gamma #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"geneg_prompteta","gen e/#gamma corrected #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"geneg_promptphi","gen e/#gamma corrected #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"geneg_vx","gen e/#gamma v_{x} / cm",50000,-250,250));
  all1dhists.push_back(new TH1F(selection+"geneg_vy","gen e/#gamma v_{y} / cm",50000,-250,250));
  all1dhists.push_back(new TH1F(selection+"geneg_vz","gen e/#gamma v_{z} / cm",50000,-500,500));
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

// Combine the different filter matches for each gen electron
void trigFilterAnalyzer::combineFiltmchForGen(vector<int> genidx, vector< pair<int,int> >filtgenmchpair, vector<bool> *decidx) {

  for( unsigned int genctr=0; genctr<genidx.size(); genctr++) {
    bool foundfiltmch = false;
    for(pair<int,int> filtgenmch : filtgenmchpair) {
      if(filtgenmch.first==genidx[genctr] && filtgenmch.second!=-1) {
	foundfiltmch = true;
	break;
      }
    }
    decidx->at(genctr) = foundfiltmch;
  }
}
