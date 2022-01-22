#include "DDM10_robustanalyzer.hh"
#include <iostream>
#include <numeric>

#include "TMath.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

using namespace std;

// Initialize and open the root file in the constructor
DDM10_robustanalyzer::DDM10_robustanalyzer(TString filename, TString outfilename, bool issimu){

  outFileName = outfilename;
  
  isMC = issimu;
  
  inputChain = new TChain("events");
  cout<<"Initializing for file: "<<filename<<endl;
  inputChain->Add(filename);

  inputChain->SetBranchAddress("bsx", &bsx);
  inputChain->SetBranchAddress("bsy", &bsy);
  inputChain->SetBranchAddress("bsz", &bsz);
  inputChain->SetBranchAddress("trig_DoubleMu33NoFiltersNoVtxDisplaced", &trigDoubleMu33);
  inputChain->SetBranchAddress("trig_DoubleL2Mu23NoVtx2Cha", &trigDoubleL2Mu23NV2Cha);
  inputChain->SetBranchAddress("trig_DoubleL2Mu23NoVtx2ChaCosmicSeed", &trigDoubleL2Mu23NV2ChaCS);
  inputChain->SetBranchAddress("trig_DoubleL3Mu10NoVtx_Displaced", &trigDoubleL3Mu10NoVtxDisplaced);
  inputChain->SetBranchAddress("trig_DoubleL2Mu10NoVtx_2Cha_PromptL3Mu0Veto", &trigDoubleL2Mu10NoVtx2ChaPromptL3Mu0Veto);
  inputChain->SetBranchAddress("l3dim33Filtn", &l3dim33FiltN);
  inputChain->SetBranchAddress("l3dim33Filt_pt", &l3dim33FiltPt);
  inputChain->SetBranchAddress("l3dim33Filt_eta", &l3dim33FiltEta);
  inputChain->SetBranchAddress("l3dim33Filt_phi", &l3dim33FiltPhi);
  inputChain->SetBranchAddress("l2dim23Filtn", &l2dim23FiltN);
  inputChain->SetBranchAddress("l2dim23Filt_pt", &l2dim23FiltPt);
  inputChain->SetBranchAddress("l2dim23Filt_eta", &l2dim23FiltEta);
  inputChain->SetBranchAddress("l2dim23Filt_phi", &l2dim23FiltPhi);
  inputChain->SetBranchAddress("l2dim23csFiltn", &l2dim23csFiltN);
  inputChain->SetBranchAddress("l2dim23csFilt_pt", &l2dim23csFiltPt);
  inputChain->SetBranchAddress("l2dim23csFilt_eta", &l2dim23csFiltEta);
  inputChain->SetBranchAddress("l2dim23csFilt_phi", &l2dim23csFiltPhi);
  inputChain->SetBranchAddress("l3ddm10Filtn", &l3ddm10FiltN);
  inputChain->SetBranchAddress("l3ddm10Filt_pt", &l3ddm10FiltPt);
  inputChain->SetBranchAddress("l3ddm10Filt_eta", &l3ddm10FiltEta);
  inputChain->SetBranchAddress("l3ddm10Filt_phi", &l3ddm10FiltPhi);
  inputChain->SetBranchAddress("l2ddm10Filtn", &l2ddm10FiltN);
  inputChain->SetBranchAddress("l2ddm10Filt_pt", &l2ddm10FiltPt);
  inputChain->SetBranchAddress("l2ddm10Filt_eta", &l2ddm10FiltEta);
  inputChain->SetBranchAddress("l2ddm10Filt_phi", &l2ddm10FiltPhi);

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
  
  outfile = new TFile(outfilename+".root","RECREATE");
}

// Fill the root file, close the root file, and handle deletions
DDM10_robustanalyzer::~DDM10_robustanalyzer() {

  inputChain->Delete();
}

// Analyzer for a single file
void DDM10_robustanalyzer::analyzer() {

  int totEntries = inputChain->GetEntries();
  cout<<"Total number of entries: "<<totEntries<<endl;

  // Count events passing certain selections
  int evtgenmumu=0, evtgenmumubasicsel=0;
  int evttrigDoubleMu33=0, evttrigDoubleL2Mu23NV2Cha=0, evttrigDoubleL2Mu23NV2ChaCS=0, evttrigDoubleL3Mu10NoVtxDisplaced=0, evttrigDoubleL2Mu10NoVtx2ChaPromptL3Mu0Veto=0;
  int nosel=0, basicsel=0, sel2=0, sel3=0, sel20=0, sel30=0, sel40=0, sel100=0;

  // Add gen hist
  addhistGenPart("nosel");
  addhistGenPart("basicsel");

  // Vector to store the gen muon index
  std::vector<int> genmuidx, genmubasicselidx;
  
  // Loop beginning on events
  for(unsigned int event=0; event<totEntries; event++) {

    inputChain->GetEntry(event);
    //if(event>20) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    // Loop on genParticles to count the number of genevents
    bool genmu=true, genmubasicsel=true;
    for(unsigned int i=0; i<genLepN; i++) {

      genmu = true;
      genmu *= abs(genLepPid[i])==13;
      genmu *= abs(genLepMomPid[i])==9000007 || abs(genLepMomPid[i])==6000113;
      if(genmu)	genmuidx.push_back(i);

      genmubasicsel = true;
      genmubasicsel *= abs(genLepPid[i])==13;
      genmubasicsel *= abs(genLepMomPid[i])==9000007 || abs(genLepMomPid[i])==6000113;
      genmubasicsel *= genLepPt[i]>10;
      genmubasicsel *= abs(genLepEta[i])<2.5;
      if(genmubasicsel)	genmubasicselidx.push_back(i);
    }
    
    if(genmuidx.size()>=2) {
      evtgenmumu++;
      fillhistpereventGenPart("nosel",genmuidx);
    }
    if(genmubasicselidx.size()>=2) {
      evtgenmumubasicsel++;
      fillhistpereventGenPart("basicsel",genmubasicselidx);
    }
    if(trigDoubleMu33>0) evttrigDoubleMu33++;
    if(trigDoubleL2Mu23NV2Cha>0) evttrigDoubleL2Mu23NV2Cha++;
    if(trigDoubleL2Mu23NV2ChaCS>0) evttrigDoubleL2Mu23NV2ChaCS++;
    if(trigDoubleL3Mu10NoVtxDisplaced) evttrigDoubleL3Mu10NoVtxDisplaced++;
    if(trigDoubleL2Mu10NoVtx2ChaPromptL3Mu0Veto) evttrigDoubleL2Mu10NoVtx2ChaPromptL3Mu0Veto++;

    // Clear all defined vectors
    genmuidx.clear();
    genmubasicselidx.clear();
  } // End of loop on events

  cout<<"Events with atleast 2 muons - "<<evtgenmumu<<endl;
  cout<<"Events with atleast 2 muons, pT>10, eta<2.5 - "<<evtgenmumubasicsel<<endl;
  cout<<"Events firing trigger : DoubleMu33NoFiltersNoVtxDisplaced - "<<evttrigDoubleMu33<<endl;
  cout<<"Events firing trigger : DoubleL2Mu23NoVtx2Cha - "<<evttrigDoubleL2Mu23NV2Cha<<endl;
  cout<<"Events firing trigger : DoubleL2Mu23NoVtx2ChaCosmicSeed - "<<evttrigDoubleL2Mu23NV2ChaCS<<endl;
  cout<<"Events firing trigger : DoubleL3Mu10NoVtxDisplaced - "<<evttrigDoubleL3Mu10NoVtxDisplaced<<endl;
  cout<<"Events firing trigger : DoubleL2Mu10NoVtx2ChaPromptL3Mu0Veto - "<<evttrigDoubleL2Mu10NoVtx2ChaPromptL3Mu0Veto<<endl;

  outfile->Write();
  outfile->Close();
}

// Function to fill a set of histograms in the event
void DDM10_robustanalyzer::fillhistpereventGenPart(TString selection, std::vector<int> genmuidx) {
  
  TH1F* mumult = (TH1F*) outfile->Get(selection+"gen_mumult");
  TH1F* mupt = (TH1F*) outfile->Get(selection+"gen_mupt");
  TH1F* mueta = (TH1F*) outfile->Get(selection+"gen_mueta");
  TH1F* muphi = (TH1F*) outfile->Get(selection+"gen_muphi");
  TH1F* mud0 = (TH1F*) outfile->Get(selection+"gen_mud0");
  TH1F* mulog10d0 = (TH1F*) outfile->Get(selection+"gen_mulog10d0");
  TH1F* ordptsubleadmupt = (TH1F*) outfile->Get(selection+"gen_ordptsubleadmupt");
  TH1F* ordd0subleadmud0 = (TH1F*) outfile->Get(selection+"gen_ordd0subleadmud0");
  TH1F* ordd0subleadmulog10d0 = (TH1F*) outfile->Get(selection+"gen_ordd0subleadmulog10d0");

  TH1F* l3dim33genmchpt = (TH1F*) outfile->Get(selection+"l3dim33genmch_mupt");
  TH1F* l3dim33genmcheta = (TH1F*) outfile->Get(selection+"l3dim33genmch_mueta");
  TH1F* l3dim33genmchphi = (TH1F*) outfile->Get(selection+"l3dim33genmch_muphi");
  TH1F* l3dim33genmchd0 = (TH1F*) outfile->Get(selection+"l3dim33genmch_mud0");
  TH1F* l3dim33genmchlog10d0 = (TH1F*) outfile->Get(selection+"l3dim33genmch_mulog10d0");
  TH1F* l3dim33genmchordptsubleadmupt = (TH1F*) outfile->Get(selection+"l3dim33genmch_ordptsubleadmupt");
  TH1F* l3dim33genmchordd0subleadmud0 = (TH1F*) outfile->Get(selection+"l3dim33genmch_ordd0subleadmud0");
  TH1F* l3dim33genmchordd0subleadmulog10d0 = (TH1F*) outfile->Get(selection+"l3dim33genmch_ordd0subleadmulog10d0");

  TH1F* l2dim23genmchpt = (TH1F*) outfile->Get(selection+"l2dim23genmch_mupt");
  TH1F* l2dim23genmcheta = (TH1F*) outfile->Get(selection+"l2dim23genmch_mueta");
  TH1F* l2dim23genmchphi = (TH1F*) outfile->Get(selection+"l2dim23genmch_muphi");
  TH1F* l2dim23genmchd0 = (TH1F*) outfile->Get(selection+"l2dim23genmch_mud0");
  TH1F* l2dim23genmchlog10d0 = (TH1F*) outfile->Get(selection+"l2dim23genmch_mulog10d0");
  TH1F* l2dim23genmchordptsubleadmupt = (TH1F*) outfile->Get(selection+"l2dim23genmch_ordptsubleadmupt");
  TH1F* l2dim23genmchordd0subleadmud0 = (TH1F*) outfile->Get(selection+"l2dim23genmch_ordd0subleadmud0");
  TH1F* l2dim23genmchordd0subleadmulog10d0 = (TH1F*) outfile->Get(selection+"l2dim23genmch_ordd0subleadmulog10d0");

  TH1F* l2dim23csgenmchpt = (TH1F*) outfile->Get(selection+"l2dim23csgenmch_mupt");
  TH1F* l2dim23csgenmcheta = (TH1F*) outfile->Get(selection+"l2dim23csgenmch_mueta");
  TH1F* l2dim23csgenmchphi = (TH1F*) outfile->Get(selection+"l2dim23csgenmch_muphi");
  TH1F* l2dim23csgenmchd0 = (TH1F*) outfile->Get(selection+"l2dim23csgenmch_mud0");
  TH1F* l2dim23csgenmchlog10d0 = (TH1F*) outfile->Get(selection+"l2dim23csgenmch_mulog10d0");
  TH1F* l2dim23csgenmchordptsubleadmupt = (TH1F*) outfile->Get(selection+"l2dim23csgenmch_ordptsubleadmupt");
  TH1F* l2dim23csgenmchordd0subleadmud0 = (TH1F*) outfile->Get(selection+"l2dim23csgenmch_ordd0subleadmud0");
  TH1F* l2dim23csgenmchordd0subleadmulog10d0 = (TH1F*) outfile->Get(selection+"l2dim23csgenmch_ordd0subleadmulog10d0");

  TH1F* l3ddm10genmchpt = (TH1F*) outfile->Get(selection+"l3ddm10genmch_mupt");
  TH1F* l3ddm10genmcheta = (TH1F*) outfile->Get(selection+"l3ddm10genmch_mueta");
  TH1F* l3ddm10genmchphi = (TH1F*) outfile->Get(selection+"l3ddm10genmch_muphi");
  TH1F* l3ddm10genmchd0 = (TH1F*) outfile->Get(selection+"l3ddm10genmch_mud0");
  TH1F* l3ddm10genmchlog10d0 = (TH1F*) outfile->Get(selection+"l3ddm10genmch_mulog10d0");
  TH1F* l3ddm10genmchordptsubleadmupt = (TH1F*) outfile->Get(selection+"l3ddm10genmch_ordptsubleadmupt");
  TH1F* l3ddm10genmchordd0subleadmud0 = (TH1F*) outfile->Get(selection+"l3ddm10genmch_ordd0subleadmud0");
  TH1F* l3ddm10genmchordd0subleadmulog10d0 = (TH1F*) outfile->Get(selection+"l3ddm10genmch_ordd0subleadmulog10d0");

  TH1F* l2ddm10genmchpt = (TH1F*) outfile->Get(selection+"l2ddm10genmch_mupt");
  TH1F* l2ddm10genmcheta = (TH1F*) outfile->Get(selection+"l2ddm10genmch_mueta");
  TH1F* l2ddm10genmchphi = (TH1F*) outfile->Get(selection+"l2ddm10genmch_muphi");
  TH1F* l2ddm10genmchd0 = (TH1F*) outfile->Get(selection+"l2ddm10genmch_mud0");
  TH1F* l2ddm10genmchlog10d0 = (TH1F*) outfile->Get(selection+"l2ddm10genmch_mulog10d0");
  TH1F* l2ddm10genmchordptsubleadmupt = (TH1F*) outfile->Get(selection+"l2ddm10genmch_ordptsubleadmupt");
  TH1F* l2ddm10genmchordd0subleadmud0 = (TH1F*) outfile->Get(selection+"l2ddm10genmch_ordd0subleadmud0");
  TH1F* l2ddm10genmchordd0subleadmulog10d0 = (TH1F*) outfile->Get(selection+"l2ddm10genmch_ordd0subleadmulog10d0");

  TH1F* mul3dim33deta = (TH1F*) outfile->Get(selection+"genmul3dim33_deta");
  TH1F* mul3dim33dphi = (TH1F*) outfile->Get(selection+"genmul3dim33_dphi");
  TH1F* mul2dim23deta = (TH1F*) outfile->Get(selection+"genmul2dim23_deta");
  TH1F* mul2dim23dphi = (TH1F*) outfile->Get(selection+"genmul2dim23_dphi");
  TH1F* mul2dim23csdeta = (TH1F*) outfile->Get(selection+"genmul2dim23cs_deta");
  TH1F* mul2dim23csdphi = (TH1F*) outfile->Get(selection+"genmul2dim23cs_dphi");
  TH1F* mul3ddm10deta = (TH1F*) outfile->Get(selection+"genmul3ddm10_deta");
  TH1F* mul3ddm10dphi = (TH1F*) outfile->Get(selection+"genmul3ddm10_dphi");
  TH1F* mul2ddm10deta = (TH1F*) outfile->Get(selection+"genmul2ddm10_deta");
  TH1F* mul2ddm10dphi = (TH1F*) outfile->Get(selection+"genmul2ddm10_dphi");

  // Make a vector for d0 values
  std::vector<double> mud0vec;
  // Make a vector for gen matching decisions
  std::vector<bool> l3dim33genmch;
  std::vector<bool> l2dim23genmch;
  std::vector<bool> l2dim23csgenmch;
  std::vector<bool> l3ddm10genmch;
  std::vector<bool> l2ddm10genmch;

  // Fill histograms
  mumult->Fill(genmuidx.size());
  for(unsigned int i=0; i<genmuidx.size(); i++) {
    mupt->Fill(genLepPt[genmuidx[i]]);
    mueta->Fill(genLepEta[genmuidx[i]]);
    muphi->Fill(genLepPhi[genmuidx[i]]);
    TLorentzVector genmu;
    genmu.SetPtEtaPhiM(genLepPt[genmuidx[i]],genLepEta[genmuidx[i]],genLepPhi[genmuidx[i]],0.106);
    double d0 = (genLepVx[genmuidx[i]]-bsx)*genmu.Py()-(genLepVy[genmuidx[i]]-bsy)*genmu.Px();
    d0 /= genLepPt[genmuidx[i]];
    mud0->Fill(d0);
    mulog10d0->Fill(TMath::Log10(TMath::Abs(d0)));
    mud0vec.push_back(d0);

    // Start of individual trigger efficiency sections
    // For HLT_DoubleMu33No....
    if(trigDoubleMu33>0) {
      double mindeta=1000, mindphi=1000;
      for(unsigned int j=0; j<l3dim33FiltN; j++) {
	if(abs(genLepEta[genmuidx[i]]-l3dim33FiltEta[j])<abs(mindeta)) mindeta = genLepEta[genmuidx[i]]-l3dim33FiltEta[j];
	if(abs(genLepPhi[genmuidx[i]]-l3dim33FiltPhi[j])<abs(mindphi)) mindphi = genLepPhi[genmuidx[i]]-l3dim33FiltPhi[j];
      }
      if(l3dim33FiltN>0) {
	mul3dim33deta->Fill(mindeta);
	mul3dim33dphi->Fill(mindphi);
	if(abs(mindeta)<0.04) {
	  l3dim33genmchpt->Fill(genLepPt[genmuidx[i]]);
	  l3dim33genmcheta->Fill(genLepEta[genmuidx[i]]);
	  l3dim33genmchphi->Fill(genLepPhi[genmuidx[i]]);
	  l3dim33genmchd0->Fill(d0);
	  l3dim33genmchlog10d0->Fill(TMath::Log10(TMath::Abs(d0)));
	  l3dim33genmch.push_back(true);
	}
	else {
	  l3dim33genmch.push_back(false);
	}
      }
      else {
	l3dim33genmch.push_back(false);
      }
    }
    else {
      l3dim33genmch.push_back(false);
    }

    // For HLT_DoubleL2Mu23No....
    if(trigDoubleL2Mu23NV2Cha>0) {
      double mindeta=1000, mindphi=1000;
      for(unsigned int j=0; j<l2dim23FiltN; j++) {
	if(abs(genLepEta[genmuidx[i]]-l2dim23FiltEta[j])<abs(mindeta)) mindeta = genLepEta[genmuidx[i]]-l2dim23FiltEta[j];
	if(abs(genLepPhi[genmuidx[i]]-l2dim23FiltPhi[j])<abs(mindphi)) mindphi = genLepPhi[genmuidx[i]]-l2dim23FiltPhi[j];
      }
      if(l2dim23FiltN>0) {
	mul2dim23deta->Fill(mindeta);
	mul2dim23dphi->Fill(mindphi);
	if(abs(mindeta)<0.04) {
	  l2dim23genmchpt->Fill(genLepPt[genmuidx[i]]);
	  l2dim23genmcheta->Fill(genLepEta[genmuidx[i]]);
	  l2dim23genmchphi->Fill(genLepPhi[genmuidx[i]]);
	  l2dim23genmchd0->Fill(d0);
	  l2dim23genmchlog10d0->Fill(TMath::Log10(TMath::Abs(d0)));
	  l2dim23genmch.push_back(true);
	}
	else {
	  l2dim23genmch.push_back(false);
	}
      }
      else {
	l2dim23genmch.push_back(false);
      }
    }
    else {
      l2dim23genmch.push_back(false);
    }

    // For HLT_DoubleL2Mu23CSNo....
    if(trigDoubleL2Mu23NV2ChaCS>0) {
      double mindeta=1000, mindphi=1000;
      for(unsigned int j=0; j<l2dim23csFiltN; j++) {
	if(abs(genLepEta[genmuidx[i]]-l2dim23csFiltEta[j])<abs(mindeta)) mindeta = genLepEta[genmuidx[i]]-l2dim23csFiltEta[j];
	if(abs(genLepPhi[genmuidx[i]]-l2dim23csFiltPhi[j])<abs(mindphi)) mindphi = genLepPhi[genmuidx[i]]-l2dim23csFiltPhi[j];
      }
      if(l2dim23csFiltN>0) {
	mul2dim23csdeta->Fill(mindeta);
	mul2dim23csdphi->Fill(mindphi);
	if(abs(mindeta)<0.04) {
	  l2dim23csgenmchpt->Fill(genLepPt[genmuidx[i]]);
	  l2dim23csgenmcheta->Fill(genLepEta[genmuidx[i]]);
	  l2dim23csgenmchphi->Fill(genLepPhi[genmuidx[i]]);
	  l2dim23csgenmchd0->Fill(d0);
	  l2dim23csgenmchlog10d0->Fill(TMath::Log10(TMath::Abs(d0)));
	  l2dim23csgenmch.push_back(true);
	}
	else {
	  l2dim23csgenmch.push_back(false);
	}
      }
      else {
	l2dim23csgenmch.push_back(false);
      }
    }
    else {
      l2dim23csgenmch.push_back(false);
    }

    // For HLT_DoubleMu10No....
    if(trigDoubleL3Mu10NoVtxDisplaced>0) {
      double mindeta=1000, mindphi=1000;
      for(unsigned int j=0; j<l3ddm10FiltN; j++) {
	if(abs(genLepEta[genmuidx[i]]-l3ddm10FiltEta[j])<abs(mindeta)) mindeta = genLepEta[genmuidx[i]]-l3ddm10FiltEta[j];
	if(abs(genLepPhi[genmuidx[i]]-l3ddm10FiltPhi[j])<abs(mindphi)) mindphi = genLepPhi[genmuidx[i]]-l3ddm10FiltPhi[j];
      }
      if(l3ddm10FiltN>0) {
	mul3ddm10deta->Fill(mindeta);
	mul3ddm10dphi->Fill(mindphi);
	if(abs(mindeta)<0.04) {
	  l3ddm10genmchpt->Fill(genLepPt[genmuidx[i]]);
	  l3ddm10genmcheta->Fill(genLepEta[genmuidx[i]]);
	  l3ddm10genmchphi->Fill(genLepPhi[genmuidx[i]]);
	  l3ddm10genmchd0->Fill(d0);
	  l3ddm10genmchlog10d0->Fill(TMath::Log10(TMath::Abs(d0)));
	  l3ddm10genmch.push_back(true);
	}
	else {
	  l3ddm10genmch.push_back(false);
	}
      }
      else {
	l3ddm10genmch.push_back(false);
      }
    }
    else {
      l3ddm10genmch.push_back(false);
    }

    // For HLT_DoubleL2Mu10No....
    if(trigDoubleL2Mu10NoVtx2ChaPromptL3Mu0Veto>0) {
      double mindeta=1000, mindphi=1000;
      for(unsigned int j=0; j<l2ddm10FiltN; j++) {
	if(abs(genLepEta[genmuidx[i]]-l2ddm10FiltEta[j])<abs(mindeta)) mindeta = genLepEta[genmuidx[i]]-l2ddm10FiltEta[j];
	if(abs(genLepPhi[genmuidx[i]]-l2ddm10FiltPhi[j])<abs(mindphi)) mindphi = genLepPhi[genmuidx[i]]-l2ddm10FiltPhi[j];
      }
      if(l2ddm10FiltN>0) {
	mul2ddm10deta->Fill(mindeta);
	mul2ddm10dphi->Fill(mindphi);
	if(abs(mindeta)<0.04) {
	  l2ddm10genmchpt->Fill(genLepPt[genmuidx[i]]);
	  l2ddm10genmcheta->Fill(genLepEta[genmuidx[i]]);
	  l2ddm10genmchphi->Fill(genLepPhi[genmuidx[i]]);
	  l2ddm10genmchd0->Fill(d0);
	  l2ddm10genmchlog10d0->Fill(TMath::Log10(TMath::Abs(d0)));
	  l2ddm10genmch.push_back(true);
	}
	else {
	  l2ddm10genmch.push_back(false);
	}
      }
      else {
	l2ddm10genmch.push_back(false);
      }
    }
    else {
      l2ddm10genmch.push_back(false);
    } // End of loop on individual trigger efficiency section
    
  } // End of loop on gen-muon indexes

  // Check that the size of vectors are equal and fill the subleading histograms
  if(mud0vec.size()==genmuidx.size() &&
     l3dim33genmch.size()==genmuidx.size() &&
     l2dim23genmch.size()==genmuidx.size() &&
     l2dim23csgenmch.size()==genmuidx.size() &&
     l3ddm10genmch.size()==genmuidx.size() &&
     l2ddm10genmch.size()==genmuidx.size() ) {

    int minpospt = findmin(&genmuidx[0], genLepPt, genmuidx.size(), true, false);
    ordptsubleadmupt->Fill(genLepPt[genmuidx[minpospt]]);
    if(l3dim33genmch[minpospt]) l3dim33genmchordptsubleadmupt->Fill(genLepPt[genmuidx[minpospt]]);
    if(l2dim23genmch[minpospt]) l2dim23genmchordptsubleadmupt->Fill(genLepPt[genmuidx[minpospt]]);
    if(l2dim23csgenmch[minpospt]) l2dim23csgenmchordptsubleadmupt->Fill(genLepPt[genmuidx[minpospt]]);
    if(l3ddm10genmch[minpospt]) l3ddm10genmchordptsubleadmupt->Fill(genLepPt[genmuidx[minpospt]]);
    if(l2ddm10genmch[minpospt]) l2ddm10genmchordptsubleadmupt->Fill(genLepPt[genmuidx[minpospt]]);

    int minposd0 = findmin(&genmuidx[0], &mud0vec[0], genmuidx.size(), false, true);
    ordd0subleadmud0->Fill(mud0vec[minposd0]);
    ordd0subleadmulog10d0->Fill(TMath::Log10(TMath::Abs(mud0vec[minposd0])));
    if(l3dim33genmch[minposd0]) {
      l3dim33genmchordd0subleadmud0->Fill(mud0vec[minposd0]);
      l3dim33genmchordd0subleadmulog10d0->Fill(TMath::Log10(TMath::Abs(mud0vec[minposd0])));
    }
    if(l2dim23genmch[minposd0]) {
      l2dim23genmchordd0subleadmud0->Fill(mud0vec[minposd0]);
      l2dim23genmchordd0subleadmulog10d0->Fill(TMath::Log10(TMath::Abs(mud0vec[minposd0])));
    }
    if(l2dim23csgenmch[minposd0]) {
      l2dim23csgenmchordd0subleadmud0->Fill(mud0vec[minposd0]);
      l2dim23csgenmchordd0subleadmulog10d0->Fill(TMath::Log10(TMath::Abs(mud0vec[minposd0])));
    }
    if(l3ddm10genmch[minposd0]) {
      l3ddm10genmchordd0subleadmud0->Fill(mud0vec[minposd0]);
      l3ddm10genmchordd0subleadmulog10d0->Fill(TMath::Log10(TMath::Abs(mud0vec[minposd0])));
    }
    if(l2ddm10genmch[minposd0]) {
      l2ddm10genmchordd0subleadmud0->Fill(mud0vec[minposd0]);
      l2ddm10genmchordd0subleadmulog10d0->Fill(TMath::Log10(TMath::Abs(mud0vec[minposd0])));
    }
  }
  else {
    cout<<"Error!! Some input not filled for a gen muon. Mismatching vector."<<endl;
  }

  // End of filling histograms

  // Clear vectors
  mud0vec.clear();
  l3dim33genmch.clear();
  l2dim23genmch.clear();
  l2dim23csgenmch.clear();
  l3ddm10genmch.clear();
  l2ddm10genmch.clear();
}

// Function to add a set of histograms for a selection
void DDM10_robustanalyzer::addhistGenPart(TString selection) {

  all1dhists.push_back(new TH1F(selection+"gen_mumult","#mu multiplicity",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"gen_mupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"gen_mueta","#mu #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"gen_muphi","#mu #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"gen_mud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"gen_mulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_ordptsubleadmupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"gen_ordd0subleadmud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"gen_ordd0subleadmulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  
  all1dhists.push_back(new TH1F(selection+"l3dim33genmch_mupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"l3dim33genmch_mueta","#mu #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"l3dim33genmch_muphi","#mu #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"l3dim33genmch_mud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"l3dim33genmch_mulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"l3dim33genmch_ordptsubleadmupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"l3dim33genmch_ordd0subleadmud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"l3dim33genmch_ordd0subleadmulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  
  all1dhists.push_back(new TH1F(selection+"l2dim23genmch_mupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"l2dim23genmch_mueta","#mu #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"l2dim23genmch_muphi","#mu #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"l2dim23genmch_mud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"l2dim23genmch_mulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"l2dim23genmch_ordptsubleadmupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"l2dim23genmch_ordd0subleadmud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"l2dim23genmch_ordd0subleadmulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));

  all1dhists.push_back(new TH1F(selection+"l2dim23csgenmch_mupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"l2dim23csgenmch_mueta","#mu #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"l2dim23csgenmch_muphi","#mu #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"l2dim23csgenmch_mud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"l2dim23csgenmch_mulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"l2dim23csgenmch_ordptsubleadmupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"l2dim23csgenmch_ordd0subleadmud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"l2dim23csgenmch_ordd0subleadmulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  
  all1dhists.push_back(new TH1F(selection+"l3ddm10genmch_mupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"l3ddm10genmch_mueta","#mu #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"l3ddm10genmch_muphi","#mu #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"l3ddm10genmch_mud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"l3ddm10genmch_mulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"l3ddm10genmch_ordptsubleadmupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"l3ddm10genmch_ordd0subleadmud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"l3ddm10genmch_ordd0subleadmulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  
  all1dhists.push_back(new TH1F(selection+"l2ddm10genmch_mupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"l2ddm10genmch_mueta","#mu #eta",52,-2.6,2.6));
  all1dhists.push_back(new TH1F(selection+"l2ddm10genmch_muphi","#mu #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"l2ddm10genmch_mud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"l2ddm10genmch_mulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"l2ddm10genmch_ordptsubleadmupt","#mu p_{T} / GeV",500,-5,495));
  all1dhists.push_back(new TH1F(selection+"l2ddm10genmch_ordd0subleadmud0","#mu d_{0} / cm",24000,-12,12));
  all1dhists.push_back(new TH1F(selection+"l2ddm10genmch_ordd0subleadmulog10d0","#mu log_{10}d_{0} / log_{10}cm",1000,-5,5));
  
  all1dhists.push_back(new TH1F(selection+"genmul3dim33_deta","min. #Delta#eta(gen#mu, L3DiMu33 trig#mu)",10000,-1,1));
  all1dhists.push_back(new TH1F(selection+"genmul3dim33_dphi","min. #Delta#phi(gen#mu, L3DiMu33 trig#mu)",10000,-1,1));
  all1dhists.push_back(new TH1F(selection+"genmul2dim23_deta","min. #Delta#eta(gen#mu, L2DiMu23 trig#mu)",10000,-1,1));
  all1dhists.push_back(new TH1F(selection+"genmul2dim23_dphi","min. #Delta#phi(gen#mu, L2DiMu23 trig#mu)",10000,-1,1));
  all1dhists.push_back(new TH1F(selection+"genmul2dim23cs_deta","min. #Delta#eta(gen#mu, L2DiMu23 CS trig#mu)",10000,-1,1));
  all1dhists.push_back(new TH1F(selection+"genmul2dim23cs_dphi","min. #Delta#phi(gen#mu, L2DiMu23 CS trig#mu)",10000,-1,1));
  all1dhists.push_back(new TH1F(selection+"genmul3ddm10_deta","min. #Delta#eta(gen#mu, L3DDMu10 trig#mu)",10000,-1,1));
  all1dhists.push_back(new TH1F(selection+"genmul3ddm10_dphi","min. #Delta#phi(gen#mu, L3DDMu10 trig#mu)",10000,-1,1));
  all1dhists.push_back(new TH1F(selection+"genmul2ddm10_deta","min. #Delta#eta(gen#mu, L2DDMu10 trig#mu)",10000,-1,1));
  all1dhists.push_back(new TH1F(selection+"genmul2ddm10_dphi","min. #Delta#phi(gen#mu, L2DDMu10 trig#mu)",10000,-1,1));
}

// Function to sort the indices based on a factor (Usually pT)
void DDM10_robustanalyzer::sort(int* idx, double* factor, int n) {
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

// Function to find the minimum of a list 
int DDM10_robustanalyzer::findmin(int* idx, double* factor, int n, bool useidx=true, bool useabs=false) {
  int minpos = 0;
  if(useidx) {
    for(unsigned int i=0; i<n; i++) {
      if((useabs && (abs(*(factor+*(idx+i)))<abs(*(factor+*(idx+minpos))))) ||
	 (!useabs && (*(factor+*(idx+i))<*(factor+*(idx+minpos))))
	 )
	minpos = i;
    }
  }
  if(!useidx) {
    for(unsigned int i=0; i<n; i++) {
      if((useabs && (abs(*(factor+i))<abs(*(factor+minpos)))) ||
	 (!useabs && (*(factor+i)<*(factor+minpos)))
	 )
	minpos = i;
    }
  }
  return minpos;
}
