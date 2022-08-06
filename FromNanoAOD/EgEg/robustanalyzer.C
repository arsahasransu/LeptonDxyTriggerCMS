/*
 * AUTHOR: Abanti Ranadhir Sahasransu - asahasra@cern.ch
 * The code was made to read NanoAOD provately created by the EXO group
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
  inputChain->SetBranchAddress("bunchCrossing", &bunch);

  inputChain->SetBranchAddress("nElectron", &eln);
  //inputChain->SetBranchAddress("Electron_deltaEtaSC", &eldetasc);
  //inputChain->SetBranchAddress("Electron_dxy", &eldxy);
  //inputChain->SetBranchAddress("Electron_dxyErr", &eldxyerr);
  //inputChain->SetBranchAddress("Electron_dz", &eldz);
  //inputChain->SetBranchAddress("Electron_dzerr", &eldzerr);
  //inputChain->SetBranchAddress("Electron_eInvMinusPInv", &elooemoop);
  //inputChain->SetBranchAddress("Electron_energyErr", &elenergyerr);
  inputChain->SetBranchAddress("Electron_eta", &eleta);
  //inputChain->SetBranchAddress("Electron_hoe", &elhoe);
  inputChain->SetBranchAddress("Electron_phi", &elphi);
  inputChain->SetBranchAddress("Electron_pt", &elpt);
  //inputChain->SetBranchAddress("Electron_r9", &elr9);
  //inputChain->SetBranchAddress("Electron_scEtOverPt", &elscetoverpt);
  //inputChain->SetBranchAddress("Electron_sieie", &elsieie);
  //inputChain->SetBranchAddress("Electron_charge", &elcharge);
  //inputChain->SetBranchAddress("Electron_cutBased", &elcutbasedid);
  //inputChain->SetBranchAddress("Electron_pdgId", &elpdgid);
  //inputChain->SetBranchAddress("Electron_photonIdx", &elphotonIdx);

  inputChain->SetBranchAddress("nLowPtElectron", &lowpteln);
  //inputChain->SetBranchAddress("LowPtElectron_ID", &lowptelid);
  //inputChain->SetBranchAddress("LowPtElectron_deltaEtaSC", &lowpteldetasc);
  //inputChain->SetBranchAddress("LowPtElectron_dxy", &lowpteldxy);
  //inputChain->SetBranchAddress("LowPtElectron_dxyErr", &lowpteldxyerr);
  //inputChain->SetBranchAddress("LowPtElectron_dz", &lowpteldz);
  //inputChain->SetBranchAddress("LowPtElectron_dzErr", &lowpteldzerr);
  //inputChain->SetBranchAddress("LowPtElectron_eInvMinusPInv", &lowptelooemoop);
  //inputChain->SetBranchAddress("LowPtElectron_energyErr", &lowptelenergyerr);
  inputChain->SetBranchAddress("LowPtElectron_eta", &lowpteleta);
  //inputChain->SetBranchAddress("LowPtElectron_hoe", &lowptelhoe);
  inputChain->SetBranchAddress("LowPtElectron_phi", &lowptelphi);
  inputChain->SetBranchAddress("LowPtElectron_pt", &lowptelpt);
  //inputChain->SetBranchAddress("LowPtElectron_ptbiased", &lowptelptbias);
  //inputChain->SetBranchAddress("LowPtElectron_scEtOverPt", &lowptelscetoverpt);
  //inputChain->SetBranchAddress("LowPtElectron_sieie", &lowptelsieie);
  //inputChain->SetBranchAddress("LowPtElectron_unbiased", &lowptelunbias);
  //inputChain->SetBranchAddress("LowPtElectron_charge", &lowptelcharge);
  //inputChain->SetBranchAddress("LowPtElectron_electronIdx", &lowptelid);
  //inputChain->SetBranchAddress("LowPtElectron_pdgId", &lowptelpdgid);

  inputChain->SetBranchAddress("nPhoton", &phn);
  inputChain->SetBranchAddress("Photon_pt", &phpt);
  inputChain->SetBranchAddress("Photon_eta", &pheta);
  inputChain->SetBranchAddress("Photon_phi", &phphi);
  
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
  int nCores = 6; // Assume parallel processing over 7 cores where
  // there is a lesser no.of events in the last core
  int beginevent = splitCnt*(totEntries/nCores);
  int endevent = (splitCnt+1)*(totEntries/nCores);
  endevent = endevent<totEntries?endevent:totEntries; // Verfied that this logic to parallelize works

  // Count events passing certain selections
  int nosel=0;

  addhist("noselelectron");
  addhist("nosellowptelectron");
  addhist("noselphoton");
  
  // Loop beginning on events
  for(unsigned int event=beginevent; event<endevent; event++) {

    vector<int> noselelidx;
    vector<int> nosellowptelidx;
    vector<int> noselphidx;
  
    inputChain->GetEntry(event);
    //if(event>1000) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;
    
    // Loop beginning on electrons
    for(unsigned int idx=0; idx<eln; idx++) {

      noselelidx.push_back(idx);
      
    } // End of loop on electrons
       
    // Loop beginning on low pt electrons
    for(unsigned int idx=0; idx<lowpteln; idx++) {

      nosellowptelidx.push_back(idx);
      
    } // End of loop on low pt electrons
       
    // Loop beginning on photons
    for(unsigned int idx=0; idx<phn; idx++) {

      noselphidx.push_back(idx);
      
    } // End of loop on photons
       
    fillhistinevent("noselelectron", noselelidx);
    fillhistinevent("nosellowptelectron", nosellowptelidx);
    fillhistinevent("noselphoton", noselphidx);

    // Clear all the vectors
    noselelidx.clear();
    nosellowptelidx.clear();
    noselphidx.clear();

  } // End of loop on events

  cout<<totEntries<<"\t"<<nosel<<endl;
}

// Function to fill a set of histograms in the event
void robustanalyzer::fillhistinevent(TString selection, vector<int> elidx) {

  // nothing here for now
  TH1F* elmult = (TH1F*) outfile->Get(selection+"electron_mult");
  TH1F* leadelpt = (TH1F*) outfile->Get(selection+"electron_lead_pt");
  TH1F* leadeleta = (TH1F*) outfile->Get(selection+"electron_lead_eta");
  TH1F* leadelphi = (TH1F*) outfile->Get(selection+"electron_lead_phi");
  TH1F* subleadelpt = (TH1F*) outfile->Get(selection+"electron_sublead_pt");
  TH1F* subleadeleta = (TH1F*) outfile->Get(selection+"electron_sublead_eta");
  TH1F* subleadelphi = (TH1F*) outfile->Get(selection+"electron_sublead_phi");
  /*
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
  */

  if(elidx.size()>0) {
    elmult->Fill(elidx.size());
    leadelpt->Fill(elpt[elidx[0]]);
    leadeleta->Fill(eleta[elidx[0]]);
    leadelphi->Fill(elphi[elidx[0]]);
    /*
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

    */
  } // End of condition requiring atleast one eg object
  
  if(elidx.size()>=2) { // Condition requiring atleast two eg object
    //TLorentzVector leadeg, subleadeg;
    //leadeg.SetPtEtaPhiM(egRecoPt[egidx[0]],egRecoEta[egidx[0]],egRecoPhi[egidx[0]],0.106);
    //subleadeg.SetPtEtaPhiM(egRecoPt[egidx[1]],egRecoEta[egidx[1]],egRecoPhi[egidx[1]],0.106);

    subleadelpt->Fill(elpt[elidx[1]]);
    subleadeleta->Fill(eleta[elidx[1]]);
    subleadelphi->Fill(elphi[elidx[1]]);
    /*
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
      }*/
  } // End of condition requiring atleast two eg object

}

// Function to add a set of histograms for a selection
void robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"electron_mult","N electron",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"electron_lead_pt","electron p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"electron_lead_eta","electron #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"electron_lead_phi","electron #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"electron_sublead_pt","electron p_{T} / GeV",550,-50,500));
  all1dhists.push_back(new TH1F(selection+"electron_sublead_eta","electron #eta",100,-5,5));
  all1dhists.push_back(new TH1F(selection+"electron_sublead_phi","electron #phi",66,-3.3,3.3));
  /*
  // barrel variables - lead pT e/gamma
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegclustershape","barrel electron clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegin5x5clusshape","barrel electron #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegin5x5noiseclnd","barrel electron #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadeghovere","barrel electron H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegscenergy","barrel electron SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegecalpfclustiso","barrel electron ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegecalpfclustisoovere","barrel electron ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadeghcalpfclustiso","barrel electron hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadeghcalpfclustisoovere","barrel electron hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegpixelmchvar_s2","barrel electron pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegchi2","barrel electron #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegdeta","barrel electron #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegdetaseed","barrel electron #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegdphi","barrel electron #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegmhits","barrel electron missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegnlayerit","barrel electron num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegooeseedoop","barrel electron 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegooesclsoop","barrel electron 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegvalhits","barrel electron valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegtrkiso","barrel electron track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadeghovereoversupcluse","barrel electron H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoeb_leadegseedclustime","barrel electron_{seed} time / ns",20000,-10,10));

  // barrel variables - sub-lead pT e/gamma
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegclustershape","barrel electron clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegin5x5clusshape","barrel electron #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegin5x5noiseclnd","barrel electron #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadeghovere","barrel electron H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegscenergy","barrel electron SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegecalpfclustiso","barrel electron ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegecalpfclustisoovere","barrel electron ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadeghcalpfclustiso","barrel electron hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadeghcalpfclustisoovere","barrel electron hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegpixelmchvar_s2","barrel electron pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegchi2","barrel electron #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegdeta","barrel electron #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegdetaseed","barrel electron #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegdphi","barrel electron #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegmhits","barrel electron missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegnlayerit","barrel electron num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegooeseedoop","barrel electron 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegooesclsoop","barrel electron 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegvalhits","barrel electron valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegtrkiso","barrel electron track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadeghovereoversupcluse","barrel electron H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoeb_subleadegseedclustime","barrel electron_{seed} time / ns",20000,-10,10));

  // barrel variables - invariant mass
  all1dhists.push_back(new TH1F(selection+"recoeb_leadsubleadM","M(electron_{1},electron_{2}) / GeV",500,0,500));

  // end-cap variables - lead pT e/gamma
  all1dhists.push_back(new TH1F(selection+"recoee_leadegclustershape","end-cap electron clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegin5x5clusshape","end-cap electron #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegin5x5noiseclnd","end-cap electron #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_leadeghovere","end-cap electron H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegscenergy","end-cap electron SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegecalpfclustiso","end-cap electron ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegecalpfclustisoovere","end-cap electron ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoee_leadeghcalpfclustiso","end-cap electron hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_leadeghcalpfclustisoovere","end-cap electron hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegpixelmchvar_s2","end-cap electron pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegchi2","end-cap electron #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegdeta","end-cap electron #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegdetaseed","end-cap electron #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegdphi","end-cap electron #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegmhits","end-cap electron missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegnlayerit","end-cap electron num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegooeseedoop","end-cap electron 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegooesclsoop","end-cap electron 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegvalhits","end-cap electron valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegtrkiso","end-cap electron track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_leadeghovereoversupcluse","end-cap electron H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoee_leadegseedclustime","end-cap electron_{seed} time / ns",20000,-10,10));

  // end-cap variables - sub-lead pT e/gamma
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegclustershape","end-cap electron clus.shape",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegin5x5clusshape","end-cap electron #sigmaE i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegin5x5noiseclnd","end-cap electron #sigmaE(noise clean) i#etai#eta 5x5",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadeghovere","end-cap electron H / GeV",1000,0,100));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegscenergy","end-cap electron SC energy",5000,0,5000));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegecalpfclustiso","end-cap electron ecal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegecalpfclustisoovere","end-cap electron ecal PF Iso./E",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadeghcalpfclustiso","end-cap electron hcal PF Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadeghcalpfclustisoovere","end-cap electron hcal PF Iso.",1000,-0.5,9.5));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegpixelmchvar_s2","end-cap electron pixelmachvar",1000,-50,950));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegchi2","end-cap electron #chi^{2}",1000,-10,90));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegdeta","end-cap electron #Delta#eta",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegdetaseed","end-cap electron #Delta#eta seed",1000,-0.03,0.07));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegdphi","end-cap electron #Delta#phi",1000,-0.3,0.7));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegmhits","end-cap electron missing hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegnlayerit","end-cap electron num IT layer",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegooeseedoop","end-cap electron 1/E_seed-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegooesclsoop","end-cap electron 1/E_sc-1/p",1000,-0.1,0.9));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegvalhits","end-cap electron valid hits",40,-5,35));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegtrkiso","end-cap electron track Iso.",1000,-5,95));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadeghovereoversupcluse","end-cap electron H/E",1000,0,10));
  all1dhists.push_back(new TH1F(selection+"recoee_subleadegseedclustime","end-cap electron_{seed} time / ns",20000,-10,10));

  // end-cap variables - invariant mass
  all1dhists.push_back(new TH1F(selection+"recoee_leadsubleadM","M(electron_{1},electron_{2}) / GeV",500,0,500));
  */
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
