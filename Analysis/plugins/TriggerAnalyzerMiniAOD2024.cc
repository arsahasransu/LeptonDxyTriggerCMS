// system include files
#include <memory>
#include <string>
#include <float.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
#include "TTree.h"

class TriggerAnalyzerMiniAOD2024 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit TriggerAnalyzerMiniAOD2024(const edm::ParameterSet&);
  ~TriggerAnalyzerMiniAOD2024();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void clearVars();

  edm::EDGetTokenT< edm::TriggerResults > trgResultsToken_;
  edm::EDGetTokenT< pat::TriggerObjectStandAloneCollection > trgObjectsToken_;
  edm::EDGetTokenT< edm::SortedCollection< EcalRecHit, edm::StrictWeakOrdering< EcalRecHit >>> rechiteb_token;
  edm::EDGetTokenT< edm::SortedCollection< EcalRecHit, edm::StrictWeakOrdering< EcalRecHit >>> rechitee_token;
  //edm::EDGetTokenT< std::vector< pat::Electron >> electron_token;
  edm::EDGetTokenT< std::vector< reco::SuperCluster >> suprclus_token;
  //edm::EDGetTokenT< std::vector< pat::Electron >> lowptelectron_token;
  edm::EDGetTokenT< std::vector< pat::Photon >> photon_token;
  //edm::EDGetTokenT< std::vector< pat::Photon >> ootphoton_token;
  edm::EDGetTokenT< reco::BeamSpot > BS_token;
  edm::EDGetTokenT< double > rho_token;
  edm::EDGetTokenT< std::vector< reco::Vertex >> PV_token;

  edm::Service< TFileService > fs;
  TTree* tree;
  int run;
  int lumSec;

  bool HLT_DiPhoton10Time1p4ns, HLT_DiPhoton10Time1ns, HLT_DiPhoton10_CaloIdL;
  bool HLTOR_METTrig;
  bool HLTOR_METTrigFull;
  bool HLTOR_JetTrigFull;

  double bs_x;
  double bs_y;
  double bs_z;
  double rho;

  int pho_n;
  vector<double> pho_e;
  vector<double> pho_pt;
  vector<double> pho_eta;
  vector<double> pho_phi;
  vector<double> pho_seedtime;
  vector<double> pho_smin;
  vector<double> pho_smaj;
  vector<double> pho_sinin_noiseclnd;
  vector<double> pho_hoe;
  vector<double> pho_chargedhadroniso;
  vector<double> pho_neutralhadroniso;
  vector<double> pho_photoniso;

  int dieg10time1ns_usfinfilt_n;
  vector<double> dieg10time1ns_usfinfilt_pt;
  vector<double> dieg10time1ns_usfinfilt_eta;
  vector<double> dieg10time1ns_usfinfilt_phi;

  int pv_n;
  vector<double> pv_x;
  vector<double> pv_xerr;
  vector<double> pv_y;
  vector<double> pv_yerr;
  vector<double> pv_z;
  vector<double> pv_zerr;
  vector<bool> pv_isvalid;
};

TriggerAnalyzerMiniAOD2024::TriggerAnalyzerMiniAOD2024(const edm::ParameterSet& iConfig) {
  trgResultsToken_= consumes< edm::TriggerResults >( edm::InputTag("TriggerResults::HLT") );
  trgObjectsToken_ = consumes< pat::TriggerObjectStandAloneCollection >( edm::InputTag("slimmedPatTrigger") );
  rechiteb_token = consumes< edm::SortedCollection< EcalRecHit, edm::StrictWeakOrdering< EcalRecHit >>>( edm::InputTag("reducedEgamma:reducedEBRecHits") );
  rechitee_token = consumes< edm::SortedCollection< EcalRecHit, edm::StrictWeakOrdering< EcalRecHit >>>( edm::InputTag("reducedEgamma:reducedEERecHits") );
  //electron_token = consumes< std::vector< pat::Electron >>( edm::InputTag("slimmedElectrons") );
  suprclus_token = consumes< std::vector< reco::SuperCluster >>( edm::InputTag("reducedEgamma:reducedSuperClusters") );
  //lowptelectron_token = consumes< std::vector< pat::Electron >>( edm::InputTag("slimmedLowPtElectrons") );
  photon_token = consumes< std::vector< pat::Photon >>( edm::InputTag("slimmedPhotons") );
  //ootphoton_token = consumes< std::vector< pat::Photon >>( edm::InputTag("slimmedOOTPhotons") );
  BS_token = consumes< reco::BeamSpot > ( edm::InputTag("offlineBeamSpot"));
  rho_token = consumes< double > ( edm::InputTag("fixedGridRhoAll"));
  PV_token = consumes< std::vector< reco::Vertex >> ( edm::InputTag("offlineSlimmedPrimaryVertices"));
  
  usesResource("TFileService");

  tree = fs->make<TTree>("tree", "tree");

  tree->Branch("run", &run, "run/I");
  tree->Branch("lumSec", &lumSec, "lumSec/I");

  tree->Branch("HLT_DiPhoton10Time1p4ns", &HLT_DiPhoton10Time1p4ns, "HLT_DiPhoton10Time1p4ns/O");
  tree->Branch("HLT_DiPhoton10Time1ns", &HLT_DiPhoton10Time1ns, "HLT_DiPhoton10Time1ns/O");
  tree->Branch("HLT_DiPhoton10_CaloIdL", &HLT_DiPhoton10_CaloIdL, "HLT_DiPhoton10_CaloIdL/O");
  tree->Branch("HLTOR_METTrig", &HLTOR_METTrig, "HLTOR_METTrig/O");
  tree->Branch("HLTOR_METTrigFull", &HLTOR_METTrigFull, "HLTOR_METTrigFull/O");
  tree->Branch("HLTOR_JetTrigFull", &HLTOR_JetTrigFull, "HLTOR_JetTrigFull/O");

  tree->Branch("bs_x", &bs_x, "bs_x/D");
  tree->Branch("bs_y", &bs_y, "bs_y/D");
  tree->Branch("bs_z", &bs_z, "bs_z/D");
  tree->Branch("rho", &rho, "rho/D");

  tree->Branch("pho_n", &pho_n, "pho_n/I");
  tree->Branch("pho_e", &pho_e);
  tree->Branch("pho_pt", &pho_pt);
  tree->Branch("pho_eta", &pho_eta);
  tree->Branch("pho_phi", &pho_phi);
  tree->Branch("pho_seedtime", &pho_seedtime);
  tree->Branch("pho_smin", &pho_smin);
  tree->Branch("pho_smaj", &pho_smaj);
  tree->Branch("pho_sinin_noiseclnd", &pho_sinin_noiseclnd);
  tree->Branch("pho_hoe", &pho_hoe);
  tree->Branch("pho_chargedhadroniso", &pho_chargedhadroniso);
  tree->Branch("pho_neutralhadroniso", &pho_neutralhadroniso);
  tree->Branch("pho_photoniso", &pho_photoniso);

  tree->Branch("dieg10time1ns_usfinfilt_n", &dieg10time1ns_usfinfilt_n, "dieg10time1ns_usfinfilt_n/i");
  tree->Branch("dieg10time1ns_usfinfilt_pt", &dieg10time1ns_usfinfilt_pt);
  tree->Branch("dieg10time1ns_usfinfilt_eta", &dieg10time1ns_usfinfilt_eta);
  tree->Branch("dieg10time1ns_usfinfilt_phi", &dieg10time1ns_usfinfilt_phi);

  tree->Branch("pv_n", &pv_n, "pv_n/I");
  tree->Branch("pv_x", &pv_x);
  tree->Branch("pv_xerr", &pv_xerr);
  tree->Branch("pv_y", &pv_y);
  tree->Branch("pv_yerr", &pv_yerr);
  tree->Branch("pv_z", &pv_z);
  tree->Branch("pv_zerr", &pv_zerr);
  tree->Branch("pv_isvalid", &pv_isvalid);
}


TriggerAnalyzerMiniAOD2024::~TriggerAnalyzerMiniAOD2024() {
}

void TriggerAnalyzerMiniAOD2024::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;

  const double def_max = 3.7e10;
  
  run = iEvent.eventAuxiliary().run();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();

  // Beam Spot
  edm::Handle<reco::BeamSpot> bsH;
  iEvent.getByToken(BS_token, bsH);
  if(bsH.isValid()) {
    bs_x = bsH->x0();
    bs_y = bsH->y0();
    bs_z = bsH->z0();
  }
  else {
    bs_x = def_max;
    bs_y = def_max;
    bs_z = def_max;
  }

  // Rho token
  edm::Handle<double> rhH;
  iEvent.getByToken(rho_token, rhH);
  if(rhH.isValid()) {
    rho = (*rhH);
  }

  // Primary Vertex
  edm::Handle<std::vector<reco::Vertex>> pvH;
  iEvent.getByToken(PV_token, pvH);
  pv_n = 0;
  if(pvH.isValid()) {
    for(auto pv_iter=pvH->begin(); pv_iter!=pvH->end(); ++pv_iter) {
      pv_x.push_back(pv_iter->x());
      pv_xerr.push_back(pv_iter->xError());
      pv_y.push_back(pv_iter->y());
      pv_yerr.push_back(pv_iter->yError());
      pv_z.push_back(pv_iter->z());
      pv_zerr.push_back(pv_iter->zError());
      pv_isvalid.push_back(pv_iter->isValid());
      pv_n++;
    }
  }

  HLT_DiPhoton10Time1p4ns = false;
  HLT_DiPhoton10Time1ns = false;
  HLT_DiPhoton10_CaloIdL = false;
  HLTOR_METTrig = false;
  HLTOR_METTrigFull = false;
  HLTOR_JetTrigFull = false;

  std::string all_HLT_Jet_paths[] = {"HLT_AK8DiPFJet250_250_MassSD30_v", "HLT_AK8DiPFJet250_250_MassSD50_v",
                                   "HLT_AK8DiPFJet260_260_MassSD30_v", "HLT_AK8DiPFJet260_260_MassSD50_v", 
                                   "HLT_AK8DiPFJet270_270_MassSD30_v", "HLT_AK8DiPFJet280_280_MassSD30_v", 
                                   "HLT_AK8DiPFJet290_290_MassSD30_v", "HLT_AK8PFJet140_v", "HLT_AK8PFJet200_v", 
                                   "HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50_v", 
                                   "HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53_v", 
                                   "HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55_v", 
                                   "HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60_v", 
                                   "HLT_AK8PFJet220_SoftDropMass40_v", "HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06_v", 
                                   "HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10_v", 
                                   "HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03_v", 
                                   "HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05_v", 
                                   "HLT_AK8PFJet230_SoftDropMass40_v", 
                                   "HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06_v", 
                                   "HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10_v", 
                                   "HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03_v", 
                                   "HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05_v", "HLT_AK8PFJet260_v", 
                                   "HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06_v", 
                                   "HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10_v", 
                                   "HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03_v", 
                                   "HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05_v", 
                                   "HLT_AK8PFJet320_v", "HLT_AK8PFJet400_MassSD30_v", "HLT_AK8PFJet400_v", 
                                   "HLT_AK8PFJet40_v", "HLT_AK8PFJet420_MassSD30_v", 
                                   "HLT_AK8PFJet425_SoftDropMass40_v", "HLT_AK8PFJet450_MassSD30_v", 
                                   "HLT_AK8PFJet450_SoftDropMass40_v", "HLT_AK8PFJet450_v", 
                                   "HLT_AK8PFJet470_MassSD30_v", "HLT_AK8PFJet500_MassSD30_v", 
                                   "HLT_AK8PFJet500_v", "HLT_AK8PFJet550_v", "HLT_AK8PFJet60_v", 
                                   "HLT_AK8PFJet80_v", "HLT_AK8PFJetFwd140_v", "HLT_AK8PFJetFwd15_v", 
                                   "HLT_AK8PFJetFwd200_v", "HLT_AK8PFJetFwd25_v", "HLT_AK8PFJetFwd260_v", 
                                   "HLT_AK8PFJetFwd320_v", "HLT_AK8PFJetFwd400_v", "HLT_AK8PFJetFwd40_v", 
                                   "HLT_AK8PFJetFwd450_v", "HLT_AK8PFJetFwd500_v", "HLT_AK8PFJetFwd60_v", 
                                   "HLT_AK8PFJetFwd80_v", "HLT_CaloJet500_NoJetID_v", "HLT_CaloJet550_NoJetID_v", 
                                   "HLT_CaloMET350_NotCleaned_v", "HLT_CaloMET90_NotCleaned_v", "HLT_CaloMHT90_v", 
                                   "HLT_DiPFJetAve100_HFJEC_v", "HLT_DiPFJetAve140_v", "HLT_DiPFJetAve160_HFJEC_v", 
                                   "HLT_DiPFJetAve200_v", "HLT_DiPFJetAve220_HFJEC_v", "HLT_DiPFJetAve260_HFJEC_v", 
                                   "HLT_DiPFJetAve260_v", "HLT_DiPFJetAve300_HFJEC_v", "HLT_DiPFJetAve320_v", 
                                   "HLT_DiPFJetAve400_v", "HLT_DiPFJetAve40_v", "HLT_DiPFJetAve500_v", 
                                   "HLT_DiPFJetAve60_HFJEC_v", "HLT_DiPFJetAve60_v", "HLT_DiPFJetAve80_HFJEC_v", 
                                   "HLT_DiPFJetAve80_v", "HLT_DoublePFJets100_PFBTagDeepJet_p71_v", 
                                   "HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71_v", 
                                   "HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71_v", 
                                   "HLT_DoublePFJets200_PFBTagDeepJet_p71_v", 
                                   "HLT_DoublePFJets350_PFBTagDeepJet_p71_v", "HLT_DoublePFJets40_PFBTagDeepJet_p71_v", 
                                   "HLT_Mu12_DoublePFJets100_PFBTagDeepJet_p71_v", 
                                   "HLT_Mu12_DoublePFJets200_PFBTagDeepJet_p71_v", 
                                   "HLT_Mu12_DoublePFJets350_PFBTagDeepJet_p71_v", 
                                   "HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepJet_p71_v", 
                                   "HLT_Mu12_DoublePFJets40_PFBTagDeepJet_p71_v", 
                                   "HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepJet_p71_v", 
                                   "HLT_Mu12eta2p3_PFJet40_v", 
                                   "HLT_PFJet110_v", "HLT_PFJet140_v", "HLT_PFJet200_v", "HLT_PFJet260_v", 
                                   "HLT_PFJet320_v", "HLT_PFJet400_v", "HLT_PFJet40_v", "HLT_PFJet450_v", 
                                   "HLT_PFJet500_v", "HLT_PFJet550_v15", "HLT_PFJet60_v", "HLT_PFJet80_v", 
                                   "HLT_PFJetFwd140_v", "HLT_PFJetFwd200_v", "HLT_PFJetFwd260_v", 
                                   "HLT_PFJetFwd320_v", "HLT_PFJetFwd400_v", "HLT_PFJetFwd40_v", 
                                   "HLT_PFJetFwd450_v", "HLT_PFJetFwd500_v", "HLT_PFJetFwd60_v", 
                                   "HLT_PFJetFwd80_v", "HLT_QuadPFJet100_88_70_30_v", 
                                   "HLT_QuadPFJet103_88_75_15_v", "HLT_QuadPFJet105_88_75_30_v", 
                                   "HLT_QuadPFJet105_88_76_15_v", "HLT_QuadPFJet111_90_80_15_v", 
                                   "HLT_QuadPFJet111_90_80_30_v"
                                  };

  std::string all_MET_paths[] = {"HLT_MET105_IsoTrk50_v", "HLT_MET120_IsoTrk50_v",
                               "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v", 
                               "HLT_PFHT500_PFMET110_PFMHT110_IDTight_v",
                               "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v",
                               "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v",
                               "HLT_PFMET105_IsoTrk50_v", "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
                               "HLT_PFMET120_PFMHT120_IDTight_v", "HLT_PFMET130_PFMHT130_IDTight_v",
                               "HLT_PFMET140_PFMHT140_IDTight_v", "HLT_PFMET200_BeamHaloCleaned_v",
                               "HLT_PFMET200_NotCleaned_v", "HLT_PFMET250_NotCleaned_v",
                               "HLT_PFMET300_NotCleaned_v",
                               "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF_v",
                               "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF_v",
                               "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",
                               "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", 
                               "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF_v", 
                               "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v", 
                               "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF_v", 
                               "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v", 
                               "HLT_PFMETTypeOne140_PFMHT140_IDTight_v", 
                               "HLT_PFMETTypeOne200_BeamHaloCleaned_v",
                              };

  // MET Trigger bits
  bool HLT_PFMET120_PFMHT120_IDTight=false, HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF=false, HLT_PFMETNoMu120_PFMHTNoMu120_IDTight=false, HLT_CaloMET80_NotCleaned=false, HLT_PFMET200_NotCleaned=false, HLT_PFMET200_BeamHaloCleaned=false, HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight=false;

  //Accessing trigger bits:
  //This works in both RAW, AOD or MINIAOD 
  //Here we access the decision provided by the HLT (i.e. original trigger step). 
  edm::Handle<edm::TriggerResults> trgResultsH;
  iEvent.getByToken(trgResultsToken_, trgResultsH);
  if( !trgResultsH.failedToGet() ) {
    int N_Triggers = trgResultsH->size();
    //cout<<"Number of triggers: "<<N_Triggers<<endl;
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trgResultsH);
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trgResultsH.product()->accept(i_Trig)) {
  	    //cout << "Path: " <<trigName.triggerName(i_Trig)<<"Results: "<<trgResultsH.product()->accept(i_Trig)<<endl;
	      TString TrigPath =trigName.triggerName(i_Trig);
	      if(TrigPath.Index("HLT_DiPhoton10Time1p4ns_v") >=0) HLT_DiPhoton10Time1p4ns = true; 
	      if(TrigPath.Index("HLT_DiPhoton10Time1ns_v") >=0) HLT_DiPhoton10Time1ns = true; 
	      if(TrigPath.Index("HLT_DiPhoton10_CaloIdL_v") >=0) HLT_DiPhoton10_CaloIdL = true; 
	      if(TrigPath.Index("HLT_PFMET120_PFMHT120_IDTight_v") >=0) HLT_PFMET120_PFMHT120_IDTight=true;
	      if(TrigPath.Index("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF_v") >=0) HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF=true;
	      if(TrigPath.Index("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") >=0) HLT_PFMETNoMu120_PFMHTNoMu120_IDTight=true;
	      if(TrigPath.Index("HLT_CaloMET80_NotCleaned_v") >=0) HLT_CaloMET80_NotCleaned=true;
	      if(TrigPath.Index("HLT_PFMET200_NotCleaned_v") >=0) HLT_PFMET200_NotCleaned=true;
	      if(TrigPath.Index("HLT_PFMET200_BeamHaloCleaned_v") >=0) HLT_PFMET200_BeamHaloCleaned=true;
	      if(TrigPath.Index("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v") >=0) HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight=true;

        // Loop over the string of HLT paths and obtain the OR of all decision bits
        for (auto const& path : all_HLT_Jet_paths) {
          if( TrigPath.Index(path.c_str()) >=0 ) {
            HLTOR_JetTrigFull = true;
            break;
          }
        }
        for (auto const& path : all_MET_paths) {
          if( TrigPath.Index(path.c_str()) >=0 ) {
            HLTOR_METTrigFull = true;
            break;
          }
        }
      }
    }
    //if(HLT_DiPhoton10sminlt0p12 || HLT_DiPhoton10Time1p4ns || HLT_DiPhoton10_CaloIdL) cout<<"Passing one of the desired triggers!"<<endl;
  } // End of loop for accessing the trigger bits
  HLTOR_METTrig = (HLT_PFMET120_PFMHT120_IDTight || 
		   HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF || 
		   HLT_PFMETNoMu120_PFMHTNoMu120_IDTight || 
		   HLT_CaloMET80_NotCleaned || 
		   HLT_PFMET200_NotCleaned || 
		   HLT_PFMET200_BeamHaloCleaned || 
		   HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);


  // Trigger Filter Objects
  edm::Handle<pat::TriggerObjectStandAloneCollection> trgObjectsH;
  iEvent.getByToken(trgObjectsToken_, trgObjectsH);
  dieg10time1ns_usfinfilt_n = 0;
  std::string filterName_ = "hltDiEG10CaloIdLTime1nsUnseededFilter";
  if( trgObjectsH.isValid() && trgResultsH.isValid() ) {
    if( trgObjectsH->empty() ) {
      edm::LogWarning("TriggerAnalyzerMiniAOD2024") << "Trigger objects collection is empty";
    }
    if( trgResultsH->size() == 0 ) {
      edm::LogWarning("TriggerAnalyzerMiniAOD2024") << "Trigger results collection is empty";
    }
    const edm::TriggerNames &names = iEvent.triggerNames(*trgResultsH);
    for (pat::TriggerObjectStandAlone obj : *trgObjectsH) {
      obj.unpackFilterLabels(iEvent, *trgResultsH);
      obj.unpackPathNames(names);
      if( obj.filterLabels().empty() ) {
        std::cout<<"No filter labels found for this trigger object"<<std::endl;
        continue;
      }
      
      for (auto& filterLabel : obj.filterLabels()) {
	      if (filterLabel == filterName_) {
          dieg10time1ns_usfinfilt_pt.push_back( obj.pt() );
          dieg10time1ns_usfinfilt_eta.push_back( obj.eta() );
          dieg10time1ns_usfinfilt_phi.push_back( obj.phi() );
          dieg10time1ns_usfinfilt_n++;
	      }
      }
    }
  }

  // ECAL rechits  
  edm::Handle< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > rechitebH;
  iEvent.getByToken(rechiteb_token, rechitebH);
  edm::Handle< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > rechiteeH;
  iEvent.getByToken(rechitee_token, rechiteeH);

  // Super Cluster collection
  edm::Handle< std::vector<reco::SuperCluster> > sclusterH;
  iEvent.getByToken(suprclus_token, sclusterH);

  // Photon
  edm::Handle<std::vector<pat::Photon> > photonH;
  iEvent.getByToken(photon_token, photonH);
  pho_n = 0;
  if(photonH.isValid()) {
    double seedtime = def_max;
    for(auto pho_iter=photonH->begin(); pho_iter!=photonH->end(); pho_iter++) {
      pho_e.push_back(pho_iter->energy());
      pho_pt.push_back(pho_iter->pt());
      pho_eta.push_back(pho_iter->eta());
      pho_phi.push_back(pho_iter->phi());
      
      seedtime = def_max;
      DetId SCseedID = pho_iter->seed()->seed();
      if(rechitebH.isValid() && seedtime==def_max) {
	auto rechitseed = rechitebH->find(SCseedID);
	if(rechitseed!=rechitebH->end()) {
	  seedtime = rechitseed->time();
	}
      }
      if(rechiteeH.isValid() && seedtime==def_max) {
	auto rechitseed = rechiteeH->find(SCseedID);
	if(rechitseed!=rechiteeH->end()) {
	  seedtime = rechitseed->time();
	}
      }
      pho_seedtime.push_back(seedtime);

      if( (sclusterH.isValid()) &&
	  ( ( (std::abs(pho_iter->eta())<1.479) && rechitebH.isValid() ) || 
	    ( (std::abs(pho_iter->eta())>=1.479) && rechiteeH.isValid() ) 
	    ) 
	  ) {

	const EcalRecHitCollection* rechits = (std::abs(pho_iter->eta()) < 1.479) ? rechitebH.product() : rechiteeH.product();
	reco::CaloClusterPtr SCseed = pho_iter->superCluster()->seed();
	Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*SCseed, *rechits);
	pho_smin.push_back(moments.sMin);
	pho_smaj.push_back(moments.sMaj);

      }
      else {

	pho_smin.push_back(def_max);
	pho_smaj.push_back(def_max);

      }
      
      pho_sinin_noiseclnd.push_back(pho_iter->full5x5_sigmaIetaIeta());
      pho_hoe.push_back(pho_iter->hadronicOverEm());
      pho_chargedhadroniso.push_back(pho_iter->chargedHadronIso());
      pho_neutralhadroniso.push_back(pho_iter->neutralHadronIso());
      pho_photoniso.push_back(pho_iter->photonIso());

      pho_n++;
    }
  } // End of photon header

  tree->Fill();
  clearVars();
}


void TriggerAnalyzerMiniAOD2024::beginJob() {
}

void TriggerAnalyzerMiniAOD2024::endJob() {
}

void TriggerAnalyzerMiniAOD2024::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TriggerAnalyzerMiniAOD2024::clearVars() {
  pho_e.clear();
  pho_pt.clear();
  pho_eta.clear();
  pho_phi.clear();
  pho_seedtime.clear();
  pho_smin.clear();
  pho_smaj.clear();
  pho_sinin_noiseclnd.clear();
  pho_hoe.clear();
  pho_chargedhadroniso.clear();
  pho_neutralhadroniso.clear();
  pho_photoniso.clear();
  
  dieg10time1ns_usfinfilt_pt.clear();
  dieg10time1ns_usfinfilt_eta.clear();
  dieg10time1ns_usfinfilt_phi.clear();

  pv_x.clear();
  pv_xerr.clear();
  pv_y.clear();
  pv_yerr.clear();
  pv_z.clear();
  pv_zerr.clear();
  pv_isvalid.clear();
};


//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzerMiniAOD2024);
