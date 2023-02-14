// system include files
#include <memory>
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

class TriggerAnalyzerMiniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit TriggerAnalyzerMiniAOD(const edm::ParameterSet&);
  ~TriggerAnalyzerMiniAOD();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void clearVars();

  edm::EDGetTokenT<edm::TriggerResults> trgResultsToken_;
  edm::EDGetTokenT< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > rechiteb_token;
  edm::EDGetTokenT< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > rechitee_token;
  edm::EDGetTokenT<std::vector<pat::Electron> > electron_token;
  edm::EDGetTokenT<std::vector<reco::SuperCluster> > suprclus_token;
  //edm::EDGetTokenT<std::vector<pat::Electron> > lowptelectron_token;
  //edm::EDGetTokenT<std::vector<pat::Photon> > photon_token;
  //edm::EDGetTokenT<std::vector<pat::Photon> > ootphoton_token;
  edm::EDGetTokenT<reco::BeamSpot> BS_token;
  edm::EDGetTokenT<double> rho_token;
  edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;

  edm::Service<TFileService> fs;
  TTree* tree;
  int run;
  int lumSec;

  bool HLT_DiPhoton10sminlt0p12, HLT_DiPhoton10Time1p4ns, HLT_DiPhoton10_CaloIdL;
  bool HLTOR_METTrig;

  double bs_x;
  double bs_y;
  double bs_z;
  double rho;

  int ele_n;
  vector<double> ele_e;
  vector<double> ele_pt;
  vector<double> ele_eta;
  vector<double> ele_phi;
  vector<double> ele_d0;
  vector<double> ele_dz;
  vector<double> ele_seedtime;
  vector<double> ele_smin;
  vector<double> ele_smaj;
  vector<double> ele_sinin_noiseclnd;
  vector<double> ele_detaseed;
  vector<double> ele_dphiin;
  vector<double> ele_hoe;
  vector<double> ele_chargedhadroniso;
  vector<double> ele_neutralhadroniso;
  vector<double> ele_photoniso;
  vector<double> ele_ooemoop;
  vector<double> ele_innerhits;
  vector<bool> ele_convveto;

  int pv_n;
  vector<double> pv_x;
  vector<double> pv_xerr;
  vector<double> pv_y;
  vector<double> pv_yerr;
  vector<double> pv_z;
  vector<double> pv_zerr;
  vector<bool> pv_isvalid;
};

TriggerAnalyzerMiniAOD::TriggerAnalyzerMiniAOD(const edm::ParameterSet& iConfig) {
  trgResultsToken_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT") );
  rechiteb_token = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> >(edm::InputTag("reducedEgamma:reducedEBRecHits") );
  rechitee_token = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> >(edm::InputTag("reducedEgamma:reducedEERecHits") );
  electron_token = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons") );
  suprclus_token = consumes<std::vector<reco::SuperCluster> >(edm::InputTag("reducedEgamma:reducedSuperClusters") );
  //lowptelectron_token = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedLowPtElectrons") );
  //photon_token = consumes<std::vector<pat::Photon> >(edm::InputTag("slimmedPhotons") );
  //ootphoton_token = consumes<std::vector<pat::Photon> >(edm::InputTag("slimmedOOTPhotons") );
  BS_token = consumes<reco::BeamSpot> (edm::InputTag("offlineBeamSpot"));
  rho_token = consumes<double> (edm::InputTag("fixedGridRhoAll"));
  PV_token = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));
  
  usesResource("TFileService");

  tree = fs->make<TTree>("tree", "tree");

  tree->Branch("run", &run, "run/I");
  tree->Branch("lumSec", &lumSec, "lumSec/I");

  tree->Branch("HLT_DiPhoton10sminlt0p12", &HLT_DiPhoton10sminlt0p12, "HLT_DiPhoton10sminlt0p12/O");
  tree->Branch("HLT_DiPhoton10Time1p4ns", &HLT_DiPhoton10Time1p4ns, "HLT_DiPhoton10Time1p4ns/O");
  tree->Branch("HLT_DiPhoton10_CaloIdL", &HLT_DiPhoton10_CaloIdL, "HLT_DiPhoton10_CaloIdL/O");
  tree->Branch("HLTOR_METTrig", &HLTOR_METTrig, "HLTOR_METTrig/O");

  tree->Branch("bs_x", &bs_x, "bs_x/D");
  tree->Branch("bs_y", &bs_y, "bs_y/D");
  tree->Branch("bs_z", &bs_z, "bs_z/D");
  tree->Branch("rho", &rho, "rho/D");

  tree->Branch("ele_n", &ele_n, "ele_n/I");
  tree->Branch("ele_e", &ele_e);
  tree->Branch("ele_pt", &ele_pt);
  tree->Branch("ele_eta", &ele_eta);
  tree->Branch("ele_phi", &ele_phi);
  tree->Branch("ele_d0", &ele_d0);
  tree->Branch("ele_dz", &ele_dz);
  tree->Branch("ele_seedtime", &ele_seedtime);
  tree->Branch("ele_smin", &ele_smin);
  tree->Branch("ele_smaj", &ele_smaj);
  tree->Branch("ele_sinin_noiseclnd", &ele_sinin_noiseclnd);
  tree->Branch("ele_detaseed", &ele_detaseed);
  tree->Branch("ele_dphiin", &ele_dphiin);
  tree->Branch("ele_hoe", &ele_hoe);
  tree->Branch("ele_chargedhadroniso", &ele_chargedhadroniso);
  tree->Branch("ele_neutralhadroniso", &ele_neutralhadroniso);
  tree->Branch("ele_photoniso", &ele_photoniso);
  tree->Branch("ele_ooemoop", &ele_ooemoop);
  tree->Branch("ele_innerhits", &ele_innerhits);
  tree->Branch("ele_convveto", &ele_convveto);

  tree->Branch("pv_n", &pv_n, "pv_n/I");
  tree->Branch("pv_x", &pv_x);
  tree->Branch("pv_xerr", &pv_xerr);
  tree->Branch("pv_y", &pv_y);
  tree->Branch("pv_yerr", &pv_yerr);
  tree->Branch("pv_z", &pv_z);
  tree->Branch("pv_zerr", &pv_zerr);
  tree->Branch("pv_isvalid", &pv_isvalid);
}


TriggerAnalyzerMiniAOD::~TriggerAnalyzerMiniAOD() {
}

void TriggerAnalyzerMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;

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
    bs_x = DBL_MAX;
    bs_y = DBL_MAX;
    bs_z = DBL_MAX;
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

  HLT_DiPhoton10sminlt0p12 = false;
  HLT_DiPhoton10Time1p4ns = false;
  HLT_DiPhoton10_CaloIdL = false;
  HLTOR_METTrig = false;

  // MET Trigger bits
  bool HLT_PFMET120_PFMHT120_IDTight=false, HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF=false, HLT_PFMETNoMu120_PFMHTNoMu120_IDTight=false, HLT_CaloMET80_NotCleaned=false, HLT_PFMET200_NotCleaned=false, HLT_PFMET200_BeamHaloCleaned=false, HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight=false;

  //Accessing trigger bits:
  //This works in both RAW, AOD or MINIAOD 
  //Here we access the decision provided by the HLT (i.e. original trigger step). 
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgResultsToken_, trigResults);
  if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    //cout<<"Number of triggers: "<<N_Triggers<<endl;
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trigResults.product()->accept(i_Trig)) {
	//cout << "Path: " <<trigName.triggerName(i_Trig)<<"Results: "<<trigResults.product()->accept(i_Trig)<<endl;
	TString TrigPath =trigName.triggerName(i_Trig);
	if(TrigPath.Index("HLT_DiPhoton10sminlt0p12_v") >=0) HLT_DiPhoton10sminlt0p12 = true; 
	if(TrigPath.Index("HLT_DiPhoton10Time1p4ns_v") >=0) HLT_DiPhoton10Time1p4ns = true; 
	if(TrigPath.Index("HLT_DiPhoton10_CaloIdL_v") >=0) HLT_DiPhoton10_CaloIdL = true; 
	if(TrigPath.Index("HLT_PFMET120_PFMHT120_IDTight_v") >=0) HLT_PFMET120_PFMHT120_IDTight=true;
	if(TrigPath.Index("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF_v") >=0) HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF=true;
	if(TrigPath.Index("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") >=0) HLT_PFMETNoMu120_PFMHTNoMu120_IDTight=true;
	if(TrigPath.Index("HLT_CaloMET80_NotCleaned_v") >=0) HLT_CaloMET80_NotCleaned=true;
	if(TrigPath.Index("HLT_PFMET200_NotCleaned_v") >=0) HLT_PFMET200_NotCleaned=true;
	if(TrigPath.Index("HLT_PFMET200_BeamHaloCleaned_v") >=0) HLT_PFMET200_BeamHaloCleaned=true;
	if(TrigPath.Index("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v") >=0) HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight=true;
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

  // ECAL rechits  
  edm::Handle< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > rechitebH;
  iEvent.getByToken(rechiteb_token, rechitebH);
  edm::Handle< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > rechiteeH;
  iEvent.getByToken(rechitee_token, rechiteeH);

  // Super Cluster collection
  edm::Handle< std::vector<reco::SuperCluster> > sclusterH;
  iEvent.getByToken(suprclus_token, sclusterH);

  // Electron
  edm::Handle<std::vector<pat::Electron> > electronH;
  iEvent.getByToken(electron_token, electronH);
  ele_n = 0;
  if(electronH.isValid()) {
    double seedtime = DBL_MAX;
    for(auto ele_iter=electronH->begin(); ele_iter!=electronH->end(); ele_iter++) {
      ele_e.push_back(ele_iter->energy());
      ele_pt.push_back(ele_iter->pt());
      ele_eta.push_back(ele_iter->eta());
      ele_phi.push_back(ele_iter->phi());
      ele_d0.push_back(ele_iter->dB(pat::Electron::PV3D));
      ele_dz.push_back(ele_iter->dB(pat::Electron::PVDZ));
      
      seedtime = DBL_MAX;
      DetId SCseedID = ele_iter->seed()->seed();
      if(rechitebH.isValid() && seedtime==DBL_MAX) {
	auto rechitseed = rechitebH->find(SCseedID);
	if(rechitseed!=rechitebH->end()) {
	  seedtime = rechitseed->time();
	}
      }
      if(rechiteeH.isValid() && seedtime==DBL_MAX) {
	auto rechitseed = rechiteeH->find(SCseedID);
	if(rechitseed!=rechiteeH->end()) {
	  seedtime = rechitseed->time();
	}
      }
      ele_seedtime.push_back(seedtime);

      if( (sclusterH.isValid()) &&
	  ( ( (std::abs(ele_iter->eta())<1.479) && rechitebH.isValid() ) || 
	    ( (std::abs(ele_iter->eta())>=1.479) && rechiteeH.isValid() ) 
	    ) 
	  ) {

	const EcalRecHitCollection* rechits = (std::abs(ele_iter->eta()) < 1.479) ? rechitebH.product() : rechiteeH.product();
	reco::CaloClusterPtr SCseed = ele_iter->superCluster()->seed();
	Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*SCseed, *rechits);
	ele_smin.push_back(moments.sMin);
	ele_smaj.push_back(moments.sMaj);

      }
      else {

	ele_smin.push_back(DBL_MAX);
	ele_smaj.push_back(DBL_MAX);

      }
      
      ele_sinin_noiseclnd.push_back(ele_iter->full5x5_sigmaIetaIeta());
      ele_detaseed.push_back(ele_iter->deltaEtaSuperClusterTrackAtVtx());
      ele_dphiin.push_back(ele_iter->deltaPhiSuperClusterTrackAtVtx());
      ele_hoe.push_back(ele_iter->hadronicOverEm());
      ele_chargedhadroniso.push_back(ele_iter->chargedHadronIso());
      ele_neutralhadroniso.push_back(ele_iter->neutralHadronIso());
      ele_photoniso.push_back(ele_iter->photonIso());
      ele_ooemoop.push_back( (1.0/ele_iter->correctedEcalEnergy()) - (1.0/ele_iter->p()) );
      ele_innerhits.push_back(ele_iter->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
      ele_convveto.push_back(ele_iter->passConversionVeto());

      ele_n++;
    }
  } // End of electron header

  tree->Fill();
  clearVars();
}


void TriggerAnalyzerMiniAOD::beginJob() {
}

void TriggerAnalyzerMiniAOD::endJob() {
}

void TriggerAnalyzerMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TriggerAnalyzerMiniAOD::clearVars() {
  ele_e.clear();
  ele_pt.clear();
  ele_eta.clear();
  ele_phi.clear();
  ele_d0.clear();
  ele_dz.clear();
  ele_seedtime.clear();
  ele_smin.clear();
  ele_smaj.clear();
  ele_sinin_noiseclnd.clear();
  ele_detaseed.clear();
  ele_dphiin.clear();
  ele_hoe.clear();
  ele_chargedhadroniso.clear();
  ele_neutralhadroniso.clear();
  ele_photoniso.clear();
  ele_ooemoop.clear();
  ele_innerhits.clear();
  ele_convveto.clear();
  
  pv_x.clear();
  pv_xerr.clear();
  pv_y.clear();
  pv_yerr.clear();
  pv_z.clear();
  pv_zerr.clear();
  pv_isvalid.clear();
};


//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzerMiniAOD);
