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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
#include "TTree.h"

class TriggerAnalyzerRAWMiniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit TriggerAnalyzerRAWMiniAOD(const edm::ParameterSet&);
  ~TriggerAnalyzerRAWMiniAOD();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void clearVars();

  edm::EDGetTokenT<edm::TriggerResults> trgResultsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjectsToken_;
  edm::EDGetTokenT< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > rechiteb_token;
  edm::EDGetTokenT< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > rechitee_token;
  //edm::EDGetTokenT<std::vector<pat::Electron> > electron_token;
  //edm::EDGetTokenT<std::vector<pat::Electron> > lowptelectron_token;
  edm::EDGetTokenT<std::vector<pat::Photon> > photon_token;
  //edm::EDGetTokenT<std::vector<pat::Photon> > ootphoton_token;
  edm::EDGetTokenT<reco::BeamSpot> BS_token;
  edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;

  edm::Service<TFileService> fs;
  TTree* tree;
  int run;
  int lumSec;

  bool HLT_DiPhoton10sminlt0p12, HLT_DiPhoton10Time1p4ns, HLT_DiPhoton10_CaloIdL;
  int dieg10sminlt0p12_usfinfilt_n;
  vector<double> dieg10sminlt0p12_usfinfilt_pt;
  vector<double> dieg10sminlt0p12_usfinfilt_eta;
  vector<double> dieg10sminlt0p12_usfinfilt_phi;
  int dieg10time1p4ns_usfinfilt_n;
  vector<double> dieg10time1p4ns_usfinfilt_pt;
  vector<double> dieg10time1p4ns_usfinfilt_eta;
  vector<double> dieg10time1p4ns_usfinfilt_phi;
  int dieg10caloidl_usfinfilt_n;
  vector<double> dieg10caloidl_usfinfilt_pt;
  vector<double> dieg10caloidl_usfinfilt_eta;
  vector<double> dieg10caloidl_usfinfilt_phi;
  int pho_n;
  vector<double> pho_pt;
  vector<double> pho_eta;
  vector<double> pho_phi;
  vector<double> pho_seedtime;
  vector<double> pho_smin;
  vector<double> pho_smax;
  int pv_n;
  vector<double> pv_x;
  vector<double> pv_xerr;
  vector<double> pv_y;
  vector<double> pv_yerr;
  vector<double> pv_z;
  vector<double> pv_zerr;
  vector<bool> pv_isvalid;
  double bs_x;
  double bs_y;
  double bs_z;
};

TriggerAnalyzerRAWMiniAOD::TriggerAnalyzerRAWMiniAOD(const edm::ParameterSet& iConfig) {
  trigObjectsToken_=consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("slimmedPatTrigger::RECO"));  
  trgResultsToken_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT") );
  rechiteb_token = consumes< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> >(edm::InputTag("reducedEgamma:reducedEBRecHits") );
  rechitee_token = consumes< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> >(edm::InputTag("reducedEgamma:reducedEERecHits") );
  //electron_token = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons") );
  //lowptelectron_token = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedLowPtElectrons") );
  photon_token = consumes<std::vector<pat::Photon> >(edm::InputTag("slimmedPhotons") );
  //ootphoton_token = consumes<std::vector<pat::Photon> >(edm::InputTag("slimmedOOTPhotons") );
  BS_token = consumes<reco::BeamSpot> (edm::InputTag("offlineBeamSpot"));
  PV_token = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));
  
  usesResource("TFileService");

  tree = fs->make<TTree>("tree", "tree");
  tree->Branch("run", &run, "run/i");
  tree->Branch("lumSec", &lumSec, "lumSec/i");
  tree->Branch("HLT_DiPhoton10sminlt0p12", &HLT_DiPhoton10sminlt0p12, "HLT_DiPhoton10sminlt0p12/b");
  tree->Branch("HLT_DiPhoton10Time1p4ns", &HLT_DiPhoton10Time1p4ns, "HLT_DiPhoton10Time1p4ns/b");
  tree->Branch("HLT_DiPhoton10_CaloIdL", &HLT_DiPhoton10_CaloIdL, "HLT_DiPhoton10_CaloIdL/b");
  tree->Branch("dieg10sminlt0p12_usfinfilt_n", &dieg10sminlt0p12_usfinfilt_n, "dieg10sminlt0p12_usfinfilt_n/i");
  tree->Branch("dieg10sminlt0p12_usfinfilt_pt", &dieg10sminlt0p12_usfinfilt_pt);
  tree->Branch("dieg10sminlt0p12_usfinfilt_eta", &dieg10sminlt0p12_usfinfilt_eta);
  tree->Branch("dieg10sminlt0p12_usfinfilt_phi", &dieg10sminlt0p12_usfinfilt_phi);
  tree->Branch("dieg10time1p4ns_usfinfilt_n", &dieg10time1p4ns_usfinfilt_n, "dieg10time1p4ns_usfinfilt_n/i");
  tree->Branch("dieg10time1p4ns_usfinfilt_pt", &dieg10time1p4ns_usfinfilt_pt);
  tree->Branch("dieg10time1p4ns_usfinfilt_eta", &dieg10time1p4ns_usfinfilt_eta);
  tree->Branch("dieg10time1p4ns_usfinfilt_phi", &dieg10time1p4ns_usfinfilt_phi);
  tree->Branch("dieg10caloidl_usfinfilt_n", &dieg10caloidl_usfinfilt_n, "dieg10caloidl_usfinfilt_n/i");
  tree->Branch("dieg10caloidl_usfinfilt_pt", &dieg10caloidl_usfinfilt_pt);
  tree->Branch("dieg10caloidl_usfinfilt_eta", &dieg10caloidl_usfinfilt_eta);
  tree->Branch("dieg10caloidl_usfinfilt_phi", &dieg10caloidl_usfinfilt_phi);
  tree->Branch("pho_n", &pho_n, "pho_n/i");
  tree->Branch("pho_pt", &pho_pt);
  tree->Branch("pho_eta", &pho_eta);
  tree->Branch("pho_phi", &pho_phi);
  tree->Branch("pho_seedtime", &pho_seedtime);
  tree->Branch("pho_smin", &pho_smin);
  tree->Branch("pho_smax", &pho_smax);
  tree->Branch("pv_n", &pv_n, "pv_n/i");
  tree->Branch("pv_x", &pv_x);
  tree->Branch("pv_xerr", &pv_xerr);
  tree->Branch("pv_y", &pv_y);
  tree->Branch("pv_yerr", &pv_yerr);
  tree->Branch("pv_z", &pv_z);
  tree->Branch("pv_zerr", &pv_zerr);
  tree->Branch("pv_isvalid", &pv_isvalid);
  tree->Branch("bs_x", &bs_x, "bs_x/d");
  tree->Branch("bs_y", &bs_y, "bs_y/d");
  tree->Branch("bs_z", &bs_z, "bs_z/d");
}


TriggerAnalyzerRAWMiniAOD::~TriggerAnalyzerRAWMiniAOD() {
}

void TriggerAnalyzerRAWMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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

  //HLT_DiPhoton10sminlt0p1 = false;
  HLT_DiPhoton10sminlt0p12 = false;
  //HLT_DiPhoton10Time2ns = false;
  //HLT_DiPhoton10Time1p8ns = false;
  //HLT_DiPhoton10Time1p6ns = false;
  HLT_DiPhoton10Time1p4ns = false;
  HLT_DiPhoton10_CaloIdL = false;

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
      }
    }
    //if(HLT_DiPhoton10sminlt0p12 || HLT_DiPhoton10Time1p4ns || HLT_DiPhoton10_CaloIdL) cout<<"Passing one of the desired triggers!"<<endl;
  } // End of loop for accessing the trigger bits

  //Accessing the trigger objects in MINIAOD
  //This recipe works for MINIAOD only
  dieg10sminlt0p12_usfinfilt_n = 0;
  dieg10time1p4ns_usfinfilt_n = 0;
  dieg10caloidl_usfinfilt_n = 0;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(trigObjectsToken_, triggerObjects);
  const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackFilterLabels(iEvent,*trigResults);
    obj.unpackPathNames(names);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
      string myfillabl=obj.filterLabels()[h];
      //cout << "Trigger object name, pt, eta, phi: "
      //<< myfillabl<<", " << obj.pt()<<", "<<obj.eta()<<", "<<obj.phi() << endl;
      if(myfillabl.compare("hltDiEG10CaloIdLsminlt0p12UnseededFilter")==0) {
	dieg10sminlt0p12_usfinfilt_pt.push_back(obj.pt());
	dieg10sminlt0p12_usfinfilt_eta.push_back(obj.eta());
	dieg10sminlt0p12_usfinfilt_phi.push_back(obj.phi());
	dieg10sminlt0p12_usfinfilt_n++;
      }
      if(myfillabl.compare("hltDiEG10CaloIdLTime1p4nsUnseededFilter")==0) {
	dieg10time1p4ns_usfinfilt_pt.push_back(obj.pt());
	dieg10time1p4ns_usfinfilt_eta.push_back(obj.eta());
	dieg10time1p4ns_usfinfilt_phi.push_back(obj.phi());
	dieg10time1p4ns_usfinfilt_n++;
      }
      if(myfillabl.compare("hltDiEG10CaloIdLClusterShapeUnseededFilter")==0) {
	dieg10caloidl_usfinfilt_pt.push_back(obj.pt());
	dieg10caloidl_usfinfilt_eta.push_back(obj.eta());
	dieg10caloidl_usfinfilt_phi.push_back(obj.phi());
	dieg10caloidl_usfinfilt_n++;
      }
    }
  } // End of loop for accessing the trigger objects

  // ECAL renhits  
  edm::Handle< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > rechitebH;
  iEvent.getByToken(rechiteb_token, rechitebH);
  edm::Handle< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > rechiteeH;
  iEvent.getByToken(rechitee_token, rechiteeH);

  // Photon
  edm::Handle<std::vector<pat::Photon> > photonH;
  iEvent.getByToken(photon_token, photonH);
  pho_n = 0;
  if(photonH.isValid()) {
    double seedtime = DBL_MAX;
    //double seedtime2 = DBL_MAX;
    for(auto pho_iter=photonH->begin(); pho_iter!=photonH->end(); pho_iter++) {
      pho_pt.push_back(pho_iter->pt());
      pho_eta.push_back(pho_iter->eta());
      pho_phi.push_back(pho_iter->phi());

      auto& ssFull5x5 = pho_iter->full5x5_showerShapeVariables();
      pho_smin.push_back(ssFull5x5.smMinor);
      pho_smax.push_back(ssFull5x5.smMajor);

      seedtime = DBL_MAX;
      //seedtime2 = DBL_MAX;
      DetId SCseedID = pho_iter->seed()->seed();
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
      /*auto rechits2 = pho_iter->recHits();
      auto rechitseed2 = rechits2->find(SCseedID);
      if(rechitseed2!=rechits2->end()) {
	seedtime2 = rechitseed2->time();
      }
      std::cout<<"Compare seedtime: "<<seedtime<<"with seedtime2:"<<seedtime2<<std::endl;*/
      pho_seedtime.push_back(seedtime);
    }
  }

  tree->Fill();
  clearVars();
}


void TriggerAnalyzerRAWMiniAOD::beginJob() {
}

void TriggerAnalyzerRAWMiniAOD::endJob() {
}

void TriggerAnalyzerRAWMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TriggerAnalyzerRAWMiniAOD::clearVars() {
  dieg10sminlt0p12_usfinfilt_pt.clear();
  dieg10sminlt0p12_usfinfilt_eta.clear();
  dieg10sminlt0p12_usfinfilt_phi.clear();
  dieg10time1p4ns_usfinfilt_pt.clear();
  dieg10time1p4ns_usfinfilt_eta.clear();
  dieg10time1p4ns_usfinfilt_phi.clear();
  dieg10caloidl_usfinfilt_pt.clear();
  dieg10caloidl_usfinfilt_eta.clear();
  dieg10caloidl_usfinfilt_phi.clear();
  pho_pt.clear();
  pho_eta.clear();
  pho_phi.clear();
  pho_seedtime.clear();
  pho_smin.clear();
  pv_x.clear();
  pv_xerr.clear();
  pv_y.clear();
  pv_yerr.clear();
  pv_z.clear();
  pv_zerr.clear();
  pv_isvalid.clear();
};


//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzerRAWMiniAOD);
