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

class QuickHistoMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit QuickHistoMaker(const edm::ParameterSet&);
  ~QuickHistoMaker();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void clearVars();

  edm::EDGetTokenT<edm::TriggerResults> trgResultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigEventToken_;

  edm::Service<TFileService> fs;
  TTree* tree;
  int run;
  int lumSec;

  bool HLT_DiPhoton10sminlt0p12, HLT_DiPhoton10Time1p4ns, HLT_DiPhoton10_CaloIdL;
  int dieg10sminlt0p12_ptfilt_n;
  vector<double> dieg10sminlt0p12_ptfilt_pt;
  vector<double> dieg10sminlt0p12_ptfilt_eta;
  vector<double> dieg10sminlt0p12_ptfilt_phi;
  int dieg10sminlt0p12_etafilt_n;
  vector<double> dieg10sminlt0p12_etafilt_pt;
  vector<double> dieg10sminlt0p12_etafilt_eta;
  vector<double> dieg10sminlt0p12_etafilt_phi;
  int dieg10sminlt0p12_usptfilt_n;
  vector<double> dieg10sminlt0p12_usptfilt_pt;
  vector<double> dieg10sminlt0p12_usptfilt_eta;
  vector<double> dieg10sminlt0p12_usptfilt_phi;
  int dieg10sminlt0p12_usetafilt_n;
  vector<double> dieg10sminlt0p12_usetafilt_pt;
  vector<double> dieg10sminlt0p12_usetafilt_eta;
  vector<double> dieg10sminlt0p12_usetafilt_phi;
  int dieg10sminlt0p12_usfinfilt_n;
  vector<double> dieg10sminlt0p12_usfinfilt_pt;
  vector<double> dieg10sminlt0p12_usfinfilt_eta;
  vector<double> dieg10sminlt0p12_usfinfilt_phi;
  int dieg10time1p4ns_usfinfilt_n;
  vector<double> dieg10time1p4ns_usfinfilt_pt;
  vector<double> dieg10time1p4ns_usfinfilt_eta;
  vector<double> dieg10time1p4ns_usfinfilt_phi;
};

QuickHistoMaker::QuickHistoMaker(const edm::ParameterSet& iConfig) {
  trigEventToken_= consumes<trigger::TriggerEvent>(edm::InputTag("hltTriggerSummaryAOD::HLT2"));  
  trgResultsToken_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT2") );
  
  usesResource("TFileService");

  tree = fs->make<TTree>("tree", "tree");
  tree->Branch("run", &run, "run/i");
  tree->Branch("lumSec", &lumSec, "lumSec/i");
  tree->Branch("HLT_DiPhoton10sminlt0p12", &HLT_DiPhoton10sminlt0p12, "HLT_DiPhoton10sminlt0p12/b");
  tree->Branch("HLT_DiPhoton10Time1p4ns", &HLT_DiPhoton10Time1p4ns, "HLT_DiPhoton10Time1p4ns/b");
  tree->Branch("HLT_DiPhoton10_CaloIdL", &HLT_DiPhoton10_CaloIdL, "HLT_DiPhoton10_CaloIdL/b");
  tree->Branch("dieg10sminlt0p12_ptfilt_n", &dieg10sminlt0p12_ptfilt_n, "dieg10sminlt0p12_ptfilt_n/i");
  tree->Branch("dieg10sminlt0p12_ptfilt_pt", &dieg10sminlt0p12_ptfilt_pt);
  tree->Branch("dieg10sminlt0p12_ptfilt_eta", &dieg10sminlt0p12_ptfilt_eta);
  tree->Branch("dieg10sminlt0p12_ptfilt_phi", &dieg10sminlt0p12_ptfilt_phi);
  tree->Branch("dieg10sminlt0p12_etafilt_n", &dieg10sminlt0p12_etafilt_n, "dieg10sminlt0p12_etafilt_n/i");
  tree->Branch("dieg10sminlt0p12_etafilt_pt", &dieg10sminlt0p12_etafilt_pt);
  tree->Branch("dieg10sminlt0p12_etafilt_eta", &dieg10sminlt0p12_etafilt_eta);
  tree->Branch("dieg10sminlt0p12_etafilt_phi", &dieg10sminlt0p12_etafilt_phi);
  tree->Branch("dieg10sminlt0p12_usptfilt_n", &dieg10sminlt0p12_usptfilt_n, "dieg10sminlt0p12_usptfilt_n/i");
  tree->Branch("dieg10sminlt0p12_usptfilt_pt", &dieg10sminlt0p12_usptfilt_pt);
  tree->Branch("dieg10sminlt0p12_usptfilt_eta", &dieg10sminlt0p12_usptfilt_eta);
  tree->Branch("dieg10sminlt0p12_usptfilt_phi", &dieg10sminlt0p12_usptfilt_phi);
  tree->Branch("dieg10sminlt0p12_usetafilt_n", &dieg10sminlt0p12_usetafilt_n, "dieg10sminlt0p12_usetafilt_n/i");
  tree->Branch("dieg10sminlt0p12_usetafilt_pt", &dieg10sminlt0p12_usetafilt_pt);
  tree->Branch("dieg10sminlt0p12_usetafilt_eta", &dieg10sminlt0p12_usetafilt_eta);
  tree->Branch("dieg10sminlt0p12_usetafilt_phi", &dieg10sminlt0p12_usetafilt_phi);
  tree->Branch("dieg10sminlt0p12_usfinfilt_n", &dieg10sminlt0p12_usfinfilt_n, "dieg10sminlt0p12_usfinfilt_n/i");
  tree->Branch("dieg10sminlt0p12_usfinfilt_pt", &dieg10sminlt0p12_usfinfilt_pt);
  tree->Branch("dieg10sminlt0p12_usfinfilt_eta", &dieg10sminlt0p12_usfinfilt_eta);
  tree->Branch("dieg10sminlt0p12_usfinfilt_phi", &dieg10sminlt0p12_usfinfilt_phi);
  tree->Branch("dieg10time1p4ns_usfinfilt_n", &dieg10time1p4ns_usfinfilt_n, "dieg10time1p4ns_usfinfilt_n/i");
  tree->Branch("dieg10time1p4ns_usfinfilt_pt", &dieg10time1p4ns_usfinfilt_pt);
  tree->Branch("dieg10time1p4ns_usfinfilt_eta", &dieg10time1p4ns_usfinfilt_eta);
  tree->Branch("dieg10time1p4ns_usfinfilt_phi", &dieg10time1p4ns_usfinfilt_phi);
}


QuickHistoMaker::~QuickHistoMaker() {
}

void QuickHistoMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;

  run = iEvent.eventAuxiliary().run();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();

  HLT_DiPhoton10sminlt0p12 = false;
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
  } // End of loop for accessing the trigger bits

  //Accessing the trigger objects in RAW
  //This recipe works for RAW or RECO only
  dieg10sminlt0p12_ptfilt_n = 0;
  dieg10sminlt0p12_etafilt_n = 0;
  dieg10sminlt0p12_usptfilt_n = 0;
  dieg10sminlt0p12_usetafilt_n = 0;
  dieg10sminlt0p12_usfinfilt_n = 0;
  dieg10time1p4ns_usfinfilt_n = 0;

  edm::Handle<trigger::TriggerEvent> triggerEvt;
  iEvent.getByToken(trigEventToken_, triggerEvt);
  if(triggerEvt.isValid()) {
    trigger::TriggerObjectCollection allTriggerObjects = triggerEvt->getObjects();

    size_t hltEG10EtIndex = (*triggerEvt).filterIndex( edm::InputTag("hltEG10EtFilter","","HLT2") );
    if(hltEG10EtIndex < (*triggerEvt).sizeFilters()) {
      const trigger::Keys &keys = (*triggerEvt).filterKeys(hltEG10EtIndex);
      for(size_t j=0; j<keys.size(); j++) {
	trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
	dieg10sminlt0p12_ptfilt_pt.push_back(foundObject.pt());
	dieg10sminlt0p12_ptfilt_eta.push_back(foundObject.eta());
	dieg10sminlt0p12_ptfilt_phi.push_back(foundObject.phi());
	dieg10sminlt0p12_ptfilt_n++;
	//cout<<"hltEG10EtFilter["<<j<<"]: "<<foundObject.pt()<<"\t"<<foundObject.eta()<<"\t"<<foundObject.phi()<<endl;
      }
    }

    size_t hltEG10EtaIndex = (*triggerEvt).filterIndex( edm::InputTag("hltEG10EtaFilter","","HLT2") );
    if(hltEG10EtaIndex < (*triggerEvt).sizeFilters()) {
      const trigger::Keys &keys = (*triggerEvt).filterKeys(hltEG10EtaIndex);
      for(size_t j=0; j<keys.size(); j++) {
	trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
	dieg10sminlt0p12_etafilt_pt.push_back(foundObject.pt());
	dieg10sminlt0p12_etafilt_eta.push_back(foundObject.eta());
	dieg10sminlt0p12_etafilt_phi.push_back(foundObject.phi());
	dieg10sminlt0p12_etafilt_n++;
	//cout<<"hltEG10EtaFilter["<<j<<"]: "<<foundObject.pt()<<"\t"<<foundObject.eta()<<"\t"<<foundObject.phi()<<endl;
      }
    }

    size_t hltDiEG10EtUnseededIndex = (*triggerEvt).filterIndex( edm::InputTag("hltDiEG10EtUnseededFilter","","HLT2") );
    if(hltDiEG10EtUnseededIndex < (*triggerEvt).sizeFilters()) {
      const trigger::Keys &keys = (*triggerEvt).filterKeys(hltDiEG10EtUnseededIndex);
      for(size_t j=0; j<keys.size(); j++) {
	trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
	dieg10sminlt0p12_usptfilt_pt.push_back(foundObject.pt());
	dieg10sminlt0p12_usptfilt_eta.push_back(foundObject.eta());
	dieg10sminlt0p12_usptfilt_phi.push_back(foundObject.phi());
	dieg10sminlt0p12_usptfilt_n++;
	//cout<<"hltEG10EtUnseededFilter["<<j<<"]: "<<foundObject.pt()<<"\t"<<foundObject.eta()<<"\t"<<foundObject.phi()<<endl;
      }
    }

    size_t hltDiEG10EtaUnseededIndex = (*triggerEvt).filterIndex( edm::InputTag("hltDiEG10EtaUnseededFilter","","HLT2") );
    if(hltDiEG10EtaUnseededIndex < (*triggerEvt).sizeFilters()) {
      const trigger::Keys &keys = (*triggerEvt).filterKeys(hltDiEG10EtaUnseededIndex);
      for(size_t j=0; j<keys.size(); j++) {
	trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
	dieg10sminlt0p12_usetafilt_pt.push_back(foundObject.pt());
	dieg10sminlt0p12_usetafilt_eta.push_back(foundObject.eta());
	dieg10sminlt0p12_usetafilt_phi.push_back(foundObject.phi());
	dieg10sminlt0p12_usetafilt_n++;
	//cout<<"hltEG10EtaUnseededFilter["<<j<<"]: "<<foundObject.pt()<<"\t"<<foundObject.eta()<<"\t"<<foundObject.phi()<<endl;
      }
    }

    size_t hltDiEG10CaloIdLsminlt0p12UnseededIndex = (*triggerEvt).filterIndex( edm::InputTag("hltDiEG10CaloIdLsminlt0p12UnseededFilter","","HLT2") );
    if(hltDiEG10CaloIdLsminlt0p12UnseededIndex < (*triggerEvt).sizeFilters()) {
      const trigger::Keys &keys = (*triggerEvt).filterKeys(hltDiEG10CaloIdLsminlt0p12UnseededIndex);
      for(size_t j=0; j<keys.size(); j++) {
	trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
	dieg10sminlt0p12_usfinfilt_pt.push_back(foundObject.pt());
	dieg10sminlt0p12_usfinfilt_eta.push_back(foundObject.eta());
	dieg10sminlt0p12_usfinfilt_phi.push_back(foundObject.phi());
	dieg10sminlt0p12_usfinfilt_n++;
	//cout<<"hltDiEG10CaloIdLsminlt0p12UnseededFilter["<<j<<"]: "<<foundObject.pt()<<"\t"<<foundObject.eta()<<"\t"<<foundObject.phi()<<endl;
      }
    }

    size_t hltDiEG10CaloIdLTime1p4nsUnseededIndex = (*triggerEvt).filterIndex( edm::InputTag("hltDiEG10CaloIdLTime1p4nsUnseededFilter","","HLT2") );
    if(hltDiEG10CaloIdLTime1p4nsUnseededIndex < (*triggerEvt).sizeFilters()) {
      const trigger::Keys &keys = (*triggerEvt).filterKeys(hltDiEG10CaloIdLTime1p4nsUnseededIndex);
      for(size_t j=0; j<keys.size(); j++) {
	trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
	dieg10time1p4ns_usfinfilt_pt.push_back(foundObject.pt());
        dieg10time1p4ns_usfinfilt_eta.push_back(foundObject.eta());
        dieg10time1p4ns_usfinfilt_phi.push_back(foundObject.phi());
        dieg10time1p4ns_usfinfilt_n++;
      }
    }

  }
  
  tree->Fill();
  clearVars();
}


void QuickHistoMaker::beginJob() {
}

void QuickHistoMaker::endJob() {
}

void QuickHistoMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void QuickHistoMaker::clearVars() {
  dieg10sminlt0p12_ptfilt_pt.clear();
  dieg10sminlt0p12_ptfilt_eta.clear();
  dieg10sminlt0p12_ptfilt_phi.clear();
  dieg10sminlt0p12_etafilt_pt.clear();
  dieg10sminlt0p12_etafilt_eta.clear();
  dieg10sminlt0p12_etafilt_phi.clear();
  dieg10sminlt0p12_ptfilt_pt.clear();
  dieg10sminlt0p12_ptfilt_eta.clear();
  dieg10sminlt0p12_ptfilt_phi.clear();
  dieg10sminlt0p12_etafilt_pt.clear();
  dieg10sminlt0p12_etafilt_eta.clear();
  dieg10sminlt0p12_etafilt_phi.clear();
  dieg10sminlt0p12_usfinfilt_pt.clear();
  dieg10sminlt0p12_usfinfilt_eta.clear();
  dieg10sminlt0p12_usfinfilt_phi.clear();
  dieg10time1p4ns_usfinfilt_pt.clear();
  dieg10time1p4ns_usfinfilt_eta.clear();
  dieg10time1p4ns_usfinfilt_phi.clear();
};


//define this as a plug-in
DEFINE_FWK_MODULE(QuickHistoMaker);
