#ifndef DATA_ROBUSTANALYZER_H
#define DATA_ROBUSTANALYZER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"

using namespace std;

class data_robustanalyzer {

 public:
  data_robustanalyzer(TString, TString, bool);
  ~data_robustanalyzer();
  
  void analyzersinglefile();
  int doGenMatchingUnseeded(vector<int>);
  void addgenhist(TString);
  void addhist(TString);
  void addhistunseeded(TString);
  void addhistcomparegenrecounseeded(TString);
  void fillgenhistinevent(TString, vector<int>);
  void fillhistinevent(TString, vector<int>);
  void fillhistineventunseeded(TString, vector<int>);
  void fillhistcomparegenrecounseeded(TString, vector<int>);
  bool isL1EgSeeded(int);
  void sort(int*, double*, int);
  
 private:
  
  TChain* inputChain;

  TFile* outfile;

  bool isMC;
  int genLepN, egRecoN, egusRecoN, l1FiltN, l1egObjN;
  double genLepPid[100], genLepPt[100], genLepEta[100], genLepPhi[100], genLepVx[100], genLepVy[100], genLepVz[100], genLepNMom[100], genLepMomPid[100], genLepTimeLight[100], genLepTimeAct[100], genLepTimeDiff[100];
  double l1FiltPt[100], l1FiltEta[100], l1FiltPhi[100];
  double l1egObjPt[100], l1egObjEta[100], l1egObjPhi[100];
  double egRecoPt[100], egRecoEta[100], egRecoPhi[100];
  double egusRecoPt[100], egusRecoEta[100], egusRecoPhi[100];
  double eghltEgammaClusterShape[100], eghltEgammaClusterShape_sigmaIEtaIEta5x5[100], eghltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[100], eghltEgammaEcalPFClusterIso[100], eghltEgammaHcalPFClusterIso[100], eghltEgammaHoverE[100], eghltEgammaSuperClusterEnergy[100], eghltEcalSeedClusterTime[100];
  double egushltEgammaClusterShape[100], egushltEgammaClusterShape_sigmaIEtaIEta5x5[100], egushltEgammaClusterShape_sigmaIEtaIEta5x5NoiseCleaned[100], egushltEgammaEcalPFClusterIso[100], egushltEgammaHcalPFClusterIso[100], egushltEgammaHoverE[100], egushltEgammaSuperClusterEnergy[100], egushltEcalSeedClusterTime[100];
  double eghltEgammaPixelMatchVars_s2[100], eghltEgammaEleGsfTrackIso[100], eghltEgammaGsfTrackVars_Chi2[100], eghltEgammaGsfTrackVars_Deta[100], eghltEgammaGsfTrackVars_DetaSeed[100], eghltEgammaGsfTrackVars_Dphi[100], eghltEgammaGsfTrackVars_MissingHits[100], eghltEgammaGsfTrackVars_NLayerIT[100], eghltEgammaGsfTrackVars_OneOESeedMinusOneOP[100], eghltEgammaGsfTrackVars_OneOESuperMinusOneOP[100], eghltEgammaGsfTrackVars_ValidHits[100];
  double egushltEgammaPixelMatchVars_s2[100], egushltEgammaEleGsfTrackIso[100], egushltEgammaGsfTrackVars_Chi2[100], egushltEgammaGsfTrackVars_Deta[100], egushltEgammaGsfTrackVars_DetaSeed[100], egushltEgammaGsfTrackVars_Dphi[100], egushltEgammaGsfTrackVars_MissingHits[100], egushltEgammaGsfTrackVars_NLayerIT[100], egushltEgammaGsfTrackVars_OneOESeedMinusOneOP[100], egushltEgammaGsfTrackVars_OneOESuperMinusOneOP[100], egushltEgammaGsfTrackVars_ValidHits[100];

  std::vector<TH1F*> all1dhists;
};

#endif
