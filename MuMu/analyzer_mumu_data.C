#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std;

void sort(double *pt, double *eta, double *phi, double *dxy, double *dxy_sig, int n) {

  // Sort according to pt
  for(unsigned int i=0; i<n; i++) {
    for(unsigned int j=0; j<n; j++) {
      if(*(pt+i)>*(pt+j)) {
	double temp = *(pt+i);
	*(pt+i) = *(pt+j);
	*(pt+j) = temp;
	temp = *(eta+i);
	*(eta+i) = *(eta+j);
	*(eta+j) = temp;
	temp = *(phi+i);
	*(phi+i) = *(phi+j);
	*(phi+j) = temp;
	temp = *(dxy+i);
	*(dxy+i) = *(dxy+j);
	*(dxy+j) = temp;
	temp = *(dxy_sig+i);
	*(dxy_sig+i) = *(dxy_sig+j);
	*(dxy_sig+j) = temp;
      }
    }
  }
}

// Functions required for computing sphericity and spherocity
TMatrixD makeMomentumTensor3D(std::vector<TLorentzVector*> lvarray) {
  TMatrixD total(3,3);
  double normaliser = 0.0;
  TVector3 vect;

  for(int ii=0; ii<3; ii++) {
    for(int jj=0; jj<3; jj++) {
      normaliser=0.0;
      total[ii][jj]=0.0;
      for(unsigned int vec=0; vec<lvarray.size(); vec++) {
	vect = lvarray[vec]->Vect();
	total[ii][jj] += vect[ii]*vect[jj];
	normaliser += vect.Mag()*vect.Mag();
      }
      if(normaliser>0) {
	total[ii][jj] /= normaliser;
      }
    }
  }
  return total;
}

// Functions required for computing sphericity and spherocity
double sphericity(std::vector<TLorentzVector*> lvarray) {
  TMatrixD total = makeMomentumTensor3D(lvarray);

  TVectorD eigenvals;
  auto eigenvectors = total.EigenVectors(eigenvals);

  double sphericity = 1.5*(eigenvals[1]+eigenvals[2]);
  return sphericity;
}

double transversespherocity(std::vector<TLorentzVector*> lvarray) {

  Double_t lowestval = 100000000.0;

  TVector3 workvec;
  double workerNumerator = 0.0;
  double workerDenominator = 0.0;
  TVector3 unitvec;

  for(unsigned int vecCtr=0; vecCtr<lvarray.size(); vecCtr++) {
    TLorentzVector* candvec = lvarray[vecCtr];

    unitvec.SetXYZ(candvec->Px()/candvec->Pt(), candvec->Py()/candvec->Pt(), 0.0);
    workvec.SetXYZ(0, 0, 0);
    workerNumerator = 0.0;
    workerDenominator = 0.0;

    for(unsigned int vecCtr2=0; vecCtr2<lvarray.size(); vecCtr2++) {
      TLorentzVector* workLVvec = lvarray[vecCtr2];

      workerDenominator += workLVvec->Pt();
      workvec.SetXYZ(workLVvec->Px(), workLVvec->Py(), 0.0);
      workerNumerator += workvec.Cross(unitvec).Mag();
    }
    lowestval = TMath::Min(lowestval, workerNumerator/workerDenominator);
  }

  double spherocity = TMath::Power(lowestval*TMath::Pi()*0.5, 2);
  return spherocity;
}

int analyzer_mumu_singlefile(TString inrootfile, TString outrootfile) {

  auto chain = new TChain("events");
  chain->Add(inrootfile);
  
  int muRecoN, muFiltN, htFiltN;
  double muFiltPt[100];
  double muFiltEta[100];
  double muFiltPhi[100];
  double muRecoPt[100];
  double muRecoEta[100];
  double muRecoPhi[100];
  double muRecoDxy[100];
  double muRecoDxySig[100];
  double htFiltPt[2];
  double htFiltEta[2];
  double htFiltPhi[2];
  
  chain->SetBranchAddress("muFiltn", &muFiltN);
  chain->SetBranchAddress("muFilt_pt", &muFiltPt);
  chain->SetBranchAddress("muFilt_eta", &muFiltEta);
  chain->SetBranchAddress("muFilt_phi", &muFiltPhi);
  chain->SetBranchAddress("mun", &muRecoN);
  chain->SetBranchAddress("mu_pt", &muRecoPt);
  chain->SetBranchAddress("mu_eta", &muRecoEta);
  chain->SetBranchAddress("mu_phi", &muRecoPhi);
  chain->SetBranchAddress("mu_dxy", &muRecoDxy);
  chain->SetBranchAddress("mu_dxy_sig", &muRecoDxySig);
  chain->SetBranchAddress("htFiltn", &htFiltN);
  chain->SetBranchAddress("htFilt_pt", &htFiltPt);
  chain->SetBranchAddress("htFilt_eta", &htFiltEta);
  chain->SetBranchAddress("htFilt_phi", &htFiltPhi);

  int totEntries = chain->GetEntries();
  std::cout<<"Events to process: "<<totEntries<<std::endl;

  auto outfile = new TFile(outrootfile,"RECREATE");

  // all objects
  auto muFiltmult = new TH1F("muFiltmult","N #mu",10,0,10);
  auto muFiltpt = new TH1F("muFiltpt","#mu p_{T} / GeV",500,0,500);
  auto muFilteta = new TH1F("muFilteta","#mu #eta",52,-2.6,2.6);
  auto muFiltphi = new TH1F("muFiltphi","#mu #phi",66,-3.3,3.3);
  auto mumult = new TH1F("mumult","N #mu",10,0,10);
  auto mupt = new TH1F("mupt","#mu p_{T} / GeV",500,0,500);
  auto mueta = new TH1F("mueta","#mu #eta",52,-2.6,2.6);
  auto muphi = new TH1F("muphi","#mu #phi",66,-3.3,3.3);
  auto mudxy = new TH1F("mudxy","#mu d_{0} / cm",20000,-20,20);
  auto mulog10dxy = new TH1F("mulog10dxy","#mu log_{10}d_{0}",500,-3,2);
  auto mudxysig = new TH1F("mudxysig","sig. #mu d_{0}",100,-1,99);
  auto ht = new TH1F("ht","HT /GeV",500,0,500);
  auto mht = new TH1F("mht","MHT /GeV",500,0,500);
  auto mhtphi = new TH1F("mhtphi","MHT #phi",66,-3.3,3.3);

  // selective objects
  auto mu2pt = new TH1F("mu2pt","sub-lead #mu p_{T} / GeV",500,0,500);
  auto mu2eta = new TH1F("mu2eta","sub-lead #mu #eta",52,-2.6,2.6);
  auto mu2phi = new TH1F("mu2phi","sub-lead #mu #phi",66,-3.3,3.3);
  auto mu2dxy = new TH1F("mu2dxy","sub-lead #mu d_{0} / cm",20000,-20,20);
  auto mu2log10dxy = new TH1F("mu2log10dxy","sub-lead #mu log_{10}d_{0} / cm",500,-3,2);
  auto mu2dxysig = new TH1F("mu2dxysig","sig. sub-lead #mu d_{0}",100,-1,99);
  auto mindxy = new TH1F("mindxy","min #mu d_{0} / cm",20000,-20,20);
  auto minlog10dxy = new TH1F("minlog10dxy","min #mu log_{10}d_{0} / cm",500,-3,2);
  auto mindxysig = new TH1F("mindxysig","sig. min #mu d_{0}",100,-1,99);

  // angles between objects
  auto mumudR = new TH1F("mumudR","#DeltaR(#mu,#mu)",140,0,7);
  auto mumudPhi = new TH1F("mumudPhi","#Delta#phi(#mu,#mu)",66,0,3.3);
  auto mumudEta = new TH1F("mumudEta","#Delta#eta(#mu,#mu)",60,0,6);
  auto mu2MhtdPhi = new TH1F("mu2MhtdPhi","#Delta#phi(#mu_{2},MHT)",66,0,3.3);
  auto mumuMhtdPhi = new TH1F("mumuMhtdPhi","#Delta#phi(#mu#mu,MHT)",66,0,3.3);
  auto h_sphericity = new TH1F("sphericity","sphericity",100,0,1);
  auto h_spherocity = new TH1F("spherocity","spherocity",100,0,1);
  auto mumudPhivdEta = new TH2F("mumudPhivdEta","#Delta#phi(#mu,#mu)v#Delta#eta(#mu,#mu)",66,0,3.3,60,0,6);

  // mumu variables
  auto mumuM = new TH1F("mumuM","M(#mu, #mu)",200, 0, 200);
  auto mumulog10M = new TH1F("mumulog10M","log10M(#mu, #mu)",300, -1, 2);
  auto muptmupt = new TH2F("muptmupt","p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);
  auto muetamueta = new TH2F("muetamueta","#eta(#mu_{1}) #sim #eta(#mu_{2})",52,-2.6,2.6,52,-2.6,2.6);
  auto muphimuphi = new TH2F("muphimuphi","#phi(#mu_{1}) #sim #phi(#mu_{2})",66,-3.3,3.3,66,-3.3,3.3);

  // Run2 displaced MuMu trigger selection
  auto orig_mumuM = new TH1F("orig_mumuM","orig M(#mu, #mu)",200, 0, 200);
  auto orig_mumult = new TH1F("orig_mumult","orig N #mu",10,0,10);
  auto orig_mu2pt = new TH1F("orig_mu2pt","orig sub-lead #mu p_{T} / GeV",500,0,500);
  auto orig_mu1pt = new TH1F("orig_mu1pt","orig lead #mu p_{T} / GeV",500,0,500);
  auto orig_mindxy = new TH1F("orig_mindxy","orig min #mu d_{0} / cm",20000,-20,20);
  auto orig_minlog10dxy = new TH1F("orig_minlog10dxy","orig min #mu log_{10}d_{0} / cm",500,-3,2);
  auto orig_mudxy = new TH1F("orig_mudxy","orig #mu d_{0} / cm",20000,-20,20);
  auto orig_mulog10dxy = new TH1F("orig_mulog10dxy","orig #mu log_{10}d_{0}",500,-3,2);
  auto orig_mu1log10dxysig = new TH1F("orig_mu1log10dxysig","orig sig. #mu_{1} d_{0} log",1000,-5,5);
  auto orig_mu2log10dxysig = new TH1F("orig_mu2log10dxysig","orig sig. #mu_{2} d_{0} log",1000,-5,5);
  auto orig_mu1dxysig = new TH1F("orig_mu1dxysig","orig sig. #mu_{1} d_{0}",1000,0,1000);
  auto orig_mu2dxysig = new TH1F("orig_mu2dxysig","orig sig. #mu_{2} d_{0}",1000,0,1000);
  auto orig_mumudR = new TH1F("orig_mumudR","orig #DeltaR(#mu,#mu)",100,0,5);
  auto orig_sphericity = new TH1F("orig_sphericity","orig Sphericity",100,0,1);
  auto orig_spherocity = new TH1F("orig_spherocity","orig Spherocity",100,0,1);
  auto orig_ht = new TH1F("orig_ht","orig HT /GeV",500,0,500);
  auto orig_mht = new TH1F("orig_mht","orig MHT /GeV",500,0,500);

  // All selected events
  auto sel_mumult = new TH1F("sel_mumult","sel N #mu",10,0,10);
  auto sel_mu2pt = new TH1F("sel_mu2pt","sel sub-lead #mu p_{T} / GeV",500,0,500);
  auto sel_mu1pt = new TH1F("sel_mu1pt","sel lead #mu p_{T} / GeV",500,0,500);
  auto sel_mindxy = new TH1F("sel_mindxy","sel min #mu d_{0} / cm",20000,-20,20);
  auto sel_minlog10dxy = new TH1F("sel_minlog10dxy","sel min #mu log_{10} d_{0} / log_{10} cm",500,-3,2);
  auto sel_mumuM = new TH1F("sel_mumuM","sel M(#mu, #mu)",200, 0, 200);
  auto sel_mudxy = new TH1F("sel_mudxy","sel #mu d_{0} / cm",20000,-20,20);
  auto sel_mulog10dxy = new TH1F("sel_mulog10dxy","sel #mu log_{10}d_{0}",500,-3,2);
  auto sel_mu1log10dxysig = new TH1F("sel_mu1log10dxysig","sel sig. #mu_{1} d_{0} log",1000,-5,5);
  auto sel_mu2log10dxysig = new TH1F("sel_mu2log10dxysig","sel sig. #mu_{2} d_{0} log",1000,-5,5);
  auto sel_mu1dxysig = new TH1F("sel_mu1dxysig","sel sig. #mu_{1} d_{0}",1000,0,1000);
  auto sel_mu2dxysig = new TH1F("sel_mu2dxysig","sel sig. #mu_{2} d_{0}",1000,0,1000);
  auto sel_mumudR = new TH1F("sel_mumudR","sel #DeltaR(#mu,#mu)",100,0,5);
  auto sel_mumudEta = new TH1F("sel_mumudEta","sel #Delta#eta(#mu,#mu)",100,0,5);
  auto sel_mumudPhi = new TH1F("sel_mumudPhi","sel #Delta#phi(#mu,#mu)",100,0,3.15);
  auto sel_sphericity = new TH1F("sel_sphericity","sel Sphericity",100,0,1);
  auto sel_spherocity = new TH1F("sel_spherocity","sel Spherocity",100,0,1);
  auto sel_ht = new TH1F("sel_ht","sel HT /GeV",500,0,500);
  auto sel_mht = new TH1F("sel_mht","sel MHT /GeV",500,0,500);
  auto sel_muptmupt = new TH2F("sel_muptmupt","sel p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);
  auto sel_muetamueta = new TH2F("sel_muetamueta","sel #eta(#mu_{1}) #sim #eta(#mu_{2})",52,-2.6,2.6,52,-2.6,2.6);
  auto sel_muphimuphi = new TH2F("sel_muphimuphi","sel #phi(#mu_{1}) #sim #phi(#mu_{2})",66,-3.3,3.3,66,-3.3,3.3);
  auto sel_mu1dxysigabsdxy = new TH2F("sel_mu1dxysigabsdxy","sig. #mu_{1} d_{0} #sim abs(#mu_{1} d_{0})",1000,0,100,200,0,2);
  auto sel_mu2dxysigabsdxy = new TH2F("sel_mu2dxysigabsdxy","sig. #mu_{2} d_{0} #sim abs(#mu_{2} d_{0})",1000,0,100,200,0,2);

  // Z peak variables 84 - 98
  auto Z_mumuM = new TH1F("Z_mumuM","Z M(#mu, #mu)",200, 0, 200);
  auto Z_mu2pt = new TH1F("Z_mu2pt","Z sub-lead #mu p_{T} / GeV",500,0,500);
  auto Z_mindxy = new TH1F("Z_mindxy","Z min #mu d_{0} / cm",20000,-20,20);
  auto Z_minlog10dxy = new TH1F("Z_minlog10dxy","Z min #mu log_{10} d_{0} / log_{10} cm",500,-3,2);
  auto Z_mudxy = new TH1F("Z_mudxy","Z #mu d_{0} / cm",20000,-20,20);
  auto Z_mulog10dxy = new TH1F("Z_mulog10dxy","Z #mu log_{10}d_{0}",500,-3,2);
  auto Z_mu1log10dxysig = new TH1F("Z_mu1log10dxysig","Z sig. #mu_{1} d_{0} log",1000,-5,5);
  auto Z_mu2log10dxysig = new TH1F("Z_mu2log10dxysig","Z sig. #mu_{2} d_{0} log",1000,-5,5);
  auto Z_mu1dxysig = new TH1F("Z_mu1dxysig","Z sig. #mu_{1} d_{0}",1000,0,1000);
  auto Z_mu2dxysig = new TH1F("Z_mu2dxysig","Z sig. #mu_{2} d_{0}",1000,0,1000);
  auto Z_mumudR = new TH1F("Z_mumudR","Z #DeltaR(#mu,#mu)",100,0,5);
  auto Z_muptmupt = new TH2F("Z_muptmupt","Z p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);
  auto Z_muetamueta = new TH2F("Z_muetamueta","Z #eta(#mu_{1}) #sim #eta(#mu_{2})",52,-2.6,2.6,52,-2.6,2.6);
  auto Z_muphimuphi = new TH2F("Z_muphimuphi","Z #phi(#mu_{1}) #sim #phi(#mu_{2})",66,-3.3,3.3,66,-3.3,3.3);
  auto Z_Mdxy = new TH2F("Z_Mdxy","Z M(#mu, #mu) #mu d_{0} / cm",200,-20,20,20,80,100);
  auto Z_dxysigdxy = new TH2F("Z_dxysigdxy","Z dxysigdxy / cm",200,-20,20,100,0,100);

  // J/psi peak variables 2 - 4
  auto Jpsi_mumuM = new TH1F("Jpsi_mumuM","Jpsi M(#mu, #mu)",200, 0, 200);
  auto Jpsi_mu2pt = new TH1F("Jpsi_mu2pt","Jpsi sub-lead #mu p_{T} / GeV",500,0,500);
  auto Jpsi_mindxy = new TH1F("Jpsi_mindxy","Jpsi min #mu d_{0} / cm",20000,-20,20);
  auto Jpsi_minlog10dxy = new TH1F("Jpsi_minlog10dxy","Jpsi min #mu log_{10} d_{0} / log_{10} cm",500,-3,2);
  auto Jpsi_mudxy = new TH1F("Jpsi_mudxy","Jpsi #mu d_{0} / cm",20000,-20,20);
  auto Jpsi_mulog10dxy = new TH1F("Jpsi_mulog10dxy","Jpsi #mu log_{10}d_{0}",500,-3,2);
  auto Jpsi_mu1log10dxysig = new TH1F("Jpsi_mu1log10dxysig","Jpsi sig. #mu_{1} d_{0} log",1000,-5,5);
  auto Jpsi_mu2log10dxysig = new TH1F("Jpsi_mu2log10dxysig","Jpsi sig. #mu_{2} d_{0} log",1000,-5,5);
  auto Jpsi_mu1dxysig = new TH1F("Jpsi_mu1dxysig","Jpsi sig. #mu_{1} d_{0}",1000,0,1000);
  auto Jpsi_mu2dxysig = new TH1F("Jpsi_mu2dxysig","Jpsi sig. #mu_{2} d_{0}",1000,0,1000);
  auto Jpsi_mumudR = new TH1F("Jpsi_mumudR","Jpsi #DeltaR(#mu,#mu)",100,0,5);
  auto Jpsi_muptmupt = new TH2F("Jpsi_muptmupt","Jpsi p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);
  auto Jpsi_muetamueta = new TH2F("Jpsi_muetamueta","Jpsi #eta(#mu_{1}) #sim #eta(#mu_{2})",52,-2.6,2.6,52,-2.6,2.6);
  auto Jpsi_muphimuphi = new TH2F("Jpsi_muphimuphi","Jpsi #phi(#mu_{1}) #sim #phi(#mu_{2})",66,-3.3,3.3,66,-3.3,3.3);

  // M plateau 35 - 84
  auto Mplt_mumuM = new TH1F("Mplt_mumuM","Mplt M(#mu, #mu)",200, 0, 200);
  auto Mplt_mu2pt = new TH1F("Mplt_mu2pt","Mplt sub-lead #mu p_{T} / GeV",500,0,500);
  auto Mplt_mindxy = new TH1F("Mplt_mindxy","Mplt min #mu d_{0} / cm",20000,-20,20);
  auto Mplt_minlog10dxy = new TH1F("Mplt_minlog10dxy","Mplt min #mu log_{10} d_{0} / log_{10} cm",500,-3,2);
  auto Mplt_mudxy = new TH1F("Mplt_mudxy","Mplt #mu d_{0} / cm",20000,-20,20);
  auto Mplt_mulog10dxy = new TH1F("Mplt_mulog10dxy","Mplt #mu log_{10}d_{0}",500,-3,2);
  auto Mplt_mu1log10dxysig = new TH1F("Mplt_mu1log10dxysig","Mplt sig. #mu_{1} d_{0} log",1000,-5,5);
  auto Mplt_mu2log10dxysig = new TH1F("Mplt_mu2log10dxysig","Mplt sig. #mu_{2} d_{0} log",1000,-5,5);
  auto Mplt_mu1dxysig = new TH1F("Mplt_mu1dxysig","Mplt sig. #mu_{1} d_{0}",1000,0,1000);
  auto Mplt_mu2dxysig = new TH1F("Mplt_mu2dxysig","Mplt sig. #mu_{2} d_{0}",1000,0,1000);
  auto Mplt_mumudR = new TH1F("Mplt_mumudR","Mplt #DeltaR(#mu,#mu)",100,0,5);
  auto Mplt_muptmupt = new TH2F("Mplt_muptmupt","Mplt p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);
  auto Mplt_muetamueta = new TH2F("Mplt_muetamueta","Mplt #eta(#mu_{1}) #sim #eta(#mu_{2})",52,-2.6,2.6,52,-2.6,2.6);
  auto Mplt_muphimuphi = new TH2F("Mplt_muphimuphi","Mplt #phi(#mu_{1}) #sim #phi(#mu_{2})",66,-3.3,3.3,66,-3.3,3.3);
    
  // others
  auto oth_mumuM = new TH1F("oth_mumuM","oth M(#mu, #mu)",200, 0, 200);
  auto oth_mu2pt = new TH1F("oth_mu2pt","oth sub-lead #mu p_{T} / GeV",500,0,500);
  auto oth_mindxy = new TH1F("oth_mindxy","oth min #mu d_{0} / cm",20000,-20,20);
  auto oth_minlog10dxy = new TH1F("oth_minlog10dxy","oth min #mu log_{10} d_{0} / log_{10} cm",500,-3,2);
  auto oth_mudxy = new TH1F("oth_mudxy","oth #mu d_{0} / cm",20000,-20,20);
  auto oth_mulog10dxy = new TH1F("oth_mulog10dxy","oth #mu log_{10}d_{0}",500,-3,2);
  auto oth_mu1log10dxysig = new TH1F("oth_mu1log10dxysig","oth sig. #mu_{1} d_{0} log",1000,-5,5);
  auto oth_mu2log10dxysig = new TH1F("oth_mu2log10dxysig","oth sig. #mu_{2} d_{0} log",1000,-5,5);
  auto oth_mu1dxysig = new TH1F("oth_mu1dxysig","oth sig. #mu_{1} d_{0}",1000,0,1000);
  auto oth_mu2dxysig = new TH1F("oth_mu2dxysig","oth sig. #mu_{2} d_{0}",1000,0,1000);
  auto oth_mumudR = new TH1F("oth_mumudR","oth #DeltaR(#mu,#mu)",100,0,5);
  auto oth_muptmupt = new TH2F("oth_muptmupt","oth p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);
  auto oth_muetamueta = new TH2F("oth_muetamueta","oth #eta(#mu_{1}) #sim #eta(#mu_{2})",52,-2.6,2.6,52,-2.6,2.6);
  auto oth_muphimuphi = new TH2F("oth_muphimuphi","oth #phi(#mu_{1}) #sim #phi(#mu_{2})",66,-3.3,3.3,66,-3.3,3.3);

  // dR>1
  auto cut1_mumuM = new TH1F("cut1_mumuM","cut1 M(#mu, #mu)",200, 0, 200);
  auto cut1_mumult = new TH1F("cut1_mumult","cut1 N #mu",10,0,10);
  auto cut1_mu2pt = new TH1F("cut1_mu2pt","cut1 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut1_mu1pt = new TH1F("cut1_mu1pt","cut1 lead #mu p_{T} / GeV",500,0,500);
  auto cut1_mindxy = new TH1F("cut1_mindxy","cut1 min #mu d_{0} / cm",20000,-20,20);
  auto cut1_minlog10dxy = new TH1F("cut1_minlog10dxy","cut1 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut1_mudxy = new TH1F("cut1_mudxy","cut1 #mu d_{0} / cm",20000,-20,20);
  auto cut1_mulog10dxy = new TH1F("cut1_mulog10dxy","cut1 #mu log_{10}d_{0}",500,-3,2);
  auto cut1_mu1log10dxysig = new TH1F("cut1_mu1log10dxysig","cut1 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut1_mu2log10dxysig = new TH1F("cut1_mu2log10dxysig","cut1 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut1_mu1dxysig = new TH1F("cut1_mu1dxysig","cut1 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut1_mu2dxysig = new TH1F("cut1_mu2dxysig","cut1 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut1_mumudR = new TH1F("cut1_mumudR","cut1 #DeltaR(#mu,#mu)",100,0,5);
  auto cut1_sphericity = new TH1F("cut1_sphericity","cut1 Sphericity",100,0,1);
  auto cut1_spherocity = new TH1F("cut1_spherocity","cut1 Spherocity",100,0,1);
  auto cut1_ht = new TH1F("cut1_ht","cut1 HT /GeV",500,0,500);
  auto cut1_mht = new TH1F("cut1_mht","cut1 MHT /GeV",500,0,500);
  auto cut1_muptmupt = new TH2F("cut1_muptmupt","cut1 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // 84<M>98
  auto cut2_mumuM = new TH1F("cut2_mumuM","cut2 M(#mu, #mu)",200, 0, 200);
  auto cut2_mumult = new TH1F("cut2_mumult","cut2 N #mu",10,0,10);
  auto cut2_mu2pt = new TH1F("cut2_mu2pt","cut2 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut2_mu1pt = new TH1F("cut2_mu1pt","cut2 lead #mu p_{T} / GeV",500,0,500);
  auto cut2_mindxy = new TH1F("cut2_mindxy","cut2 min #mu d_{0} / cm",20000,-20,20);
  auto cut2_minlog10dxy = new TH1F("cut2_minlog10dxy","cut2 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut2_mudxy = new TH1F("cut2_mudxy","cut2 #mu d_{0} / cm",20000,-20,20);
  auto cut2_mulog10dxy = new TH1F("cut2_mulog10dxy","cut2 #mu log_{10}d_{0}",500,-3,2);
  auto cut2_mu1log10dxysig = new TH1F("cut2_mu1log10dxysig","cut2 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut2_mu2log10dxysig = new TH1F("cut2_mu2log10dxysig","cut2 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut2_mu1dxysig = new TH1F("cut2_mu1dxysig","cut2 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut2_mu2dxysig = new TH1F("cut2_mu2dxysig","cut2 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut2_mumudR = new TH1F("cut2_mumudR","cut2 #DeltaR(#mu,#mu)",100,0,5);
  auto cut2_sphericity = new TH1F("cut2_sphericity","cut2 Sphericity",100,0,1);
  auto cut2_spherocity = new TH1F("cut2_spherocity","cut2 Spherocity",100,0,1);
  auto cut2_ht = new TH1F("cut2_ht","cut2 HT /GeV",500,0,500);
  auto cut2_mht = new TH1F("cut2_mht","cut2 MHT /GeV",500,0,500);
  auto cut2_muptmupt = new TH2F("cut2_muptmupt","cut2 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // 84<M>98 and dR>1
  auto cut3_mumuM = new TH1F("cut3_mumuM","cut3 M(#mu, #mu)",200, 0, 200);
  auto cut3_mumult = new TH1F("cut3_mumult","cut3 N #mu",10,0,10);
  auto cut3_mu2pt = new TH1F("cut3_mu2pt","cut3 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut3_mu1pt = new TH1F("cut3_mu1pt","cut3 lead #mu p_{T} / GeV",500,0,500);
  auto cut3_mindxy = new TH1F("cut3_mindxy","cut3 min #mu d_{0} / cm",20000,-20,20);
  auto cut3_minlog10dxy = new TH1F("cut3_minlog10dxy","cut3 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut3_mudxy = new TH1F("cut3_mudxy","cut3 #mu d_{0} / cm",20000,-20,20);
  auto cut3_mulog10dxy = new TH1F("cut3_mulog10dxy","cut3 #mu log_{10}d_{0}",500,-3,2);
  auto cut3_minlog10dxysig = new TH1F("cut3_minlog10dxysig","cut3 min sig. #mu d_{0} log",1000,-5,5);
  auto cut3_maxlog10dxysig = new TH1F("cut3_maxlog10dxysig","cut3 max sig. #mu d_{0} log",1000,-5,5);
  auto cut3_mu1log10dxysig = new TH1F("cut3_mu1log10dxysig","cut3 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut3_mu2log10dxysig = new TH1F("cut3_mu2log10dxysig","cut3 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut3_mu1dxysig = new TH1F("cut3_mu1dxysig","cut3 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut3_mu2dxysig = new TH1F("cut3_mu2dxysig","cut3 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut3_mumudR = new TH1F("cut3_mumudR","cut3 #DeltaR(#mu,#mu)",100,0,5);
  auto cut3_sphericity = new TH1F("cut3_sphericity","cut3 Sphericity",100,0,1);
  auto cut3_spherocity = new TH1F("cut3_spherocity","cut3 Spherocity",100,0,1);
  auto cut3_ht = new TH1F("cut3_ht","cut3 HT /GeV",500,0,500);
  auto cut3_mht = new TH1F("cut3_mht","cut3 MHT /GeV",500,0,500);
  auto cut3_muptmupt = new TH2F("cut3_muptmupt","cut3 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // d0>0.01
  auto cut4_mumuM = new TH1F("cut4_mumuM","cut4 M(#mu, #mu)",200, 0, 200);
  auto cut4_mumult = new TH1F("cut4_mumult","cut4 N #mu",10,0,10);
  auto cut4_mu2pt = new TH1F("cut4_mu2pt","cut4 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut4_mu1pt = new TH1F("cut4_mu1pt","cut4 lead #mu p_{T} / GeV",500,0,500);
  auto cut4_mindxy = new TH1F("cut4_mindxy","cut4 min #mu d_{0} / cm",20000,-20,20);
  auto cut4_minlog10dxy = new TH1F("cut4_minlog10dxy","cut4 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut4_mudxy = new TH1F("cut4_mudxy","cut4 #mu d_{0} / cm",20000,-20,20);
  auto cut4_mulog10dxy = new TH1F("cut4_mulog10dxy","cut4 #mu log_{10}d_{0}",500,-3,2);
  auto cut4_mu1log10dxysig = new TH1F("cut4_mu1log10dxysig","cut4 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut4_mu2log10dxysig = new TH1F("cut4_mu2log10dxysig","cut4 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut4_mu1dxysig = new TH1F("cut4_mu1dxysig","cut4 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut4_mu2dxysig = new TH1F("cut4_mu2dxysig","cut4 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut4_mumudR = new TH1F("cut4_mumudR","cut4 #DeltaR(#mu,#mu)",100,0,5);
  auto cut4_sphericity = new TH1F("cut4_sphericity","cut4 Sphericity",100,0,1);
  auto cut4_spherocity = new TH1F("cut4_spherocity","cut4 Spherocity",100,0,1);
  auto cut4_ht = new TH1F("cut4_ht","cut4 HT /GeV",500,0,500);
  auto cut4_mht = new TH1F("cut4_mht","cut4 MHT /GeV",500,0,500);
  auto cut4_muptmupt = new TH2F("cut4_muptmupt","cut4 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // dxysigofmuon>0.1, dR>1, 84<M>98
  auto cut5_mumuM = new TH1F("cut5_mumuM","cut5 M(#mu, #mu)",200, 0, 200);
  auto cut5_mumult = new TH1F("cut5_mumult","cut5 N #mu",10,0,10);
  auto cut5_mu1pt = new TH1F("cut5_mu1pt","cut5 lead #mu p_{T} / GeV",500,0,500);
  auto cut5_mu2pt = new TH1F("cut5_mu2pt","cut5 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut5_mindxy = new TH1F("cut5_mindxy","cut5 min #mu d_{0} / cm",20000,-20,20);
  auto cut5_minlog10dxy = new TH1F("cut5_minlog10dxy","cut5 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut5_mudxy = new TH1F("cut5_mudxy","cut5 #mu d_{0} / cm",20000,-20,20);
  auto cut5_mulog10dxy = new TH1F("cut5_mulog10dxy","cut5 #mu log_{10}d_{0}",500,-3,2);
  auto cut5_mu1log10dxysig = new TH1F("cut5_mu1log10dxysig","cut5 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut5_mu2log10dxysig = new TH1F("cut5_mu2log10dxysig","cut5 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut5_minlog10dxysig = new TH1F("cut5_minlog10dxysig","cut5 min sig. #mu d_{0} log",1000,-5,5);
  auto cut5_maxlog10dxysig = new TH1F("cut5_maxlog10dxysig","cut5 max sig. #mu d_{0} log",1000,-5,5);
  auto cut5_mu1dxysig = new TH1F("cut5_mu1dxysig","cut5 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut5_mu2dxysig = new TH1F("cut5_mu2dxysig","cut5 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut5_mumudR = new TH1F("cut5_mumudR","cut5 #DeltaR(#mu,#mu)",100,0,5);
  auto cut5_sphericity = new TH1F("cut5_sphericity","cut5 Sphericity",100,0,1);
  auto cut5_spherocity = new TH1F("cut5_spherocity","cut5 Spherocity",100,0,1);
  auto cut5_ht = new TH1F("cut5_ht","cut5 HT /GeV",500,0,500);
  auto cut5_mht = new TH1F("cut5_mht","cut5 MHT /GeV",500,0,500);
  auto cut5_muptmupt = new TH2F("cut5_muptmupt","cut5 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // dxysigofmuon>1, dR>1, 84<M>98
  auto cut6_mumuM = new TH1F("cut6_mumuM","cut6 M(#mu, #mu)",200, 0, 200);
  auto cut6_mumult = new TH1F("cut6_mumult","cut6 N #mu",10,0,10);
  auto cut6_mu2pt = new TH1F("cut6_mu2pt","cut6 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut6_mu1pt = new TH1F("cut6_mu1pt","cut6 lead #mu p_{T} / GeV",500,0,500);
  auto cut6_mindxy = new TH1F("cut6_mindxy","cut6 min #mu d_{0} / cm",20000,-20,20);
  auto cut6_minlog10dxy = new TH1F("cut6_minlog10dxy","cut6 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut6_mudxy = new TH1F("cut6_mudxy","cut6 #mu d_{0} / cm",20000,-20,20);
  auto cut6_mulog10dxy = new TH1F("cut6_mulog10dxy","cut6 #mu log_{10}d_{0}",500,-3,2);
  auto cut6_mu1log10dxysig = new TH1F("cut6_mu1log10dxysig","cut6 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut6_mu2log10dxysig = new TH1F("cut6_mu2log10dxysig","cut6 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut6_minlog10dxysig = new TH1F("cut6_minlog10dxysig","cut6 min sig. #mu d_{0} log",1000,-5,5);
  auto cut6_maxlog10dxysig = new TH1F("cut6_maxlog10dxysig","cut6 max sig. #mu d_{0} log",1000,-5,5);
  auto cut6_mu1dxysig = new TH1F("cut6_mu1dxysig","cut6 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut6_mu2dxysig = new TH1F("cut6_mu2dxysig","cut6 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut6_mumudR = new TH1F("cut6_mumudR","cut6 #DeltaR(#mu,#mu)",100,0,5);
  auto cut6_sphericity = new TH1F("cut6_sphericity","cut6 Sphericity",100,0,1);
  auto cut6_spherocity = new TH1F("cut6_spherocity","cut6 Spherocity",100,0,1);
  auto cut6_ht = new TH1F("cut6_ht","cut6 HT /GeV",500,0,500);
  auto cut6_mht = new TH1F("cut6_mht","cut6 MHT /GeV",500,0,500);
  auto cut6_muptmupt = new TH2F("cut6_muptmupt","cut6 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // dxysigofmuon>10, dR>1, 84<M>98
  auto cut7_mumuM = new TH1F("cut7_mumuM","cut7 M(#mu, #mu)",200, 0, 200);
  auto cut7_mumult = new TH1F("cut7_mumult","cut7 N #mu",10,0,10);
  auto cut7_mu2pt = new TH1F("cut7_mu2pt","cut7 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut7_mu1pt = new TH1F("cut7_mu1pt","cut7 lead #mu p_{T} / GeV",500,0,500);
  auto cut7_mindxy = new TH1F("cut7_mindxy","cut7 min #mu d_{0} / cm",20000,-20,20);
  auto cut7_minlog10dxy = new TH1F("cut7_minlog10dxy","cut7 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut7_mudxy = new TH1F("cut7_mudxy","cut7 #mu d_{0} / cm",20000,-20,20);
  auto cut7_mulog10dxy = new TH1F("cut7_mulog10dxy","cut7 #mu log_{10}d_{0}",500,-3,2);
  auto cut7_mu1log10dxysig = new TH1F("cut7_mu1log10dxysig","cut7 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut7_mu2log10dxysig = new TH1F("cut7_mu2log10dxysig","cut7 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut7_minlog10dxysig = new TH1F("cut7_minlog10dxysig","cut7 min sig. #mu d_{0} log",1000,-5,5);
  auto cut7_maxlog10dxysig = new TH1F("cut7_maxlog10dxysig","cut7 max sig. #mu d_{0} log",1000,-5,5);
  auto cut7_mu1dxysig = new TH1F("cut7_mu1dxysig","cut7 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut7_mu2dxysig = new TH1F("cut7_mu2dxysig","cut7 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut7_mumudR = new TH1F("cut7_mumudR","cut7 #DeltaR(#mu,#mu)",100,0,5);
  auto cut7_sphericity = new TH1F("cut7_sphericity","cut7 Sphericity",100,0,1);
  auto cut7_spherocity = new TH1F("cut7_spherocity","cut7 Spherocity",100,0,1);
  auto cut7_ht = new TH1F("cut7_ht","cut7 HT /GeV",500,0,500);
  auto cut7_mht = new TH1F("cut7_mht","cut7 MHT /GeV",500,0,500);
  auto cut7_muptmupt = new TH2F("cut7_muptmupt","cut7 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // dxy>0.01, dR>1, 84<M>98
  auto cut8_mumuM = new TH1F("cut8_mumuM","cut8 M(#mu, #mu)",200, 0, 200);
  auto cut8_mumult = new TH1F("cut8_mumult","cut8 N #mu",10,0,10);
  auto cut8_mu2pt = new TH1F("cut8_mu2pt","cut8 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut8_mu1pt = new TH1F("cut8_mu1pt","cut8 lead #mu p_{T} / GeV",500,0,500);
  auto cut8_mindxy = new TH1F("cut8_mindxy","cut8 min #mu d_{0} / cm",20000,-20,20);
  auto cut8_minlog10dxy = new TH1F("cut8_minlog10dxy","cut8 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut8_mudxy = new TH1F("cut8_mudxy","cut8 #mu d_{0} / cm",20000,-20,20);
  auto cut8_mulog10dxy = new TH1F("cut8_mulog10dxy","cut8 #mu log_{10}d_{0}",500,-3,2);
  auto cut8_mu1log10dxysig = new TH1F("cut8_mu1log10dxysig","cut8 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut8_mu2log10dxysig = new TH1F("cut8_mu2log10dxysig","cut8 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut8_mu1dxysig = new TH1F("cut8_mu1dxysig","cut8 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut8_mu2dxysig = new TH1F("cut8_mu2dxysig","cut8 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut8_mumudR = new TH1F("cut8_mumudR","cut8 #DeltaR(#mu,#mu)",100,0,5);
  auto cut8_sphericity = new TH1F("cut8_sphericity","cut8 Sphericity",100,0,1);
  auto cut8_spherocity = new TH1F("cut8_spherocity","cut8 Spherocity",100,0,1);
  auto cut8_ht = new TH1F("cut8_ht","cut8 HT /GeV",500,0,500);
  auto cut8_mht = new TH1F("cut8_mht","cut8 MHT /GeV",500,0,500);
  auto cut8_muptmupt = new TH2F("cut8_muptmupt","cut8 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // dxysigofmuon>1, dxy>0.01, dR>1, 84<M>98
  auto cut9_mumuM = new TH1F("cut9_mumuM","cut9 M(#mu, #mu)",200, 0, 200);
  auto cut9_mumult = new TH1F("cut9_mumult","cut9 N #mu",10,0,10);
  auto cut9_mu2pt = new TH1F("cut9_mu2pt","cut9 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut9_mu1pt = new TH1F("cut9_mu1pt","cut9 lead #mu p_{T} / GeV",500,0,500);
  auto cut9_mindxy = new TH1F("cut9_mindxy","cut9 min #mu d_{0} / cm",20000,-20,20);
  auto cut9_maxdxy = new TH1F("cut9_maxdxy","cut9 max #mu d_{0} / cm",20000,-20,20);
  auto cut9_minlog10dxy = new TH1F("cut9_minlog10dxy","cut9 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut9_maxlog10dxy = new TH1F("cut9_maxlog10dxy","cut9 max #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut9_mu1log10dxysig = new TH1F("cut9_mu1log10dxysig","cut9 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut9_mu2log10dxysig = new TH1F("cut9_mu2log10dxysig","cut9 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut9_mu1dxysig = new TH1F("cut9_mu1dxysig","cut9 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut9_mu2dxysig = new TH1F("cut9_mu2dxysig","cut9 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut9_minlog10dxysig = new TH1F("cut9_minlog10dxysig","cut9 min sig. #mu d_{0} log",1000,-5,5);
  auto cut9_maxlog10dxysig = new TH1F("cut9_maxlog10dxysig","cut9 max sig. #mu d_{0} log",1000,-5,5);
  auto cut9_mumudR = new TH1F("cut9_mumudR","cut9 #DeltaR(#mu,#mu)",100,0,5);
  auto cut9_ht = new TH1F("cut9_ht","cut9 HT /GeV",500,0,500);
  auto cut9_mht = new TH1F("cut9_mht","cut9 MHT /GeV",500,0,500);
  auto cut9_muptmupt = new TH2F("cut9_muptmupt","cut9 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // mu2pt>20
  auto cut10_mumuM = new TH1F("cut10_mumuM","cut10 M(#mu, #mu)",200, 0, 200);
  auto cut10_mumult = new TH1F("cut10_mumult","cut10 N #mu",10,0,10);
  auto cut10_mu2pt = new TH1F("cut10_mu2pt","cut10 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut10_mu1pt = new TH1F("cut10_mu1pt","cut10 lead #mu p_{T} / GeV",500,0,500);
  auto cut10_mindxy = new TH1F("cut10_mindxy","cut10 min #mu d_{0} / cm",20000,-20,20);
  auto cut10_maxdxy = new TH1F("cut10_maxdxy","cut10 max #mu d_{0} / cm",20000,-20,20);
  auto cut10_minlog10dxy = new TH1F("cut10_minlog10dxy","cut10 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut10_maxlog10dxy = new TH1F("cut10_maxlog10dxy","cut10 max #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut10_mu1log10dxysig = new TH1F("cut10_mu1log10dxysig","cut10 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut10_mu2log10dxysig = new TH1F("cut10_mu2log10dxysig","cut10 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut10_mu1dxysig = new TH1F("cut10_mu1dxysig","cut10 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut10_mu2dxysig = new TH1F("cut10_mu2dxysig","cut10 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut10_minlog10dxysig = new TH1F("cut10_minlog10dxysig","cut10 min sig. #mu d_{0} log",1000,-5,5);
  auto cut10_maxlog10dxysig = new TH1F("cut10_maxlog10dxysig","cut10 max sig. #mu d_{0} log",1000,-5,5);
  auto cut10_mumudR = new TH1F("cut10_mumudR","cut10 #DeltaR(#mu,#mu)",100,0,5);
  auto cut10_ht = new TH1F("cut10_ht","cut10 HT /GeV",500,0,500);
  auto cut10_mht = new TH1F("cut10_mht","cut10 MHT /GeV",500,0,500);
  auto cut10_muptmupt = new TH2F("cut10_muptmupt","cut10 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // dxysigofmuon>1, mu2pt>20, dxy>0.01, dR>1, 84<M>98
  auto cut11_mumuM = new TH1F("cut11_mumuM","cut11 M(#mu, #mu)",200, 0, 200);
  auto cut11_mumult = new TH1F("cut11_mumult","cut11 N #mu",10,0,10);
  auto cut11_mu2pt = new TH1F("cut11_mu2pt","cut11 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut11_mu1pt = new TH1F("cut11_mu1pt","cut11 lead #mu p_{T} / GeV",500,0,500);
  auto cut11_mindxy = new TH1F("cut11_mindxy","cut11 min #mu d_{0} / cm",20000,-20,20);
  auto cut11_maxdxy = new TH1F("cut11_maxdxy","cut11 max #mu d_{0} / cm",20000,-20,20);
  auto cut11_minlog10dxy = new TH1F("cut11_minlog10dxy","cut11 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut11_maxlog10dxy = new TH1F("cut11_maxlog10dxy","cut11 max #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut11_mu1log10dxysig = new TH1F("cut11_mu1log10dxysig","cut11 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut11_mu2log10dxysig = new TH1F("cut11_mu2log10dxysig","cut11 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut11_mu1dxysig = new TH1F("cut11_mu1dxysig","cut11 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut11_mu2dxysig = new TH1F("cut11_mu2dxysig","cut11 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut11_minlog10dxysig = new TH1F("cut11_minlog10dxysig","cut11 min sig. #mu d_{0} log",1000,-5,5);
  auto cut11_maxlog10dxysig = new TH1F("cut11_maxlog10dxysig","cut11 max sig. #mu d_{0} log",1000,-5,5);
  auto cut11_mumudR = new TH1F("cut11_mumudR","cut11 #DeltaR(#mu,#mu)",100,0,5);
  auto cut11_ht = new TH1F("cut11_ht","cut11 HT /GeV",500,0,500);
  auto cut11_mht = new TH1F("cut11_mht","cut11 MHT /GeV",500,0,500);
  auto cut11_muptmupt = new TH2F("cut11_muptmupt","cut11 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // dxymuon>0.01, mu2pt>20, dxy>0.01, dR>1, 84<M>98
  auto cut12_mumuM = new TH1F("cut12_mumuM","cut12 M(#mu, #mu)",200, 0, 200);
  auto cut12_mumult = new TH1F("cut12_mumult","cut12 N #mu",10,0,10);
  auto cut12_mu2pt = new TH1F("cut12_mu2pt","cut12 sub-lead #mu p_{T} / GeV",500,0,500);
  auto cut12_mu1pt = new TH1F("cut12_mu1pt","cut12 lead #mu p_{T} / GeV",500,0,500);
  auto cut12_mindxy = new TH1F("cut12_mindxy","cut12 min #mu d_{0} / cm",20000,-20,20);
  auto cut12_maxdxy = new TH1F("cut12_maxdxy","cut12 max #mu d_{0} / cm",20000,-20,20);
  auto cut12_minlog10dxy = new TH1F("cut12_minlog10dxy","cut12 min #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut12_maxlog10dxy = new TH1F("cut12_maxlog10dxy","cut12 max #mu log_{10}d_{0} / cm",500,-3,2);
  auto cut12_mu1log10dxysig = new TH1F("cut12_mu1log10dxysig","cut12 sig. #mu_{1} d_{0} log",1000,-5,5);
  auto cut12_mu2log10dxysig = new TH1F("cut12_mu2log10dxysig","cut12 sig. #mu_{2} d_{0} log",1000,-5,5);
  auto cut12_mu1dxysig = new TH1F("cut12_mu1dxysig","cut12 sig. #mu_{1} d_{0}",1000,0,1000);
  auto cut12_mu2dxysig = new TH1F("cut12_mu2dxysig","cut12 sig. #mu_{2} d_{0}",1000,0,1000);
  auto cut12_minlog10dxysig = new TH1F("cut12_minlog10dxysig","cut12 min sig. #mu d_{0} log",1000,-5,5);
  auto cut12_maxlog10dxysig = new TH1F("cut12_maxlog10dxysig","cut12 max sig. #mu d_{0} log",1000,-5,5);
  auto cut12_mumudR = new TH1F("cut12_mumudR","cut12 #DeltaR(#mu,#mu)",100,0,5);
  auto cut12_ht = new TH1F("cut12_ht","cut12 HT /GeV",500,0,500);
  auto cut12_mht = new TH1F("cut12_mht","cut12 MHT /GeV",500,0,500);
  auto cut12_muptmupt = new TH2F("cut12_muptmupt","cut12 p_{T}(#mu_{1}) #sim p_{T}(#mu_{2})",200,0,200,200,0,200);

  // cut flow tests
  int cutflow1 = 0;
  int cutflow2 = 0;
  int cutflow3 = 0;
  int cutflow4 = 0;
  int cutflowdatarate = 0;
  int cutflowcut1 = 0;
  int cutflowcut2 = 0;
  int cutflowcut3 = 0;
  int cutflowcut4 = 0;
  int cutflowcut5 = 0;
  int cutflowcut6 = 0;
  int cutflowcut7 = 0;
  int cutflowcut8 = 0;
  int cutflowcut9 = 0;
  int cutflowcut10 = 0;
  int cutflowcut11 = 0;
  int cutflowcut12 = 0;

  double datarate = 0.5;
  
  for(unsigned int event=0; event<totEntries; event++) {

    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;
    chain->GetEntry(event);

    cutflow1++;
    
    if(muRecoN<=0) continue;
    cutflow2++;
    
    muFiltmult->Fill(muFiltN);
    mumult->Fill(muRecoN);

    for(unsigned int mucnt=0; mucnt<muFiltN; mucnt++) {
      muFiltpt->Fill(muFiltPt[mucnt]);
      muFilteta->Fill(muFiltEta[mucnt]);
      muFiltphi->Fill(muFiltPhi[mucnt]);
    }
    for(unsigned int mucnt=0; mucnt<muRecoN; mucnt++) {
      mupt->Fill(muRecoPt[mucnt]);
      mueta->Fill(muRecoEta[mucnt]);
      muphi->Fill(muRecoPhi[mucnt]);
      mudxy->Fill(muRecoDxy[mucnt]);
      mulog10dxy->Fill(TMath::Log10(std::abs(muRecoDxy[mucnt])));
      mudxysig->Fill(muRecoDxySig[mucnt]);
    }
    if(htFiltN>0) {
      ht->Fill(htFiltPt[0]);
      mht->Fill(htFiltPt[1]);
      mhtphi->Fill(htFiltPhi[1]);
    }

    if(muRecoN<=1) continue;
    cutflow3++;
    
    sort(muRecoPt, muRecoEta, muRecoPhi, muRecoDxy, muRecoDxySig, muRecoN);

    // Fill the events for Run 2 DisplacedDoubleMu33 trigger
    bool foundmu1 = false;
    int pos1 = -1;
    bool foundmu2 = false;
    int pos2 = -1;
    for(int ctr=0; ctr<muRecoN; ctr++) {
      if(muRecoPt[ctr]<33) continue;
      if(TMath::Abs(muRecoEta[ctr])>2.5) continue;
      if(TMath::Abs(muRecoDxy[ctr])<0.01) continue;

      if(!foundmu1) {
	pos1 = ctr;
	foundmu1 = true;
	continue;
      }
      if(!foundmu2) {
	pos2 = ctr;
	foundmu2 = true;
	continue;
      }
    }
    if(foundmu1 && foundmu2) {
      cutflowdatarate++;
      TLorentzVector *mu0 = new TLorentzVector();
      TLorentzVector *mu1 = new TLorentzVector();
      TLorentzVector *mht = new TLorentzVector();
      mu0->SetPtEtaPhiM(muRecoPt[pos1],muRecoEta[pos1],muRecoPhi[pos1],0.1057);
      mu1->SetPtEtaPhiM(muRecoPt[pos2],muRecoEta[pos2],muRecoPhi[pos2],0.1057);
      mht->SetPtEtaPhiM(htFiltPt[1],0,htFiltPhi[1],0);
      std::vector<TLorentzVector*> lvarr;
      lvarr.push_back(mu0);
      lvarr.push_back(mu1);
      orig_mu1pt->Fill(muRecoPt[pos1]);
      orig_mu2pt->Fill(muRecoPt[pos2]);
      double minDXY = TMath::Abs(muRecoDxy[pos1])<TMath::Abs(muRecoDxy[pos2])?muRecoDxy[pos1]:muRecoDxy[pos2];
      orig_mindxy->Fill(minDXY);
      orig_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      orig_mumuM->Fill(((*mu0)+(*mu1)).M());
      orig_mumudR->Fill(mu0->DeltaR((*mu1)));
      orig_mu1dxysig->Fill(muRecoDxySig[pos1]);
      orig_mu2dxysig->Fill(muRecoDxySig[pos2]);
      orig_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1]));
      orig_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2]));
      orig_ht->Fill(htFiltPt[0]);
      orig_mht->Fill(htFiltPt[1]);
      orig_sphericity->Fill(sphericity(lvarr));
      orig_spherocity->Fill(transversespherocity(lvarr));
      /*  auto orig_mumult = new TH1F("orig_mumult","orig N #mu",10,0,10);
	  auto orig_mudxy = new TH1F("orig_mudxy","orig #mu d_{0} / cm",20000,-20,20);
	  auto orig_mulog10dxy = new TH1F("orig_mulog10dxy","orig #mu log_{10}d_{0}",500,-3,2);
      */
    }

    // Define variables for the muons with highest pT
    TLorentzVector *mu0 = new TLorentzVector();
    TLorentzVector *mu1 = new TLorentzVector();
    TLorentzVector *mht = new TLorentzVector();
    mu0->SetPtEtaPhiM(muRecoPt[0],muRecoEta[0],muRecoPhi[0],0.1057);
    mu1->SetPtEtaPhiM(muRecoPt[1],muRecoEta[1],muRecoPhi[1],0.1057);
    mht->SetPtEtaPhiM(htFiltPt[1],0,htFiltPhi[1],0);
    std::vector<TLorentzVector*> lvarr;
    lvarr.push_back(mu0);
    lvarr.push_back(mu1);

    double dEta = std::abs(mu0->Eta()-mu1->Eta());
    double dPhi = std::abs(mu0->DeltaPhi((*mu1)));
    double dR = mu0->DeltaR((*mu1));
    double mu2mhtdPhi = std::abs(mht->DeltaPhi((*mu1)));
    double mumumhtdphi = std::abs(mht->DeltaPhi((*mu1)+(*mu0)));
    double mumuinvm = ((*mu0)+(*mu1)).M();
    double minDxy = TMath::Abs(muRecoDxy[0])<TMath::Abs(muRecoDxy[1])?TMath::Abs(muRecoDxy[0]):TMath::Abs(muRecoDxy[1]);

    // Fill all the muons for the selected events
    mu2pt->Fill(muRecoPt[1]);
    mu2eta->Fill(muRecoEta[1]);
    mu2phi->Fill(muRecoPhi[1]);
    mu2dxy->Fill(muRecoDxy[1]);
    mu2log10dxy->Fill(TMath::Log10(std::abs(muRecoDxy[1])));
    mu2dxysig->Fill(muRecoDxySig[1]);
    mindxy->Fill(*std::min_element(muRecoDxy,muRecoDxy+muRecoN));
    minlog10dxy->Fill(TMath::Log10(std::abs(*std::min_element(muRecoDxy,muRecoDxy+muRecoN))));
    

    // Fill variables for selected events
    std::vector<TLorentzVector*> all_muarr;
    for(unsigned int mu1=0; mu1<muRecoN; mu1++) {
      TLorentzVector *mutemp1 = new TLorentzVector();
      mutemp1->SetPtEtaPhiM(muRecoPt[mu1],muRecoEta[mu1],muRecoPhi[mu1],0.1057);
      all_muarr.push_back(mutemp1);
      for(unsigned int mu2=mu1+1; mu2<muRecoN; mu2++) {
	TLorentzVector *mutemp2 = new TLorentzVector();
	mutemp2->SetPtEtaPhiM(muRecoPt[mu2],muRecoEta[mu2],muRecoPhi[mu2],0.1057);
	mumuM->Fill(((*mutemp1)+(*mutemp2)).M());
	mumulog10M->Fill(TMath::Log10(TMath::Abs(((*mutemp1)+(*mutemp2)).M())));
	mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
	mumudEta->Fill(TMath::Abs(mutemp1->Eta()-mutemp2->Eta()));
	mumudPhi->Fill(mutemp1->DeltaPhi((*mutemp2)));
	mumudPhivdEta->Fill(mutemp1->DeltaPhi((*mutemp2)),TMath::Abs(mutemp1->Eta()-mutemp2->Eta()));
	muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
	muetamueta->Fill(mutemp1->Eta(),mutemp2->Eta());
	muphimuphi->Fill(mutemp1->Phi(),mutemp2->Phi());
      }
    }
    h_sphericity->Fill(sphericity(all_muarr));
    h_spherocity->Fill(transversespherocity(all_muarr));
    mu2MhtdPhi->Fill(mu2mhtdPhi);
    mumuMhtdPhi->Fill(mumumhtdphi);

    sel_mu1pt->Fill(muRecoPt[0]);
    sel_mu2pt->Fill(muRecoPt[1]);
    sel_muptmupt->Fill(muRecoPt[0],muRecoPt[1]);
    sel_mindxy->Fill(minDxy);
    sel_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDxy)));
    sel_mumuM->Fill(mumuinvm);
    sel_mudxy->Fill(muRecoDxy[0]);
    sel_mudxy->Fill(muRecoDxy[1]);
    sel_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[0])));
    sel_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[1])));
    sel_mu1dxysig->Fill(muRecoDxySig[0]);
    sel_mu2dxysig->Fill(muRecoDxySig[1]);
    sel_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[0]));
    sel_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[1]));
    sel_mumudR->Fill(dR);
    sel_mumudEta->Fill(dEta);
    sel_mumudPhi->Fill(dPhi);
    sel_sphericity->Fill(sphericity(lvarr));
    sel_spherocity->Fill(transversespherocity(lvarr));
    sel_ht->Fill(htFiltPt[0]);
    sel_mht->Fill(htFiltPt[1]);
    sel_mu1dxysigabsdxy->Fill(muRecoDxySig[0], TMath::Abs(muRecoDxy[0]));
    sel_mu2dxysigabsdxy->Fill(muRecoDxySig[1], TMath::Abs(muRecoDxy[1]));
    
    if(mumuinvm>=84 && mumuinvm<=98) {
      Z_mu2pt->Fill(muRecoPt[1]);
      Z_mindxy->Fill(minDxy);
      Z_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDxy)));
      Z_mumuM->Fill(mumuinvm);
      Z_mudxy->Fill(muRecoDxy[0]);
      Z_mudxy->Fill(muRecoDxy[1]);
      Z_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[0])));
      Z_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[1])));
      Z_mu1dxysig->Fill(muRecoDxySig[0]);
      Z_mu2dxysig->Fill(muRecoDxySig[1]);
      Z_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[0]));
      Z_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[1]));
      Z_mumudR->Fill(mu0->DeltaR((*mu1)));
      Z_muptmupt->Fill(mu0->Pt(),mu1->Pt());
      Z_muetamueta->Fill(mu0->Eta(),mu1->Eta());
      Z_muphimuphi->Fill(mu0->Phi(),mu1->Phi());
      if(abs(muRecoDxy[0])>abs(muRecoDxy[1])) {
	Z_Mdxy->Fill(muRecoDxy[0],mumuinvm);
	Z_dxysigdxy->Fill(muRecoDxy[0],muRecoDxySig[0]);
      }
      else {
	Z_Mdxy->Fill(muRecoDxy[1],mumuinvm);
	Z_dxysigdxy->Fill(muRecoDxy[1],muRecoDxySig[1]);
      }
    }
    else if(mumuinvm>=2 && mumuinvm<=4) {
      Jpsi_mu2pt->Fill(muRecoPt[1]);
      Jpsi_mindxy->Fill(minDxy);
      Jpsi_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDxy)));
      Jpsi_mumuM->Fill(mumuinvm);
      Jpsi_mudxy->Fill(muRecoDxy[0]);
      Jpsi_mudxy->Fill(muRecoDxy[1]);
      Jpsi_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[0])));
      Jpsi_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[1])));
      Jpsi_mu1dxysig->Fill(muRecoDxySig[0]);
      Jpsi_mu2dxysig->Fill(muRecoDxySig[1]);
      Jpsi_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[0]));
      Jpsi_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[1]));
      Jpsi_mumudR->Fill(mu0->DeltaR((*mu1)));
      Jpsi_muptmupt->Fill(mu0->Pt(),mu1->Pt());
      Jpsi_muetamueta->Fill(mu0->Eta(),mu1->Eta());
      Jpsi_muphimuphi->Fill(mu0->Phi(),mu1->Phi());
    }
    else if(mumuinvm>=35 && mumuinvm<=84) {
      Mplt_mu2pt->Fill(muRecoPt[1]);
      Mplt_mindxy->Fill(minDxy);
      Mplt_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDxy)));
      Mplt_mumuM->Fill(mumuinvm);
      Mplt_mudxy->Fill(muRecoDxy[0]);
      Mplt_mudxy->Fill(muRecoDxy[1]);
      Mplt_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[0])));
      Mplt_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[1])));
      Mplt_mu1dxysig->Fill(muRecoDxySig[0]);
      Mplt_mu2dxysig->Fill(muRecoDxySig[1]);
      Mplt_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[0]));
      Mplt_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[1]));
      Mplt_mumudR->Fill(mu0->DeltaR((*mu1)));
      Mplt_muptmupt->Fill(mu0->Pt(),mu1->Pt());
      Mplt_muetamueta->Fill(mu0->Eta(),mu1->Eta());
      Mplt_muphimuphi->Fill(mu0->Phi(),mu1->Phi());
    }
    else {
      oth_mu2pt->Fill(muRecoPt[1]);
      oth_mindxy->Fill(minDxy);
      oth_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDxy)));
      oth_mumuM->Fill(mumuinvm);
      oth_mudxy->Fill(muRecoDxy[0]);
      oth_mudxy->Fill(muRecoDxy[1]);
      oth_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[0])));
      oth_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[1])));
      oth_mu1dxysig->Fill(muRecoDxySig[0]);
      oth_mu2dxysig->Fill(muRecoDxySig[1]);
      oth_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[0]));
      oth_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[1]));
      oth_mumudR->Fill(mu0->DeltaR((*mu1)));
      oth_muptmupt->Fill(mu0->Pt(),mu1->Pt());
      oth_muetamueta->Fill(mu0->Eta(),mu1->Eta());
      oth_muphimuphi->Fill(mu0->Phi(),mu1->Phi());
    }
    
    // Fill histograms for cut 1
    std::vector<TLorentzVector*> cut1_muarr;
    bool cut1cond = false;
    TLorentzVector *mutemp1 = new TLorentzVector();
    int pos1_cut1 = -1;
    TLorentzVector *mutemp2 = new TLorentzVector();
    int pos2_cut1 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut1 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut1 = cnt2;
	double deltaR = mutemp1->DeltaR((*mutemp2));
	if(deltaR>1) {
	  cut1cond = true;
	  break;
	}
      }
      if(cut1cond) break;
    }
    if(cut1cond) {
      cutflowcut1++;
      cut1_muarr.push_back(mutemp1);
      cut1_muarr.push_back(mutemp2);
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut1])<TMath::Abs(muRecoDxy[pos2_cut1])?muRecoDxy[pos1_cut1]:muRecoDxy[pos2_cut1];
      cut1_mumuM->Fill(invM);
      cut1_mumult->Fill(muRecoN);
      cut1_mu1pt->Fill(mutemp1->Pt());
      cut1_mu2pt->Fill(mutemp2->Pt());
      cut1_mindxy->Fill(minDXY);
      cut1_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut1_mudxy->Fill(muRecoDxy[pos1_cut1]);
      cut1_mudxy->Fill(muRecoDxy[pos2_cut1]);
      cut1_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos1_cut1])));
      cut1_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos2_cut1])));
      cut1_mu1dxysig->Fill(muRecoDxySig[pos1_cut1]);
      cut1_mu2dxysig->Fill(muRecoDxySig[pos2_cut1]);
      cut1_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut1]));
      cut1_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut1]));
      cut1_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut1_sphericity->Fill(sphericity(cut1_muarr));
      cut1_spherocity->Fill(transversespherocity(cut1_muarr));
      cut1_ht->Fill(htFiltPt[0]);
      cut1_mht->Fill(htFiltPt[1]);
      cut1_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }
        
    // Fill histograms for cut 2
    std::vector<TLorentzVector*> cut2_muarr;
    bool cut2cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut2 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut2 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut2 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut2 = cnt2;
	double invM = ((*mutemp1)+(*mutemp2)).M();
	if(invM>98 || invM<84) {
	  cut2cond = true;
	  break;
	}
      }
      if(cut2cond) break;
    }
    if(cut2cond) {
      cutflowcut2++;
      cut2_muarr.push_back(mutemp1);
      cut2_muarr.push_back(mutemp2);
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut2])<TMath::Abs(muRecoDxy[pos2_cut2])?muRecoDxy[pos1_cut2]:muRecoDxy[pos2_cut2];
      cut2_mumuM->Fill(invM);
      cut2_mumult->Fill(muRecoN);
      cut2_mu1pt->Fill(mutemp1->Pt());
      cut2_mu2pt->Fill(mutemp2->Pt());
      cut2_mindxy->Fill(minDXY);
      cut2_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut2_mudxy->Fill(muRecoDxy[pos1_cut2]);
      cut2_mudxy->Fill(muRecoDxy[pos2_cut2]);
      cut2_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos1_cut2])));
      cut2_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos2_cut2])));
      cut2_mu1dxysig->Fill(muRecoDxySig[pos1_cut2]);
      cut2_mu2dxysig->Fill(muRecoDxySig[pos2_cut2]);
      cut2_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut2]));
      cut2_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut2]));
      cut2_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut2_sphericity->Fill(sphericity(cut2_muarr));
      cut2_spherocity->Fill(transversespherocity(cut2_muarr));
      cut2_ht->Fill(htFiltPt[0]);
      cut2_mht->Fill(htFiltPt[1]);
      cut2_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

    // Fill histograms for cut 3
    std::vector<TLorentzVector*> cut3_muarr;
    bool cut3cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut3 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut3 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut3 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut3 = cnt2;
	double deltaR = mutemp1->DeltaR((*mutemp2));
	double invM = ((*mutemp1)+(*mutemp2)).M();
	if((invM>98 || invM<84)&&(deltaR>1)) {
	  cut3cond = true;
	  break;
	}
      }
      if(cut3cond) break;
    }
    if(cut3cond) {
      cutflowcut3++;
      cut3_muarr.push_back(mutemp1);
      cut3_muarr.push_back(mutemp2);
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut3])<TMath::Abs(muRecoDxy[pos2_cut3])?muRecoDxy[pos1_cut3]:muRecoDxy[pos2_cut3];
      cut3_mumuM->Fill(invM);
      cut3_mumult->Fill(muRecoN);
      cut3_mu1pt->Fill(mutemp1->Pt());
      cut3_mu2pt->Fill(mutemp2->Pt());
      cut3_mindxy->Fill(minDXY);
      cut3_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      double mindxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut3])<TMath::Log10(muRecoDxySig[pos2_cut3])?TMath::Log10(muRecoDxySig[pos1_cut3]):TMath::Log10(muRecoDxySig[pos2_cut3]);
      double maxdxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut3])>TMath::Log10(muRecoDxySig[pos2_cut3])?TMath::Log10(muRecoDxySig[pos1_cut3]):TMath::Log10(muRecoDxySig[pos2_cut3]);
      cut3_minlog10dxysig->Fill(mindxysigl10);
      cut3_maxlog10dxysig->Fill(maxdxysigl10);
      cut3_mudxy->Fill(muRecoDxy[pos1_cut3]);
      cut3_mudxy->Fill(muRecoDxy[pos2_cut3]);
      cut3_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos1_cut3])));
      cut3_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos2_cut3])));
      cut3_mu1dxysig->Fill(muRecoDxySig[pos1_cut3]);
      cut3_mu2dxysig->Fill(muRecoDxySig[pos2_cut3]);
      cut3_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut3]));
      cut3_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut3]));
      cut3_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut3_sphericity->Fill(sphericity(cut3_muarr));
      cut3_spherocity->Fill(transversespherocity(cut3_muarr));
      cut3_ht->Fill(htFiltPt[0]);
      cut3_mht->Fill(htFiltPt[1]);
      cut3_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

    // Fill histograms for cut 4
    std::vector<TLorentzVector*> cut4_muarr;
    bool cut4cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut4 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut4 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      if(TMath::Abs(muRecoDxy[cnt1])<0.01) continue;
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut4 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	if(TMath::Abs(muRecoDxy[cnt2])<0.01) continue;

	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut4 = cnt2;
	cut4cond = true;
	break;
      }
      if(cut4cond) break;
    }
    if(cut4cond) {
      cutflowcut4++;
      cut4_muarr.push_back(mutemp1);
      cut4_muarr.push_back(mutemp2);
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut4])<TMath::Abs(muRecoDxy[pos2_cut4])?muRecoDxy[pos1_cut4]:muRecoDxy[pos2_cut4];
      cut4_mumuM->Fill(invM);
      cut4_mumult->Fill(muRecoN);
      cut4_mu1pt->Fill(mutemp1->Pt());
      cut4_mu2pt->Fill(mutemp2->Pt());
      cut4_mindxy->Fill(minDXY);
      cut4_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut4_mudxy->Fill(muRecoDxy[pos1_cut4]);
      cut4_mudxy->Fill(muRecoDxy[pos2_cut4]);
      cut4_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos1_cut4])));
      cut4_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos2_cut4])));
      cut4_mu1dxysig->Fill(muRecoDxySig[pos1_cut4]);
      cut4_mu2dxysig->Fill(muRecoDxySig[pos2_cut4]);
      cut4_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut4]));
      cut4_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut4]));
      cut4_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut4_sphericity->Fill(sphericity(cut4_muarr));
      cut4_spherocity->Fill(transversespherocity(cut4_muarr));
      cut4_ht->Fill(htFiltPt[0]);
      cut4_mht->Fill(htFiltPt[1]);
      cut4_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

    // Fill histograms for cut 5
    std::vector<TLorentzVector*> cut5_muarr;
    bool cut5cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut5 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut5 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      if(TMath::Abs(muRecoDxySig[cnt1])<0.1) continue;
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut5 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	if(TMath::Abs(muRecoDxySig[cnt2])<0.1) continue;

	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut5 = cnt2;
	double deltaR = mutemp1->DeltaR((*mutemp2));
	double invM = ((*mutemp1)+(*mutemp2)).M();
	if((invM>98 || invM<84)&&(deltaR>1)) {
	  cut5cond = true;
	  break;
	}
      }
      if(cut5cond) break;
    }
    if(cut5cond) {
      cutflowcut5++;
      cut5_muarr.push_back(mutemp1);
      cut5_muarr.push_back(mutemp2);
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut5])<TMath::Abs(muRecoDxy[pos2_cut5])?muRecoDxy[pos1_cut5]:muRecoDxy[pos2_cut5];
      cut5_mumuM->Fill(invM);
      cut5_mumult->Fill(muRecoN);
      cut5_mu1pt->Fill(mutemp1->Pt());
      cut5_mu2pt->Fill(mutemp2->Pt());
      cut5_mindxy->Fill(minDXY);
      cut5_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut5_mudxy->Fill(muRecoDxy[pos1_cut5]);
      cut5_mudxy->Fill(muRecoDxy[pos2_cut5]);
      cut5_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos1_cut5])));
      cut5_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos2_cut5])));
      cut5_mu1dxysig->Fill(muRecoDxySig[pos1_cut5]);
      cut5_mu2dxysig->Fill(muRecoDxySig[pos2_cut5]);
      cut5_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut5]));
      cut5_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut5]));
      double mindxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut5])<TMath::Log10(muRecoDxySig[pos2_cut5])?TMath::Log10(muRecoDxySig[pos1_cut5]):TMath::Log10(muRecoDxySig[pos2_cut5]);
      double maxdxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut5])>TMath::Log10(muRecoDxySig[pos2_cut5])?TMath::Log10(muRecoDxySig[pos1_cut5]):TMath::Log10(muRecoDxySig[pos2_cut5]);
      cut5_minlog10dxysig->Fill(mindxysigl10);
      cut5_maxlog10dxysig->Fill(maxdxysigl10);
      cut5_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut5_sphericity->Fill(sphericity(cut5_muarr));
      cut5_spherocity->Fill(transversespherocity(cut5_muarr));
      cut5_ht->Fill(htFiltPt[0]);
      cut5_mht->Fill(htFiltPt[1]);
      cut5_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

    // Fill histograms for cut 6
    std::vector<TLorentzVector*> cut6_muarr;
    bool cut6cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut6 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut6 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      if(TMath::Abs(muRecoDxySig[cnt1])<1) continue;
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut6 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	if(TMath::Abs(muRecoDxySig[cnt2])<1) continue;

	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut6 = cnt2;
	double deltaR = mutemp1->DeltaR((*mutemp2));
	double invM = ((*mutemp1)+(*mutemp2)).M();
	if((invM>98 || invM<84)&&(deltaR>1)) {
	  cut6cond = true;
	  break;
	}
      }
      if(cut6cond) break;
    }
    if(cut6cond) {
      cutflowcut6++;
      cut6_muarr.push_back(mutemp1);
      cut6_muarr.push_back(mutemp2);
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut6])<TMath::Abs(muRecoDxy[pos2_cut6])?muRecoDxy[pos1_cut6]:muRecoDxy[pos2_cut6];
      cut6_mumuM->Fill(invM);
      cut6_mumult->Fill(muRecoN);
      cut6_mu1pt->Fill(mutemp1->Pt());
      cut6_mu2pt->Fill(mutemp2->Pt());
      cut6_mindxy->Fill(minDXY);
      cut6_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut6_mudxy->Fill(muRecoDxy[pos1_cut6]);
      cut6_mudxy->Fill(muRecoDxy[pos2_cut6]);
      cut6_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos1_cut6])));
      cut6_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos2_cut6])));
      cut6_mu1dxysig->Fill(muRecoDxySig[pos1_cut6]);
      cut6_mu2dxysig->Fill(muRecoDxySig[pos2_cut6]);
      cut6_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut6]));
      cut6_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut6]));
      double mindxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut6])<TMath::Log10(muRecoDxySig[pos2_cut6])?TMath::Log10(muRecoDxySig[pos1_cut6]):TMath::Log10(muRecoDxySig[pos2_cut6]);
      double maxdxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut6])>TMath::Log10(muRecoDxySig[pos2_cut6])?TMath::Log10(muRecoDxySig[pos1_cut6]):TMath::Log10(muRecoDxySig[pos2_cut6]);
      cut6_minlog10dxysig->Fill(mindxysigl10);
      cut6_maxlog10dxysig->Fill(maxdxysigl10);
      cut6_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut6_sphericity->Fill(sphericity(cut6_muarr));
      cut6_spherocity->Fill(transversespherocity(cut6_muarr));
      cut6_ht->Fill(htFiltPt[0]);
      cut6_mht->Fill(htFiltPt[1]);
      cut6_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

    // Fill histograms for cut 7
    std::vector<TLorentzVector*> cut7_muarr;
    bool cut7cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut7 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut7 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      if(TMath::Abs(muRecoDxySig[cnt1])<10) continue;
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut7 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	if(TMath::Abs(muRecoDxySig[cnt2])<10) continue;

	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut7 = cnt2;
	double deltaR = mutemp1->DeltaR((*mutemp2));
	double invM = ((*mutemp1)+(*mutemp2)).M();
	if((invM>98 || invM<84)&&(deltaR>1)) {
	  cut7cond = true;
	  break;
	}
      }
      if(cut7cond) break;
    }
    if(cut7cond) {
      cutflowcut7++;
      cut7_muarr.push_back(mutemp1);
      cut7_muarr.push_back(mutemp2);
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut7])<TMath::Abs(muRecoDxy[pos2_cut7])?muRecoDxy[pos1_cut7]:muRecoDxy[pos2_cut7];
      cut7_mumuM->Fill(invM);
      cut7_mumult->Fill(muRecoN);
      cut7_mu1pt->Fill(mutemp1->Pt());
      cut7_mu2pt->Fill(mutemp2->Pt());
      cut7_mindxy->Fill(minDXY);
      cut7_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut7_mudxy->Fill(muRecoDxy[pos1_cut7]);
      cut7_mudxy->Fill(muRecoDxy[pos2_cut7]);
      cut7_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos1_cut7])));
      cut7_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos2_cut7])));
      cut7_mu1dxysig->Fill(muRecoDxySig[pos1_cut7]);
      cut7_mu2dxysig->Fill(muRecoDxySig[pos2_cut7]);
      cut7_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut7]));
      cut7_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut7]));
      double mindxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut7])<TMath::Log10(muRecoDxySig[pos2_cut7])?TMath::Log10(muRecoDxySig[pos1_cut7]):TMath::Log10(muRecoDxySig[pos2_cut7]);
      double maxdxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut7])>TMath::Log10(muRecoDxySig[pos2_cut7])?TMath::Log10(muRecoDxySig[pos1_cut7]):TMath::Log10(muRecoDxySig[pos2_cut7]);
      cut7_minlog10dxysig->Fill(mindxysigl10);
      cut7_maxlog10dxysig->Fill(maxdxysigl10);
      cut7_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut7_sphericity->Fill(sphericity(cut7_muarr));
      cut7_spherocity->Fill(transversespherocity(cut7_muarr));
      cut7_ht->Fill(htFiltPt[0]);
      cut7_mht->Fill(htFiltPt[1]);
      cut7_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

    // Fill histograms for cut 8
    std::vector<TLorentzVector*> cut8_muarr;
    bool cut8cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut8 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut8 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      if(TMath::Abs(muRecoDxy[cnt1])<0.01) continue;
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut8 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	if(TMath::Abs(muRecoDxy[cnt2])<0.01) continue;

	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut8 = cnt2;
	double deltaR = mutemp1->DeltaR((*mutemp2));
	double invM = ((*mutemp1)+(*mutemp2)).M();
	if((invM>98 || invM<84)&&(deltaR>1)) {
	  cut8cond = true;
	  break;
	}
      }
      if(cut8cond) break;
    }
    if(cut8cond) {
      cutflowcut8++;
      cut8_muarr.push_back(mutemp1);
      cut8_muarr.push_back(mutemp2);
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut8])<TMath::Abs(muRecoDxy[pos2_cut8])?muRecoDxy[pos1_cut8]:muRecoDxy[pos2_cut8];
      cut8_mumuM->Fill(invM);
      cut8_mumult->Fill(muRecoN);
      cut8_mu1pt->Fill(mutemp1->Pt());
      cut8_mu2pt->Fill(mutemp2->Pt());
      cut8_mindxy->Fill(minDXY);
      cut8_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut8_mudxy->Fill(muRecoDxy[pos1_cut8]);
      cut8_mudxy->Fill(muRecoDxy[pos2_cut8]);
      cut8_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos1_cut8])));
      cut8_mulog10dxy->Fill(TMath::Log10(TMath::Abs(muRecoDxy[pos2_cut8])));
      cut8_mu1dxysig->Fill(muRecoDxySig[pos1_cut8]);
      cut8_mu2dxysig->Fill(muRecoDxySig[pos2_cut8]);
      cut8_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut8]));
      cut8_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut8]));
      cut8_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut8_sphericity->Fill(sphericity(cut8_muarr));
      cut8_spherocity->Fill(transversespherocity(cut8_muarr));
      cut8_ht->Fill(htFiltPt[0]);
      cut8_mht->Fill(htFiltPt[1]);
      cut8_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

    // Fill histograms for cut 9
    bool cut9cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut9 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut9 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      if(TMath::Abs(muRecoDxy[cnt1])<0.01 || TMath::Abs(muRecoDxySig[cnt1])<1) continue;
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut9 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	if(TMath::Abs(muRecoDxy[cnt2])<0.01 || TMath::Abs(muRecoDxySig[cnt2])<1) continue;

	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut9 = cnt2;
	double deltaR = mutemp1->DeltaR((*mutemp2));
	double invM = ((*mutemp1)+(*mutemp2)).M();
	if((invM>98 || invM<84)&&(deltaR>1)) {
	  cut9cond = true;
	  break;
	}
      }
      if(cut9cond) break;
    }
    if(cut9cond) {
      cutflowcut9++;
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut9])<TMath::Abs(muRecoDxy[pos2_cut9])?muRecoDxy[pos1_cut9]:muRecoDxy[pos2_cut9];
      double maxDXY = TMath::Abs(muRecoDxy[pos1_cut9])>TMath::Abs(muRecoDxy[pos2_cut9])?muRecoDxy[pos1_cut9]:muRecoDxy[pos2_cut9];
      cut9_mumuM->Fill(invM);
      cut9_mumult->Fill(muRecoN);
      cut9_mu1pt->Fill(mutemp1->Pt());
      cut9_mu2pt->Fill(mutemp2->Pt());
      cut9_mindxy->Fill(minDXY);
      cut9_maxdxy->Fill(maxDXY);
      cut9_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut9_maxlog10dxy->Fill(TMath::Log10(TMath::Abs(maxDXY)));
      cut9_mu1dxysig->Fill(muRecoDxySig[pos1_cut9]);
      cut9_mu2dxysig->Fill(muRecoDxySig[pos2_cut9]);
      cut9_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut9]));
      cut9_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut9]));
      double mindxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut9])<TMath::Log10(muRecoDxySig[pos2_cut9])?TMath::Log10(muRecoDxySig[pos1_cut9]):TMath::Log10(muRecoDxySig[pos2_cut9]);
      double maxdxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut9])>TMath::Log10(muRecoDxySig[pos2_cut9])?TMath::Log10(muRecoDxySig[pos1_cut9]):TMath::Log10(muRecoDxySig[pos2_cut9]);
      cut9_minlog10dxysig->Fill(mindxysigl10);
      cut9_maxlog10dxysig->Fill(maxdxysigl10);
      cut9_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut9_ht->Fill(htFiltPt[0]);
      cut9_mht->Fill(htFiltPt[1]);
      cut9_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

    // Fill histograms for cut 10
    bool cut10cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut10 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut10 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      if(muRecoPt[cnt1]<20) continue;
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut10 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	if(muRecoPt[cnt2]<20) continue;

	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut10 = cnt2;
	double deltaR = mutemp1->DeltaR((*mutemp2));
	double invM = ((*mutemp1)+(*mutemp2)).M();
	if((invM>98 || invM<84)&&(deltaR>1)) {
	  cut10cond = true;
	  break;
	}
      }
      if(cut10cond) break;
    }
    if(cut10cond) {
      cutflowcut10++;
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut10])<TMath::Abs(muRecoDxy[pos2_cut10])?muRecoDxy[pos1_cut10]:muRecoDxy[pos2_cut10];
      double maxDXY = TMath::Abs(muRecoDxy[pos1_cut10])>TMath::Abs(muRecoDxy[pos2_cut10])?muRecoDxy[pos1_cut10]:muRecoDxy[pos2_cut10];
      cut10_mumuM->Fill(invM);
      cut10_mumult->Fill(muRecoN);
      cut10_mu1pt->Fill(mutemp1->Pt());
      cut10_mu2pt->Fill(mutemp2->Pt());
      cut10_mindxy->Fill(minDXY);
      cut10_maxdxy->Fill(maxDXY);
      cut10_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut10_maxlog10dxy->Fill(TMath::Log10(TMath::Abs(maxDXY)));
      cut10_mu1dxysig->Fill(muRecoDxySig[pos1_cut10]);
      cut10_mu2dxysig->Fill(muRecoDxySig[pos2_cut10]);
      cut10_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut10]));
      cut10_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut10]));
      double mindxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut10])<TMath::Log10(muRecoDxySig[pos2_cut10])?TMath::Log10(muRecoDxySig[pos1_cut10]):TMath::Log10(muRecoDxySig[pos2_cut10]);
      double maxdxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut10])>TMath::Log10(muRecoDxySig[pos2_cut10])?TMath::Log10(muRecoDxySig[pos1_cut10]):TMath::Log10(muRecoDxySig[pos2_cut10]);
      cut10_minlog10dxysig->Fill(mindxysigl10);
      cut10_maxlog10dxysig->Fill(maxdxysigl10);
      cut10_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut10_ht->Fill(htFiltPt[0]);
      cut10_mht->Fill(htFiltPt[1]);
      cut10_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

    // Fill histograms for cut 11
    bool cut11cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut11 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut11 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      if(TMath::Abs(muRecoDxySig[cnt1])<1 || muRecoPt[cnt1]<20) continue;
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut11 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	if(TMath::Abs(muRecoDxySig[cnt2])<1 || muRecoPt[cnt2]<20) continue;

	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut11 = cnt2;
	double deltaR = mutemp1->DeltaR((*mutemp2));
	double invM = ((*mutemp1)+(*mutemp2)).M();
	if((invM>98 || invM<84)&&(deltaR>1)) {
	  cut11cond = true;
	  break;
	}
      }
      if(cut11cond) break;
    }
    if(cut11cond) {
      cutflowcut11++;
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut11])<TMath::Abs(muRecoDxy[pos2_cut11])?muRecoDxy[pos1_cut11]:muRecoDxy[pos2_cut11];
      double maxDXY = TMath::Abs(muRecoDxy[pos1_cut11])>TMath::Abs(muRecoDxy[pos2_cut11])?muRecoDxy[pos1_cut11]:muRecoDxy[pos2_cut11];
      cut11_mumuM->Fill(invM);
      cut11_mumult->Fill(muRecoN);
      cut11_mu1pt->Fill(mutemp1->Pt());
      cut11_mu2pt->Fill(mutemp2->Pt());
      cut11_mindxy->Fill(minDXY);
      cut11_maxdxy->Fill(maxDXY);
      cut11_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut11_maxlog10dxy->Fill(TMath::Log10(TMath::Abs(maxDXY)));
      cut11_mu1dxysig->Fill(muRecoDxySig[pos1_cut11]);
      cut11_mu2dxysig->Fill(muRecoDxySig[pos2_cut11]);
      cut11_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut11]));
      cut11_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut11]));
      double mindxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut11])<TMath::Log10(muRecoDxySig[pos2_cut11])?TMath::Log10(muRecoDxySig[pos1_cut11]):TMath::Log10(muRecoDxySig[pos2_cut11]);
      double maxdxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut11])>TMath::Log10(muRecoDxySig[pos2_cut11])?TMath::Log10(muRecoDxySig[pos1_cut11]):TMath::Log10(muRecoDxySig[pos2_cut11]);
      cut11_minlog10dxysig->Fill(mindxysigl10);
      cut11_maxlog10dxysig->Fill(maxdxysigl10);
      cut11_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut11_ht->Fill(htFiltPt[0]);
      cut11_mht->Fill(htFiltPt[1]);
      cut11_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

    // Fill histograms for cut 12
    bool cut12cond = false;
    mutemp1 = new TLorentzVector();
    int pos1_cut12 = -1;
    mutemp2 = new TLorentzVector();
    int pos2_cut12 = -1;
    for(unsigned int cnt1=0; cnt1<muRecoN; cnt1++) {
      if(TMath::Abs(muRecoDxy[cnt1])<0.01 || muRecoPt[cnt1]<20) continue;
      mutemp1->SetPtEtaPhiM(muRecoPt[cnt1],muRecoEta[cnt1],muRecoPhi[cnt1],0.1057);
      pos1_cut12 = cnt1;
      for(unsigned int cnt2=cnt1+1; cnt2<muRecoN; cnt2++) {
	if(TMath::Abs(muRecoDxy[cnt2])<0.01 || muRecoPt[cnt2]<20) continue;

	mutemp2->SetPtEtaPhiM(muRecoPt[cnt2],muRecoEta[cnt2],muRecoPhi[cnt2],0.1057);
	pos2_cut12 = cnt2;
	double deltaR = mutemp1->DeltaR((*mutemp2));
	double invM = ((*mutemp1)+(*mutemp2)).M();
	if((invM>98 || invM<84)&&(deltaR>1)) {
	  cut12cond = true;
	  break;
	}
      }
      if(cut12cond) break;
    }
    if(cut12cond) {
      cutflowcut12++;
      double invM = ((*mutemp1)+(*mutemp2)).M();
      double minDXY = TMath::Abs(muRecoDxy[pos1_cut12])<TMath::Abs(muRecoDxy[pos2_cut12])?muRecoDxy[pos1_cut12]:muRecoDxy[pos2_cut12];
      double maxDXY = TMath::Abs(muRecoDxy[pos1_cut12])>TMath::Abs(muRecoDxy[pos2_cut12])?muRecoDxy[pos1_cut12]:muRecoDxy[pos2_cut12];
      cut12_mumuM->Fill(invM);
      cut12_mumult->Fill(muRecoN);
      cut12_mu1pt->Fill(mutemp1->Pt());
      cut12_mu2pt->Fill(mutemp2->Pt());
      cut12_mindxy->Fill(minDXY);
      cut12_maxdxy->Fill(maxDXY);
      cut12_minlog10dxy->Fill(TMath::Log10(TMath::Abs(minDXY)));
      cut12_maxlog10dxy->Fill(TMath::Log10(TMath::Abs(maxDXY)));
      cut12_mu1dxysig->Fill(muRecoDxySig[pos1_cut12]);
      cut12_mu2dxysig->Fill(muRecoDxySig[pos2_cut12]);
      cut12_mu1log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos1_cut12]));
      cut12_mu2log10dxysig->Fill(TMath::Log10(muRecoDxySig[pos2_cut12]));
      double mindxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut12])<TMath::Log10(muRecoDxySig[pos2_cut12])?TMath::Log10(muRecoDxySig[pos1_cut12]):TMath::Log10(muRecoDxySig[pos2_cut12]);
      double maxdxysigl10 = TMath::Log10(muRecoDxySig[pos1_cut12])>TMath::Log10(muRecoDxySig[pos2_cut12])?TMath::Log10(muRecoDxySig[pos1_cut12]):TMath::Log10(muRecoDxySig[pos2_cut12]);
      cut12_minlog10dxysig->Fill(mindxysigl10);
      cut12_maxlog10dxysig->Fill(maxdxysigl10);
      cut12_mumudR->Fill(mutemp1->DeltaR((*mutemp2)));
      cut12_ht->Fill(htFiltPt[0]);
      cut12_mht->Fill(htFiltPt[1]);
      cut12_muptmupt->Fill(mutemp1->Pt(),mutemp2->Pt());
    }

  }
  
  outfile->Write();
  outfile->Close();

  std::cout<<cutflow1<<"\t"<<cutflow2<<"\t"<<cutflow3<<"\t"<<cutflowdatarate<<"\t"<<cutflowcut1<<"\t"<<cutflowcut2<<"\t"<<cutflowcut3<<"\t"<<cutflowcut4<<"\t"<<cutflowcut5<<"\t"<<cutflowcut6<<"\t"<<cutflowcut7<<"\t"<<cutflowcut8<<"\t"<<cutflowcut9<<"\t"<<cutflowcut10<<"\t"<<cutflowcut11<<"\t"<<cutflowcut12<<std::endl;
  
  chain->Delete();
  return -1;
}

int main() {

  TString infile = "./data/Efmrl_MuMu16DisplacedSkim1.root";
  TString outfile = "hists_efmrl_MuMu16DisplacedSkim1.root";
  analyzer_mumu_singlefile(infile, outfile);
  
  return -1;
}
