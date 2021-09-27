#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 1.1/213;
double sig3cmsf = 100.0/330;
double sig30cmsf = 100.0/330;
double sig3msf = 100.0/330;
double sig1msf = 100.0/330;
TString cutdeets = "Cut details";
TFile* datafile = TFile::Open("hists_efmrl_MuMu16DisplacedSkim1.root","READ");
TFile* sig3cmfile = TFile::Open("hists_STHDM3cm_doublemu.root","READ");
TFile* sig30cmfile = TFile::Open("hists_STHDM30cm_doublemu.root","READ");
TFile* sig1mfile = TFile::Open("hists_STHDM1m_doublemu.root","READ");
TFile* sig3mfile = TFile::Open("hists_STHDM3m_doublemu.root","READ");

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
int comparemultihistnormalized(TString cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool isMConly=false, bool overflow=false, float legPos[]=(float []){0.7,0.75,0.95,1}) {

  std::vector<TString> legNam;
  legNam.push_back("c#tau = 3 cm");
  legNam.push_back("c#tau = 30 cm");
  legNam.push_back("c#tau = 1 m");
  legNam.push_back("c#tau = 3 m");

  TString histname = cutname+"_"+var;
  TH1F* datahist;
  if(!isMConly) datahist = (TH1F*) datafile->Get(histname);
  TH1F* sig3cmhist = (TH1F*) sig3cmfile->Get(histname);
  TH1F* sig30cmhist = (TH1F*) sig30cmfile->Get(histname);
  TH1F* sig1mhist = (TH1F*) sig1mfile->Get(histname);
  TH1F* sig3mhist = (TH1F*) sig3mfile->Get(histname);

  // Get the title from histogram title
  TString xtitle = sig3cmhist->GetTitle();

  std::vector<TH1F*> allhists;
  if(!isMConly) allhists.push_back(datahist);
  allhists.push_back(sig3cmhist);
  allhists.push_back(sig30cmhist);
  allhists.push_back(sig1mhist);
  allhists.push_back(sig3mhist);

  if(rebin==-1) rebin = 1;
  xbinlow = xbinlow==-1?0:xbinlow;
  xbinhigh = xbinhigh==-1?(overflow?allhists[0]->GetNbinsX()+1:allhists[0]->GetNbinsX()):xbinhigh;
  xbinlow = xbinlow/rebin;
  xbinhigh = xbinhigh/rebin;

  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {
    if(rebin!=1) allhists[histctr]->Rebin(rebin);

    // Make changes to sig and bkg to enable good basic plotting
    double err = 0.0;
    //allhists[histctr]->SetBinContent(xbinlow,allhists[histctr]->IntegralAndError(0,xbinlow,err));
    //allhists[histctr]->SetBinError(xbinlow,err);
    err = 0.0;
    if(overflow) {
      allhists[histctr]->SetBinContent(xbinhigh,allhists[histctr]->IntegralAndError(xbinhigh,allhists[histctr]->GetNbinsX()+1,err));
      allhists[histctr]->SetBinError(xbinhigh,err);
    }

    allhists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allhists[histctr]->GetXaxis()->SetTitle(xtitle);
    allhists[histctr]->GetYaxis()->SetTitle("normalized number of events (a.u.)");

    allhists[histctr]->SetTitle("");

    allhists[histctr]->SetLineWidth(2);
    if(isMConly) allhists[histctr]->SetLineColor(coloropt[histctr+1]);
    else allhists[histctr]->SetLineColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  //gStyle->SetOptStat(0);
  //allhists[0]->DrawNormalized("hist e1");
  //for(unsigned int histctr=1; histctr<allhists.size(); histctr++) {
  //allhists[histctr]->DrawNormalized("hist e1 same");
  //}
  //c1->SetLogy(logY);
  c1 = enhance_plotter(allhists, legNam, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY);
  c1->SaveAs("./dirplots/"+cutname+"/"+cutname+"_"+var+".png");
  
  return -1;
}

int plotter() {
  
  comparemultihistnormalized("genmumu", "log10d0", 20, 80, -1, false, true, false, (float []){0.15,0.65,0.35,0.95});
  comparemultihistnormalized("genmumu", "pt", 0, 60, -1, true, true, false, (float []){0.7,0.65,0.9,0.95});
  comparemultihistnormalized("genmumu", "eta", 2, -1, 2, true, true, false, (float []){0.4,0.2,0.6,0.6});
  
  return -1;
}
