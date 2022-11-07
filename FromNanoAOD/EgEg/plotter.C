#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

TString cutdeets = "Cut details";
TFile* metdatafile = TFile::Open("hists_metdata_55608.root","READ");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};

int efficiency(std::vector<TFile*> file, std::vector<TString> cutnames, double legpos[2], int nbins=0, double *rebin=0, TString xaxistitle="p_{T} / GeV") {

  if(cutnames.size()/2!=file.size()) {
    cout<<"Inconsistent entries for file size and cutnames: Suggested cutnameSize = 2* fileSize"<<endl;
    return -1;
  }
  
  std::vector<TEfficiency*> pEff;
  TH1F* demo;

  for(unsigned int filenum=0; filenum<file.size(); filenum++) {
    
    TH1F* nosel = (TH1F*) file[filenum]->Get(cutnames[filenum*2]);
    TH1F* sel = (TH1F*) file[filenum]->Get(cutnames[filenum*2+1]);
    
    nosel = (TH1F*) nosel->Rebin(nbins,"newx",rebin);
    sel = (TH1F*) sel->Rebin(nbins,"newx",rebin);
    if(filenum==0) demo = (TH1F*) sel->Clone();
    
    pEff.push_back(0);
    if(TEfficiency::CheckConsistency((*sel),(*nosel))) {
      pEff[filenum] = new TEfficiency((*sel),(*nosel));
    }
    
    sel->SetLineColor(kRed);
    gStyle->SetOptStat(0);
    pEff[filenum]->SetLineWidth(3);
    pEff[filenum]->SetLineColor(coloropt[filenum]);
  }
  
  std::vector<TH1F*> allhists;
  demo->SetTitle("");
  demo->GetXaxis()->SetTitle(xaxistitle);
  demo->GetYaxis()->SetTitle("RECO eff.");
  demo->SetLineColorAlpha(kWhite,1);
  demo->SetFillColorAlpha(kWhite,1);
  allhists.push_back(demo);
  
  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(), (float []){0.7,0.65,0.9,0.95},false,(float []){0.0,1.1},false);
  auto pad = c1->GetPad(3);
  for(unsigned int filenum=0; filenum<file.size(); filenum++) {
    pEff[filenum]->Draw("same");
  }
  TLatex latex;
  latex.DrawLatex((*legpos),(*(legpos+1)),seltext[0]);
  latex.DrawLatex((*legpos),(*(legpos+1))-0.1,seltext[1]);
  c1->SaveAs("./dirplots/"+((TString)file[0]->GetName()).ReplaceAll(".root","")+"/"+cutnames[0]+"_eff.png");
  
  return -1;
}

int comparesamevariable(std::vector<TFile*> file, std::vector<TString> cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool underflow=false, bool overflow=false, float yrange[]=(float []){0.1,100}, float legPos[]=(float []){0.7,0.75,0.95,1}, bool normalize=true, TString xaxistitle="xaxis") {

  // Check if the size of file vector is same as the cutnames
  if(file.size()!=cutname.size()) {
    cout<<"Error! Mismatching size of vectors"<<endl;
    return -1;
  }

  TString foldername = "";
  foldername += ((TString)file[0]->GetName()).ReplaceAll(".root","").ReplaceAll("hists_","")+"_"+cutname[0]+"_";

  std::vector<TH1F*> allhists;
  for(unsigned int histctr=0; histctr<cutname.size(); histctr++) {
    TString cutwithvar = cutname[histctr]+"_"+var;
    allhists.push_back((TH1F*)file[histctr]->Get(cutwithvar));
    cutname[histctr] += "_"+((TString)file[histctr]->GetName()).ReplaceAll(".root","").ReplaceAll("hists_","");
  }
  
  // Get the title from histogram title
  TString dummytitle = "xaxis";
  TString histtitle = allhists[0]->GetTitle();
  TString xtitle = xaxistitle.CompareTo(dummytitle)?xaxistitle:histtitle;

  if(rebin==-1) rebin = 1;
  xbinlow = xbinlow==-1?(underflow?0:1):xbinlow;
  xbinhigh = xbinhigh==-1?(overflow?allhists[0]->GetNbinsX()+1:allhists[0]->GetNbinsX()):xbinhigh;
  xbinlow = xbinlow/rebin;
  xbinhigh = xbinhigh/rebin;

  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {
    if(rebin!=1) allhists[histctr]->Rebin(rebin);

    // Make changes to sig and bkg to enable good basic plotting
    double err = 0.0;
    if(underflow) {
      allhists[histctr]->SetBinContent(xbinlow,allhists[histctr]->IntegralAndError(0,xbinlow,err));
      allhists[histctr]->SetBinError(xbinlow,err);
    }
    err = 0.0;
    if(overflow) {
      allhists[histctr]->SetBinContent(xbinhigh,allhists[histctr]->IntegralAndError(xbinhigh,allhists[histctr]->GetNbinsX()+1,err));
      allhists[histctr]->SetBinError(xbinhigh,err);
    }

    allhists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allhists[histctr]->GetXaxis()->SetTitle(xtitle);
    if(normalize) allhists[histctr]->GetYaxis()->SetTitle("normalized number of events (a.u.)");
    else allhists[histctr]->GetYaxis()->SetTitle("number of events");

    allhists[histctr]->SetTitle("");

    allhists[histctr]->SetLineWidth(2);
    allhists[histctr]->SetLineColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,normalize);
  c1->SaveAs("./dirplots/"+foldername+"/"+var+".png");

  return -1;
}

int plotter() {

  std::vector<int> coloroptrate{1, kRed+2, kRed-3, kRed-7, kRed-9};

  std::vector<TFile*> file;
  std::vector<TString> name;
  std::vector<TString> legend;

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(metdatafile);
  name.push_back("ec_electron");
  legend.push_back("|#eta|>1.479");
  coloropt.push_back(kRed);
  file.push_back(metdatafile);
  name.push_back("bar_electron");
  legend.push_back("|#eta|<1.479");
  coloropt.push_back(kBlue);

  legendEntries = legend;  
  //comparesamevariable(file, name, "event_trigdec", 40, 55, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.5,0.8,0.75,0.99}, false, "Trigger Decisions");
  //comparesamevariable(file, name, "pt", 50, 100, 1, true, true, true, (float []){0.8,1e4}, (float []){0.6,0.6,0.89,0.95}, false, "photon filter p_{T} / GeV");
  //comparesamevariable(file, name, "eta", -1, -1, 2, false, true, true, (float []){0,700}, (float []){0.11,0.7,0.36,0.99}, false, "photon filter #eta");
  //comparesamevariable(file, name, "phi", -1, -1, 2, false, true, true, (float []){0,350}, (float []){0.25,0.3,0.49,0.55}, false, "photon filter #phi");

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();

  file.push_back(metdatafile);
  name.push_back("met_mid_gt2_bar_el_el_sublead_pt");
  name.push_back("sm12_mid_gt2_bar_el_el_sublead_pt");
  legend.push_back("bla");
  legendEntries = legend;
  vector<double> binspt{2,3,4,7,10,12,14,16,18,20};
  efficiency(file, name, (double []){11,0.7}, binspt.size()-1, &binspt[0], "p_{T} [GeV]");
  name.clear();
  legend.clear();
  name.push_back("met_mid_gt2_bar_el_el_lead_log10d0");
  name.push_back("sm12_mid_gt2_bar_el_el_lead_log10d0");
  legend.push_back("bla");
  legendEntries = legend;
  vector<double> binslog10d0{-3,-2,-1,0,1,2};
  efficiency(file, name, (double []){0,0.7}, binslog10d0.size()-1, &binslog10d0[0], "log_{10}d_{0} [log_{10}cm]");
  
  return -1;
}
