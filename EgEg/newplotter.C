#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 1.09/19; // Data rate - 1.1 Hz, 22 events from parent trigger selection
double sig3cmsf = 1.0/15;
double sig30cmsf = 1.0/28;
double sig1msf = 1.0/19;
double sig3msf = 1.0/7;
TString cutdeets = "Cut details";
TFile* datafile = TFile::Open("hists_Efmrl.root","READ");
TFile* dyfile = TFile::Open("hists_DY.root","READ");
TFile* sig3cmfile = TFile::Open("hists_M200dM20ctau3cm.root","READ");
TFile* sig30cmfile = TFile::Open("hists_M200dM20ctau30cm.root","READ");
TFile* sig1mfile = TFile::Open("hists_M200dM20ctau1m.root","READ");
TFile* sig3mfile = TFile::Open("hists_M200dM20ctau3m.root","READ");

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

int newplotter() {

  std::vector<int> coloroptrate{1, kRed+2, kRed-3, kRed-7, kRed-9};

  std::vector<TFile*> file;
  std::vector<TString> name;
  std::vector<TString> legend;

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig3cmfile);
  file.push_back(sig30cmfile);
  file.push_back(sig1mfile);
  name.push_back("genptgt10geneg_eta");
  name.push_back("genptgt10Anoselusrecomchgenel_eta");
  name.push_back("genptgt10geneg_eta");
  name.push_back("genptgt10Anoselusrecomchgenel_eta");
  name.push_back("genptgt10geneg_eta");
  name.push_back("genptgt10Anoselusrecomchgenel_eta");
  coloropt.push_back(kRed-2);
  coloropt.push_back(kBlue-2);
  coloropt.push_back(kGreen-2);
  legend.push_back("");
  legendEntries = legend;
  vector<double> binseta{-3,-2.9,-2.8,-2.7,-2.6,-2.5,-2.3,-2.2,-1.9,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.9,2.2,2.3,2.5,2.6,2.7,2.8,2.9,3};
  //seltext[0] = "2 gen e, mom Z, p_{T}>10 GeV, |#eta|<1.4 or 1.6|#eta|<2.4";
  seltext[0] = "2 gen e, mom #chi^{#pm}, p_{T}>10 GeV, |#eta|<1.2";
  seltext[1] = "1 unseeded e/#gamma";
  vector<double> legpos{-2,1};
  //efficiency(file, name, &legpos[0], 40, &binseta[0], "#eta");
  name.clear();
  name.push_back("gennoselgeneg_pt");
  name.push_back("gennoselAnoselusrecomchgenel_pt");
  name.push_back("gennoselgeneg_pt");
  name.push_back("gennoselAnoselusrecomchgenel_pt");
  name.push_back("gennoselgeneg_pt");
  name.push_back("gennoselAnoselusrecomchgenel_pt");
  vector<double> binspt{1,2,3,4,5,6,7,8,9,10,13,15,19,23,27,34,42,50,80,150};
  legpos = {80, 0.8};
  //efficiency(file, name, &legpos[0], 19, &binspt[0], "p_{T} / GeV");
  name.clear();
  name.push_back("genptgt10etalt12geneg_log10d0");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_log10d0");
  name.push_back("genptgt10etalt12geneg_log10d0");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_log10d0");
  name.push_back("genptgt10etalt12geneg_log10d0");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_log10d0");
  vector<double> binslog10d0{-5,-3.5,-3,-2.5,-2,-1.5,-1,-0.6,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2,2.2,2.4,3};
  legpos = {-3, 0.4};
  //efficiency(file, name, &legpos[0], 22, &binslog10d0[0], "log_{10}d_{0} / log_{10}cm");
  name.clear();
  name.push_back("genptgt10etalt12geneg_d0");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_d0");
  name.push_back("genptgt10etalt12geneg_d0");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_d0");
  name.push_back("genptgt10etalt12geneg_d0");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_d0");
  vector<double> binsd0{-11,-8,-6,-5,-4,-3,-2,-1,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,1,2,3,4,5,6,7,8,9,10,11};
  legpos = {0,1};
  //efficiency(file, name, &legpos[0], 22, &binsd0[0], "d_{0} / cm");
  name.clear();
  name.push_back("genptgt10etalt12geneg_lxy");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_lxy");
  name.push_back("genptgt10etalt12geneg_lxy");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_lxy");
  name.push_back("genptgt10etalt12geneg_lxy");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_lxy");
  vector<double> binslxy{0,1,2,3,4,5,6,7,8,9,10,15,25,40,60,100};
  legpos = {10,1};
  //efficiency(file, name, &legpos[0], 15, &binslxy[0], "l_{xy} / cm");
  name.clear();
  name.push_back("genptgt10etalt12geneg_log10lxy");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_log10lxy");
  name.push_back("genptgt10etalt12geneg_log10lxy");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_log10lxy");
  name.push_back("genptgt10etalt12geneg_log10lxy");
  name.push_back("genptgt10etalt12Anoselusrecomchgenel_log10lxy");
  vector<double> binslog10lxy{-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.5,2,3,3.5};
  legpos = {10,1};
  //efficiency(file, name, &legpos[0], 14, &binslog10lxy[0], "l_{xy} / cm");

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig3cmfile);
  name.push_back("genptgt10etalt12geneg");
  legend.push_back("signal 3 cm");
  coloropt.push_back(kRed);
  file.push_back(sig3cmfile);
  name.push_back("genptgt10etalt12Anoselusrecomchgenel");
  legend.push_back("3 cm trig. ch");
  coloropt.push_back(kRed-2);
  file.push_back(sig30cmfile);
  name.push_back("genptgt10etalt12geneg");
  legend.push_back("signal 30 cm");
  coloropt.push_back(kBlue);
  file.push_back(sig30cmfile);
  name.push_back("genptgt10etalt12Anoselusrecomchgenel");
  legend.push_back("30 cm trig. mch");
  coloropt.push_back(kBlue-2);
  file.push_back(sig1mfile);
  name.push_back("genptgt10etalt12geneg");
  legend.push_back("signal 1m");
  coloropt.push_back(kGreen);
  file.push_back(sig1mfile);
  name.push_back("genptgt10etalt12Anoselusrecomchgenel");
  legend.push_back("1 m trig. mch");
  coloropt.push_back(kGreen-2);
  legendEntries = legend;  
  //comparesamevariable(file, name, "egmompid", -1, -1, 1, true, true, true, (float []){8e-1,1e7}, (float []){0.5,0.6,0.75,0.89}, false, "electron mother pdg id");
  //comparesamevariable(file, name, "egmult", 5, 15, 1, true, true, true, (float []){8e-1,1e7}, (float []){0.5,0.6,0.75,0.89}, false, "electron multiplicity");
  //comparesamevariable(file, name, "pt", 50, 100, 1, true, true, true, (float []){0.8,1e5}, (float []){0.6,0.6,0.89,0.95}, false, "e/#gamma p_{T} / GeV");
  //comparesamevariable(file, name, "eta", -1, -1, 2, false, true, true, (float []){0,600}, (float []){0.11,0.7,0.36,0.99}, false, "e/#gamma #eta");
  //comparesamevariable(file, name, "phi", -1, -1, 2, false, true, true, (float []){0,1500}, (float []){0.25,0.3,0.49,0.55}, false, "e/#gamma #phi");
  //comparesamevariable(file, name, "log10d0", 50, -1, 20, true, true, true, (float []){0.9,1e6}, (float []){0.7,0.6,0.95,0.89}, false, "e/#gamma log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(file, name, "d0", 9000, 11000, 20, true, true, true, (float []){0.9,1e6}, (float []){0.7,0.6,0.95,0.89}, false, "e/#gamma d_{0} / cm");
  //comparesamevariable(file, name, "log10lxy", -1, -1, 10, true, true, true, (float []){0.9,1e6}, (float []){0.6,0.6,0.85,0.96}, false, "e/#gamma log_{10}l_{xy} / log_{10}cm");
  //comparesamevariable(file, name, "lxy", 1000, 11500, 400, true, true, true, (float []){0.9,1e6}, (float []){0.5,0.6,0.75,0.96}, false, "e/#gamma l_{xy} / cm");
  
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig3cmfile);
  name.push_back("genbasicselptgt15Abasicselusgenmchrecoebus");
  legend.push_back("3 cm gen mch");
  coloropt.push_back(kRed+2);
  file.push_back(sig3cmfile);
  name.push_back("genbasicselptgt15Abasicselusgenmchrecoebus_pxlmch22");
  legend.push_back("pxlmch22 gen mch");
  coloropt.push_back(kBlue+2);
  file.push_back(sig3cmfile);
  name.push_back("basicselusrecoebus");
  legend.push_back("signal 3 cm");
  coloropt.push_back(kRed);
  file.push_back(sig3cmfile);
  name.push_back("basicselusrecoebus_pxlmch22");
  legend.push_back("pxlmch22");
  coloropt.push_back(kBlue);
  legendEntries = legend;  
  //comparesamevariable(file, name, "egpixelmchvar_s2", 45, 250, 10, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma loose pixel match var s2");
  //comparesamevariable(file, name, "leadegooeseedoop", -1, -1, 10, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma_{1} track 1/E-1/p / GeV^{-1}");
  //comparesamevariable(file, name, "subleadegooeseedoop", -1, -1, 10, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma_{2} track 1/E-1/p / GeV^{-1}");
  //comparesamevariable(file, name, "leadegchi2", -1, -1, 1, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma track #chi^2");
  //comparesamevariable(file, name, "leadegdeta", -1, -1, 1, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma track #Delta#eta_{in}");
  //comparesamevariable(file, name, "leadegdetaseed", -1, -1, 1, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma track #Delta#eta_{in}^{seed}");
  //comparesamevariable(file, name, "leadegdphi", -1, -1, 1, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma track #Delta#phi_{in}");
  //comparesamevariable(file, name, "leadegmhits", -1, -1, 1, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma track missing hits");
  //comparesamevariable(file, name, "leadegnlayerit", -1, -1, 1, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma track layer num");
  //comparesamevariable(file, name, "leadegooesclsoop", 50, 200, 5, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma_{1} 1/E_{sclus}-1/p / GeV^{-1}");
  //comparesamevariable(file, name, "leadegvalhits", -1, -1, 1, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma valid hits");
  //comparesamevariable(file, name, "leadegtrkiso", 10, 200, 5, true, true, true, (float []){0.5,1e5}, (float []){0.5,0.74,0.75,0.99}, false, "e/#gamma_{1} trk. iso. / GeV");
  
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(dyfile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("dy, veto id, Z");
  coloropt.push_back(kBlue);
  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("data, veto id, Z");
  coloropt.push_back(kBlack);
  file.push_back(sig3cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("3 cm gen mch");
  coloropt.push_back(kRed+3);
  file.push_back(sig30cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("30 cm gen mch");
  coloropt.push_back(kRed-2);
  file.push_back(sig1mfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("1 m gen mch");
  coloropt.push_back(kRed-6);
  file.push_back(sig3mfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("3 m gen mch");
  coloropt.push_back(kRed-9);
  /*
  coloropt.push_back(kBlue);
  file.push_back(datafile);
  name.push_back("basicselusrecoebus");
  legend.push_back("data, no cut");
  coloropt.push_back(kGray);
  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("data, veto id, Z");
  */
  legendEntries = legend;
  comparesamevariable(file, name, "egseedclustime", 7500, 13500, 250, true, true, true, (float []){1e-3,1}, (float []){0.65,0.6,0.8,0.95}, true, "e/#gamma ecal seed cluster time / ns");

  return -1;
}
