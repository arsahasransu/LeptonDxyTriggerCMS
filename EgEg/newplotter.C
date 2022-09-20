#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 8.3/771; // Data rate - 8.3 Hz, 771 events from prescaled dig33_caloidl selection
// Data rate equivalent to 6.8 Hz, 76 events from unprescaled dieg70 selection
double sig3cmsf = 1.0/28;
double sig30cmsf = 1.0/26;
double sig1msf = 1.0/11;
double sig3msf = 1.0/3;
TString cutdeets = "Cut details";
TFile* datafile = TFile::Open("hists_Efmrl.root","READ");
TFile* dyfile = TFile::Open("hists_DY.root","READ");
TFile* sig3cmfile = TFile::Open("hists_M200dM20ctau3cm.root","READ");
TFile* sig30cmfile = TFile::Open("hists_M200dM20ctau30cm.root","READ");
TFile* sig1mfile = TFile::Open("hists_M200dM20ctau1m.root","READ");
TFile* sig3mfile = TFile::Open("hists_M200dM20ctau3m.root","READ");
TFile* sigfile = TFile::Open("hists_M200dM20.root","READ");

TString seltext[4] = {"line1", "line2", "line3", "line4"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};

int makeratehist(TString cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool overflow=false, float legPos[]=(float []){0.7,0.75,0.95,1}, float seltextpos[]=(float []){0.1,1}, float yrange[]=(float []){0.1,1}, float drawsignalline=-1.0, TString xaxistitle="p_{T} / GeV") {

  std::vector<TString> legNam;
  legNam.push_back("2018 data");
  legNam.push_back("c#tau = 3 cm");
  legNam.push_back("c#tau = 30 cm");
  legNam.push_back("c#tau = 1 m");
  legNam.push_back("c#tau = 3 m");

  TString histname = cutname+"_"+var;
  TH1F* datahist;
  datahist = (TH1F*) datafile->Get(histname);
  datahist->Scale(datasf);
  TH1F* sig3cmhist = (TH1F*) sig3cmfile->Get(histname);
  sig3cmhist->Scale(sig3cmsf);
  TH1F* sig30cmhist = (TH1F*) sig30cmfile->Get(histname);
  sig30cmhist->Scale(sig30cmsf);
  TH1F* sig1mhist = (TH1F*) sig1mfile->Get(histname);
  sig1mhist->Scale(sig1msf);
  TH1F* sig3mhist = (TH1F*) sig3mfile->Get(histname);
  sig3mhist->Scale(sig3msf);

  std::vector<TH1F*> allhists;
  std::vector<TH1F*> allratehists;
  allhists.push_back(datahist);
  allhists.push_back(sig3cmhist);
  allhists.push_back(sig30cmhist);
  allhists.push_back(sig1mhist);
  allhists.push_back(sig3mhist);

  double axisscale = 0.0;
  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {

    allratehists.push_back((TH1F*) allhists[histctr]->Clone());
    for(unsigned int bincnt=0; bincnt<allhists[histctr]->GetNbinsX(); bincnt++) {
      double err = 0;
      allratehists[histctr]->SetBinContent(bincnt, allhists[histctr]->IntegralAndError(bincnt,allhists[histctr]->GetNbinsX()+1,err));
      allratehists[histctr]->SetBinError(bincnt, err);
    }
    allratehists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allratehists[histctr]->GetXaxis()->SetTitle(xaxistitle);
    allratehists[histctr]->GetYaxis()->SetTitle("rate / Hz");

    allratehists[histctr]->SetTitle("");

    allratehists[histctr]->SetLineWidth(2);
    allratehists[histctr]->SetLineColor(coloropt[histctr]);
    if(histctr==1) {
      axisscale = allratehists[0]->GetBinContent(0)*0.8/allratehists[1]->GetBinContent(0);
      allratehists[histctr]->Scale(axisscale);
    }
    if(histctr>1) {
      allratehists[histctr]->Scale(axisscale);
    }
  }

  cout<<allratehists[1]->GetMaximum()<<endl;
  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter_rate(allratehists, legNam, allratehists[0]->GetXaxis()->GetTitle(),allratehists[0]->GetYaxis()->GetTitle(),legPos,yrange,logY,false);

  TPad* pad = (TPad*) c1->FindObject("pad3");
  pad->cd();
  TGaxis *axis;
  axis = new TGaxis(allratehists[0]->GetBinLowEdge(xbinhigh+1),yrange[0],allratehists[0]->GetBinLowEdge(xbinhigh+1),yrange[1],yrange[0]/axisscale,yrange[1]/axisscale,510,"-L");
  axis->SetLineColor(coloropt[1]);
  axis->SetLabelColor(coloropt[1]);
  axis->SetLabelFont(132);
  axis->SetLabelSize(0.04);
  axis->SetLabelOffset(-0.035);
  axis->Draw();

  if(drawsignalline!=-1.0) {
    drawsignalline = allratehists[1]->GetMaximum();
    TLine *signalline = new TLine(allratehists[0]->GetBinLowEdge(xbinlow),drawsignalline,allratehists[0]->GetBinLowEdge(xbinhigh+1),drawsignalline);
    signalline->SetLineWidth(2);
    signalline->SetLineColor(coloropt[1]);
    signalline->SetLineStyle(9);
    signalline->Draw();
  }

  TLatex sel;
  sel.SetTextFont(132);
  sel.SetTextSize(0.065);
  sel.DrawLatex(seltextpos[0], seltextpos[1]+0.3*(yrange[1]-yrange[0]), seltext[0]);
  sel.DrawLatex(seltextpos[0], seltextpos[1]+0.2*(yrange[1]-yrange[0]), seltext[1]);
  sel.DrawLatex(seltextpos[0], seltextpos[1]+0.1*(yrange[1]-yrange[0]), seltext[2]);
  sel.DrawLatex(seltextpos[0], seltextpos[1], seltext[3]);
  
  c1->SaveAs("./dirplots/"+cutname+"/"+cutname+"_"+var+"_ratehist.png");
  
  return -1;
}

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

  // Plots for gen
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig3cmfile);
  name.push_back("genbasicselbargeneg");
  legend.push_back("3 cm ");
  coloropt.push_back(2);
  file.push_back(sig30cmfile);
  name.push_back("genbasicselbargeneg");
  legend.push_back("30 cm ");
  coloropt.push_back(210);
  file.push_back(sig1mfile);
  name.push_back("genbasicselbargeneg");
  legend.push_back("1 m ");
  coloropt.push_back(218);
  file.push_back(sig3mfile);
  name.push_back("genbasicselbargeneg");
  legend.push_back("3 m ");
  coloropt.push_back(4);
  legendEntries = legend;
  //comparesamevariable(file, name, "egmult", 5, 10, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "electron multiplicity");
  //comparesamevariable(file, name, "pt", 50, 120, 1, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "electron p_{T} / GeV");
  //comparesamevariable(file, name, "eta", -1, -1, 1, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "electron #eta");
  //comparesamevariable(file, name, "phi", -1, -1, 1, true, true, true, (float []){8e-1,1e3}, (float []){0.3,0.3,0.65,0.75}, false, "electron #phi");
  //comparesamevariable(file, name, "egmompid", -1, -1, 1, true, true, true, (float []){8e-1,1e6}, (float []){0.7,0.6,0.85,0.95}, false, "electron mother PDG");
  //comparesamevariable(file, name, "deltaetamom", -1, -1, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "#Delta#eta(mom, gen) / GeV");
  //comparesamevariable(file, name, "deltaphimom", -1, -1, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "#Delta#phi(mom, gen) / GeV");
  //comparesamevariable(file, name, "deltaRmom", -1, -1, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "#DeltaR(mom, gen) / GeV");
  //comparesamevariable(file, name, "d0", 9000, 11000, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "electron d_{0} / cm");
  //comparesamevariable(file, name, "lxy", 1000, 11500, 100, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "electron l_{xy} / cm");
  //comparesamevariable(file, name, "log10d0", 50, -1, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "electron log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(file, name, "log10lxy", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "electron log_{10}l_{xy} / log_{10}cm");
  //comparesamevariable(file, name, "vx", 20000, 30000, 40, true, true, true, (float []){8e-1,1e4}, (float []){0.7,0.6,0.9,0.95}, false, "electron v_{x} / cm");
  //comparesamevariable(file, name, "vy", 20000, 30000, 40, true, true, true, (float []){8e-1,1e4}, (float []){0.7,0.6,0.9,0.95}, false, "electron v_{y} / cm");
  //comparesamevariable(file, name, "vz", 20000, 30000, 80, true, true, true, (float []){8e-1,1e4}, (float []){0.7,0.6,0.9,0.95}, false, "electron v_{z} / cm");
  //comparesamevariable(file, name, "pvx", 24990, 25010, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "mother v_{x} / cm");
  //comparesamevariable(file, name, "pvy", 24990, 25010, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "mother v_{y} / cm");
  //comparesamevariable(file, name, "pvz", 24000, 26000, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "mother v_{z} / cm");
  //comparesamevariable(file, name, "prompteta", -1, -1, 1, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "prompt equivalent electron #eta");
  //comparesamevariable(file, name, "promptphi", -1, -1, 1, true, true, true, (float []){8e-1,1e3}, (float []){0.3,0.3,0.65,0.75}, false, "prompt equivalent electron #phi");

  // Plots for comparing gen matching
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig3cmfile);
  name.push_back("gennoselAnoselusgeneltrigebus");
  legend.push_back("3 cm ");
  coloropt.push_back(2);
  file.push_back(sig3cmfile);
  name.push_back("gennoselAnoselusgenmchgeneltrigebus");
  legend.push_back("3 cm ");
  coloropt.push_back(50);
  file.push_back(sig30cmfile);
  name.push_back("gennoselAnoselusgeneltrigebus");
  legend.push_back("30 cm ");
  coloropt.push_back(210);
  file.push_back(sig30cmfile);
  name.push_back("gennoselAnoselusgenmchgeneltrigebus");
  legend.push_back("30 cm ");
  coloropt.push_back(212);
  file.push_back(sig1mfile);
  name.push_back("gennoselAnoselusgeneltrigebus");
  legend.push_back("1 m ");
  coloropt.push_back(218);
  file.push_back(sig1mfile);
  name.push_back("gennoselAnoselusgenmchgeneltrigebus");
  legend.push_back("1 m ");
  coloropt.push_back(219);
  file.push_back(sig3mfile);
  name.push_back("gennoselAnoselusgeneltrigebus");
  legend.push_back("3 m ");
  coloropt.push_back(4);
  file.push_back(sig3mfile);
  name.push_back("gennoselAnoselusgenmchgeneltrigebus");
  legend.push_back("3 m ");
  coloropt.push_back(38);
  legendEntries = legend;
  //comparesamevariable(file, name, "dE", -1, -1, 200, true, true, true, (float []){8e-1,1e4}, (float []){-1,0.55,0.85,0.99}, false, "#DeltaE(gen, trig. us) / GeV");
  //comparesamevariable(file, name, "dPt", -1, -1, 200, true, true, true, (float []){8e-1,1e4}, (float []){-1,0.55,0.85,0.99}, false, "#Deltap_{T}(gen, trig. us) / GeV");
  //comparesamevariable(file, name, "dEta", -1, -1, 40, true, true, true, (float []){8e-1,1e4}, (float []){-1,0.55,0.85,0.99}, false, "#Delta#eta(gen, trig. us)");
  //comparesamevariable(file, name, "qdPhi", -1, -1, 40, true, true, true, (float []){8e-1,1e7}, (float []){-1,0.55,0.85,0.99}, false, "Q_{gen}#Delta#phi(gen, trig. us)");

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig3cmfile);
  name.push_back("genbasicselbarAnoselusgeneltrigebus");
  legend.push_back("3 cm ");
  coloropt.push_back(2);
  file.push_back(sig3cmfile);
  name.push_back("genbasicselbarAnoselusgenmchgeneltrigebus");
  legend.push_back("3 cm ");
  coloropt.push_back(50);
  file.push_back(sig30cmfile);
  name.push_back("genbasicselbarAnoselusgeneltrigebus");
  legend.push_back("30 cm ");
  coloropt.push_back(210);
  file.push_back(sig30cmfile);
  name.push_back("genbasicselbarAnoselusgenmchgeneltrigebus");
  legend.push_back("30 cm ");
  coloropt.push_back(212);
  file.push_back(sig1mfile);
  name.push_back("genbasicselbarAnoselusgeneltrigebus");
  legend.push_back("1 m ");
  coloropt.push_back(218);
  file.push_back(sig1mfile);
  name.push_back("genbasicselbarAnoselusgenmchgeneltrigebus");
  legend.push_back("1 m ");
  coloropt.push_back(219);
  file.push_back(sig3mfile);
  name.push_back("genbasicselbarAnoselusgeneltrigebus");
  legend.push_back("3 m ");
  coloropt.push_back(4);
  file.push_back(sig3mfile);
  name.push_back("genbasicselbarAnoselusgenmchgeneltrigebus");
  legend.push_back("3 m ");
  coloropt.push_back(38);
  legendEntries = legend;
  //comparesamevariable(file, name, "dE", -1, -1, 200, true, true, true, (float []){8e-1,1e4}, (float []){-1,0.55,0.85,0.99}, false, "#DeltaE(gen, trig. us) / GeV");
  //comparesamevariable(file, name, "dPt", -1, -1, 200, true, true, true, (float []){8e-1,1e4}, (float []){-1,0.55,0.85,0.99}, false, "#Deltap_{T}(gen, trig. us) / GeV");
  //comparesamevariable(file, name, "dEta", 3000, 5000, 10, true, true, true, (float []){8e-1,1e4}, (float []){-1,0.55,0.85,0.99}, false, "#Delta#eta(gen, trig. us)");
  //comparesamevariable(file, name, "qdPhi", 2000, 5000, 10, true, true, true, (float []){8e-1,1e7}, (float []){-1,0.55,0.85,0.99}, false, "Q_{gen}#Delta#phi(gen, trig. us)");
  //comparesamevariable(file, name, "dPromptEta", 3000, 5000, 10, true, true, true, (float []){8e-1,1e4}, (float []){-1,0.55,0.85,0.99}, false, "#Delta#eta(prompt eq. gen, trig. us)");
  //comparesamevariable(file, name, "qdPromptPhi", 2000, 5000, 10, true, true, true, (float []){8e-1,1e7}, (float []){-1,0.55,0.85,0.99}, false, "Q_{gen}#Delta#phi(prompt eq. gen, trig. us)");
  
  // Plots for comparing gen matching
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbarAcuttimedelayonlyusrecomchgenel");
  legend.push_back("3m HLT_DiPhoton10Time1p4ns");
  coloropt.push_back(4);
  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbarAdieg33caloidlusrecomchgenel");
  legend.push_back("3m HLT_DoublePhoton33");
  coloropt.push_back(38);
  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAcuttimedelayonlyusrecomchgenel");
  legend.push_back("1m HLT_DiPhoton10Time1p4ns");
  coloropt.push_back(218);
  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAdieg33caloidlusrecomchgenel");
  legend.push_back("1m HLT_DoublePhoton33");
  coloropt.push_back(219);
  legendEntries = legend;
  //comparesamevariable(file, name, "pt", 50, 120, 1, true, true, true, (float []){8e-1,2e2}, (float []){0.4,0.6,0.55,0.95}, false, "electron p_{T} / GeV");
  //comparesamevariable(file, name, "d0", 9000, 11000, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "electron d_{0} / cm");
  //comparesamevariable(file, name, "lxy", 1000, 11500, 100, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "electron l_{xy} / cm");
  //comparesamevariable(file, name, "log10d0", 250, -1, 20, true, true, true, (float []){8e-1,2e2}, (float []){-1,0.6,0.75,0.95}, false, "electron log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(file, name, "log10lxy", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "electron log_{10}l_{xy} / log_{10}cm");

  // Plots for comparing gen matching
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selbarAcuttimedelaysminusrecomchgenel");
  legend.push_back("30cm HLT_DiPhoton10Timelt2ns_sminlt0p1");
  coloropt.push_back(210);
  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selbarAdieg33caloidlusrecomchgenel");
  legend.push_back("30cm HLT_DoublePhoton33");
  coloropt.push_back(212);
  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAcuttimedelaysminusrecomchgenel");
  legend.push_back("3cm HLT_DiPhoton10Timelt2ns_sminlt0p1");
  coloropt.push_back(2);
  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAdieg33caloidlusrecomchgenel");
  legend.push_back("3cm HLT_DoublePhoton33");
  coloropt.push_back(50);
  legendEntries = legend;
  //comparesamevariable(file, name, "pt", 50, 120, 1, true, true, true, (float []){8e-1,2e2}, (float []){-1,0.6,-1,0.95}, false, "electron p_{T} / GeV");
  //comparesamevariable(file, name, "d0", 9000, 11000, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "electron d_{0} / cm");
  //comparesamevariable(file, name, "lxy", 1000, 11500, 100, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "electron l_{xy} / cm");
  //comparesamevariable(file, name, "log10d0", 50, -1, 20, true, true, true, (float []){8e-1,2e2}, (float []){-1,0.6,-1,0.95}, false, "electron log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(file, name, "log10lxy", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "electron log_{10}l_{xy} / log_{10}cm");

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig3cmfile);
  file.push_back(sig30cmfile);
  file.push_back(sig1mfile);
  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbargeneg_eta");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_eta");
  name.push_back("genbasicptgt10selbargeneg_eta");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_eta");
  name.push_back("genbasicptgt10selbargeneg_eta");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_eta");
  name.push_back("genbasicptgt10selbargeneg_eta");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_eta");
  coloropt.push_back(2);
  coloropt.push_back(210);
  coloropt.push_back(218);
  coloropt.push_back(4);
  legend.push_back("");
  legendEntries = legend;
  vector<double> binseta{-3,-2.9,-2.8,-2.7,-2.6,-2.5,-2.3,-2.2,-1.9,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.9,2.2,2.3,2.5,2.6,2.7,2.8,2.9,3};
  //seltext[0] = "2 gen e, mom Z, p_{T}>10 GeV, |#eta|<1.4 or 1.6|#eta|<2.4";
  seltext[0] = "2 gen e, mom #chi^{#pm}, p_{T}>10 GeV, |#eta|<1.479, ecal fiducial";
  seltext[1] = "1 unseeded e/#gamma";
  vector<double> legpos{-2,0.2};
  //efficiency(file, name, &legpos[0], 40, &binseta[0], "#eta");
  name.clear();
  name.push_back("genbasicselbargeneg_pt");
  name.push_back("genbasicselbarAnoselusrecomchgenel_pt");
  name.push_back("genbasicselbargeneg_pt");
  name.push_back("genbasicselbarAnoselusrecomchgenel_pt");
  name.push_back("genbasicselbargeneg_pt");
  name.push_back("genbasicselbarAnoselusrecomchgenel_pt");
  name.push_back("genbasicselbargeneg_pt");
  name.push_back("genbasicselbarAnoselusrecomchgenel_pt");
  vector<double> binspt{1,2,3,4,5,6,7,8,9,10,13,15,19,23,27,34,42,50,80,150};
  legpos = {40, 0.5};
  //efficiency(file, name, &legpos[0], 19, &binspt[0], "p_{T} / GeV");
  name.clear();
  name.push_back("genbasicptgt10selbargeneg_log10d0");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_log10d0");
  name.push_back("genbasicptgt10selbargeneg_log10d0");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_log10d0");
  name.push_back("genbasicptgt10selbargeneg_log10d0");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_log10d0");
  name.push_back("genbasicptgt10selbargeneg_log10d0");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_log10d0");
  vector<double> binslog10d0{-5,-3.5,-3,-2.5,-2,-1.5,-1,-0.6,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2,2.2,2.4,3};
  legpos = {-3, 0.4};
  //efficiency(file, name, &legpos[0], 22, &binslog10d0[0], "log_{10}d_{0} / log_{10}cm");
  name.clear();
  name.push_back("genbasicptgt10selbargeneg_d0");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_d0");
  name.push_back("genbasicptgt10selbargeneg_d0");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_d0");
  name.push_back("genbasicptgt10selbargeneg_d0");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_d0");
  name.push_back("genbasicptgt10selbargeneg_d0");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_d0");
  vector<double> binsd0{-11,-8,-6,-5,-4,-3,-2,-1,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,1,2,3,4,5,6,7,8,9,10,11};
  legpos = {0,1};
  //efficiency(file, name, &legpos[0], 22, &binsd0[0], "d_{0} / cm");
  name.clear();
  name.push_back("genbasicptgt10selbargeneg_lxy");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_lxy");
  name.push_back("genbasicptgt10selbargeneg_lxy");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_lxy");
  name.push_back("genbasicptgt10selbargeneg_lxy");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_lxy");
  name.push_back("genbasicptgt10selbargeneg_lxy");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_lxy");
  vector<double> binslxy{0,1,2,3,4,5,6,7,8,9,10,15,25,40,60,100};
  legpos = {10,1};
  //efficiency(file, name, &legpos[0], 15, &binslxy[0], "l_{xy} / cm");
  name.clear();
  name.push_back("genbasicptgt10selbargeneg_log10lxy");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_log10lxy");
  name.push_back("genbasicptgt10selbargeneg_log10lxy");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_log10lxy");
  name.push_back("genbasicptgt10selbargeneg_log10lxy");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_log10lxy");
  name.push_back("genbasicptgt10selbargeneg_log10lxy");
  name.push_back("genbasicptgt10selbarAnoselusrecomchgenel_log10lxy");
  vector<double> binslog10lxy{-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.82,1.84,1.86,1.88,1.9,1.92,1.94,1.96,1.98,2,2.02,2.04,2.06,2.08,2.1,2.12,2.14,2.16,2.18,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.5};
  legpos = {10,1};
  //efficiency(file, name, &legpos[0], 47, &binslog10lxy[0], "log_{10} l_{xy} / log_{10}cm");

  // Study the unseeded egamma object properties
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("data, Z wind, veto id");
  coloropt.push_back(kBlack);
  file.push_back(datafile);
  name.push_back("noselusrecoebus");
  legend.push_back("data all");
  coloropt.push_back(16);
  file.push_back(dyfile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("DY");
  coloropt.push_back(4);
  file.push_back(sig3cmfile);
  name.push_back("noselusrecoebus");
  legend.push_back("3cm all");
  coloropt.push_back(50);
  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("3cm gen mch");
  coloropt.push_back(2);
  file.push_back(sig30cmfile);
  name.push_back("noselusrecoebus");
  legend.push_back("30cm all");
  coloropt.push_back(212);
  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("30cm gen mch");
  coloropt.push_back(210);
  file.push_back(sig1mfile);
  name.push_back("noselusrecoebus");
  legend.push_back("1m all");
  coloropt.push_back(219);
  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("1m gen mch");
  coloropt.push_back(218);
  file.push_back(sig3mfile);
  name.push_back("noselusrecoebus");
  legend.push_back("3m all");
  coloropt.push_back(224);
  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("3m gen mch");
  coloropt.push_back(222);
  legendEntries = legend;    
  //comparesamevariable(file, name, "egseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma ecal seed cluster time / ns");
  //comparesamevariable(file, name, "leadegseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{1} ecal seed cluster time / ns");
  //comparesamevariable(file, name, "leadegclustershape", 5, 300, 10, true, true, true, (float []){1e-5,1}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{1} #sigmai#etai#eta_{5x5} (noise clnd.)");
  //comparesamevariable(file, name, "leadegsmin", -1, -1, 10, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} s_{min}");
  //comparesamevariable(file, name, "leadegsmaj", -1, -1, 10, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} s_{maj}");
  //comparesamevariable(file, name, "subleadegseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{2} ecal seed cluster time / ns");
  //comparesamevariable(file, name, "subleadegclustershape", 5, 300, 10, true, true, true, (float []){1e-5,1}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{2} #sigmai#etai#eta_{5x5} (noise clnd.)");
  //comparesamevariable(file, name, "subleadegsmin", -1, -1, 10, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} s_{min}");
  //comparesamevariable(file, name, "subleadegsmaj", -1, -1, 10, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} s_{maj}");

  // Plotting for rate estimation
  coloropt.clear();
  coloropt.push_back(kBlack);
  coloropt.push_back(kRed);
  coloropt.push_back(kRed-7);
  coloropt.push_back(kRed-5);
  coloropt.push_back(kRed-2);

  seltext[0] = "1 L1 seeded e/#gamma: p_{T}>33 GeV, CaloIdL";
  seltext[1] = "2 unseeded e/#gamma: p_{T}>33 GeV, CaloIdL";  
  //makeratehist("dieg33caloidlidusrecoegus", "subleadegpt", 50, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,11}, 0.875);
  //makeratehist("dieg33caloidlidusrecoegus", "leadegpt", 50, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,11}, 0.875);

  seltext[0] = "2 e/#gamma: p_{T}>15 GeV, |#eta|<2.5";
  seltext[1] = "ECAL Tight Id: ";
  seltext[2] = "#sigmai#etai#eta < 0.0104(0.0353)";
  seltext[3] = "H/E < 0.026(0.0188)";
  //makeratehist("basicselusrecoegus", "leadegpt", 50, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,1000}, (float []){0,1e4}, 0.875);
  //makeratehist("basicselusrecoegus", "subleadegpt", 50, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,1000}, (float []){0,1e4}, 0.875);
  //makeratehist("seleletightcaloidusrecoegus", "subleadegpt", 50, 110, 1, false, false, (float []){0.65,0.65,0.85,0.975}, (float []){28,6}, (float []){0,28}, 0.875, "e/#gamma p_{T} [GeV]");
  
  // Study the unseeded egamma object properties
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("data, Z wind, veto id");
  coloropt.push_back(kBlack);
  file.push_back(datafile);
  name.push_back("basicselusrecoebus");
  legend.push_back("data all");
  coloropt.push_back(16);
  file.push_back(dyfile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("DY");
  coloropt.push_back(4);
  file.push_back(sig3cmfile);
  name.push_back("basicselusrecoebus");
  legend.push_back("3cm all");
  coloropt.push_back(50);
  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAbasicselusgenmchrecoebus");
  legend.push_back("3cm gen mch");
  coloropt.push_back(2);
  file.push_back(sig30cmfile);
  name.push_back("basicselusrecoebus");
  legend.push_back("30cm all");
  coloropt.push_back(212);
  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selbarAbasicselusgenmchrecoebus");
  legend.push_back("30cm gen mch");
  coloropt.push_back(210);
  file.push_back(sig1mfile);
  name.push_back("basicselusrecoebus");
  legend.push_back("1m all");
  coloropt.push_back(219);
  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAbasicselusgenmchrecoebus");
  legend.push_back("1m gen mch");
  coloropt.push_back(218);
  file.push_back(sig3mfile);
  name.push_back("basicselusrecoebus");
  legend.push_back("3m all");
  coloropt.push_back(224);
  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbarAbasicselusgenmchrecoebus");
  legend.push_back("3m gen mch");
  coloropt.push_back(222);
  legendEntries = legend;    
  //comparesamevariable(file, name, "egseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma ecal seed cluster time / ns");
  //comparesamevariable(file, name, "leadegseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{1} ecal seed cluster time / ns");
  //comparesamevariable(file, name, "leadegin5x5noiseclnd", 5, 300, 10, true, true, true, (float []){1e-3,1}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{1} #sigmai#etai#eta_{5x5} (noise clnd.)");
  //comparesamevariable(file, name, "leadegsmin", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} s_{min}");
  //comparesamevariable(file, name, "leadegsmaj", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} s_{maj}");
  //comparesamevariable(file, name, "leadegecalpfclustisoovere", -1, -1, 1, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} ecal iso./E");
  //comparesamevariable(file, name, "leadeghovereoversupcluse", -1, 150, 4, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} H/E");
  //comparesamevariable(file, name, "subleadegseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{2} ecal seed cluster time / ns");
  //comparesamevariable(file, name, "subleadegin5x5noiseclnd", 5, 300, 10, true, true, true, (float []){1e-3,1}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{2} #sigmai#etai#eta_{5x5} (noise clnd.)");
  //comparesamevariable(file, name, "subleadegsmin", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} s_{min}");
  //comparesamevariable(file, name, "subleadegsmaj", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} s_{maj}");
  //comparesamevariable(file, name, "subleadegecalpfclustisoovere", -1, -1, 1, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} ecal iso./E");
  //comparesamevariable(file, name, "subleadeghovereoversupcluse", -1, 150, 4, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} H/E");
  
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoeeus");
  legend.push_back("data, Z wind, veto id");
  coloropt.push_back(kBlack);
  file.push_back(datafile);
  name.push_back("basicselusrecoeeus");
  legend.push_back("data all");
  coloropt.push_back(16);
  file.push_back(dyfile);
  name.push_back("selelevetozwindidusrecoeeus");
  legend.push_back("DY");
  coloropt.push_back(4);
  file.push_back(sig3cmfile);
  name.push_back("basicselusrecoeeus");
  legend.push_back("3cm all");
  coloropt.push_back(50);
  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selecAbasicselusgenmchrecoeeus");
  legend.push_back("3cm gen mch");
  coloropt.push_back(2);
  file.push_back(sig30cmfile);
  name.push_back("basicselusrecoeeus");
  legend.push_back("30cm all");
  coloropt.push_back(212);
  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selecAbasicselusgenmchrecoeeus");
  legend.push_back("30cm gen mch");
  coloropt.push_back(210);
  file.push_back(sig1mfile);
  name.push_back("basicselusrecoeeus");
  legend.push_back("1m all");
  coloropt.push_back(219);
  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selecAbasicselusgenmchrecoeeus");
  legend.push_back("1m gen mch");
  coloropt.push_back(218);
  file.push_back(sig3mfile);
  name.push_back("basicselusrecoeeus");
  legend.push_back("3m all");
  coloropt.push_back(224);
  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selecAbasicselusgenmchrecoeeus");
  legend.push_back("3m gen mch");
  coloropt.push_back(222);
  legendEntries = legend;    
  //comparesamevariable(file, name, "egseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma ecal seed cluster time / ns");
  //comparesamevariable(file, name, "leadegseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{1} ecal seed cluster time / ns");
  //comparesamevariable(file, name, "leadegin5x5noiseclnd", 5, 600, 20, true, true, true, (float []){1e-3,1}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{1} #sigmai#etai#eta_{5x5} (noise clnd.)");
  //comparesamevariable(file, name, "leadegsmin", -1, -1, 40, true, true, true, (float []){5e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} s_{min}");
  //comparesamevariable(file, name, "leadegsmaj", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} s_{maj}");
  //comparesamevariable(file, name, "leadegecalpfclustisoovere", -1, -1, 1, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} ecal iso./E");
  //comparesamevariable(file, name, "leadeghovereoversupcluse", -1, 150, 4, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} H/E");
  //comparesamevariable(file, name, "subleadegseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{2} ecal seed cluster time / ns");
  //comparesamevariable(file, name, "subleadegin5x5noiseclnd", 5, 600, 20, true, true, true, (float []){1e-3,1}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{2} #sigmai#etai#eta_{5x5} (noise clnd.)");
  //comparesamevariable(file, name, "subleadegsmin", -1, -1, 40, true, true, true, (float []){5e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} s_{min}");
  //comparesamevariable(file, name, "subleadegsmaj", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} s_{maj}");
  //comparesamevariable(file, name, "subleadegecalpfclustisoovere", -1, -1, 1, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} ecal iso./E");
  //comparesamevariable(file, name, "subleadeghovereoversupcluse", -1, 150, 4, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} H/E");

  // Plotting for rate estimation
  coloropt.clear();
  coloropt.push_back(kBlack);
  coloropt.push_back(2);
  coloropt.push_back(210);
  coloropt.push_back(218);
  coloropt.push_back(222);

  seltext[0] = "2 unseeded e/#gamma: p_{T}>10 GeV, |#eta|<2.5";
  seltext[1] = "#sigmai#etai#eta<0.016(0.04), H/E<0.2";  
  //makeratehist("cut1usrecoegus", "subleadegpt", 50, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,100}, (float []){0,1000}, 0.875);
  //makeratehist("cut1usrecoegus", "leadegpt", 50, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,1000}, 0.875);
  //makeratehist("cut1usrecoebus", "leadegseedclustime", 10000, 17000, 250, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,100}, 0.875);

  // Study the unseeded egamma object properties
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("data, Z wind, veto id");
  coloropt.push_back(kBlack);
  file.push_back(datafile);
  name.push_back("cut1usrecoebus");
  legend.push_back("data all");
  coloropt.push_back(16);
  file.push_back(dyfile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("DY");
  coloropt.push_back(4);
  file.push_back(sig3cmfile);
  name.push_back("cut1usrecoebus");
  legend.push_back("3cm all");
  coloropt.push_back(50);
  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAcut1usgenmchrecoebus");
  legend.push_back("3cm gen mch");
  coloropt.push_back(2);
  file.push_back(sig30cmfile);
  name.push_back("cut1usrecoebus");
  legend.push_back("30cm all");
  coloropt.push_back(212);
  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selbarAcut1usgenmchrecoebus");
  legend.push_back("30cm gen mch");
  coloropt.push_back(210);
  file.push_back(sig1mfile);
  name.push_back("cut1usrecoebus");
  legend.push_back("1m all");
  coloropt.push_back(219);
  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAcut1usgenmchrecoebus");
  legend.push_back("1m gen mch");
  coloropt.push_back(218);
  file.push_back(sig3mfile);
  name.push_back("cut1usrecoebus");
  legend.push_back("3m all");
  coloropt.push_back(224);
  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbarAcut1usgenmchrecoebus");
  legend.push_back("3m gen mch");
  coloropt.push_back(222);
  legendEntries = legend;    
  //comparesamevariable(file, name, "egseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma ecal seed cluster time / ns");
  //comparesamevariable(file, name, "leadegseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{1} ecal seed cluster time / ns");
  //comparesamevariable(file, name, "leadegin5x5noiseclnd", 5, 300, 10, true, true, true, (float []){1e-3,1}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{1} #sigmai#etai#eta_{5x5} (noise clnd.)");
  //comparesamevariable(file, name, "leadegsmin", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} s_{min}");
  //comparesamevariable(file, name, "leadegsmaj", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} s_{maj}");
  //comparesamevariable(file, name, "leadegecalpfclustisoovere", -1, -1, 1, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} ecal iso./E");
  //comparesamevariable(file, name, "leadeghovereoversupcluse", -1, 150, 4, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} H/E");
  //comparesamevariable(file, name, "subleadegseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{2} ecal seed cluster time / ns");
  //comparesamevariable(file, name, "subleadegin5x5noiseclnd", 5, 300, 10, true, true, true, (float []){1e-3,1}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{2} #sigmai#etai#eta_{5x5} (noise clnd.)");
  //comparesamevariable(file, name, "subleadegsmin", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} s_{min}");
  //comparesamevariable(file, name, "subleadegsmaj", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} s_{maj}");
  //comparesamevariable(file, name, "subleadegecalpfclustisoovere", -1, -1, 1, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} ecal iso./E");
  //comparesamevariable(file, name, "subleadeghovereoversupcluse", -1, 150, 4, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} H/E");
  
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoeeus");
  legend.push_back("data, Z wind, veto id");
  coloropt.push_back(kBlack);
  file.push_back(datafile);
  name.push_back("cut1usrecoeeus");
  legend.push_back("data all");
  coloropt.push_back(16);
  file.push_back(dyfile);
  name.push_back("selelevetozwindidusrecoeeus");
  legend.push_back("DY");
  coloropt.push_back(4);
  file.push_back(sig3cmfile);
  name.push_back("cut1usrecoeeus");
  legend.push_back("3cm all");
  coloropt.push_back(50);
  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selecAcut1usgenmchrecoeeus");
  legend.push_back("3cm gen mch");
  coloropt.push_back(2);
  file.push_back(sig30cmfile);
  name.push_back("cut1usrecoeeus");
  legend.push_back("30cm all");
  coloropt.push_back(212);
  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selecAcut1usgenmchrecoeeus");
  legend.push_back("30cm gen mch");
  coloropt.push_back(210);
  file.push_back(sig1mfile);
  name.push_back("cut1usrecoeeus");
  legend.push_back("1m all");
  coloropt.push_back(219);
  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selecAcut1usgenmchrecoeeus");
  legend.push_back("1m gen mch");
  coloropt.push_back(218);
  file.push_back(sig3mfile);
  name.push_back("cut1usrecoeeus");
  legend.push_back("3m all");
  coloropt.push_back(224);
  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selecAcut1usgenmchrecoeeus");
  legend.push_back("3m gen mch");
  coloropt.push_back(222);
  legendEntries = legend;    
  //comparesamevariable(file, name, "egseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma ecal seed cluster time / ns");
  //comparesamevariable(file, name, "leadegseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{1} ecal seed cluster time / ns");
  //comparesamevariable(file, name, "leadegin5x5noiseclnd", 5, 600, 20, true, true, true, (float []){1e-3,1}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{1} #sigmai#etai#eta_{5x5} (noise clnd.)");
  //comparesamevariable(file, name, "leadegsmin", -1, -1, 40, true, true, true, (float []){5e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} s_{min}");
  //comparesamevariable(file, name, "leadegsmaj", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} s_{maj}");
  //comparesamevariable(file, name, "leadegecalpfclustisoovere", -1, -1, 1, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} ecal iso./E");
  //comparesamevariable(file, name, "leadeghovereoversupcluse", -1, 150, 4, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{1} H/E");
  //comparesamevariable(file, name, "subleadegseedclustime", 7500, 17000, 250, true, true, true, (float []){1e-3,10}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{2} ecal seed cluster time / ns");
  //comparesamevariable(file, name, "subleadegin5x5noiseclnd", 5, 600, 20, true, true, true, (float []){1e-3,1}, (float []){0.7,0.4,0.8,0.95}, true, "e/#gamma_{2} #sigmai#etai#eta_{5x5} (noise clnd.)");
  //comparesamevariable(file, name, "subleadegsmin", -1, -1, 40, true, true, true, (float []){5e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} s_{min}");
  //comparesamevariable(file, name, "subleadegsmaj", -1, -1, 40, true, true, true, (float []){5e-3,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} s_{maj}");
  //comparesamevariable(file, name, "subleadegecalpfclustisoovere", -1, -1, 1, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} ecal iso./E");
  //comparesamevariable(file, name, "subleadeghovereoversupcluse", -1, 150, 4, true, true, true, (float []){1e-5,1}, (float []){-1,0.4,0.8,0.95}, true, "e/#gamma_{2} H/E");

  // Plotting for rate estimation
  coloropt.clear();
  coloropt.push_back(kBlack);
  coloropt.push_back(2);
  coloropt.push_back(210);
  coloropt.push_back(218);
  coloropt.push_back(222);

  seltext[0] = "2 unseeded e/#gamma: p_{T}>10 GeV, |#eta|<2.5, smin<0.4";
  seltext[1] = "#sigmai#etai#eta<0.016(0.04), H/E<0.5(0.6), time<2ns";  
  //makeratehist("cut2usrecoegus", "subleadegpt", 50, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,1000}, 0.875);
  //makeratehist("cut2usrecoegus", "leadegpt", 50, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,1000}, 0.875);
  //makeratehist("cut2usrecoebus", "subleadegsmin", -1, 500, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,1000}, 0.875);
  //makeratehist("cut2usrecoebus", "leadegsmin", -1, 500, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,1000}, 0.875);
  //makeratehist("cut2usrecoebus", "subleadegsmaj", -1, 500, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,1000}, 0.875);
  //makeratehist("cut2usrecoebus", "leadegsmaj", -1, 500, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,1000}, 0.875);

  seltext[0] = "2 unseeded e/#gamma: p_{T}>10 GeV, |#eta|<2.5, smin<0.16";
  seltext[1] = "#sigmai#etai#eta<0.016(0.04), H/E<0.5(0.6), time<2ns";  
  //makeratehist("cut3usrecoegus", "subleadegpt", 50, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){55,4}, (float []){0,10}, 0.875);

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sigfile);
  name.push_back("genptgt10barselgeneg");
  legend.push_back("3 cm ");
  coloropt.push_back(2);
  file.push_back(sigfile);
  name.push_back("genptgt10barsel_dieg33mchgeneg");
  legend.push_back("dieg33");
  coloropt.push_back(45);
  file.push_back(sigfile);
  name.push_back("genptgt10barsel_time1nsAsminlt0p16mchgeneg");
  legend.push_back("time1nsAsminlt0p16");
  coloropt.push_back(46);
  legendEntries = legend;
  //comparesamevariable(file, name, "egmult", 5, 10, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "electron multiplicity");
  //comparesamevariable(file, name, "pt", 50, 120, 1, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "electron p_{T} / GeV");
  //comparesamevariable(file, name, "eta", -1, -1, 1, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "electron #eta");
  //comparesamevariable(file, name, "phi", -1, -1, 1, true, true, true, (float []){8e-1,1e3}, (float []){0.3,0.3,0.65,0.75}, false, "electron #phi");
  //comparesamevariable(file, name, "d0", 9000, 11000, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "electron d_{0} / cm");
  //comparesamevariable(file, name, "lxy", -1, -1, 100, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "electron l_{xy} / cm");
  //comparesamevariable(file, name, "log10d0", 50, -1, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "electron log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(file, name, "log10lxy", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "electron log_{10}l_{xy} / log_{10}cm");
  //comparesamevariable(file, name, "vx", 20000, 30000, 40, true, true, true, (float []){8e-1,1e4}, (float []){0.7,0.6,0.9,0.95}, false, "electron v_{x} / cm");
  //comparesamevariable(file, name, "vy", 20000, 30000, 40, true, true, true, (float []){8e-1,1e4}, (float []){0.7,0.6,0.9,0.95}, false, "electron v_{y} / cm");
  //comparesamevariable(file, name, "vz", -1, -1, 80, true, true, true, (float []){8e-1,1e4}, (float []){0.7,0.6,0.9,0.95}, false, "electron v_{z} / cm");
  //comparesamevariable(file, name, "pvx", 24990, 25010, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "mother v_{x} / cm");
  //comparesamevariable(file, name, "pvy", 24990, 25010, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.8,0.6,0.95,0.95}, false, "mother v_{y} / cm");
  //comparesamevariable(file, name, "pvz", 24000, 26000, 20, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "mother v_{z} / cm");
  //comparesamevariable(file, name, "prompteta", -1, -1, 1, true, true, true, (float []){8e-1,1e4}, (float []){0.8,0.6,0.95,0.95}, false, "prompt equivalent electron #eta");
  //comparesamevariable(file, name, "promptphi", -1, -1, 1, true, true, true, (float []){8e-1,1e3}, (float []){0.3,0.3,0.65,0.75}, false, "prompt equivalent electron #phi");

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAdieg33caloidlusrecomchgenel");
  legend.push_back("dieg33");
  coloropt.push_back(kBlue);
  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAcuttimedelayonlyusrecomchgenel");
  legend.push_back(" c#tau = 1 m, t > 1.4 ns");
  coloropt.push_back(kRed-2);
  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAcuttimegt1nsonlyusrecomchgenel");
  legend.push_back(" c#tau = 1 m, t > 1 ns");
  coloropt.push_back(kRed);
  legendEntries = legend;
  comparesamevariable(file, name, "log10d0", 350, -1, 20, true, true, true, (float []){8e-1,5e2}, (float []){0.6,0.6,0.75,0.95}, false, "electron log_{10}d_{0} / log_{10}cm");

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAdieg33caloidlusrecomchgenel");
  legend.push_back("dieg33");
  coloropt.push_back(kBlue);
  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAcuttimedelaysminusrecomchgenel");
  legend.push_back(" c#tau = 3 cm, t > 1.4 ns");
  coloropt.push_back(kRed-2);
  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAcutsminlt0p16onlyusrecomchgenel");
  legend.push_back(" c#tau = 3 cm, t > 1 ns");
  coloropt.push_back(kRed);
  legendEntries = legend;
  comparesamevariable(file, name, "log10d0", 150, -1, 20, true, true, true, (float []){8e-1,5e2}, (float []){0.6,0.6,0.75,0.95}, false, "electron log_{10}d_{0} / log_{10}cm");

  return -1;
}
