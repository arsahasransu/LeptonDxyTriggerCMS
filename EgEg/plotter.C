#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 8.3/5967; // Data rate - 8.3 Hz, 5967 new, 771 old events from prescaled dig33_caloidl selection
//double datasf = 6.8/76; // Data rate equivalent to 6.8 Hz, 76 events from unprescaled dieg70 selection
double sig3cmsf = 1.0/28;
double sig30cmsf = 1.0/26;
double sig1msf = 1.0/11;
double sig3msf = 1.0/3;
TString cutdeets = "Cut details";
TFile* datafile = TFile::Open("hists_Efmrl.root","READ");
TFile* dyfile = TFile::Open("hists_DY.root","READ");
TFile* sig3cmfile = TFile::Open("hists_M200dM20ctau3cm__2.root","READ");
TFile* sig30cmfile = TFile::Open("hists_M200dM20ctau30cm__2.root","READ");
TFile* sig1mfile = TFile::Open("hists_M200dM20ctau1m__2.root","READ");
TFile* sig3mfile = TFile::Open("hists_M200dM20ctau3m__2.root","READ");
TFile* sigfile = TFile::Open("hists_M200dM20.root","READ");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};
std::vector<TString> histtype{"p e1", "hist same"};
std::vector<int> markerstyle{20, 24};
std::vector<int> markersize{10, 10};
std::vector<TString> legendmarkerstyle{"lep", "l"};
std::vector<double> scale{1, 1};

int comparesamevariable(std::vector<TFile*> file, std::vector<TString> cutname, TString var, int nbins=0, double *rebin=0, bool logY=false, bool underflow=false, bool overflow=false, float yrange[]=(float []){0.1,100}, float legPos[]=(float []){0.7,0.75,0.95,1}, float seltextpos[]=(float []){0.1,1}, bool normalize=true, TString xaxistitle="xaxis", TString yaxistitle="yaxis", bool dividebybinw=true, TString savefilename="default") {

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

  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {
    if(nbins>0) {
      allhists[histctr] = (TH1F*) allhists[histctr]->Rebin(nbins,"newx",rebin);
    }
  
    // Make changes to sig and bkg to enable good basic plotting
    double err = 0.0;
    if(underflow) {
      allhists[histctr]->GetXaxis()->SetRange(0, allhists[histctr]->GetNbinsX());
    }
    err = 0.0;
    if(overflow) {
      allhists[histctr]->GetXaxis()->SetRange(1, allhists[histctr]->GetNbinsX()+1);
    }
    if(underflow && overflow) {
      allhists[histctr]->GetXaxis()->SetRange(0, allhists[histctr]->GetNbinsX()+1);
    }

    if(dividebybinw) {
      for(unsigned int bincnt=1; bincnt<=allhists[histctr]->GetNbinsX(); bincnt++) {
	allhists[histctr]->SetBinContent(bincnt, allhists[histctr]->GetBinContent(bincnt)/allhists[histctr]->GetBinWidth(bincnt));
      }
    }
    
    allhists[histctr]->GetXaxis()->SetTitle(xtitle);
    if(normalize) allhists[histctr]->GetYaxis()->SetTitle(yaxistitle);
    else allhists[histctr]->GetYaxis()->SetTitle(yaxistitle);

    allhists[histctr]->SetTitle("");

    allhists[histctr]->SetLineWidth(2);
    allhists[histctr]->SetLineColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  c1 = fancy_enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,normalize, histtype, markerstyle, markersize, legendmarkerstyle, scale);

  TPad* pad = (TPad*) c1->FindObject("pad3");
  pad->cd();

  auto linehor = new TLine(*rebin,yrange[1],*(rebin+nbins),yrange[1]);
  linehor->SetLineColor(kBlack);
  linehor->SetLineWidth(4);
  linehor->SetLineStyle(1);
  linehor->Draw();
  
  auto linevert = new TLine(*(rebin+nbins),yrange[0],*(rebin+nbins),yrange[1]);
  linevert->SetLineColor(kBlack);
  linevert->SetLineWidth(4);
  linevert->SetLineStyle(1);
  linevert->Draw();
  
  TLatex sel;
  sel.SetTextFont(42);
  sel.SetTextSize(0.045);
  if(!logY) sel.DrawLatex(seltextpos[0], seltextpos[1]+0.1*(yrange[1]-yrange[0]), seltext[0]);
  else sel.DrawLatex(seltextpos[0], seltextpos[1]+0.02*(yrange[1]-yrange[0]), seltext[0]);
  sel.DrawLatex(seltextpos[0], seltextpos[1], seltext[1]);
  
  c1->SaveAs("./HQPlots/"+foldername+"/"+savefilename+".png");
  c1->SaveAs("./HQPlots/"+foldername+"/"+savefilename+".C");

  return -1;
}

int makeratehist(std::vector<TString> cutname, TString var, int nbins=0, double *rebin=0, bool logY=false, bool overflow=false, float legPos[]=(float []){0.7,0.75,0.95,1}, float seltextpos[]=(float []){0.1,1}, float yrange[]=(float []){0.1,1}, float drawsignalline=-1.0, TString xaxistitle="p_{T} [GeV]", TString yaxistitle="rate [Hz]") {

  std::vector<TString> legNam;
  legNam.push_back("2018 data");
  legNam.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  legNam.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  legNam.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  legNam.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");

  TString histname0 = cutname[0]+"_"+var;
  TString histname1 = cutname[1]+"_"+var;
  TString histname2 = cutname[2]+"_"+var;
  TString histname3 = cutname[3]+"_"+var;
  TString histname4 = cutname[4]+"_"+var;
  TH1F* datahist;
  datahist = (TH1F*) datafile->Get(histname0);
  datahist->Scale(datasf);
  TH1F* sig3cmhist = (TH1F*) sig3cmfile->Get(histname1);
  sig3cmhist->Scale(sig3cmsf);
  TH1F* sig30cmhist = (TH1F*) sig30cmfile->Get(histname2);
  sig30cmhist->Scale(sig30cmsf);
  TH1F* sig1mhist = (TH1F*) sig1mfile->Get(histname3);
  sig1mhist->Scale(sig1msf);
  TH1F* sig3mhist = (TH1F*) sig3mfile->Get(histname4);
  sig3mhist->Scale(sig3msf);

  std::vector<TH1F*> allhists;
  std::vector<TH1F*> allratehists;
  
  if(nbins>0) {
    datahist = (TH1F*) datahist->Rebin(nbins,"newx",rebin);
    sig3cmhist = (TH1F*) sig3cmhist->Rebin(nbins,"newx",rebin);
    sig30cmhist = (TH1F*) sig30cmhist->Rebin(nbins,"newx",rebin);
    sig1mhist = (TH1F*) sig1mhist->Rebin(nbins,"newx",rebin);
    sig3mhist = (TH1F*) sig3mhist->Rebin(nbins,"newx",rebin);
  }
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
    allratehists[histctr]->GetXaxis()->SetTitle(xaxistitle);
    allratehists[histctr]->GetYaxis()->SetTitle(yaxistitle);

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
  axis = new TGaxis(*(rebin+nbins),yrange[0],*(rebin+nbins),yrange[1],yrange[0]/axisscale,yrange[1]/axisscale,510,"-L");
  axis->SetLineColor(coloropt[1]);
  axis->SetLineWidth(4);
  axis->SetLabelColor(coloropt[1]);
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.06);
  axis->SetLabelOffset(-0.035);
  axis->Draw();

  if(drawsignalline!=-1.0) {
    drawsignalline = allratehists[1]->GetMaximum();
    TLine *signalline = new TLine(*rebin,drawsignalline,*(rebin+nbins),drawsignalline);
    signalline->SetLineWidth(2);
    signalline->SetLineColor(coloropt[1]);
    signalline->SetLineStyle(9);
    signalline->Draw();
  }

  auto line = new TLine(*rebin,yrange[1],*(rebin+nbins),yrange[1]);
  line->SetLineColor(kBlack);
  line->SetLineWidth(4);
  line->SetLineStyle(1);
  line->Draw();
  
  TLatex sel;
  sel.SetTextFont(42);
  sel.SetTextSize(0.045);
  sel.DrawLatex(seltextpos[0], seltextpos[1]+0.1*(yrange[1]-yrange[0]), seltext[0]);
  sel.DrawLatex(seltextpos[0], seltextpos[1], seltext[1]);
  
  c1->SaveAs("./HQPlots/"+cutname[0]+"/"+cutname[0]+"_"+var+"_ratehist.png");
  c1->SaveAs("./HQPlots/"+cutname[0]+"/"+cutname[0]+"_"+var+"_ratehist.C");
  
  TCanvas* c2;
  c2 = new TCanvas();
  c2 = enhance_plotter_rate(allratehists, legNam, allratehists[0]->GetXaxis()->GetTitle(),allratehists[0]->GetYaxis()->GetTitle(),legPos,yrange,logY,false);
  c2->GetCanvas()->SetGrayscale();

  TPad* pad3 = (TPad*) c2->FindObject("pad3");
  pad3->cd();
  axis->Draw();

  if(drawsignalline!=-1.0) {
    drawsignalline = allratehists[1]->GetMaximum();
    TLine *signalline = new TLine(*rebin,drawsignalline,*(rebin+nbins),drawsignalline);
    signalline->SetLineWidth(2);
    signalline->SetLineColor(coloropt[1]);
    signalline->SetLineStyle(9);
    signalline->Draw();
  }

  line->Draw();
  
  sel.DrawLatex(seltextpos[0], seltextpos[1]+0.1*(yrange[1]-yrange[0]), seltext[0]);
  sel.DrawLatex(seltextpos[0], seltextpos[1], seltext[1]);
  
  c2->SaveAs("./HQPlots/"+cutname[0]+"/"+cutname[0]+"_"+var+"_ratehist_grayscale.png");

  return -1;
}

int plotter() {

  std::vector<int> coloroptrate{1, kRed+2, kRed-3, kRed-7, kRed-9};

  std::vector<TFile*> file;
  std::vector<TString> name;
  std::vector<TString> legend;

  // REFERENCE TRIGGER PROPERTIES ON WHICH RATE IS SHOWN
  vector<double> dieg33caloidlidusrecoegus_subleadegpt_binspt{30,33,36,40,44,48,52,56,60,80,100};
  int dieg33caloidlidusrecoegus_subleadegpt_nbinspt = 10;
  seltext[0] = "#geq1 seeded, #geq2 unseeded e/#gamma";
  seltext[1] = "p_{T}>33 GeV, Calo ID";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(datafile);
  name.push_back("dieg33caloidlidusrecoegus");
  legend.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);

  file.push_back(sig3cmfile);
  name.push_back("dieg33caloidlidusrecoegus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("dieg33caloidlidusrecoegus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("dieg33caloidlidusrecoegus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("dieg33caloidlidusrecoegus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  legendEntries = legend;    
  comparesamevariable(file, name, "subleadegpt", dieg33caloidlidusrecoegus_subleadegpt_nbinspt, &dieg33caloidlidusrecoegus_subleadegpt_binspt[0], true, false, false, (float []){3e-2,3}, (float []){0.525,0.62,0.775,0.97}, (float []){64,0.15}, true, "e/#gamma_{2} p_{T} [GeV]", "normalized events / GeV", true, "dieg33caloidlidusrecoegus_subleadegpt");
  
  makeratehist({"dieg33caloidlidusrecoegus","dieg33caloidlidusrecoegus","dieg33caloidlidusrecoegus","dieg33caloidlidusrecoegus","dieg33caloidlidusrecoegus"}, "subleadegpt", dieg33caloidlidusrecoegus_subleadegpt_nbinspt, &dieg33caloidlidusrecoegus_subleadegpt_binspt[0], false, false, (float []){0.48,0.65,0.68,0.975}, (float []){64,4.45}, (float []){0,11}, 0.875, "e/#gamma_{2} p_{T} [GeV]");

  vector<double> cut1usrecoegus_subleadegpt_binspt{9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,33,36,40};
  int cut1usrecoegus_subleadegpt_nbinspt = cut1usrecoegus_subleadegpt_binspt.size()-1;
  seltext[0] = "#geq1 seeded, #geq2 unseeded e/#gamma";
  seltext[1] = "p_{T}>10 GeV, Calo ID";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(datafile);
  name.push_back("cut1usrecoegus");
  legend.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);

  file.push_back(sig3cmfile);
  name.push_back("cut1usrecoegus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("cut1usrecoegus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("cut1usrecoegus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("cut1usrecoegus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "subleadegpt", cut1usrecoegus_subleadegpt_nbinspt, &cut1usrecoegus_subleadegpt_binspt[0], true, false, false, (float []){3e-4,3}, (float []){0.525,0.62,0.775,0.97}, (float []){64,0.15}, true, "e/#gamma_{2} p_{T} [GeV]", "normalized events / GeV", true, "genptgt10Acut1usgenmchrecoegus_subleadegpt");
  
  //makeratehist({"cut1usrecoegus","cut1usrecoegus","cut1usrecoegus","cut1usrecoegus","cut1usrecoegus"}, "subleadegpt", cut1usrecoegus_subleadegpt_nbinspt, &cut1usrecoegus_subleadegpt_binspt[0], false, false, (float []){0.48,0.675,0.68,0.975}, (float []){24,60}, (float []){0,150}, 0.875, "e/#gamma_{2} p_{T} [GeV]");

  vector<double> genbasicptgt10selbarAnoselus_binsseedclustime{-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0};
  int genbasicptgt10selbarAnoselus_nbinsseedclustime = genbasicptgt10selbarAnoselus_binsseedclustime.size()-1;
  vector<double> genbasicptgt10selbarAnoselus_binsmin{0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6};
  int genbasicptgt10selbarAnoselus_nbinsmin = genbasicptgt10selbarAnoselus_binsmin.size()-1;
  seltext[0] = "gen. p_{T} > 10 GeV, |#eta|<1.479";
  seltext[1] = "#geq1 unseeded e/#gamma";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(datafile);
  name.push_back("noselusrecoebus");
  legend.push_back("2018 data, all e");
  coloropt.push_back(kBlack);
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);

  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("2018 data, prompt e");
  coloropt.push_back(kBlack);
  histtype.push_back("p e1 same");
  markerstyle.push_back(24);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dyfile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(4);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "egseedclustime", genbasicptgt10selbarAnoselus_nbinsseedclustime, &genbasicptgt10selbarAnoselus_binsseedclustime[0], true, false, false, (float []){3e-3,0.7}, (float []){0.525,0.6,0.775,0.99}, (float []){-1.1,0.4}, true, "e/#gamma ecal seed cluster time [ns]","normalized events / 0.2 ns", false, "genbasicptgt10selbarAnoselusgenmchrecoebus_egseedclustime");
  //comparesamevariable(file, name, "leadegsmin", genbasicptgt10selbarAnoselus_nbinsmin, &genbasicptgt10selbarAnoselus_binsmin[0], true, false, false, (float []){3e-3,0.7}, (float []){0.525,0.6,0.775,0.99}, (float []){0.02,0.4}, true, "e/#gamma smin","normalized events / 0.01 units", false, "genbasicptgt10selbarAnoselusgenmchrecoebus_leadsmin");
  //comparesamevariable(file, name, "leadegsmin", 10, 600, 10, true, false, false, (float []){3e-3,0.9}, (float []){0.525,0.52,0.775,0.97}, true, "e/#gamma s_{min}");
    
  return -1;
}
