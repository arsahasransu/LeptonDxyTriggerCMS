#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 8.3/771; // Data rate - 8.3 Hz, 771 events from prescaled dig33_caloidl selection
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

int makeratehist(TString cutname, TString var, int nbins=0, double *rebin=0, bool logY=false, bool overflow=false, float legPos[]=(float []){0.7,0.75,0.95,1}, float seltextpos[]=(float []){0.1,1}, float yrange[]=(float []){0.1,1}, float drawsignalline=-1.0, TString xaxistitle="p_{T} [GeV]", TString yaxistitle="rate [Hz]") {

  std::vector<TString> legNam;
  legNam.push_back("2018 data");
  legNam.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  legNam.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  legNam.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  legNam.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");

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
  
  c1->SaveAs("./HQPlots/"+cutname+"/"+cutname+"_"+var+"_ratehist.png");
  c1->SaveAs("./HQPlots/"+cutname+"/"+cutname+"_"+var+"_ratehist.C");
  
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
  
  c2->SaveAs("./HQPlots/"+cutname+"/"+cutname+"_"+var+"_ratehist_grayscale.png");

  return -1;
}

int plotter() {

  std::vector<int> coloroptrate{1, kRed+2, kRed-3, kRed-7, kRed-9};

  std::vector<TFile*> file;
  std::vector<TString> name;
  std::vector<TString> legend;

  // GENERATOR LEPTON PROPERTIES
  vector<double> gen_binslog10d0{-4.0,-3.6,-3.2,-2.8,-2.4,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0};
  int gen_nbinslog10d0 = gen_binslog10d0.size()-1;
  vector<double> gen_binspt{10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,34,38,42,46,50};
  int gen_nbinspt = gen_binspt.size()-1;
  seltext[0] = " ";
  seltext[1] = " ";
  
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();
  
  file.push_back(sig3cmfile);
  name.push_back("gennoselgeneg");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("gennoselgeneg");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig1mfile);
  name.push_back("gennoselgeneg");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig3mfile);
  name.push_back("gennoselgeneg");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "log10d0", gen_nbinslog10d0, &gen_binslog10d0[0], true, false, false, (float []){6e-5,6}, (float []){0.525,0.7,0.775,0.99}, (float []){0.015,0.015}, true, "gen. lepton log_{10}d_{0} [log_{10}cm]", "normalized events / log_{10}cm", true, "genleptond0");
  //comparesamevariable(file, name, "subleadpt", gen_nbinspt, &gen_binspt[0], true, false, false, (float []){6e-5,6}, (float []){0.525,0.7,0.775,0.99}, (float []){0.015,0.015}, true, "gen. lepton_{2} p_{T} [GeV]", "normalized events / GeV", true, "genleptoneg2pt");
    
  // ANGULAR MATCHING - BEFORE PROMPT CORRECTION
  vector<double> bar_deta{-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.45,-0.4,-0.35,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1,1.4};
  int bar_deta_nbins = 49;
  vector<double> bar_qdphi{-0.8,-0.7,-0.6,-0.5,-0.45,-0.4,-0.35,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.4,0.45,0.5};
  int bar_qdphi_nbins = bar_qdphi.size()-1;
  seltext[0] = "gen e: p_{T}>10 GeV, |#eta_{prompt}|<1.479";
  seltext[1] = ">1 unseeded e/#gamma";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kRed-5);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kRed-2);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "dEta", bar_deta_nbins, &bar_deta[0], true, false, false, (float []){8e-5,0.9}, (float []){0.525,0.72,0.775,0.97}, (float []){-0.9,0.1}, true, "#Delta#eta(gen,SC)", "normalized events / unit", true);
  //comparesamevariable(file, name, "qdPhi", bar_qdphi_nbins, &bar_qdphi[0], true, false, false, (float []){2e-4,0.9}, (float []){0.15,0.72,0.4,0.97}, (float []){-0.75,0.02}, true, "charge#times#Delta#phi(gen,SC)", "normalized events / unit", true);
  
  seltext[0] = "gen e: p_{T}>10 GeV, 1.479<|#eta_{prompt}|<2.4";
  seltext[1] = ">1 unseeded e/#gamma";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kRed-5);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kRed-2);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "dEta", bar_deta_nbins, &bar_deta[0], true, false, false, (float []){8e-5,0.9}, (float []){0.525,0.72,0.775,0.97}, (float []){-0.9,0.1}, true, "#Delta#eta(gen,SC)", "normalized events / unit", true);
  //comparesamevariable(file, name, "qdPhi", bar_qdphi_nbins, &bar_qdphi[0], true, false, false, (float []){2e-4,0.9}, (float []){0.15,0.72,0.4,0.97}, (float []){-0.75,0.02}, true, "charge#times#Delta#phi(gen,SC)", "normalized events / unit", true);
  
  // ANGULAR MATCHING - AFTER PROMPT CORRECTION
  seltext[0] = "gen e: p_{T}>10 GeV, |#eta_{prompt}|<1.479";
  seltext[1] = ">1 unseeded e/#gamma";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kRed-5);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kRed-2);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "dPromptEta", bar_deta_nbins, &bar_deta[0], true, false, false, (float []){8e-5,0.9}, (float []){0.525,0.72,0.775,0.97}, (float []){0.2,0.01}, true, "#Delta#eta_{prompt}(gen,SC)", "normalized events / unit", true);
  //comparesamevariable(file, name, "qdPromptPhi", bar_qdphi_nbins, &bar_qdphi[0], true, false, false, (float []){2e-4,0.9}, (float []){0.15,0.72,0.4,0.97}, (float []){-0.75,0.02}, true, "charge#times#Delta#phi_{prompt}(gen,SC)", "normalized events / unit", true);
  
  seltext[0] = "gen e: p_{T}>10 GeV, 1.479<|#eta_{prompt}|<2.4";
  seltext[1] = ">1 unseeded e/#gamma";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kRed-5);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kRed-2);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "dPromptEta", bar_deta_nbins, &bar_deta[0], true, false, false, (float []){8e-5,0.9}, (float []){0.525,0.72,0.775,0.97}, (float []){0.2,0.01}, true, "#Delta#eta_{prompt}(gen,SC)", "normalized events / unit", true);
  //comparesamevariable(file, name, "qdPromptPhi", bar_qdphi_nbins, &bar_qdphi[0], true, false, false, (float []){2e-4,0.9}, (float []){0.15,0.72,0.4,0.97}, (float []){-0.75,0.02}, true, "charge#times#Delta#phi_{prompt}(gen,SC)", "normalized events / unit", true);

  // KINEMATIC PROPERTIES OF THE SIGNAL LEPTON BEFORE AND AFTERGEN MCH COMPARED TO THE DATA AND PROMPT MC
  vector<double> genbasicselptgt10_leadebsinin_binsleadebsinin{0.0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03};
  int genbasicselptgt10_leadebsinin_nbinsleadebsinin = genbasicselptgt10_leadebsinin_binsleadebsinin.size()-1;
  seltext[0] = ">1 unseeded e/#gamma";
  seltext[1] = "p_{T}>10 GeV, |#eta|<1.48, No ID";  

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
  name.push_back("basicselusrecoebus");
  legend.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("2018 data, prompt e");
  coloropt.push_back(kGray+1);
  histtype.push_back("p same");
  markerstyle.push_back(24);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dyfile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  file.push_back(sig3cmfile);
  name.push_back("basicselusrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig3cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("c#tau = 3 cm, gen. matched");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "subleadegin5x5noiseclnd", genbasicselptgt10_leadebsinin_nbinsleadebsinin, &genbasicselptgt10_leadebsinin_binsleadebsinin[0], true, false, false, (float []){3e-3,1.4}, (float []){0.52,0.65,0.77,0.99}, (float []){0.015,0.015}, true, "e/#gamma_{2} #sigmai#etai#eta", "normalized events / unit", true, "compare_genmchToData_EB_subleadegin5x5noiseclnd");
  
  vector<double> genbasicselptgt10_leadeesinin_binsleadeesinin{0.0,0.004,0.008,0.012,0.016,0.02,0.024,0.028,0.032,0.036,0.04,0.044,0.048,0.052,0.056,0.06,0.064,0.068,0.072,0.076,0.08};
  int genbasicselptgt10_leadeesinin_nbinsleadeesinin = genbasicselptgt10_leadeesinin_binsleadeesinin.size()-1;
  seltext[0] = ">1 unseeded e/#gamma";
  seltext[1] = "p_{T}>10 GeV, 1.479<|#eta|<2.5, No ID";  

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
  name.push_back("basicselusrecoeeus");
  legend.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoeeus");
  legend.push_back("2018 data, prompt e");
  coloropt.push_back(kGray+1);
  histtype.push_back("p same");
  markerstyle.push_back(24);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dyfile);
  name.push_back("selelevetozwindidusrecoeeus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  file.push_back(sig3cmfile);
  name.push_back("basicselusrecoeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig3cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("c#tau = 3 cm, gen. matched");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "subleadegin5x5noiseclnd", genbasicselptgt10_leadeesinin_nbinsleadeesinin, &genbasicselptgt10_leadeesinin_binsleadeesinin[0], true, false, false, (float []){3e-3,9}, (float []){0.52,0.65,0.77,0.99}, (float []){0.05,0.015}, true, "e/#gamma_{2} #sigmai#etai#eta", "normalized events / unit", true, "compare_genmchToData_EE_subleadegin5x5noiseclnd");
  
  seltext[0] = ">1 unseeded e/#gamma";
  seltext[1] = "p_{T}>10 GeV, |#eta|<1.48, No ID";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();
  
  file.push_back(sig3cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig1mfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig3mfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "subleadegin5x5noiseclnd", genbasicselptgt10_leadebsinin_nbinsleadebsinin, &genbasicselptgt10_leadebsinin_binsleadebsinin[0], true, false, false, (float []){3e-3,1.4}, (float []){0.525,0.65,0.775,0.99}, (float []){0.015,0.015}, true, "e/#gamma_{2} #sigmai#etai#eta", "normalized events / unit", true, "comapre_var_EB");
  
  seltext[0] = ">1 unseeded e/#gamma";
  seltext[1] = "p_{T}>10 GeV, 1.479<|#eta|<2.5, No ID";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();
  
  file.push_back(sig3cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig1mfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig3mfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange+1);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "subleadegin5x5noiseclnd", genbasicselptgt10_leadeesinin_nbinsleadeesinin, &genbasicselptgt10_leadeesinin_binsleadeesinin[0], true, false, false, (float []){3e-3,1.4}, (float []){0.525,0.65,0.775,0.99}, (float []){0.05,0.015}, true, "e/#gamma_{2} #sigmai#etai#eta", "normalized events / unit", true, "compare_var");

  // REFERENCE TRIGGER PROPERTIES ON WHICH RATE IS SHOWN
  vector<double> dieg33caloidlidusrecoegus_subleadegpt_binspt{30,33,36,40,44,48,52,56,60,80,100};
  int dieg33caloidlidusrecoegus_subleadegpt_nbinspt = 10;
  seltext[0] = "1 seeded, 2 unseeded e/#gamma";
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
  //comparesamevariable(file, name, "subleadegpt", dieg33caloidlidusrecoegus_subleadegpt_nbinspt, &dieg33caloidlidusrecoegus_subleadegpt_binspt[0], true, false, false, (float []){3e-2,3}, (float []){0.525,0.62,0.775,0.97}, (float []){64,0.15}, true, "e/#gamma_{2} p_{T} [GeV]", "normalized events / GeV", true, "subleadegpt");
  
  //makeratehist("dieg33caloidlidusrecoegus", "subleadegpt", dieg33caloidlidusrecoegus_subleadegpt_nbinspt, &dieg33caloidlidusrecoegus_subleadegpt_binspt[0], false, false, (float []){0.48,0.65,0.68,0.975}, (float []){64,4.45}, (float []){0,11}, 0.875, "e/#gamma_{2} p_{T} [GeV]");

  seltext[0] = "1 seeded, 2 unseeded e/#gamma";
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
  //comparesamevariable(file, name, "subleadegpt", dieg33caloidlidusrecoegus_subleadegpt_nbinspt, &dieg33caloidlidusrecoegus_subleadegpt_binspt[0], true, false, false, (float []){3e-2,3}, (float []){0.525,0.62,0.775,0.97}, (float []){64,0.15}, true, "e/#gamma_{2} p_{T} [GeV]", "normalized events / GeV", true, "subleadegpt");
  
  //makeratehist("dieg33caloidlidusrecoegus", "subleadegpt", dieg33caloidlidusrecoegus_subleadegpt_nbinspt, &dieg33caloidlidusrecoegus_subleadegpt_binspt[0], false, false, (float []){0.48,0.65,0.68,0.975}, (float []){64,4.45}, (float []){0,11}, 0.875, "e/#gamma_{2} p_{T} [GeV]");

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
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kRed-5);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbarAnoselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kRed-2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "egseedclustime", 8500, 16000, 250, true, false, false, (float []){3e-3,0.9}, (float []){0.525,0.52,0.775,0.97}, true, "e/#gamma ecal seed cluster time [ns]");
  //comparesamevariable(file, name, "leadegsmin", 10, 600, 10, true, false, false, (float []){3e-3,0.9}, (float []){0.525,0.52,0.775,0.97}, true, "e/#gamma s_{min}");
    
  return -1;
}
