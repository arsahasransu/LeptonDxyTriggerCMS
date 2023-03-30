#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 0.5/(24+22+20+21); // Data rate - 8.3 Hz, 5967 new, 771 old events from prescaled dig33_caloidl selection
//double datasf = 6.8/76; // Data rate equivalent to 6.8 Hz, 76 events from unprescaled dieg70 selection
double sig3cmsf = 1.0/29;
double sig30cmsf = 1.0/10;
double sig1msf = 1.0/19;
double sig3msf = 1.0/3;
TString cutdeets = "Cut details";
TFile* datafile = TFile::Open("hists_Efmrl.root","READ");
TFile* sig3cmfile = TFile::Open("hists_M200dM20ctau3cm.root","READ");
TFile* sig30cmfile = TFile::Open("hists_M200dM20ctau30cm.root","READ");
TFile* sig1mfile = TFile::Open("hists_M200dM20ctau1m.root","READ");
TFile* sig3mfile = TFile::Open("hists_M200dM20ctau3m.root","READ");
TFile* h24mu12file = TFile::Open("hists_H2L4MuMff12.root","READ");
TFile* h24mu25file = TFile::Open("hists_H2L4MuMff25.root","READ");
TFile* h24mu50file = TFile::Open("hists_H2L4MuMff50.root","READ");

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
  else sel.DrawLatex(seltextpos[0], seltextpos[1]+0.02*seltextpos[0], seltext[0]);
  sel.DrawLatex(seltextpos[0], seltextpos[1], seltext[1]);
  
  c1->SaveAs("./dirplots/"+foldername+"/"+savefilename+".png");
  c1->SaveAs("./dirplots/"+foldername+"/"+savefilename+".C");

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

  double axisscale = 0.0, axisscale2, axisscale3, axisscale4;
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
    /*if(histctr==1) {
      axisscale = allratehists[0]->GetBinContent(0)*0.8/allratehists[1]->GetBinContent(0);
      allratehists[histctr]->Scale(axisscale);
    }
    if(histctr>1) {
      allratehists[histctr]->Scale(axisscale);
      }*/
  }
  axisscale = allratehists[0]->GetBinContent(0)*0.8/allratehists[1]->GetBinContent(0);
  axisscale2 = allratehists[0]->GetBinContent(0)*0.8/allratehists[2]->GetBinContent(0);
  axisscale = axisscale2<axisscale?axisscale2:axisscale;
  axisscale3 = allratehists[0]->GetBinContent(0)*0.8/allratehists[3]->GetBinContent(0);
  axisscale = axisscale3<axisscale?axisscale3:axisscale;
  axisscale4 = allratehists[0]->GetBinContent(0)*0.8/allratehists[4]->GetBinContent(0);
  axisscale = axisscale4<axisscale?axisscale4:axisscale;
  cout<<axisscale<<"\t"<<axisscale2<<"\t"<<axisscale3<<"\t"<<axisscale4<<endl;
  allratehists[1]->Scale(axisscale);
  allratehists[2]->Scale(axisscale);
  allratehists[3]->Scale(axisscale);
  allratehists[4]->Scale(axisscale);
  cout<<allratehists[1]->GetMaximum()<<endl;
  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter_rate(allratehists, legNam, allratehists[0]->GetXaxis()->GetTitle(),allratehists[0]->GetYaxis()->GetTitle(),legPos,yrange,logY,false);

  TPad* pad = (TPad*) c1->FindObject("pad3");
  pad->cd();
  TGaxis *axis;
  if(!logY) {
    axis = new TGaxis(*(rebin+nbins),yrange[0],*(rebin+nbins),yrange[1],yrange[0]/axisscale,yrange[1]/axisscale,510,"-L");
    axis->SetLabelOffset(-0.035);
  }
  else {
    axis = new TGaxis(*(rebin+nbins),yrange[0],*(rebin+nbins),yrange[1],yrange[0]/axisscale,yrange[1]/axisscale,510,"-LG");
    axis->SetLabelOffset(-0.035);
  }
  axis->SetLineColor(coloropt[1]);
  axis->SetLineWidth(4);
  axis->SetLabelColor(coloropt[1]);
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.06);
  axis->Draw();

  if(drawsignalline!=-1.0) {
    drawsignalline = allratehists[1]->GetMaximum();
    TLine *signalline = new TLine(*rebin,drawsignalline,*(rebin+nbins),drawsignalline);
    signalline->SetLineWidth(3);
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
  
  c1->SaveAs("./dirplots/"+cutname[0]+"/"+cutname[0]+"_"+var+"_ratehist.png");
  c1->SaveAs("./dirplots/"+cutname[0]+"/"+cutname[0]+"_"+var+"_ratehist.C");
  
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
  
  c2->SaveAs("./dirplots/"+cutname[0]+"/"+cutname[0]+"_"+var+"_ratehist_grayscale.png");

  return -1;
}

// Comparison of the type cross-check between two histogram - filled v hollow
int crossChecktwohist(std::vector<TFile*> file, std::vector<TString> cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool underflow=false, bool overflow=false, float yrange[]=(float []){0.1,100}, float legPos[]=(float []){0.7,0.75,0.95,1}, bool normalize=true, TString xaxistitle="xaxis") {
  
  if(file.size()<2 || cutname.size()<2) {
    cout<<"Need atleast 2 histograms to cross-check"<<endl;
    return -1;
  }
  
  TH1F* histfilled = (TH1F*) file[0]->Get(cutname[0]+"_"+var);
  histfilled = (TH1F*) histfilled->Clone(cutname[0]+"_"+var);
  TH1F* histhollow = (TH1F*) file[1]->Get(cutname[1]+"_"+var);
  histhollow = (TH1F*) histhollow->Clone(cutname[1]+"_"+var);
  
  std::vector<TH1F*> allhists;
  allhists.push_back(histfilled);
  allhists.push_back(histhollow);
  
  // Get the title from histogram title
  TString xtitle = xaxistitle;
  
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
    allhists[histctr]->GetYaxis()->SetTitle("number of events");
    
    allhists[histctr]->SetTitle("");
  }
  
  allhists[0]->SetFillColor(coloropt[0]);
  allhists[1]->SetLineColor(coloropt[1]);
  
  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,normalize);
  c1->SaveAs("./dirplots/"+cutname[0]+"_"+cutname[1]+"/"+cutname[0]+"_"+cutname[1]+"_"+var+".png");
  c1->SaveAs("./dirplots/"+cutname[0]+"_"+cutname[1]+"/"+cutname[0]+"_"+cutname[1]+"_"+var+".C");
  
  return -1;
}

int newplotter() {

  std::vector<int> coloroptrate{kBlack, kRed+3, kRed, kOrange+2, kOrange};
  coloropt = coloroptrate;

  std::vector<TFile*> file;
  std::vector<TString> name;
  std::vector<TString> leg;

  file.clear();
  name.clear();
  leg.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  vector<double> subleadmupt_binsptgt33{30,33,36,39,42,45,48,51,54,57,60};
  int subleadmupt_nbinsptgt33 = subleadmupt_binsptgt33.size()-1;
  seltext[0] = "N#mu#geq2, p_{T}#geq33 GeV";
  seltext[1] = "#mu |#eta|<2.5, d_{0}>0.01 cm";

  file.push_back(datafile);
  name.push_back("sel3recomu");
  leg.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(sig3cmfile);
  name.push_back("sel3recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("sel3recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("sel3recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("sel3recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
 
  legendEntries = leg;  
  comparesamevariable(file, name, "subleadpt", subleadmupt_nbinsptgt33, &subleadmupt_binsptgt33[0], true, false, false, (float []){5e-2,5}, (float []){0.575,0.6,0.825,0.99}, (float []){33,2}, true, "#mu_{2} p_{T} [GeV]", "normalized events / 3 GeV", false, "sel3recomu_subleadpt");

  seltext[0] = "N#mu#geq2, p_{T}#geq33 GeV";
  seltext[1] = "#mu |#eta|<2.5, d_{0}>0.01 cm";
  makeratehist({"sel3recomu","sel3recomu","sel3recomu","sel3recomu","sel3recomu"}, "subleadpt", subleadmupt_nbinsptgt33, &subleadmupt_binsptgt33[0], false, false, (float []){0.5,0.65,0.7,0.975}, (float []){31,0.7}, (float []){0,0.85}, 0.875, "#mu_{2} p_{T} [GeV]");

  vector<double> subleadmupt_binspt{14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,33,36,39,42,45,48,51,54,57,60};
  int subleadmupt_nbinspt = subleadmupt_binspt.size()-1;
  seltext[0] = "N#mu#geq2, p_{T}#geq16 GeV";
  seltext[1] = "#mu |#eta|<2.5, d_{0}>0.01 cm";
  makeratehist({"sel31recomu","sel31recomu","sel31recomu","sel31recomu","sel31recomu"}, "subleadpt", subleadmupt_nbinspt, &subleadmupt_binspt[0], false, false, (float []){0.5,0.65,0.7,0.975}, (float []){35,1.5}, (float []){0,5}, 0.875, "#mu_{2} p_{T} [GeV]");

  seltext[0] = "N#mu#geq2, p_{T}#geq16 GeV";
  seltext[1] = "#mu |#eta|<2.5";
  makeratehist({"sel2recomu","sel2recomu","sel2recomu","sel2recomu","sel2recomu"}, "subleadpt", subleadmupt_nbinspt, &subleadmupt_binspt[0], false, false, (float []){0.5,0.625,0.7,0.95}, (float []){40,12}, (float []){0,33}, 0.875, "#mu_{2} p_{T} [GeV]");

  file.clear();
  name.clear();
  leg.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(datafile);
  name.push_back("sel2recomumu");
  leg.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(sig3cmfile);
  name.push_back("sel2recomumu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist  same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("sel2recomumu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist  same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("sel2recomumu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist  same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("sel2recomumu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist  same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
 
  legendEntries = leg;  
  vector<double> leadsubleadmumu_binsM{0, 4, 8, 16, 28, 34, 42, 48, 56, 64, 72, 80, 84, 88, 90, 91, 92, 94, 96, 100, 108, 116, 125};
  int leadsubleadmumu_nbinsM = leadsubleadmumu_binsM.size()-1;
  comparesamevariable(file, name, "leadsubleadM", leadsubleadmumu_nbinsM, &leadsubleadmumu_binsM[0], true, false, false, (float []){7e-3,5}, (float []){0.575,0.65,0.825,0.99}, (float []){4,0.7}, true, "M(#mu_{1},#mu_{2}) [GeV]", "normalized events / GeV", true, "sel2recomumu_leadsublead_M");

  vector<double> leadsubleadmumu_binsdR{0.0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 4};
  int leadsubleadmumu_nbinsdR = leadsubleadmumu_binsdR.size()-1;
  comparesamevariable(file, name, "leadsubleaddR", leadsubleadmumu_nbinsdR, &leadsubleadmumu_binsdR[0], true, false, false, (float []){2e-3,5}, (float []){0.125,0.65,0.375,0.99}, (float []){2.5,0.7}, true, "#DeltaR(#mu_{1},#mu_{2})", "normalized events / unit", true, "sel2recomumu_leadsublead_dR");

  seltext[0] = "N#mu#geq2, p_{T}#geq16 GeV, |#eta|<2.5";
  seltext[1] = "M<80 OR M>100, dR>1";
  makeratehist({"sel100recomu","sel100recomu","sel100recomu","sel100recomu","sel100recomu"}, "subleadpt", subleadmupt_nbinspt, &subleadmupt_binspt[0], false, false, (float []){0.5,0.66,0.7,0.99}, (float []){37,0.4}, (float []){0,0.99}, 0.875, "#mu_{2} p_{T} [GeV]");

  seltext[0] = "N#mu#geq2, p_{T}#geq16 GeV";
  seltext[1] = "#mu |#eta|<2.5";

  file.clear();
  name.clear();
  leg.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(datafile);
  name.push_back("sel2recomu");
  leg.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(sig3cmfile);
  name.push_back("sel2recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("sel2recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("sel2recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("sel2recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
 
  legendEntries = leg;  

  vector<double> mu_binslog10dxy{-4.0, -3.6, -3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0};
  int mu_nbinslog10dxy = mu_binslog10dxy.size()-1;
  comparesamevariable(file, name, "subleadlog10dxy", mu_nbinslog10dxy, &mu_binslog10dxy[0], false, false, false, (float []){0,0.52}, (float []){0.55,0.68,0.825,0.99}, (float []){-3.6, 0.35}, true, "#mu_{2} log_{10}d_{0} [log_{10} cm]", "normalized events / 0.4 log_{10} cm", false, "sel2recomu_subleadmu_log10dxy");

  seltext[0] = "N#mu#geq2, p_{T}#geq16 GeV";
  seltext[1] = "#mu |#eta|<2.5, |d_{0}|>0.01 cm";

  file.clear();
  name.clear();
  leg.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(datafile);
  name.push_back("sel31recomu");
  leg.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(sig3cmfile);
  name.push_back("sel31recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("sel31recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("sel31recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("sel31recomu");
  leg.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
 
  legendEntries = leg;  

  vector<double> mu_binslog10dxysig{-0.6,-0.2,0.2,0.6,1,1.4,1.8,2.2,2.6,3,3.4,3.8,4.2};
  int mu_nbinslog10dxysig = mu_binslog10dxysig.size()-1;
  comparesamevariable(file, name, "subleadlog10dxysig", mu_nbinslog10dxysig, &mu_binslog10dxysig[0], false, false, false, (float []){0,0.52}, (float []){0.125,0.65,0.375,0.99}, (float []){2.2,0.425}, true, "#mu_{2} log_{10}d_{0} sig. [TODO]", "normalized events / 0.4 units", false, "sel31recomu_subleadmu_log10dxysig");

  seltext[0] = "N#mu#geq2, p_{T}#geq16 GeV, |#eta|<2.5";
  seltext[1] = "d_{0}>0.01 cm, d_{0} sig.>1.4";
  makeratehist({"sel30recomu","sel30recomu","sel30recomu","sel30recomu","sel30recomu"}, "subleadpt", subleadmupt_nbinspt, &subleadmupt_binspt[0], false, false, (float []){0.5,0.66,0.7,0.99}, (float []){37,0.4}, (float []){0,1.39}, 0.875, "#mu_{2} p_{T} [GeV]");

  return -1;
}
