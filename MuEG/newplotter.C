#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 1.92/(19+11+10+5); // Data rate - 1.92 Hz, 22 events from parent trigger selection
double sig3cmsf = 1.0/39;
//double sig30cmsf = 1.0/28;
double sig1msf = 1.0/12;
//double sig3msf = 1.0/7;
TString cutdeets = "Cut details";
TFile* datafile = TFile::Open("hists_Efmrl.root","READ");
//TFile* dyfile = TFile::Open("hists_DY.root","READ");
TFile* sig3cmfile = TFile::Open("hists_M200dM20ctau3cm.root","READ");
//TFile* sig30cmfile = TFile::Open("hists_M200dM20ctau30cm.root","READ");
TFile* sig1mfile = TFile::Open("hists_M200dM20ctau1m.root","READ");
//TFile* sig3mfile = TFile::Open("hists_M200dM20ctau3m.root","READ");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};

int comparemultihist(TString cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool isMConly=false, bool overflow=false, float legPos[]=(float []){0.7,0.75,0.95,1}, float yrange[]=(float []){0.1,100}, bool normalized=false) {

  std::vector<TString> legNam;
  if(!isMConly) legNam.push_back("2018 data");
  legNam.push_back("c#tau = 3 cm");
  //legNam.push_back("c#tau = 30 cm");
  //legNam.push_back("c#tau = 1 m");
  //legNam.push_back("c#tau = 3 m");

  TString histname = cutname+"_"+var;
  TH1F* datahist;
  if(!isMConly) datahist = (TH1F*) datafile->Get(histname);
  TH1F* sig3cmhist = (TH1F*) sig3cmfile->Get(histname);
  //TH1F* sig30cmhist = (TH1F*) sig30cmfile->Get(histname);
  //TH1F* sig1mhist = (TH1F*) sig1mfile->Get(histname);
  //TH1F* sig3mhist = (TH1F*) sig3mfile->Get(histname);

  // Get the title from histogram title
  TString xtitle = sig3cmhist->GetTitle();

  std::vector<TH1F*> allhists;
  if(!isMConly) allhists.push_back(datahist);
  allhists.push_back(sig3cmhist);
  //allhists.push_back(sig30cmhist);
  //allhists.push_back(sig1mhist);
  //allhists.push_back(sig3mhist);

  if(rebin==-1) rebin = 1;
  xbinlow = xbinlow==-1?1:xbinlow;
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
  c1 = enhance_plotter(allhists, legNam, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,normalized);
  c1->SaveAs("./dirplots/"+cutname+"/"+cutname+"_"+var+".png");
  
  return -1;
}

// Comparison of the type cross-check between two histogram - filled v hollow
int crossChecktwohist(TFile* file, vector<TString> cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, float legPos[]=(float []){0.7,0.75,0.95,1}, float yrange[]=(float []){0.1,100}) {

  std::vector<TString> legNam;
  legNam.push_back("L3#mu filter - DoubleMu33Displaced");
  legNam.push_back("L3#mu filter, p_{T}>16, |#eta|<2.5, d_{0}>0.1mm ");

  TH1F* histfilled = (TH1F*) file->Get(cutname[0]+"_"+var);
  TH1F* histhollow = (TH1F*) file->Get(cutname[1]+"_"+var);

  std::vector<TH1F*> allhists;
  allhists.push_back(histfilled);
  allhists.push_back(histhollow);

  // Get the title from histogram title
  TString xtitle = histfilled->GetTitle();

  if(rebin==-1) rebin = 1;
  xbinlow = xbinlow==-1?0:xbinlow;
  xbinhigh = xbinhigh==-1?allhists[0]->GetNbinsX()+1:xbinhigh;
  xbinlow = xbinlow/rebin;
  xbinhigh = xbinhigh/rebin;

  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {
    if(rebin!=1) allhists[histctr]->Rebin(rebin);

    // Make changes to sig and bkg to enable good basic plotting
    double err = 0.0;
    allhists[histctr]->SetBinContent(xbinlow,allhists[histctr]->IntegralAndError(0,xbinlow,err));
    allhists[histctr]->SetBinError(xbinlow,err);
    err = 0.0;
    allhists[histctr]->SetBinContent(xbinhigh,allhists[histctr]->IntegralAndError(xbinhigh,allhists[histctr]->GetNbinsX()+1,err));
    allhists[histctr]->SetBinError(xbinhigh,err);

    allhists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allhists[histctr]->GetXaxis()->SetTitle(xtitle);
    allhists[histctr]->GetYaxis()->SetTitle("number of events");
    
    allhists[histctr]->SetTitle("");
  }

  allhists[0]->SetFillColor(kBlue);
  allhists[1]->SetLineColor(kRed);
  
  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legNam, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,false);
  c1->SaveAs("./dirplots/"+cutname[0]+"_"+cutname[1]+"/"+cutname[0]+"_"+cutname[1]+"_"+var+".png");

  return -1;
}

int makeratehist(TString cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool overflow=false, float legPos[]=(float []){0.7,0.75,0.95,1}, float seltextpos[]=(float []){0.1,1}, float yrange[]=(float []){0.1,1}, float drawsignalline=-1.0) {

  std::vector<TString> legNam;
  legNam.push_back("2018 data");
  legNam.push_back("c#tau = 3 cm");
  //legNam.push_back("c#tau = 30 cm");
  //legNam.push_back("c#tau = 1 m");
  //legNam.push_back("c#tau = 3 m");

  TString histname = cutname+"_"+var;
  TH1F* datahist;
  datahist = (TH1F*) datafile->Get(histname);
  datahist->Scale(datasf);
  TH1F* sig3cmhist = (TH1F*) sig3cmfile->Get(histname);
  sig3cmhist->Scale(sig3cmsf);
  //TH1F* sig30cmhist = (TH1F*) sig30cmfile->Get(histname);
  //sig30cmhist->Scale(sig30cmsf);
  //TH1F* sig1mhist = (TH1F*) sig1mfile->Get(histname);
  //sig1mhist->Scale(sig1msf);
  //TH1F* sig3mhist = (TH1F*) sig3mfile->Get(histname);
  //sig3mhist->Scale(sig3msf);

  // Get the title from histogram title
  TString xtitle = sig3cmhist->GetTitle();

  std::vector<TH1F*> allhists;
  std::vector<TH1F*> allratehists;
  allhists.push_back(datahist);
  allhists.push_back(sig3cmhist);
  //allhists.push_back(sig30cmhist);
  //allhists.push_back(sig1mhist);
  //allhists.push_back(sig3mhist);

  double axisscale = 0.0;
  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {

    allratehists.push_back((TH1F*) allhists[histctr]->Clone());
    for(unsigned int bincnt=0; bincnt<allhists[histctr]->GetNbinsX(); bincnt++) {
      double err = 0;
      allratehists[histctr]->SetBinContent(bincnt, allhists[histctr]->IntegralAndError(bincnt,allhists[histctr]->GetNbinsX()+1,err));
      allratehists[histctr]->SetBinError(bincnt, err);
    }
    allratehists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allratehists[histctr]->GetXaxis()->SetTitle(xtitle);
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
  axis->SetLabelSize(0.06);
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
  sel.SetTextSize(0.05);
  sel.DrawLatex(seltextpos[0], seltextpos[1]+0.1*(yrange[1]-yrange[0]), seltext[0]);
  sel.DrawLatex(seltextpos[0], seltextpos[1], seltext[1]);
  
  c1->SaveAs("./dirplots/"+cutname+"/"+cutname+"_"+var+"_ratehist.png");
  
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
  coloropt = coloroptrate;
  seltext[0] = "N#mu#geq1, p_{T}#geq38 GeV, |#eta|<2.5, d_{0}>0.01 cm";
  seltext[1] = "Ne/#gamma#geq1, p_{T}>38 GeV, |#eta|<2.65";
  //makeratehist("sel3recoeg", "egpt", 66, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){60,0.8}, (float []){0.1,2.05}, 0.875);

  seltext[0] = "N#mu#geq1, p_{T}#geq16 GeV, |#eta|<2.5";
  seltext[1] = "Ne/#gamma#geq1, p_{T}>15 GeV, |#eta|<2.65";
  //makeratehist("sel10recoeg", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){45,5}, (float []){0,220}, 0.875);

  seltext[0] = "N#mu#geq1, p_{T}#geq16 GeV, |#eta|<2.5";
  seltext[1] = "Unseeded: Ne/#gamma#geq1, p_{T}>15 GeV, |#eta|<2.65";
  //makeratehist("sel10usrecoegus", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){45,5}, (float []){0,400}, 0.875);

  seltext[0] = "N#mu#geq1, p_{T}>20 GeV, |#eta|<2.5, d_{0} sig.>1";
  seltext[1] = "Ne/#gamma#geq1, p_{T}>20 GeV, |#eta|<2.5";
  //makeratehist("sel20recoegus", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){45,5}, (float []){0,220}, 0.875);

  seltext[0] = "N#mu#geq1, p_{T}>15 GeV, |#eta|<2.5, d_{0} sig.>4";
  seltext[1] = "Ne/#gamma#geq1, CaloIdL";
  //makeratehist("sel80recoeg", "egpt", 50, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){45,5}, (float []){0,50}, 0.875);

  seltext[0] = "N#mu#geq1, p_{T}>15 GeV, |#eta|<2.5, d_{0} sig.>4";
  seltext[1] = "Unseeded: Ne/#gamma#geq1, CaloIdL";
  //makeratehist("sel80usrecoegus", "egpt", 50, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){45,5}, (float []){0,250}, 0.875);

  std::vector<int> coloroptschemegenmch{kGreen+2, kGreen-6};
  coloropt = coloroptschemegenmch;
  std::vector<TFile*> filesel;
  std::vector<TString> namesel;
  std::vector<TString> legsel;
  //filesel.push_back(datafile);
  //namesel.push_back("sel10usrecoebus");
  //legsel.push_back("data");
  //filesel.push_back(dyfile);
  //namesel.push_back("sel10recoeb");
  //legsel.push_back("DY");
  //filesel.push_back(sig3cmfile);
  //namesel.push_back("sel10usrecoebus");
  //legsel.push_back("signal 3 cm");
  //filesel.push_back(sig30cmfile);
  //namesel.push_back("sel20recoeb");
  //legsel.push_back("signal 30 cm");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel10usrecoebusgenmch");
  legsel.push_back("3cm gen matched");
  filesel.push_back(sig1mfile);
  namesel.push_back("sel10usrecoebusgenmch");
  legsel.push_back("1m gen matched");
  //filesel.push_back(sig3mfile);
  //namesel.push_back("sel20recoeb");
  //legsel.push_back("signal 3 m");
  legendEntries = legsel;  
  //comparesamevariable(filesel, namesel, "egclustershape", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma cluster shape");
  //comparesamevariable(filesel, namesel, "in5x5clusshape", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5");
  //comparesamevariable(filesel, namesel, "in5x5noiseclnd", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5 (noise cleaned)");
  //comparesamevariable(filesel, namesel, "hovere", 1, 1000, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H / GeV");
  //comparesamevariable(filesel, namesel, "hovereoverpt", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/p_{T}");
  //comparesamevariable(filesel, namesel, "hovereoversupcluse", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/E");
  //comparesamevariable(filesel, namesel, "seedclustime", 8000, 13000, 200, true, true, true, (float []){1e-5,10}, (float []){0.11,0.7,0.36,0.99}, true, "barrel ecal seed time / ns");
  //comparesamevariable(filesel, namesel, "pixelmchvar_s2", 40, 200, 2, true, true, true, (float []){1e-5,10}, (float []){0.6,0.7,0.85,0.99}, true, "pixel match");

  filesel.clear();
  namesel.clear();
  legsel.clear();
  filesel.push_back(datafile);
  namesel.push_back("sel10usrecoeeus");
  legsel.push_back("data");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel10usrecoeeus");
  legsel.push_back("signal 3 cm");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel10usrecoeeusgenmch");
  legsel.push_back("MC gen matched");
  legendEntries = legsel;  
  //comparesamevariable(filesel, namesel, "egclustershape", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma cluster shape");
  //comparesamevariable(filesel, namesel, "in5x5clusshape", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5");
  //comparesamevariable(filesel, namesel, "in5x5noiseclnd", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5 (noise cleaned)");
  //comparesamevariable(filesel, namesel, "hovere", 1, 1000, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H / GeV");
  //comparesamevariable(filesel, namesel, "hovereoverpt", -1, 100, 4, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/p_{T}");
  //comparesamevariable(filesel, namesel, "hovereoversupcluse", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/E");
  //comparesamevariable(filesel, namesel, "seedclustime", 8000, 13000, 200, true, true, true, (float []){1e-5,10}, (float []){0.11,0.7,0.36,0.99}, true, "barrel ecal seed time / ns");

  // Unseeded sel 11
  filesel.clear();
  namesel.clear();
  legsel.clear();
  filesel.push_back(datafile);
  namesel.push_back("sel11usrecoebus");
  legsel.push_back("data");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel11usrecoebus");
  legsel.push_back("signal 3 cm");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel11usrecoebusgenmch");
  legsel.push_back("MC gen matched");
  legendEntries = legsel;  
  //comparesamevariable(filesel, namesel, "egclustershape", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma cluster shape");
  //comparesamevariable(filesel, namesel, "in5x5clusshape", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5");
  //comparesamevariable(filesel, namesel, "in5x5noiseclnd", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5 (noise cleaned)");
  //comparesamevariable(filesel, namesel, "hovere", 1, 1000, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H / GeV");
  //comparesamevariable(filesel, namesel, "hovereoverpt", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/p_{T}");
  //comparesamevariable(filesel, namesel, "hovereoversupcluse", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/E");
  //comparesamevariable(filesel, namesel, "seedclustime", 8000, 13000, 200, true, true, true, (float []){1e-5,10}, (float []){0.11,0.7,0.36,0.99}, true, "barrel ecal seed time / ns");

  filesel.clear();
  namesel.clear();
  legsel.clear();
  filesel.push_back(datafile);
  namesel.push_back("sel11usrecoeeus");
  legsel.push_back("data");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel11usrecoeeus");
  legsel.push_back("signal 3 cm");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel11usrecoeeusgenmch");
  legsel.push_back("MC gen matched");
  legendEntries = legsel;  
  //comparesamevariable(filesel, namesel, "egclustershape", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma cluster shape");
  //comparesamevariable(filesel, namesel, "in5x5clusshape", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5");
  //comparesamevariable(filesel, namesel, "in5x5noiseclnd", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5 (noise cleaned)");
  //comparesamevariable(filesel, namesel, "hovere", 1, 1000, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H / GeV");
  //comparesamevariable(filesel, namesel, "hovereoverpt", -1, 100, 4, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/p_{T}");
  //comparesamevariable(filesel, namesel, "hovereoversupcluse", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/E");
  //comparesamevariable(filesel, namesel, "seedclustime", 8000, 13000, 200, true, true, true, (float []){1e-5,10}, (float []){0.11,0.7,0.36,0.99}, true, "barrel ecal seed time / ns");

  seltext[0] = " |#eta|<2.5, #mu: p_{T}#geq16 GeV,";
  seltext[1] = "Unseeded e/#gamma: p_{T}#geq15 GeV, i#eta<0.012(0.03)";
  //makeratehist("sel11usrecoegus", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){45,5}, (float []){0,400}, 0.875);

  // Unseeded sel 12
  filesel.clear();
  namesel.clear();
  legsel.clear();
  filesel.push_back(datafile);
  namesel.push_back("sel12usrecoebus");
  legsel.push_back("data");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel12usrecoebus");
  legsel.push_back("signal 3 cm");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel12usrecoebusgenmch");
  legsel.push_back("MC gen matched");
  legendEntries = legsel;  
  //comparesamevariable(filesel, namesel, "egclustershape", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma cluster shape");
  //comparesamevariable(filesel, namesel, "in5x5clusshape", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5");
  //comparesamevariable(filesel, namesel, "in5x5noiseclnd", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5 (noise cleaned)");
  //comparesamevariable(filesel, namesel, "hovere", 1, 1000, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H / GeV");
  //comparesamevariable(filesel, namesel, "hovereoverpt", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/p_{T}");
  //comparesamevariable(filesel, namesel, "hovereoversupcluse", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/E");
  //comparesamevariable(filesel, namesel, "seedclustime", 8000, 13000, 200, true, true, true, (float []){1e-5,10}, (float []){0.11,0.7,0.36,0.99}, true, "barrel ecal seed time / ns");

  filesel.clear();
  namesel.clear();
  legsel.clear();
  filesel.push_back(datafile);
  namesel.push_back("sel12usrecoeeus");
  legsel.push_back("data");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel12usrecoeeus");
  legsel.push_back("signal 3 cm");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel12usrecoeeusgenmch");
  legsel.push_back("MC gen matched");
  legendEntries = legsel;  
  //comparesamevariable(filesel, namesel, "egclustershape", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma cluster shape");
  //comparesamevariable(filesel, namesel, "in5x5clusshape", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5");
  //comparesamevariable(filesel, namesel, "in5x5noiseclnd", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5 (noise cleaned)");
  //comparesamevariable(filesel, namesel, "hovere", 1, 1000, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H / GeV");
  //comparesamevariable(filesel, namesel, "hovereoverpt", -1, 100, 4, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/p_{T}");
  //comparesamevariable(filesel, namesel, "hovereoversupcluse", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/E");
  //comparesamevariable(filesel, namesel, "seedclustime", 8000, 13000, 200, true, true, true, (float []){1e-5,10}, (float []){0.11,0.7,0.36,0.99}, true, "barrel ecal seed time / ns");

  seltext[0] = " |#eta|<2.5, #mu: p_{T}#geq16 GeV,";
  seltext[1] = "Unseeded e/#gamma: p_{T}#geq15 GeV, i#eta<0.012(0.03), H/E<0.1";
  //makeratehist("sel12usrecoegus", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){45,5}, (float []){0,100}, 0.875);

  // Unseeded sel 13
  filesel.clear();
  namesel.clear();
  legsel.clear();
  filesel.push_back(datafile);
  namesel.push_back("sel13usrecoebus");
  legsel.push_back("data");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel13usrecoebus");
  legsel.push_back("signal 3 cm");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel13usrecoebusgenmch");
  legsel.push_back("MC gen matched");
  legendEntries = legsel;  
  //comparesamevariable(filesel, namesel, "egclustershape", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma cluster shape");
  //comparesamevariable(filesel, namesel, "in5x5clusshape", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5");
  //comparesamevariable(filesel, namesel, "in5x5noiseclnd", 1, 300, 5, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5 (noise cleaned)");
  //comparesamevariable(filesel, namesel, "hovere", 1, 1000, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H / GeV");
  //comparesamevariable(filesel, namesel, "hovereoverpt", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/p_{T}");
  //comparesamevariable(filesel, namesel, "hovereoversupcluse", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/E");
  //comparesamevariable(filesel, namesel, "ecalpfclustiso", 50, 800, 10, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "ecal iso. / GeV");
  //comparesamevariable(filesel, namesel, "ecalpfclustisoovere", 50, 550, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "ecal iso./E");
  //comparesamevariable(filesel, namesel, "seedclustime", 8000, 13000, 200, true, true, true, (float []){1e-5,10}, (float []){0.11,0.7,0.36,0.99}, true, "barrel ecal seed time / ns");

  filesel.clear();
  namesel.clear();
  legsel.clear();
  filesel.push_back(datafile);
  namesel.push_back("sel13usrecoeeus");
  legsel.push_back("data");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel13usrecoeeus");
  legsel.push_back("signal 3 cm");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel13usrecoeeusgenmch");
  legsel.push_back("MC gen matched");
  legendEntries = legsel;  
  //comparesamevariable(filesel, namesel, "egclustershape", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma cluster shape");
  //comparesamevariable(filesel, namesel, "in5x5clusshape", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5");
  //comparesamevariable(filesel, namesel, "in5x5noiseclnd", 1, 800, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma i#eta 5x5 (noise cleaned)");
  //comparesamevariable(filesel, namesel, "hovere", 1, 1000, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H / GeV");
  //comparesamevariable(filesel, namesel, "hovereoverpt", -1, 100, 4, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/p_{T}");
  //comparesamevariable(filesel, namesel, "hovereoversupcluse", -1, 40, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "e/#gamma H/E");
  //comparesamevariable(filesel, namesel, "ecalpfclustiso", 50, 800, 10, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "ecal iso. / GeV");
  //comparesamevariable(filesel, namesel, "ecalpfclustisoovere", 50, 550, 1, true, true, true, (float []){1e-3,1}, (float []){0.64,0.7,0.89,0.99}, true, "ecal iso./E");
  //comparesamevariable(filesel, namesel, "seedclustime", 8000, 13000, 200, true, true, true, (float []){1e-5,10}, (float []){0.11,0.7,0.36,0.99}, true, "barrel ecal seed time / ns");

  seltext[0] = " |#eta|<2.5, #mu: p_{T}#geq16 GeV,";
  seltext[1] = "Unseeded e/#gamma: p_{T}#geq15 GeV, i#eta<0.012(0.03), H/E<0.1";
  //makeratehist("sel13usrecoegus", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){45,5}, (float []){0,100}, 0.875);

  seltext[0] = " |#eta|<2.5, #mu: p_{T}#geq16 GeV,";
  seltext[1] = "L1Seeded e/#gamma: p_{T}#geq15 GeV, i#eta<0.012(0.03), H/E<0.1";
  //makeratehist("sel13recoeg", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){45,5}, (float []){0,100}, 0.875);

  seltext[0] = " |#eta|<2.5, #mu: p_{T}#geq16 GeV,";
  seltext[1] = "L1Seeded e/#gamma: #splitline{p_{T}#geq15 GeV, i#eta<0.012(0.03),}{H/E<0.1, ecal iso./E<0.15(0.1)}";
  //makeratehist("sel14recoeg", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){30,15}, (float []){0,30}, 0.875);

  seltext[0] = " |#eta|<2.5, #mu: p_{T}#geq16 GeV,";
  seltext[1] = "Unseeded e/#gamma: #splitline{p_{T}#geq15 GeV, i#eta<0.012(0.03),}{H/E<0.1, ecal iso./E<0.15(0.1)}";
  //makeratehist("sel14usrecoegus", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){30,15}, (float []){0,30}, 0.875);

  seltext[0] = " |#eta|<2.5, #mu: p_{T}#geq19 GeV,";
  seltext[1] = "L1Seeded e/#gamma: #splitline{p_{T}#geq19 GeV, i#eta<0.012(0.03),}{H/E<0.1, ecal iso./E<0.15(0.1)}";
  //makeratehist("sel15recoeg", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){30,15}, (float []){0,30}, 0.875);

  seltext[0] = " |#eta|<2.5, #mu: p_{T}#geq19 GeV,";
  seltext[1] = "Unseeded e/#gamma: #splitline{p_{T}#geq19 GeV, i#eta<0.012(0.03),}{H/E<0.1, ecal iso./E<0.15(0.1)}";
  //makeratehist("sel15usrecoegus", "egpt", 56, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){30,15}, (float []){0,30}, 0.875);

  // Unseeded sel 14
  filesel.clear();
  namesel.clear();
  legsel.clear();
  filesel.push_back(datafile);
  namesel.push_back("sel14usrecoebus");
  legsel.push_back("data");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel14usrecoebus");
  legsel.push_back("signal 3 cm");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel14usrecoebusgenmch");
  legsel.push_back("MC gen matched");
  legendEntries = legsel;  
  //comparesamevariable(filesel, namesel, "mhits", -1, -1, 1, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "missing hits");
  //comparesamevariable(filesel, namesel, "hcalpfclustiso", -1, -1, 1, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "hcal iso. / GeV");

  filesel.clear();
  namesel.clear();
  legsel.clear();
  filesel.push_back(datafile);
  namesel.push_back("sel14usrecoeeus");
  legsel.push_back("data");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel14usrecoeeus");
  legsel.push_back("signal 3 cm");
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel14usrecoeeusgenmch");
  legsel.push_back("MC gen matched");
  legendEntries = legsel;  
  //comparesamevariable(filesel, namesel, "mhits", -1, -1, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "missing hits");
  //comparesamevariable(filesel, namesel, "hcalpfclustiso", -1, -1, 10, true, true, true, (float []){1e-3,10}, (float []){0.64,0.7,0.89,0.99}, true, "hcal iso. / GeV");

  filesel.clear();
  namesel.clear();
  legsel.clear();
  filesel.push_back(sig3cmfile);
  namesel.push_back("sel15gen");
  legsel.push_back("MC gen");
  filesel.push_back(sig3cmfile);
  namesel.push_back("gen");
  legsel.push_back("trigger selected");
  legendEntries = legsel;  
  comparesamevariable(filesel, namesel, "recomchelpt", 50, 120, 2, true, true, true, (float []){1,4e2}, (float []){0.64,0.7,0.89,0.99}, false, "e p_{T} / GeV");
  comparesamevariable(filesel, namesel, "recomchellog10d0", 250, 750, 10, true, true, true, (float []){1,9e2}, (float []){0.64,0.7,0.89,0.99}, false, "e log_{10}d_{0} / log_{10}cm");

  //std::vector<int> coloroptschemetim{kBlack, /*kBlue, */kRed+2, kRed-3, /*kRed-7, kRed-9*/};
  return -1;
}
