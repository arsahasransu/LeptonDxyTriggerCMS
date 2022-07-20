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

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};
std::vector<TString> histtype{"p e1", "hist same"};
std::vector<int> markerstyle{20, 24};
std::vector<int> markersize{10, 10};
std::vector<TString> legendmarkerstyle{"lep", "l"};
std::vector<double> scale{1, 1};

int comparesamevariable(std::vector<TFile*> file, std::vector<TString> cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool underflow=false, bool overflow=false, float yrange[]=(float []){0.1,100}, float legPos[]=(float []){0.7,0.75,0.95,1}, bool normalize=true, TString xaxistitle="xaxis", TString yaxistitle="yaxis") {

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
    if(normalize) allhists[histctr]->GetYaxis()->SetTitle("normalized number of events");
    else allhists[histctr]->GetYaxis()->SetTitle(yaxistitle);

    allhists[histctr]->SetTitle("");

    allhists[histctr]->SetLineWidth(2);
    allhists[histctr]->SetLineColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  c1 = fancy_enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,normalize, histtype, markerstyle, markersize, legendmarkerstyle, scale);
  c1->SaveAs("./dirplots/"+foldername+"/"+var+".png");

  return -1;
}

int plotter() {

  std::vector<int> coloroptrate{1, kRed+2, kRed-3, kRed-7, kRed-9};

  std::vector<TFile*> file;
  std::vector<TString> name;
  std::vector<TString> legend;

  // Study the unseeded egamma object properties
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
  
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(sigfile);
  name.push_back("genptgt10barsel_dieg33mchgeneg");
  legend.push_back("Run 2");
  coloropt.push_back(kRed);
  histtype.push_back("p e1");
  markerstyle.push_back(24);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1.0/139.0);

  file.push_back(sigfile);
  name.push_back("genptgt10barsel_time1nsAsminlt0p16mchgeneg");
  legend.push_back("Run 3");
  coloropt.push_back(kGreen+2);
  histtype.push_back("p e1 same");
  markerstyle.push_back(26);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(6.0/139.0);

  legendEntries = legend;
  comparesamevariable(file, name, "log10d0", 250, -1, 20, true, true, true, (float []){0.009,20}, (float []){0.7,0.6,0.94,0.95}, false, "electron log_{10} d_{0} [log_{10} cm]","number of events / Run 2 integral");
  
  return -1;
}
