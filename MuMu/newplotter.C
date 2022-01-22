#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 0.5/(24+22+20+21); // Data rate - 0.5 Hz, events from parent trigger selection
double sig3cmsf = 1.0/29;
double sig30cmsf = 1.0/10;
double sig1msf = 1.0/19;
double sig3msf = 1.0/7;
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
int crossChecktwohist(std::vector<TFile*> file, std::vector<TString> cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool underflow=false, bool overflow=false, float yrange[]=(float []){0.1,100}, float legPos[]=(float []){0.7,0.75,0.95,1}, bool normalize=true, TString xaxistitle="xaxis") {
  
  if(file.size()<2 || cutname.size()<2) {
    cout<<"Need atleast 2 histograms to cross-check"<<endl;
    return -1;
  }
  
  TH1F* histfilled = (TH1F*) file[0]->Get(cutname[0]+"_"+var);
  TH1F* histhollow = (TH1F*) file[1]->Get(cutname[1]+"_"+var);
  
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
  
  return -1;
}

int makeratehist(TString cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool overflow=false, float legPos[]=(float []){0.7,0.75,0.95,1}, float seltextpos[]=(float []){0.1,1}, float yrange[]=(float []){0.1,1}, float drawsignalline=-1.0) {

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

  // Get the title from histogram title
  TString xtitle = sig3cmhist->GetTitle();

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

int comparesamevariable(std::vector<TFile*> file, std::vector<TString> cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool underflow=false, bool overflow=false, float yrange[]=(float []){0.1,100}, float legPos[]=(float []){0.7,0.75,0.95,1}, bool normalize=true, TString xaxistitle="xaxis", TString outputpathprefix="") {

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
  c1->SaveAs("./dirplots/"+outputpathprefix+foldername+"/"+var+".png");

  return -1;
}

int newplotter() {

  std::vector<int> coloroptrate{1, kRed+2, kRed-3, kRed-7, kRed-9};
  coloropt = coloroptrate;
  seltext[0] = "N#mu#geq2, p_{T}#geq33 GeV";
  seltext[1] = "#mu |#eta|<2.5, d_{0}>0.01 cm";
  //makeratehist("sel3recomu", "subleadpt", 66, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,0.35}, (float []){0,0.8}, 0.875);

  seltext[0] = "N#mu#geq2, p_{T}#geq16 GeV";
  seltext[1] = "#mu |#eta|<2.5";
  //makeratehist("sel2recomu", "subleadpt", 66, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,45}, 0.875);

  seltext[0] = "N#mu#geq2, p_{T}#geq20 GeV, |#eta|<2.5";
  seltext[1] = "#mu |d_{0}|>0.1 mm, |d_{0}| sig.>1";
  //makeratehist("sel100recomu", "subleadpt", 66, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,0.3}, (float []){0,0.96}, 0.875);

  seltext[0] = "N#mu#geq2, p_{T}#geq15 GeV, |#eta|<2.5";
  seltext[1] = "#mu d_{0} sig.>3.98";
  //makeratehist("sel40recomu", "subleadpt", 60, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,1.5}, (float []){0,5}, 0.875);

  std::vector<int> coloroptschemerecoxcheckfilt16{kBlue, kRed};
  coloropt = coloroptschemerecoxcheckfilt16;
  std::vector<TFile*> filerecoxcheckfilt16;
  std::vector<TString> namerecoxcheckfilt16;
  std::vector<TString> legrecoxcheckfilt16;
  filerecoxcheckfilt16.push_back(sig3cmfile);
  namerecoxcheckfilt16.push_back("filt16mu");
  legrecoxcheckfilt16.push_back("signal 3cm");
  filerecoxcheckfilt16.push_back(sig3cmfile);
  namerecoxcheckfilt16.push_back("sel2recomu");
  legrecoxcheckfilt16.push_back("signal 3cm");
  legendEntries = legrecoxcheckfilt16;  
  //crossChecktwohist(filerecoxcheckfilt16, namerecoxcheckfilt16, "allpt", 60, 210, 1, true, true, true, (float []){1e-1,1e3}, (float []){0.6,0.7,0.85,0.99}, false, "#mu p_{T} / GeV");
  //crossChecktwohist(filerecoxcheckfilt16, namerecoxcheckfilt16, "alleta", -1, -1, -1, false, true, true, (float []){0,130}, (float []){0.6,0.7,0.85,0.99}, false, "#mu #eta");
  //crossChecktwohist(filerecoxcheckfilt16, namerecoxcheckfilt16, "allphi", -1, -1, -1, false, true, true, (float []){0,70}, (float []){0.6,0.7,0.85,0.99}, false, "#mu #phi");

  std::vector<int> coloroptschemerecoxcheckfilt33{kBlue, kRed};
  coloropt = coloroptschemerecoxcheckfilt33;
  std::vector<TFile*> filerecoxcheckfilt33;
  std::vector<TString> namerecoxcheckfilt33;
  std::vector<TString> legrecoxcheckfilt33;
  filerecoxcheckfilt33.push_back(sig3cmfile);
  namerecoxcheckfilt33.push_back("filt33mu");
  legrecoxcheckfilt33.push_back("signal 3cm");
  filerecoxcheckfilt33.push_back(sig3cmfile);
  namerecoxcheckfilt33.push_back("sel3recomu");
  legrecoxcheckfilt33.push_back("signal 3cm");
  legendEntries = legrecoxcheckfilt33;  
  //crossChecktwohist(filerecoxcheckfilt33, namerecoxcheckfilt33, "allpt", 60, 210, 1, true, true, true, (float []){5e-1,1e2}, (float []){0.6,0.7,0.85,0.99}, false, "#mu p_{T} / GeV");
  //crossChecktwohist(filerecoxcheckfilt33, namerecoxcheckfilt33, "alleta", -1, -1, -1, false, true, true, (float []){0,40}, (float []){0.6,0.7,0.85,0.99}, false, "#mu #eta");
  //crossChecktwohist(filerecoxcheckfilt33, namerecoxcheckfilt33, "allphi", -1, -1, -1, false, true, true, (float []){0,20}, (float []){0.6,0.7,0.85,0.99}, false, "#mu #phi");

  std::vector<int> coloroptschemesel2{kRed-9, kRed-7, kRed-3, kRed+2, kBlue};
  coloropt = coloroptschemesel2;
  std::vector<TFile*> filesel2;
  std::vector<TString> namesel2;
  std::vector<TString> legsel2;
  filesel2.push_back(sig3mfile);
  namesel2.push_back("sel2recomu");
  legsel2.push_back("signal 3m");
  filesel2.push_back(sig1mfile);
  namesel2.push_back("sel2recomu");
  legsel2.push_back("signal 1m");
  filesel2.push_back(sig30cmfile);
  namesel2.push_back("sel2recomu");
  legsel2.push_back("signal 30cm");
  filesel2.push_back(sig3cmfile);
  namesel2.push_back("sel2recomu");
  legsel2.push_back("signal 3cm");
  filesel2.push_back(datafile);
  namesel2.push_back("sel2recomu");
  legsel2.push_back("data");
  legendEntries = legsel2;  
  //comparesamevariable(filesel2, namesel2, "allpt", 60, 130, 2, true, true, true, (float []){1e-3,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu p_{T} / GeV");
  //comparesamevariable(filesel2, namesel2, "alleta", -1, -1, -1, true, true, true, (float []){1e-3,1e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #eta");
  //comparesamevariable(filesel2, namesel2, "allphi", -1, -1, 3, true, true, true, (float []){1e-2,2e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #phi");
  //comparesamevariable(filesel2, namesel2, "alldxy", 11000, 13000, 100, true, true, true, (float []){9e-4,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu d_{0} / cm");
  //comparesamevariable(filesel2, namesel2, "alllog10dxy", -1, 700, 20, true, true, true, (float []){1e-3,1}, (float []){0.11,0.7,0.36,0.99}, true, "#mu log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(filesel2, namesel2, "alldxysig", -1, 40000, 1, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu sigificance d_{0}");
  //comparesamevariable(filesel2, namesel2, "alllog10dxysig", 180, 900, 20, true, true, true, (float []){1e-2,3e-1}, (float []){0.7,0.7,0.85,0.99}, true, "#mu log_{10}(significance d_{0})");
  //comparesamevariable(filesel2, namesel2, "leadpt", 60, 130, 2, true, true, true, (float []){1e-3,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu p_{T} / GeV");
  //comparesamevariable(filesel2, namesel2, "leadeta", -1, -1, -1, true, true, true, (float []){1e-3,1e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #eta");
  //comparesamevariable(filesel2, namesel2, "leadphi", -1, -1, 3, true, true, true, (float []){1e-2,2e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #phi");
  //comparesamevariable(filesel2, namesel2, "leaddxy", 11000, 13000, 100, true, true, true, (float []){9e-4,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu d_{0} / cm");
  //comparesamevariable(filesel2, namesel2, "leadlog10dxy", -1, 700, 20, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(filesel2, namesel2, "leaddxysig", -1, 24000, 100, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu sigificance d_{0}");
  //comparesamevariable(filesel2, namesel2, "leadlog10dxysig", 180, 800, 20, true, true, true, (float []){1e-2,3e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu log_{10}(significance d_{0})");
  //comparesamevariable(filesel2, namesel2, "subleadpt", 60, 130, 2, true, true, true, (float []){1e-3,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu p_{T} / GeV");
  //comparesamevariable(filesel2, namesel2, "subleadeta", -1, -1, -1, true, true, true, (float []){1e-3,1e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #eta");
  //comparesamevariable(filesel2, namesel2, "subleadphi", -1, -1, 3, true, true, true, (float []){1e-2,2e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #phi");
  //comparesamevariable(filesel2, namesel2, "subleaddxy", 11000, 13000, 100, true, true, true, (float []){9e-4,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu d_{0} / cm");
  //comparesamevariable(filesel2, namesel2, "subleadlog10dxy", -1, 700, 20, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(filesel2, namesel2, "subleaddxysig", -1, -1, 2000, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu sigificance d_{0}");
  //comparesamevariable(filesel2, namesel2, "subleadlog10dxysig", 180, 800, 20, true, true, true, (float []){1e-2,3e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu log_{10}(significance d_{0})");
  filesel2.clear();
  namesel2.clear();
  legsel2.clear();
  filesel2.push_back(sig3mfile);
  namesel2.push_back("sel2recomumu");
  legsel2.push_back("signal 3m");
  filesel2.push_back(sig1mfile);
  namesel2.push_back("sel2recomumu");
  legsel2.push_back("signal 1m");
  filesel2.push_back(sig30cmfile);
  namesel2.push_back("sel2recomumu");
  legsel2.push_back("signal 30cm");
  filesel2.push_back(sig3cmfile);
  namesel2.push_back("sel2recomumu");
  legsel2.push_back("signal 3cm");
  filesel2.push_back(datafile);
  namesel2.push_back("sel2recomumu");
  legsel2.push_back("data");
  legendEntries = legsel2;  
  //comparesamevariable(filesel2, namesel2, "leadsubleaddR", 100, 550, 20, true, true, true, (float []){1e-3,1}, (float []){-1,0.7,0.85,0.99}, true, "#DeltaR(#mu_{1},#mu_{2})");
  //comparesamevariable(filesel2, namesel2, "leadsubleadM", 100, 350, 5, true, true, true, (float []){1e-4,1}, (float []){-1,0.7,0.85,0.99}, true, "M(#mu_{1},#mu_{2})");

  std::vector<int> coloroptschemesel20{kRed-9, kRed-7, kRed-3, kRed+2, kBlue};
  coloropt = coloroptschemesel20;
  std::vector<TFile*> filesel20;
  std::vector<TString> namesel20;
  std::vector<TString> legsel20;
  filesel20.push_back(sig3mfile);
  namesel20.push_back("sel20recomu");
  legsel20.push_back("signal 3m");
  filesel20.push_back(sig1mfile);
  namesel20.push_back("sel20recomu");
  legsel20.push_back("signal 1m");
  filesel20.push_back(sig30cmfile);
  namesel20.push_back("sel20recomu");
  legsel20.push_back("signal 30cm");
  filesel20.push_back(sig3cmfile);
  namesel20.push_back("sel20recomu");
  legsel20.push_back("signal 3cm");
  filesel20.push_back(datafile);
  namesel20.push_back("sel20recomu");
  legsel20.push_back("data");
  legendEntries = legsel20;  
  //comparesamevariable(filesel20, namesel20, "allpt", 60, 130, 2, true, true, true, (float []){1e-3,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu p_{T} / GeV");
  //comparesamevariable(filesel20, namesel20, "alleta", -1, -1, -1, true, true, true, (float []){1e-3,1e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #eta");
  //comparesamevariable(filesel20, namesel20, "allphi", -1, -1, 3, true, true, true, (float []){1e-2,2e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #phi");
  //comparesamevariable(filesel20, namesel20, "alldxy", 11000, 13000, 100, true, true, true, (float []){9e-4,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu d_{0} / cm");
  //comparesamevariable(filesel20, namesel20, "alllog10dxy", -1, 700, 20, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(filesel20, namesel20, "alldxysig", -1, 24000, 100, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu sigificance d_{0}");
  //comparesamevariable(filesel20, namesel20, "alllog10dxysig", 180, 800, 20, true, true, true, (float []){1e-2,3e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu log_{10}(significance d_{0})");
  //comparesamevariable(filesel20, namesel20, "leadpt", 60, 130, 2, true, true, true, (float []){1e-3,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu p_{T} / GeV");
  //comparesamevariable(filesel20, namesel20, "leadeta", -1, -1, -1, true, true, true, (float []){1e-3,1e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #eta");
  //comparesamevariable(filesel20, namesel20, "leadphi", -1, -1, 3, true, true, true, (float []){1e-2,2e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #phi");
  //comparesamevariable(filesel20, namesel20, "leaddxy", 11000, 13000, 100, true, true, true, (float []){9e-4,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu d_{0} / cm");
  //comparesamevariable(filesel20, namesel20, "leadlog10dxy", -1, 700, 20, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(filesel20, namesel20, "leaddxysig", -1, 24000, 100, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu sigificance d_{0}");
  //comparesamevariable(filesel20, namesel20, "leadlog10dxysig", 180, 800, 20, true, true, true, (float []){1e-2,3e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu log_{10}(significance d_{0})");
  //comparesamevariable(filesel20, namesel20, "subleadpt", 60, 130, 2, true, true, true, (float []){1e-3,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu p_{T} / GeV");
  //comparesamevariable(filesel20, namesel20, "subleadeta", -1, -1, -1, true, true, true, (float []){1e-3,1e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #eta");
  //comparesamevariable(filesel20, namesel20, "subleadphi", -1, -1, 3, true, true, true, (float []){1e-2,2e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu #phi");
  //comparesamevariable(filesel20, namesel20, "subleaddxy", 11000, 13000, 100, true, true, true, (float []){9e-4,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu d_{0} / cm");
  //comparesamevariable(filesel20, namesel20, "subleadlog10dxy", -1, 700, 20, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu log_{10}d_{0} / log_{10}cm");
  //comparesamevariable(filesel20, namesel20, "subleaddxysig", -1, 24000, 100, true, true, true, (float []){1e-5,1}, (float []){-1,0.7,0.85,0.99}, true, "#mu sigificance d_{0}");
  //comparesamevariable(filesel20, namesel20, "subleadlog10dxysig", 180, 800, 20, true, true, true, (float []){1e-2,3e-1}, (float []){-1,0.7,0.85,0.99}, true, "#mu log_{10}(significance d_{0})");
  filesel20.clear();
  namesel20.clear();
  legsel20.clear();
  filesel20.push_back(sig3mfile);
  namesel20.push_back("sel20recomumu");
  legsel20.push_back("signal 3m");
  filesel20.push_back(sig1mfile);
  namesel20.push_back("sel20recomumu");
  legsel20.push_back("signal 1m");
  filesel20.push_back(sig30cmfile);
  namesel20.push_back("sel20recomumu");
  legsel20.push_back("signal 30cm");
  filesel20.push_back(sig3cmfile);
  namesel20.push_back("sel20recomumu");
  legsel20.push_back("signal 3cm");
  filesel20.push_back(datafile);
  namesel20.push_back("sel20recomumu");
  legsel20.push_back("data");
  legendEntries = legsel20;  
  //comparesamevariable(filesel20, namesel20, "leadsubleaddR", 100, 550, 20, true, true, true, (float []){1e-3,1}, (float []){-1,0.7,0.85,0.99}, true, "#DeltaR(#mu_{1},#mu_{2})");
  //comparesamevariable(filesel20, namesel20, "leadsubleadM", 100, 350, 5, true, true, true, (float []){1e-4,1}, (float []){-1,0.7,0.85,0.99}, true, "M(#mu_{1},#mu_{2})");
  
  coloropt = coloroptrate;
  seltext[0] = "N#mu#geq2, p_{T}#geq16 GeV, |#eta|<2.5";
  seltext[1] = "#DeltaR>1, M<85 or M>95";
  //makeratehist("sel20recomu", "subleadpt", 66, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,20}, 0.875);
  //makeratehist("sel20recomu", "subleadlog10dxy", 1, 700, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,20}, 0.875);
  //makeratehist("sel20recomu", "subleadlog10dxysig", 180, 800, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,20}, 0.875);

  coloropt = coloroptrate;
  seltext[0] = "N#mu#geq2, p_{T}#geq16 GeV, |#eta|<2.5";
  seltext[1] = "#DeltaR>1, M<85 or M>95";
  //makeratehist("sel30recomu", "subleadpt", 66, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,3}, 0.875);
  //makeratehist("sel30recomu", "subleadlog10dxy", 1, 700, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,3}, 0.875);
  //makeratehist("sel30recomu", "subleadlog10dxysig", 180, 800, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,3}, 0.875);
  //makeratehist("sel30recomumu", "leadsubleaddR", 100, 550, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,3}, 0.875);
  //makeratehist("sel30recomumu", "leadsubleadM", 100, 350, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,3}, 0.875);

  std::vector<int> coloroptschemegentrig{46, 48, 31, 30, 8};
  coloropt = coloroptschemegentrig;
  std::vector<TFile*> fileanglegentrig;
  std::vector<TString> nameanglegentrig;
  std::vector<TString> leganglegentrig;
  fileanglegentrig.push_back(sig1mfile);
  nameanglegentrig.push_back("noselgenmul3ddm10");
  leganglegentrig.push_back("L3Mu10");
  fileanglegentrig.push_back(sig1mfile);
  nameanglegentrig.push_back("noselgenmul2dim23cs");
  leganglegentrig.push_back("L2Mu23 Cosmic");
  fileanglegentrig.push_back(sig1mfile);
  nameanglegentrig.push_back("noselgenmul2dim23");
  leganglegentrig.push_back("L2Mu23");
  fileanglegentrig.push_back(sig1mfile);
  nameanglegentrig.push_back("noselgenmul2ddm10");
  leganglegentrig.push_back("L2Mu10");
  legendEntries = leganglegentrig;  
  //comparesamevariable(fileanglegentrig, nameanglegentrig, "deta", 4700, 5300, 10, true, true, true, (float []){1e-3,1}, (float []){0.6,0.7,0.85,0.99}, true, "#Delta#eta(gen, trig.)","DDM/");
  //comparesamevariable(fileanglegentrig, nameanglegentrig, "dphi", 3000, 7000, 80, true, true, true, (float []){1e-3,1}, (float []){0.6,0.7,0.85,0.99}, true, "#Delta#phi(gen, trig.)","DDM/");

  fileanglegentrig.clear();
  nameanglegentrig.clear();
  leganglegentrig.clear();
  fileanglegentrig.push_back(sig3cmfile);
  nameanglegentrig.push_back("noselgenmul3ddm10");
  leganglegentrig.push_back("3 cm");
  fileanglegentrig.push_back(sig30cmfile);
  nameanglegentrig.push_back("noselgenmul3ddm10");
  leganglegentrig.push_back("30 cm");
  fileanglegentrig.push_back(sig1mfile);
  nameanglegentrig.push_back("noselgenmul3ddm10");
  leganglegentrig.push_back("1 m");
  fileanglegentrig.push_back(sig3mfile);
  nameanglegentrig.push_back("noselgenmul3ddm10");
  leganglegentrig.push_back("3 m");
  legendEntries = leganglegentrig;  
  //comparesamevariable(fileanglegentrig, nameanglegentrig, "deta", 4700, 5300, 10, true, true, true, (float []){1e-3,1}, (float []){0.6,0.7,0.85,0.99}, true, "#Delta#eta(gen, trig.)","DDM/genmul3ddm10_");
  //comparesamevariable(fileanglegentrig, nameanglegentrig, "dphi", 3000, 7000, 80, true, true, true, (float []){1e-3,1}, (float []){0.6,0.7,0.85,0.99}, true, "#Delta#phi(gen, trig.)","DDM/genmul3ddm10_");

  fileanglegentrig.clear();
  nameanglegentrig.clear();
  leganglegentrig.clear();
  fileanglegentrig.push_back(sig3cmfile);
  nameanglegentrig.push_back("noselgenmul2ddm10");
  leganglegentrig.push_back("3 cm");
  fileanglegentrig.push_back(sig30cmfile);
  nameanglegentrig.push_back("noselgenmul2ddm10");
  leganglegentrig.push_back("30 cm");
  fileanglegentrig.push_back(sig1mfile);
  nameanglegentrig.push_back("noselgenmul2ddm10");
  leganglegentrig.push_back("1 m");
  fileanglegentrig.push_back(sig3mfile);
  nameanglegentrig.push_back("noselgenmul2ddm10");
  leganglegentrig.push_back("3 m");
  legendEntries = leganglegentrig;  
  //comparesamevariable(fileanglegentrig, nameanglegentrig, "deta", 4700, 5300, 10, true, true, true, (float []){1e-3,1}, (float []){0.6,0.7,0.85,0.99}, true, "#Delta#eta(gen, trig.)","DDM/genmul2ddm10_");
  //comparesamevariable(fileanglegentrig, nameanglegentrig, "dphi", 3000, 7000, 80, true, true, true, (float []){1e-3,1}, (float []){0.6,0.7,0.85,0.99}, true, "#Delta#phi(gen, trig.)","DDM/genmul2ddm10_");

  std::vector<int> coloroptschemegeneff{1, 48, 31, 30, 46, 8};
  coloropt = coloroptschemegeneff;

  std::vector<TFile*> filegeneff3cm;
  std::vector<TString> namegeneff3cm;
  std::vector<TString> leggeneff3cm;
  filegeneff3cm.push_back(sig3cmfile);
  namegeneff3cm.push_back("basicselgen");
  leggeneff3cm.push_back("gen");
  filegeneff3cm.push_back(sig3cmfile);
  namegeneff3cm.push_back("basicsell3dim33genmch");
  leggeneff3cm.push_back("DoubleL3Mu33");
  filegeneff3cm.push_back(sig3cmfile);
  namegeneff3cm.push_back("basicsell2dim23genmch");
  leggeneff3cm.push_back("DoubleL2Mu23");
  filegeneff3cm.push_back(sig3cmfile);
  namegeneff3cm.push_back("basicsell2dim23csgenmch");
  leggeneff3cm.push_back("DoubleL2Mu23CS");
  filegeneff3cm.push_back(sig3cmfile);
  namegeneff3cm.push_back("basicsell3ddm10genmch");
  leggeneff3cm.push_back("DoubleL3Mu10");
  filegeneff3cm.push_back(sig3cmfile);
  namegeneff3cm.push_back("basicsell2ddm10genmch");
  leggeneff3cm.push_back("DoubleL2Mu10");
  legendEntries = leggeneff3cm;  
  //comparesamevariable(filegeneff3cm, namegeneff3cm, "mupt", 14, 55, 1, true, true, true, (float []){9e-1,1e3}, (float []){0.6,0.6,0.85,0.99}, false, "#mu p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneff3cm, namegeneff3cm, "mulog10d0", 100, 700, 10, false, true, true, (float []){0,3.49e2}, (float []){0.11,0.6,0.36,0.99}, false, "#mu log_{10}d_{0} / log_{10}cm","DDM/");
  //comparesamevariable(filegeneff3cm, namegeneff3cm, "ordptsubleadmupt", 14, 55, 1, true, true, true, (float []){9e-1,1e3}, (float []){0.6,0.6,0.85,0.99}, false, "p_{T} ordered, #mu_{2} p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneff3cm, namegeneff3cm, "ordd0subleadmulog10d0", 100, 700, 10, false, true, true, (float []){0,3.49e2}, (float []){0.11,0.6,0.36,0.99}, false, "|d_{0}| ordered, #mu_{2} log_{10}d_{0} / log_{10}cm","DDM/");

  std::vector<TFile*> filegeneff30cm;
  std::vector<TString> namegeneff30cm;
  std::vector<TString> leggeneff30cm;
  filegeneff30cm.push_back(sig30cmfile);
  namegeneff30cm.push_back("basicselgen");
  leggeneff30cm.push_back("gen");
  filegeneff30cm.push_back(sig30cmfile);
  namegeneff30cm.push_back("basicsell3dim33genmch");
  leggeneff30cm.push_back("DoubleL3Mu33");
  filegeneff30cm.push_back(sig30cmfile);
  namegeneff30cm.push_back("basicsell2dim23genmch");
  leggeneff30cm.push_back("DoubleL2Mu23");
  filegeneff30cm.push_back(sig30cmfile);
  namegeneff30cm.push_back("basicsell2dim23csgenmch");
  leggeneff30cm.push_back("DoubleL2Mu23CS");
  filegeneff30cm.push_back(sig30cmfile);
  namegeneff30cm.push_back("basicsell3ddm10genmch");
  leggeneff30cm.push_back("DoubleL3Mu10");
  filegeneff30cm.push_back(sig30cmfile);
  namegeneff30cm.push_back("basicsell2ddm10genmch");
  leggeneff30cm.push_back("DoubleL2Mu10");
  legendEntries = leggeneff30cm;  
  //comparesamevariable(filegeneff30cm, namegeneff30cm, "mupt", 14, 55, 1, true, true, true, (float []){9e-1,1e3}, (float []){0.6,0.6,0.85,0.99}, false, "#mu p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneff30cm, namegeneff30cm, "mulog10d0", 300, 750, 10, false, true, true, (float []){0,3.49e2}, (float []){0.11,0.6,0.36,0.99}, false, "#mu log_{10}d_{0} / log_{10}cm","DDM/");
  //comparesamevariable(filegeneff30cm, namegeneff30cm, "ordptsubleadmupt", 14, 55, 1, true, true, true, (float []){9e-1,1e3}, (float []){0.6,0.6,0.85,0.99}, false, "p_{T} ordered, #mu_{2} p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneff30cm, namegeneff30cm, "ordd0subleadmulog10d0", 300, 750, 10, false, true, true, (float []){0,3.49e2}, (float []){0.11,0.6,0.36,0.99}, false, "|d_{0}| ordered, #mu_{2} log_{10}d_{0} / log_{10}cm","DDM/");

  std::vector<TFile*> filegeneff1m;
  std::vector<TString> namegeneff1m;
  std::vector<TString> leggeneff1m;
  filegeneff1m.push_back(sig1mfile);
  namegeneff1m.push_back("basicselgen");
  leggeneff1m.push_back("gen");
  filegeneff1m.push_back(sig1mfile);
  namegeneff1m.push_back("basicsell3dim33genmch");
  leggeneff1m.push_back("DoubleL3Mu33");
  filegeneff1m.push_back(sig1mfile);
  namegeneff1m.push_back("basicsell2dim23genmch");
  leggeneff1m.push_back("DoubleL2Mu23");
  filegeneff1m.push_back(sig1mfile);
  namegeneff1m.push_back("basicsell2dim23csgenmch");
  leggeneff1m.push_back("DoubleL2Mu23CS");
  filegeneff1m.push_back(sig1mfile);
  namegeneff1m.push_back("basicsell3ddm10genmch");
  leggeneff1m.push_back("DoubleL3Mu10");
  filegeneff1m.push_back(sig1mfile);
  namegeneff1m.push_back("basicsell2ddm10genmch");
  leggeneff1m.push_back("DoubleL2Mu10");
  legendEntries = leggeneff1m;  
  //comparesamevariable(filegeneff1m, namegeneff1m, "mupt", 14, 55, 1, true, true, true, (float []){9e-1,1e3}, (float []){0.6,0.6,0.85,0.99}, false, "#mu p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneff1m, namegeneff1m, "mulog10d0", 300, 800, 10, false, true, true, (float []){0,3.49e2}, (float []){0.11,0.6,0.36,0.99}, false, "#mu log_{10}d_{0} / log_{10}cm","DDM/");
  //comparesamevariable(filegeneff1m, namegeneff1m, "ordptsubleadmupt", 14, 55, 1, true, true, true, (float []){9e-1,1e3}, (float []){0.6,0.6,0.85,0.99}, false, "p_{T} ordered, #mu_{2} p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneff1m, namegeneff1m, "ordd0subleadmulog10d0", 300, 800, 10, false, true, true, (float []){0,3.49e2}, (float []){0.11,0.6,0.36,0.99}, false, "|d_{0}| ordered, #mu_{2} log_{10}d_{0} / log_{10}cm","DDM/");

  std::vector<TFile*> filegeneffmff12;
  std::vector<TString> namegeneffmff12;
  std::vector<TString> leggeneffmff12;
  filegeneffmff12.push_back(h24mu12file);
  namegeneffmff12.push_back("basicselgen");
  leggeneffmff12.push_back("gen");
  filegeneffmff12.push_back(h24mu12file);
  namegeneffmff12.push_back("basicsell3dim33genmch");
  leggeneffmff12.push_back("DoubleL3Mu33");
  filegeneffmff12.push_back(h24mu12file);
  namegeneffmff12.push_back("basicsell2dim23genmch");
  leggeneffmff12.push_back("DoubleL2Mu23");
  filegeneffmff12.push_back(h24mu12file);
  namegeneffmff12.push_back("basicsell2dim23csgenmch");
  leggeneffmff12.push_back("DoubleL2Mu23CS");
  filegeneffmff12.push_back(h24mu12file);
  namegeneffmff12.push_back("basicsell3ddm10genmch");
  leggeneffmff12.push_back("DoubleL3Mu10");
  filegeneffmff12.push_back(h24mu12file);
  namegeneffmff12.push_back("basicsell2ddm10genmch");
  leggeneffmff12.push_back("DoubleL2Mu10");
  legendEntries = leggeneffmff12;  
  //comparesamevariable(filegeneffmff12, namegeneffmff12, "mupt", 14, 100, 1, true, true, true, (float []){1,2e5}, (float []){0.6,0.6,0.85,0.99}, false, "#mu p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneffmff12, namegeneffmff12, "mulog10d0", 300, 800, 10, false, true, true, (float []){0,1.2e4}, (float []){0.11,0.6,0.36,0.99}, false, "#mu log_{10}d_{0} / log_{10}cm","DDM/");
  //comparesamevariable(filegeneffmff12, namegeneffmff12, "ordptsubleadmupt", 14, 100, 1, true, true, true, (float []){1,2e5}, (float []){0.6,0.6,0.85,0.99}, false, "p_{T} ordered, #mu_{2} p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneffmff12, namegeneffmff12, "ordd0subleadmulog10d0", 300, 800, 10, false, true, true, (float []){0,4e3}, (float []){0.11,0.6,0.36,0.99}, false, "|d_{0}| ordered, #mu_{2} log_{10}d_{0} / log_{10}cm","DDM/");

  std::vector<TFile*> filegeneffmff25;
  std::vector<TString> namegeneffmff25;
  std::vector<TString> leggeneffmff25;
  filegeneffmff25.push_back(h24mu25file);
  namegeneffmff25.push_back("basicselgen");
  leggeneffmff25.push_back("gen");
  filegeneffmff25.push_back(h24mu25file);
  namegeneffmff25.push_back("basicsell3dim33genmch");
  leggeneffmff25.push_back("DoubleL3Mu33");
  filegeneffmff25.push_back(h24mu25file);
  namegeneffmff25.push_back("basicsell2dim23genmch");
  leggeneffmff25.push_back("DoubleL2Mu23");
  filegeneffmff25.push_back(h24mu25file);
  namegeneffmff25.push_back("basicsell2dim23csgenmch");
  leggeneffmff25.push_back("DoubleL2Mu23CS");
  filegeneffmff25.push_back(h24mu25file);
  namegeneffmff25.push_back("basicsell3ddm10genmch");
  leggeneffmff25.push_back("DoubleL3Mu10");
  filegeneffmff25.push_back(h24mu25file);
  namegeneffmff25.push_back("basicsell2ddm10genmch");
  leggeneffmff25.push_back("DoubleL2Mu10");
  legendEntries = leggeneffmff25;  
  //comparesamevariable(filegeneffmff25, namegeneffmff25, "mupt", 14, 100, 1, true, true, true, (float []){1,2e5}, (float []){0.6,0.6,0.85,0.99}, false, "#mu p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneffmff25, namegeneffmff25, "mulog10d0", 300, 800, 10, false, true, true, (float []){0,1.2e4}, (float []){0.11,0.6,0.36,0.99}, false, "#mu log_{10}d_{0} / log_{10}cm","DDM/");
  //comparesamevariable(filegeneffmff25, namegeneffmff25, "ordptsubleadmupt", 14, 100, 1, true, true, true, (float []){1,2e5}, (float []){0.6,0.6,0.85,0.99}, false, "p_{T} ordered, #mu_{2} p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneffmff25, namegeneffmff25, "ordd0subleadmulog10d0", 300, 800, 10, false, true, true, (float []){0,4e3}, (float []){0.11,0.6,0.36,0.99}, false, "|d_{0}| ordered, #mu_{2} log_{10}d_{0} / log_{10}cm","DDM/");

  std::vector<TFile*> filegeneffmff50;
  std::vector<TString> namegeneffmff50;
  std::vector<TString> leggeneffmff50;
  filegeneffmff50.push_back(h24mu50file);
  namegeneffmff50.push_back("basicselgen");
  leggeneffmff50.push_back("gen");
  filegeneffmff50.push_back(h24mu50file);
  namegeneffmff50.push_back("basicsell3dim33genmch");
  leggeneffmff50.push_back("DoubleL3Mu33");
  filegeneffmff50.push_back(h24mu50file);
  namegeneffmff50.push_back("basicsell2dim23genmch");
  leggeneffmff50.push_back("DoubleL2Mu23");
  filegeneffmff50.push_back(h24mu50file);
  namegeneffmff50.push_back("basicsell2dim23csgenmch");
  leggeneffmff50.push_back("DoubleL2Mu23CS");
  filegeneffmff50.push_back(h24mu50file);
  namegeneffmff50.push_back("basicsell3ddm10genmch");
  leggeneffmff50.push_back("DoubleL3Mu10");
  filegeneffmff50.push_back(h24mu50file);
  namegeneffmff50.push_back("basicsell2ddm10genmch");
  leggeneffmff50.push_back("DoubleL2Mu10");
  legendEntries = leggeneffmff50;  
  //comparesamevariable(filegeneffmff50, namegeneffmff50, "mupt", 14, 100, 1, true, true, true, (float []){1,2e5}, (float []){0.6,0.6,0.85,0.99}, false, "#mu p_{T} / GeV","DDM/");
  //comparesamevariable(filegeneffmff50, namegeneffmff50, "mulog10d0", 300, 800, 10, false, true, true, (float []){0,1.2e4}, (float []){0.11,0.6,0.36,0.99}, false, "#mu log_{10}d_{0} / log_{10}cm","DDM/");
  comparesamevariable(filegeneffmff50, namegeneffmff50, "ordptsubleadmupt", 14, 100, 1, true, true, true, (float []){1,2e5}, (float []){0.6,0.6,0.85,0.99}, false, "p_{T} ordered, #mu_{2} p_{T} / GeV","DDM/");
  comparesamevariable(filegeneffmff50, namegeneffmff50, "ordd0subleadmulog10d0", 300, 800, 10, false, true, true, (float []){0,4e3}, (float []){0.11,0.6,0.36,0.99}, false, "|d_{0}| ordered, #mu_{2} log_{10}d_{0} / log_{10}cm","DDM/");

  return -1;
}
