#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 0.5/(2+2+3); // Data rate - 0.5 Hz, events from parent trigger selection
double sig3cmsf = 1.0/35;
double sig30cmsf = 1.0/28;
double sig1msf = 1.0/19;
double sig3msf = 1.0/7;
TString cutdeets = "Cut details";
TFile* datafile = TFile::Open("hists_Efmrl.root","READ");
TFile* sig3cmfile = TFile::Open("hists_M200dM20ctau3cm.root","READ");
TFile* sig30cmfile = TFile::Open("hists_M200dM20ctau30cm.root","READ");
TFile* sig1mfile = TFile::Open("hists_M200dM20ctau1m.root","READ");
TFile* sig3mfile = TFile::Open("hists_M200dM20ctau3m.root","READ");

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
  seltext[0] = "N#mu#geq2, p_{T}#geq33 GeV";
  seltext[1] = "#mu |#eta|<2.5, d_{0}>0.01 cm";
  //makeratehist("sel3recomu", "subleadpt", 66, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,0.35}, (float []){0,0.8}, 0.875);

  seltext[0] = "N#mu#geq2, p_{T}#geq16 GeV";
  seltext[1] = "#mu |#eta|<2.5";
  //makeratehist("sel2recomu", "subleadpt", 66, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,45}, 0.875);

  seltext[0] = "N#mu#geq2, p_{T}#geq20 GeV, |#eta|<2.5";
  seltext[1] = "#mu |d_{0}|>0.1 mm, |d_{0}| sig.>1";
  makeratehist("sel100recomu", "subleadpt", 66, 150, 1, false, false, (float []){0.55,0.65,0.75,0.975}, (float []){60,15}, (float []){0,2}, 0.875);

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

  std::vector<int> coloroptschemesel2{kBlue, kRed+2, kRed-3, kRed-7, kRed-9};
  coloropt = coloroptschemesel2;
  std::vector<TFile*> filesel2;
  std::vector<TString> namesel2;
  std::vector<TString> legsel2;
  filesel2.push_back(datafile);
  namesel2.push_back("sel2recomu");
  legsel2.push_back("data");
  filesel2.push_back(sig3cmfile);
  namesel2.push_back("sel2recomu");
  legsel2.push_back("signal 3cm");
  filesel2.push_back(sig30cmfile);
  namesel2.push_back("sel2recomu");
  legsel2.push_back("signal 30cm");
  filesel2.push_back(sig1mfile);
  namesel2.push_back("sel2recomu");
  legsel2.push_back("signal 1m");
  filesel2.push_back(sig3mfile);
  namesel2.push_back("sel2recomu");
  legsel2.push_back("signal 3m");
  legendEntries = legsel2;  
  //comparesamevariable(filesel2, namesel2, "allpt", 60, 210, 1, true, true, true, (float []){1e-1,1e3}, (float []){0.6,0.7,0.85,0.99}, false, "#mu p_{T} / GeV");
  //comparesamevariable(filesel2, namesel2, "alleta", -1, -1, -1, false, true, true, (float []){0,130}, (float []){0.6,0.7,0.85,0.99}, false, "#mu #eta");
  //comparesamevariable(filesel2, namesel2, "allphi", -1, -1, -1, false, true, true, (float []){0,70}, (float []){0.6,0.7,0.85,0.99}, false, "#mu #phi");

  return -1;
}
