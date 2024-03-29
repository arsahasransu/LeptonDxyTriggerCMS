#include <iostream>

using namespace std;

TCanvas* enhance_plotter(vector<TH1F*> histvec, vector<TString> legNam, TString xtitle, TString ytitle, float legPos[]=(float []){0.7,0.75,0.95,1}, bool logY=false, float yrange[]=(float []){0.1,1}, bool normalize=false) {

  TCanvas *c1 = new TCanvas("c1","c1",1500,1125);
  TPad *pad1 = new TPad("pad1","pad1",0,0,0.075,0.9);
  TPad *pad2 = new TPad("pad2","pad2",0.1,0,0.95,0.1);
  TPad *pad3 = new TPad("pad3","pad3",0.075,0.1,1,0.9);
  TPad *pad4 = new TPad("pad4","pad4",0.1,0.9,0.95,1);
  TPad *pad5 = new TPad("pad5","pad5",0.95,0,1,1);
  //pad1->SetFillColor(5);
  //pad2->SetFillColor(5);
  //pad3->SetFillColor(7);
  //pad4->SetFillColor(16);
  //pad5->SetFillColor(38);
  //pad5->Draw();
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();

  pad1->cd();
  TLatex laty;
  laty.SetTextSize(0.5);
  laty.SetTextFont(42);
  laty.SetTextAlign(32);
  laty.SetTextAngle(90);
  laty.DrawLatex(0.5,1,ytitle);

  pad2->cd();
  TLatex latx;
  latx.SetTextSize(0.5);
  latx.SetTextFont(42);
  latx.SetTextAlign(32);
  latx.DrawLatex(1,0.5,xtitle);

  pad4->cd();
  TLatex lat4_l;
  lat4_l.SetTextSize(0.5);
  lat4_l.SetTextFont(42);
  lat4_l.SetTextAlign(12);
  lat4_l.DrawLatex(0.07,0.5," ");
  TLatex lat4_r;
  lat4_r.SetTextSize(0.5);
  lat4_r.SetTextFont(42);
  lat4_r.SetTextAlign(32);
  lat4_r.DrawLatex(1,0.5,"2021, #sqrt{s} = 14 TeV");

  pad3->cd();
  gStyle->SetOptStat(0);
  pad3->SetLogy(logY);
  pad3->SetLeftMargin(0.1);
  pad3->SetRightMargin(0.05);
  pad3->SetTopMargin(0);
  pad3->SetBottomMargin(0.1);

  histvec[0]->GetXaxis()->SetTicks("-");
  histvec[0]->GetXaxis()->SetTickSize(0.05);
  histvec[0]->GetXaxis()->SetLabelFont(42);
  histvec[0]->GetXaxis()->SetLabelSize(0.06);
  histvec[0]->GetXaxis()->SetLabelOffset(-0.125);
  histvec[0]->GetYaxis()->SetTicks("+");
  histvec[0]->GetYaxis()->SetLabelFont(42);
  histvec[0]->GetYaxis()->SetLabelSize(0.06);
  histvec[0]->GetYaxis()->SetLabelOffset(-0.04);

  for(unsigned int ctr=0; ctr<histvec.size(); ctr++) {
    histvec[ctr]->SetLineWidth(5);
    histvec[ctr]->GetXaxis()->SetTitle("");
    histvec[ctr]->GetYaxis()->SetTitle("");
    if(normalize) histvec[ctr]->Scale(1.0/histvec[ctr]->Integral());
    histvec[ctr]->SetMinimum(yrange[0]);
    histvec[ctr]->SetMaximum(yrange[1]);
    histvec[ctr]->Draw("hist same e1");
  }

  TLegend* leg = new TLegend(legPos[0],legPos[1],legPos[2],legPos[3]);
  leg->SetTextFont(42);
  leg->SetTextSize(0.065);
  leg->SetBorderSize(0);
  for(unsigned int ctr=0; ctr<histvec.size(); ctr++) {
    leg->AddEntry(histvec[ctr],legNam[ctr],"l");
  }
  leg->Draw();
  
  return c1;
}

TCanvas* fancy_enhance_plotter(vector<TH1F*> histvec, vector<TString> legNam, TString xtitle, TString ytitle, float legPos[]=(float []){0.7,0.75,0.95,1}, bool logY=false, float yrange[]=(float []){0.1,1}, bool normalize=false, vector<TString> histtype={"p", "hist"}, vector<int> markerstyle={20,24}, vector<int> markersize={5, 5}, vector<TString> legendmarkerstyle={"lep", "l"}, vector<double> scale={1, 1}) {

  TCanvas *c1 = new TCanvas("c1","c1",1500,1125);
  TPad *pad1 = new TPad("pad1","pad1",0,0.1,0.075,0.98);
  TPad *pad2 = new TPad("pad2","pad2",0.1,0,0.99,0.1);
  TPad *pad3 = new TPad("pad3","pad3",0.075,0.1,1,0.98);
  //pad1->SetFillColor(5);
  //pad2->SetFillColor(5);
  //pad3->SetFillColor(7);
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();

  pad1->cd();
  TLatex laty;
  laty.SetTextSize(0.5);
  laty.SetTextFont(42);
  laty.SetTextAlign(32);
  laty.SetTextAngle(90);
  laty.DrawLatex(0.5,1,ytitle);

  pad2->cd();
  TLatex latx;
  latx.SetTextSize(0.5);
  latx.SetTextFont(42);
  latx.SetTextAlign(32);
  latx.DrawLatex(0.965,0.5,xtitle);

  pad3->cd();
  gStyle->SetOptStat(0);
  gStyle->SetLineWidth(4);
  pad3->SetLogy(logY);
  pad3->SetLeftMargin(0.11);
  pad3->SetRightMargin(0.05);
  pad3->SetTopMargin(0);
  pad3->SetBottomMargin(0.11);

  histvec[0]->GetXaxis()->SetTicks("-");
  histvec[0]->GetXaxis()->SetTickSize(0.05);
  histvec[0]->GetXaxis()->SetLabelFont(42);
  histvec[0]->GetXaxis()->SetLabelSize(0.06);
  histvec[0]->GetXaxis()->SetLabelOffset(-0.125);
  histvec[0]->GetYaxis()->SetTicks("+");
  histvec[0]->GetYaxis()->SetLabelFont(42);
  histvec[0]->GetYaxis()->SetLabelSize(0.06);
  histvec[0]->GetYaxis()->SetLabelOffset(-0.05);

  for(unsigned int ctr=0; ctr<histvec.size(); ctr++) {
    histvec[ctr]->SetLineWidth(5);
    histvec[ctr]->GetXaxis()->SetTitle("");
    histvec[ctr]->GetYaxis()->SetTitle("");
    if(normalize) histvec[ctr]->Scale(1.0/histvec[ctr]->Integral());
    if(scale[ctr]!=-1.0) histvec[ctr]->Scale(scale[ctr]);
    histvec[ctr]->SetMinimum(yrange[0]);
    histvec[ctr]->SetMaximum(yrange[1]);
    histvec[ctr]->SetMarkerStyle(markerstyle[ctr]);
    histvec[ctr]->SetMarkerSize(markersize[ctr]);
    histvec[ctr]->Draw(histtype[ctr]);
  }
  histvec[0]->Draw(histtype[0]+" same");  

  TLegend* leg = new TLegend(legPos[0],legPos[1],legPos[2],legPos[3]);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  for(unsigned int ctr=0; ctr<histvec.size(); ctr++) {
    leg->AddEntry(histvec[ctr],legNam[ctr],legendmarkerstyle[ctr]);
  }
  leg->Draw();
  
  return c1;
}

TCanvas* enhance_plotter_rate(vector<TH1F*> histvec, vector<TString> legNam, TString xtitle, TString ytitle, float legPos[]=(float []){0.7,0.75,0.95,1}, float yrange[]=(float []){0.1,1}, bool logY=false, bool normalize=false) {

  TCanvas *c1 = new TCanvas("c1","c1",1500,1125);
  TPad *pad1 = new TPad("pad1","pad1",0,0,0.075,0.99);
  TPad *pad2 = new TPad("pad2","pad2",0.1,0,0.925,0.1);
  TPad *pad3 = new TPad("pad3","pad3",0.075,0.1,0.925,0.99);
  TPad *pad4 = new TPad("pad4","pad4",0.925,0,1,0.99);
  //pad1->SetFillColor(kOrange); // Along y-axis
  //pad2->SetFillColor(kSpring+10); // Along x-axis
  //pad3->SetFillColor(kCyan);
  //pad4->SetFillColor(kGreen);
  
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();
  
  pad1->cd();
  TLatex laty;
  laty.SetTextSize(0.5);
  laty.SetTextFont(42);
  laty.SetTextAlign(32);
  laty.SetTextAngle(90);
  laty.SetTextColor(kBlack);
  laty.DrawLatex(0.5,0.995,ytitle);

  pad2->cd();
  TLatex latx;
  latx.SetTextSize(0.5);
  latx.SetTextFont(42);
  latx.SetTextAlign(32);
  latx.SetTextColor(kBlack);
  latx.DrawLatex(0.9,0.5,xtitle);

  pad4->cd();
  TLatex laty2;
  laty2.SetTextSize(0.5);
  laty2.SetTextFont(42);
  laty2.SetTextAlign(32);
  laty2.SetTextAngle(90);
  laty2.SetTextColor(histvec[1]->GetLineColor());
  laty2.DrawLatex(0.5,0.995,"signal acceptance / run 2 HLT");
  
  pad3->cd();
  gStyle->SetOptStat(0);
  gStyle->SetLineWidth(4);
  pad3->SetLogy(logY);
  pad3->SetLeftMargin(0.1);
  pad3->SetRightMargin(0.1);
  pad3->SetTopMargin(0.01);
  pad3->SetBottomMargin(0.1);

  histvec[0]->GetXaxis()->SetTicks("-");
  histvec[0]->GetXaxis()->SetTickSize(0.05);
  histvec[0]->GetXaxis()->SetLabelFont(42);
  histvec[0]->GetXaxis()->SetLabelSize(0.06);
  histvec[0]->GetXaxis()->SetLabelOffset(-0.12);
  histvec[0]->GetYaxis()->SetTicks("+");
  histvec[0]->GetYaxis()->SetLabelFont(42);
  histvec[0]->GetYaxis()->SetLabelSize(0.06);
  histvec[0]->GetYaxis()->SetLabelOffset(-0.035);
  histvec[0]->SetMinimum(yrange[0]);
  histvec[0]->SetMaximum(yrange[1]);
  
  histvec[0]->GetXaxis()->SetTitle("");
  histvec[0]->GetYaxis()->SetTitle("");
  if(normalize) histvec[0]->DrawNormalized("p e");
  else histvec[0]->Draw("p e");

  for(unsigned int ctr=histvec.size()-1; ctr>0; ctr--) {
    histvec[ctr]->SetLineWidth(5);
    histvec[ctr]->GetXaxis()->SetTitle("");
    histvec[ctr]->GetYaxis()->SetTitle("");
    if(normalize)
      histvec[ctr]->DrawNormalized("hist same");
    else
      histvec[ctr]->Draw("hist e same");
  }

  histvec[0]->SetMarkerStyle(20);
  histvec[0]->SetMarkerSize(2);
  if(normalize) histvec[0]->DrawNormalized("p same e");
  else histvec[0]->Draw("p same e");

  TLegend* leg = new TLegend(legPos[0],legPos[1],legPos[2],legPos[3]);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->AddEntry(histvec[0],legNam[0],"lep");
  for(unsigned int ctr=1; ctr<histvec.size(); ctr++) {
    leg->AddEntry(histvec[ctr],legNam[ctr],"l");
  }
  leg->Draw();

  return c1;
}
