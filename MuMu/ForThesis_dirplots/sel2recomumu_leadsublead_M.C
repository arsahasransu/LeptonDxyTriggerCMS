void sel2recomumu_leadsublead_M()
{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Apr 18 16:36:37 2023) by ROOT version 6.18/04
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,1500,1125);
   gStyle->SetOptStat(0);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: pad1
   TPad *pad1 = new TPad("pad1", "pad1",0,0.1,0.075,0.98);
   pad1->Draw();
   pad1->cd();
   pad1->Range(0,0,1,1);
   pad1->SetFillColor(0);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(2);
   pad1->SetFrameBorderMode(0);
   TLatex *   tex = new TLatex(0.5,1,"normalized events / GeV");
   tex->SetTextAlign(32);
   tex->SetTextFont(42);
   tex->SetTextSize(0.5);
   tex->SetTextAngle(90);
   tex->SetLineWidth(2);
   tex->Draw();
   pad1->Modified();
   c1->cd();
  
// ------------>Primitives in pad: pad2
   TPad *pad2 = new TPad("pad2", "pad2",0.1,0,0.99,0.1);
   pad2->Draw();
   pad2->cd();
   pad2->Range(0,0,1,1);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetFrameBorderMode(0);
      tex = new TLatex(0.965,0.5,"M(#mu_{1}, #mu_{2}) [GeV]");
   tex->SetTextAlign(32);
   tex->SetTextFont(42);
   tex->SetTextSize(0.5);
   tex->SetLineWidth(2);
   tex->Draw();
   pad2->Modified();
   c1->cd();
  
// ------------>Primitives in pad: pad3
   TPad *pad3 = new TPad("pad3", "pad3",0.075,0.1,1,0.98);
   pad3->Draw();
   pad3->cd();
   pad3->Range(-16.36905,-2.507628,132.4405,0.69897);
   pad3->SetFillColor(0);
   pad3->SetBorderMode(0);
   pad3->SetBorderSize(2);
   pad3->SetLogy();
   pad3->SetRightMargin(0.05);
   pad3->SetLeftMargin(0.11);
   pad3->SetTopMargin(0);
   pad3->SetFrameBorderMode(0);
   pad3->SetFrameBorderMode(0);
   gStyle->SetLineWidth(4);
   Double_t xAxis25[23] = {0, 4, 8, 16, 28, 34, 42, 48, 56, 64, 72, 80, 84, 88, 90, 91, 92, 94, 96, 100, 108, 116, 125}; 
   
   TH1F *newx__25 = new TH1F("newx__25","",22, xAxis25);
   newx__25->SetBinContent(1,0.08011824);
   newx__25->SetBinContent(2,0.02629363);
   newx__25->SetBinContent(3,0.01384282);
   newx__25->SetBinContent(4,0.005928956);
   newx__25->SetBinContent(5,0.004846277);
   newx__25->SetBinContent(6,0.01283748);
   newx__25->SetBinContent(7,0.01505439);
   newx__25->SetBinContent(8,0.01523484);
   newx__25->SetBinContent(9,0.01422949);
   newx__25->SetBinContent(10,0.01585351);
   newx__25->SetBinContent(11,0.01577618);
   newx__25->SetBinContent(12,0.02505628);
   newx__25->SetBinContent(13,0.05057657);
   newx__25->SetBinContent(14,0.1426043);
   newx__25->SetBinContent(15,0.1843647);
   newx__25->SetBinContent(16,0.178178);
   newx__25->SetBinContent(17,0.1132173);
   newx__25->SetBinContent(18,0.04918455);
   newx__25->SetBinContent(19,0.01933355);
   newx__25->SetBinContent(20,0.00920277);
   newx__25->SetBinContent(21,0.00386671);
   newx__25->SetBinContent(22,0.004399457);
   newx__25->SetBinContent(23,0.1552871);
   newx__25->SetBinError(1,0.007040386);
   newx__25->SetBinError(2,0.004033259);
   newx__25->SetBinError(3,0.002926464);
   newx__25->SetBinError(4,0.001915225);
   newx__25->SetBinError(5,0.00173155);
   newx__25->SetBinError(6,0.002818193);
   newx__25->SetBinError(7,0.003051845);
   newx__25->SetBinError(8,0.00307008);
   newx__25->SetBinError(9,0.002967055);
   newx__25->SetBinError(10,0.003131797);
   newx__25->SetBinError(11,0.003124149);
   newx__25->SetBinError(12,0.003937215);
   newx__25->SetBinError(13,0.005593781);
   newx__25->SetBinError(14,0.009392843);
   newx__25->SetBinError(15,0.01067996);
   newx__25->SetBinError(16,0.01049924);
   newx__25->SetBinError(17,0.008369262);
   newx__25->SetBinError(18,0.005516266);
   newx__25->SetBinError(19,0.003458491);
   newx__25->SetBinError(20,0.002386108);
   newx__25->SetBinError(21,0.001546684);
   newx__25->SetBinError(22,0.001649796);
   newx__25->SetBinError(23,0.009801634);
   newx__25->SetMinimum(0.007);
   newx__25->SetMaximum(5);
   newx__25->SetEntries(4823);
   newx__25->SetLineWidth(5);
   newx__25->SetMarkerStyle(20);
   newx__25->SetMarkerSize(2);
   newx__25->GetXaxis()->SetTicks("-");
   newx__25->GetYaxis()->SetTicks("+");
   newx__25->GetXaxis()->SetRange(1,22);
   newx__25->GetXaxis()->SetLabelFont(42);
   newx__25->GetXaxis()->SetLabelOffset(-0.12);
   newx__25->GetXaxis()->SetLabelSize(0.06);
   newx__25->GetXaxis()->SetTitleSize(0.035);
   newx__25->GetXaxis()->SetTickLength(0.05);
   newx__25->GetXaxis()->SetTitleOffset(1);
   newx__25->GetXaxis()->SetTitleFont(42);
   newx__25->GetYaxis()->SetLabelFont(42);
   newx__25->GetYaxis()->SetLabelOffset(-0.05);
   newx__25->GetYaxis()->SetLabelSize(0.06);
   newx__25->GetYaxis()->SetTitleSize(0.035);
   newx__25->GetYaxis()->SetTitleFont(42);
   newx__25->GetZaxis()->SetLabelFont(42);
   newx__25->GetZaxis()->SetLabelSize(0.035);
   newx__25->GetZaxis()->SetTitleSize(0.035);
   newx__25->GetZaxis()->SetTitleOffset(1);
   newx__25->GetZaxis()->SetTitleFont(42);
   newx__25->Draw("p e1");
   Double_t xAxis26[23] = {0, 4, 8, 16, 28, 34, 42, 48, 56, 64, 72, 80, 84, 88, 90, 91, 92, 94, 96, 100, 108, 116, 125}; 
   
   TH1F *newx__26 = new TH1F("newx__26","",22, xAxis26);
   newx__26->SetBinContent(1,0.006060606);
   newx__26->SetBinContent(3,0.003030303);
   newx__26->SetBinContent(4,0.008080808);
   newx__26->SetBinContent(5,0.02828283);
   newx__26->SetBinContent(6,0.07575758);
   newx__26->SetBinContent(7,0.1292929);
   newx__26->SetBinContent(8,0.1636364);
   newx__26->SetBinContent(9,0.1060606);
   newx__26->SetBinContent(10,0.06363636);
   newx__26->SetBinContent(11,0.06969697);
   newx__26->SetBinContent(12,0.04848485);
   newx__26->SetBinContent(13,0.04848485);
   newx__26->SetBinContent(14,0.02424242);
   newx__26->SetBinContent(15,0.02424242);
   newx__26->SetBinContent(16,0.07272727);
   newx__26->SetBinContent(17,0.02424242);
   newx__26->SetBinContent(18,0.04848485);
   newx__26->SetBinContent(19,0.01212121);
   newx__26->SetBinContent(20,0.02121212);
   newx__26->SetBinContent(21,0.006060606);
   newx__26->SetBinContent(22,0.01616162);
   newx__26->SetBinContent(23,0.2424242);
   newx__26->SetBinError(1,0.01212121);
   newx__26->SetBinError(3,0.008570991);
   newx__26->SetBinError(4,0.01399637);
   newx__26->SetBinError(5,0.02618481);
   newx__26->SetBinError(6,0.04285496);
   newx__26->SetBinError(7,0.05598548);
   newx__26->SetBinError(8,0.06298367);
   newx__26->SetBinError(9,0.05070667);
   newx__26->SetBinError(10,0.03927722);
   newx__26->SetBinError(11,0.04110503);
   newx__26->SetBinError(12,0.03428397);
   newx__26->SetBinError(13,0.03428397);
   newx__26->SetBinError(14,0.02424242);
   newx__26->SetBinError(15,0.02424242);
   newx__26->SetBinError(16,0.04198911);
   newx__26->SetBinError(17,0.02424242);
   newx__26->SetBinError(18,0.03428397);
   newx__26->SetBinError(19,0.01714198);
   newx__26->SetBinError(20,0.02267671);
   newx__26->SetBinError(21,0.01212121);
   newx__26->SetBinError(22,0.01979386);
   newx__26->SetBinError(23,0.07666128);
   newx__26->SetMinimum(0.007);
   newx__26->SetMaximum(5);
   newx__26->SetEntries(280);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#660000");
   newx__26->SetLineColor(ci);
   newx__26->SetLineWidth(5);
   newx__26->SetMarkerSize(0);
   newx__26->GetXaxis()->SetRange(1,1200);
   newx__26->GetXaxis()->SetLabelFont(42);
   newx__26->GetXaxis()->SetLabelSize(0.035);
   newx__26->GetXaxis()->SetTitleSize(0.035);
   newx__26->GetXaxis()->SetTitleOffset(1);
   newx__26->GetXaxis()->SetTitleFont(42);
   newx__26->GetYaxis()->SetLabelFont(42);
   newx__26->GetYaxis()->SetLabelSize(0.035);
   newx__26->GetYaxis()->SetTitleSize(0.035);
   newx__26->GetYaxis()->SetTitleFont(42);
   newx__26->GetZaxis()->SetLabelFont(42);
   newx__26->GetZaxis()->SetLabelSize(0.035);
   newx__26->GetZaxis()->SetTitleSize(0.035);
   newx__26->GetZaxis()->SetTitleOffset(1);
   newx__26->GetZaxis()->SetTitleFont(42);
   newx__26->Draw("hist  same");
   Double_t xAxis27[23] = {0, 4, 8, 16, 28, 34, 42, 48, 56, 64, 72, 80, 84, 88, 90, 91, 92, 94, 96, 100, 108, 116, 125}; 
   
   TH1F *newx__27 = new TH1F("newx__27","",22, xAxis27);
   newx__27->SetBinContent(6,0.09156977);
   newx__27->SetBinContent(7,0.1395349);
   newx__27->SetBinContent(8,0.1046512);
   newx__27->SetBinContent(9,0.1438953);
   newx__27->SetBinContent(10,0.1046512);
   newx__27->SetBinContent(11,0.06540698);
   newx__27->SetBinContent(12,0.05232558);
   newx__27->SetBinContent(13,0.07848837);
   newx__27->SetBinContent(15,0.1046512);
   newx__27->SetBinContent(19,0.02616279);
   newx__27->SetBinContent(20,0.03924419);
   newx__27->SetBinContent(21,0.02616279);
   newx__27->SetBinContent(22,0.02325581);
   newx__27->SetBinContent(23,0.627907);
   newx__27->SetBinError(6,0.0978922);
   newx__27->SetBinError(7,0.1208408);
   newx__27->SetBinError(8,0.1046512);
   newx__27->SetBinError(9,0.1227144);
   newx__27->SetBinError(10,0.1046512);
   newx__27->SetBinError(11,0.08273401);
   newx__27->SetBinError(12,0.07399955);
   newx__27->SetBinError(13,0.09063057);
   newx__27->SetBinError(15,0.1046512);
   newx__27->SetBinError(19,0.05232558);
   newx__27->SetBinError(20,0.06408549);
   newx__27->SetBinError(21,0.05232558);
   newx__27->SetBinError(22,0.04933303);
   newx__27->SetBinError(23,0.2563419);
   newx__27->SetMinimum(0.007);
   newx__27->SetMaximum(5);
   newx__27->SetEntries(89);

   ci = TColor::GetColor("#ff0000");
   newx__27->SetLineColor(ci);
   newx__27->SetLineWidth(5);
   newx__27->SetMarkerSize(0);
   newx__27->GetXaxis()->SetRange(1,1200);
   newx__27->GetXaxis()->SetLabelFont(42);
   newx__27->GetXaxis()->SetLabelSize(0.035);
   newx__27->GetXaxis()->SetTitleSize(0.035);
   newx__27->GetXaxis()->SetTitleOffset(1);
   newx__27->GetXaxis()->SetTitleFont(42);
   newx__27->GetYaxis()->SetLabelFont(42);
   newx__27->GetYaxis()->SetLabelSize(0.035);
   newx__27->GetYaxis()->SetTitleSize(0.035);
   newx__27->GetYaxis()->SetTitleFont(42);
   newx__27->GetZaxis()->SetLabelFont(42);
   newx__27->GetZaxis()->SetLabelSize(0.035);
   newx__27->GetZaxis()->SetTitleSize(0.035);
   newx__27->GetZaxis()->SetTitleOffset(1);
   newx__27->GetZaxis()->SetTitleFont(42);
   newx__27->Draw("hist  same");
   Double_t xAxis28[23] = {0, 4, 8, 16, 28, 34, 42, 48, 56, 64, 72, 80, 84, 88, 90, 91, 92, 94, 96, 100, 108, 116, 125}; 
   
   TH1F *newx__28 = new TH1F("newx__28","",22, xAxis28);
   newx__28->SetBinContent(1,0.01455133);
   newx__28->SetBinContent(3,0.021827);
   newx__28->SetBinContent(4,0.004850445);
   newx__28->SetBinContent(5,0.04850445);
   newx__28->SetBinContent(6,0.087308);
   newx__28->SetBinContent(7,0.1261116);
   newx__28->SetBinContent(8,0.1600647);
   newx__28->SetBinContent(9,0.07275667);
   newx__28->SetBinContent(10,0.065481);
   newx__28->SetBinContent(11,0.08003233);
   newx__28->SetBinContent(12,0.05820534);
   newx__28->SetBinContent(13,0.043654);
   newx__28->SetBinContent(16,0.05820534);
   newx__28->SetBinContent(17,0.02910267);
   newx__28->SetBinContent(18,0.05820534);
   newx__28->SetBinContent(19,0.043654);
   newx__28->SetBinContent(20,0.01455133);
   newx__28->SetBinContent(22,0.01293452);
   newx__28->SetBinContent(23,1.280517);
   newx__28->SetBinError(1,0.02910267);
   newx__28->SetBinError(3,0.03564334);
   newx__28->SetBinError(4,0.01680243);
   newx__28->SetBinError(5,0.05313396);
   newx__28->SetBinError(6,0.07128669);
   newx__28->SetBinError(7,0.08567594);
   newx__28->SetBinError(8,0.09652263);
   newx__28->SetBinError(9,0.06507554);
   newx__28->SetBinError(10,0.06173608);
   newx__28->SetBinError(11,0.06825181);
   newx__28->SetBinError(12,0.05820534);
   newx__28->SetBinError(13,0.0504073);
   newx__28->SetBinError(16,0.05820534);
   newx__28->SetBinError(17,0.04115739);
   newx__28->SetBinError(18,0.05820534);
   newx__28->SetBinError(19,0.0504073);
   newx__28->SetBinError(20,0.02910267);
   newx__28->SetBinError(22,0.02743826);
   newx__28->SetBinError(23,0.2730072);
   newx__28->SetMinimum(0.007);
   newx__28->SetMaximum(5);
   newx__28->SetEntries(149);

   ci = TColor::GetColor("#cc6600");
   newx__28->SetLineColor(ci);
   newx__28->SetLineWidth(5);
   newx__28->SetMarkerSize(0);
   newx__28->GetXaxis()->SetRange(1,1200);
   newx__28->GetXaxis()->SetLabelFont(42);
   newx__28->GetXaxis()->SetLabelSize(0.035);
   newx__28->GetXaxis()->SetTitleSize(0.035);
   newx__28->GetXaxis()->SetTitleOffset(1);
   newx__28->GetXaxis()->SetTitleFont(42);
   newx__28->GetYaxis()->SetLabelFont(42);
   newx__28->GetYaxis()->SetLabelSize(0.035);
   newx__28->GetYaxis()->SetTitleSize(0.035);
   newx__28->GetYaxis()->SetTitleFont(42);
   newx__28->GetZaxis()->SetLabelFont(42);
   newx__28->GetZaxis()->SetLabelSize(0.035);
   newx__28->GetZaxis()->SetTitleSize(0.035);
   newx__28->GetZaxis()->SetTitleOffset(1);
   newx__28->GetZaxis()->SetTitleFont(42);
   newx__28->Draw("hist  same");
   Double_t xAxis29[23] = {0, 4, 8, 16, 28, 34, 42, 48, 56, 64, 72, 80, 84, 88, 90, 91, 92, 94, 96, 100, 108, 116, 125}; 
   
   TH1F *newx__29 = new TH1F("newx__29","",22, xAxis29);
   newx__29->SetBinContent(1,0.1395349);
   newx__29->SetBinContent(2,0.06976745);
   newx__29->SetBinContent(3,0.03488372);
   newx__29->SetBinContent(6,0.1744186);
   newx__29->SetBinContent(7,0.09302326);
   newx__29->SetBinContent(8,0.03488372);
   newx__29->SetBinContent(9,0.06976745);
   newx__29->SetBinContent(11,0.03488372);
   newx__29->SetBinContent(12,0.06976745);
   newx__29->SetBinContent(13,0.06976745);
   newx__29->SetBinContent(14,0.1395349);
   newx__29->SetBinContent(21,0.06976745);
   newx__29->SetBinContent(23,1.953488);
   newx__29->SetBinError(1,0.1973321);
   newx__29->SetBinError(2,0.1395349);
   newx__29->SetBinError(3,0.09866606);
   newx__29->SetBinError(6,0.220624);
   newx__29->SetBinError(7,0.161121);
   newx__29->SetBinError(8,0.09866606);
   newx__29->SetBinError(9,0.1395349);
   newx__29->SetBinError(11,0.09866606);
   newx__29->SetBinError(12,0.1395349);
   newx__29->SetBinError(13,0.1395349);
   newx__29->SetBinError(14,0.1973321);
   newx__29->SetBinError(21,0.1395349);
   newx__29->SetBinError(23,0.7383492);
   newx__29->SetMinimum(0.007);
   newx__29->SetMaximum(5);
   newx__29->SetEntries(49);

   ci = TColor::GetColor("#ffcc00");
   newx__29->SetLineColor(ci);
   newx__29->SetLineWidth(5);
   newx__29->SetMarkerSize(0);
   newx__29->GetXaxis()->SetRange(1,1200);
   newx__29->GetXaxis()->SetLabelFont(42);
   newx__29->GetXaxis()->SetLabelSize(0.035);
   newx__29->GetXaxis()->SetTitleSize(0.035);
   newx__29->GetXaxis()->SetTitleOffset(1);
   newx__29->GetXaxis()->SetTitleFont(42);
   newx__29->GetYaxis()->SetLabelFont(42);
   newx__29->GetYaxis()->SetLabelSize(0.035);
   newx__29->GetYaxis()->SetTitleSize(0.035);
   newx__29->GetYaxis()->SetTitleFont(42);
   newx__29->GetZaxis()->SetLabelFont(42);
   newx__29->GetZaxis()->SetLabelSize(0.035);
   newx__29->GetZaxis()->SetTitleSize(0.035);
   newx__29->GetZaxis()->SetTitleOffset(1);
   newx__29->GetZaxis()->SetTitleFont(42);
   newx__29->Draw("hist  same");
   Double_t xAxis30[23] = {0, 4, 8, 16, 28, 34, 42, 48, 56, 64, 72, 80, 84, 88, 90, 91, 92, 94, 96, 100, 108, 116, 125}; 
   
   TH1F *newx__30 = new TH1F("newx__30","",22, xAxis30);
   newx__30->SetBinContent(1,0.08011824);
   newx__30->SetBinContent(2,0.02629363);
   newx__30->SetBinContent(3,0.01384282);
   newx__30->SetBinContent(4,0.005928956);
   newx__30->SetBinContent(5,0.004846277);
   newx__30->SetBinContent(6,0.01283748);
   newx__30->SetBinContent(7,0.01505439);
   newx__30->SetBinContent(8,0.01523484);
   newx__30->SetBinContent(9,0.01422949);
   newx__30->SetBinContent(10,0.01585351);
   newx__30->SetBinContent(11,0.01577618);
   newx__30->SetBinContent(12,0.02505628);
   newx__30->SetBinContent(13,0.05057657);
   newx__30->SetBinContent(14,0.1426043);
   newx__30->SetBinContent(15,0.1843647);
   newx__30->SetBinContent(16,0.178178);
   newx__30->SetBinContent(17,0.1132173);
   newx__30->SetBinContent(18,0.04918455);
   newx__30->SetBinContent(19,0.01933355);
   newx__30->SetBinContent(20,0.00920277);
   newx__30->SetBinContent(21,0.00386671);
   newx__30->SetBinContent(22,0.004399457);
   newx__30->SetBinContent(23,0.1552871);
   newx__30->SetBinError(1,0.007040386);
   newx__30->SetBinError(2,0.004033259);
   newx__30->SetBinError(3,0.002926464);
   newx__30->SetBinError(4,0.001915225);
   newx__30->SetBinError(5,0.00173155);
   newx__30->SetBinError(6,0.002818193);
   newx__30->SetBinError(7,0.003051845);
   newx__30->SetBinError(8,0.00307008);
   newx__30->SetBinError(9,0.002967055);
   newx__30->SetBinError(10,0.003131797);
   newx__30->SetBinError(11,0.003124149);
   newx__30->SetBinError(12,0.003937215);
   newx__30->SetBinError(13,0.005593781);
   newx__30->SetBinError(14,0.009392843);
   newx__30->SetBinError(15,0.01067996);
   newx__30->SetBinError(16,0.01049924);
   newx__30->SetBinError(17,0.008369262);
   newx__30->SetBinError(18,0.005516266);
   newx__30->SetBinError(19,0.003458491);
   newx__30->SetBinError(20,0.002386108);
   newx__30->SetBinError(21,0.001546684);
   newx__30->SetBinError(22,0.001649796);
   newx__30->SetBinError(23,0.009801634);
   newx__30->SetMinimum(0.007);
   newx__30->SetMaximum(5);
   newx__30->SetEntries(4823);
   newx__30->SetLineWidth(5);
   newx__30->SetMarkerStyle(20);
   newx__30->SetMarkerSize(2);
   newx__30->SetLineColor(kBlack);
   newx__30->GetXaxis()->SetRange(1,1200);
   newx__30->GetXaxis()->SetLabelFont(42);
   newx__30->GetXaxis()->SetLabelOffset(-0.125);
   newx__30->GetXaxis()->SetLabelSize(0.06);
   newx__30->GetXaxis()->SetTitleSize(0.035);
   newx__30->GetXaxis()->SetTickLength(0.05);
   newx__30->GetXaxis()->SetTitleOffset(1);
   newx__30->GetXaxis()->SetTitleFont(42);
   newx__30->GetYaxis()->SetLabelFont(42);
   newx__30->GetYaxis()->SetLabelOffset(-0.05);
   newx__30->GetYaxis()->SetLabelSize(0.06);
   newx__30->GetYaxis()->SetTitleSize(0.035);
   newx__30->GetYaxis()->SetTitleFont(42);
   newx__30->GetZaxis()->SetLabelFont(42);
   newx__30->GetZaxis()->SetLabelSize(0.035);
   newx__30->GetZaxis()->SetTitleSize(0.035);
   newx__30->GetZaxis()->SetTitleOffset(1);
   newx__30->GetZaxis()->SetTitleFont(42);
   newx__30->Draw("p e1 same");
   
   TLegend *leg = new TLegend(0.55,0.65,0.8,0.99,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(4);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("newx","2018 data","lep");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("newx","#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm","l");

   ci = TColor::GetColor("#660000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("newx","#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm","l");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("newx","#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m","l");

   ci = TColor::GetColor("#cc6600");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("newx","#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m","l");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   TLine *line = new TLine(0,5,125,5);
   line->SetLineWidth(4);
   line->Draw();
   line = new TLine(125,0.007,125,5);
   line->SetLineWidth(4);
   line->Draw();
      tex = new TLatex(4,2,"N #mu #geq 2, p_{T} #geq 16 GeV");
   tex->SetTextFont(42);
   tex->SetTextSize(0.045);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(4,1,"#mu |#eta| < 2.5");
   tex->SetTextFont(42);
   tex->SetTextSize(0.045);
   tex->SetLineWidth(2);
   tex->Draw();
   pad3->Modified();
   c1->cd();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   c1->SaveAs("sel2recomumu_leadsublead_M.png");
}
