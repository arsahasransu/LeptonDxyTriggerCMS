int objectAnalyzer() {

  auto chain = new TChain("demo/TrigObjects");
  chain->Add("../out_1.root");
  //chain->Add("/pnfs/iihe/cms/store/user/asahasra/EphemeralHLTPhysics1/DisplacedLeptonHLT_EphemeralHLTSkim1_asahasra/201216_220732/0000/out_*.root");

  std::cout<<chain->GetEntries()<<std::endl;

  int mu_n, eg_n, ht_n, met_n, jetak4_n, jetak8_n;
  chain->SetBranchAddress("mu_n", &mu_n);
  chain->SetBranchAddress("eg_n", &eg_n);
  chain->SetBranchAddress("ht_n", &ht_n);
  chain->SetBranchAddress("met_n", &met_n);
  chain->SetBranchAddress("jetak4_n", &jetak4_n);
  chain->SetBranchAddress("jetak8_n", &jetak8_n);
  
  auto outfile = new TFile("hists.root","RECREATE");
  auto hist_mu_n = new TH1F("mu_n","",10,0,10);
  auto hist_eg_n = new TH1F("eg_n","",10,0,10);
  auto hist_ht_n = new TH1F("ht_n","",10,0,10);
  auto hist_met_n = new TH1F("met_n","",10,0,10);
  auto hist_jetak4_n = new TH1F("jetak4_n","",10,0,10);
  auto hist_jetak8_n = new TH1F("jetak8_n","",10,0,10);
  
  for(int event=0; event<chain->GetEntries(); event++) {

    chain->GetEntry(event);
    
    hist_mu_n->Fill(mu_n);
    hist_eg_n->Fill(eg_n);
    hist_ht_n->Fill(ht_n);
    hist_met_n->Fill(met_n);
    hist_jetak4_n->Fill(jetak4_n);
    hist_jetak8_n->Fill(jetak8_n);
  }

  outfile->Write();
  outfile->Close();
  return -1;
}
