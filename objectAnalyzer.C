int objectAnalyzer() {

  auto chain = new TChain("demo/TrigObjects");
  //chain->Add("../out_1.root");
  chain->Add("/pnfs/iihe/cms/store/user/asahasra/EphemeralHLTPhysics1/DisplacedLeptonHLT_EphemeralHLTSkim1_asahasra/201216_220732/0000/out_*.root");

  std::cout<<chain->GetEntries()<<std::endl;

  int mu_n, eg_n, ht_n, met_n, jetak4_n, jetak8_n;
  std::vector<double> *mu_pt = 0;
  std::vector<double> *eg_pt = 0;
  std::vector<double> *ht_pt = 0;
  std::vector<double> *met_pt = 0;
  std::vector<double> *jetak4_pt = 0;
  std::vector<double> *jetak8_pt = 0;
  chain->SetBranchAddress("mu_n", &mu_n);
  chain->SetBranchAddress("mu_pt", &mu_pt);
  chain->SetBranchAddress("eg_n", &eg_n);
  chain->SetBranchAddress("eg_pt", &eg_pt);
  chain->SetBranchAddress("ht_n", &ht_n);
  chain->SetBranchAddress("ht_pt", &ht_pt);
  chain->SetBranchAddress("met_n", &met_n);
  chain->SetBranchAddress("met_pt", &met_pt);
  chain->SetBranchAddress("jetak4_n", &jetak4_n);
  chain->SetBranchAddress("jetak4_pt", &jetak4_pt);
  chain->SetBranchAddress("jetak8_n", &jetak8_n);
  chain->SetBranchAddress("jetak8_pt", &jetak8_pt);
  
  auto outfile = new TFile("hists.root","RECREATE");
  auto hist_mu_n = new TH1F("mu_n","",10,0,10);
  auto hist_mu_pt = new TH1F("mu_pt","",500,0,500);
  auto hist_eg_n = new TH1F("eg_n","",10,0,10);
  auto hist_eg_pt = new TH1F("eg_pt","",500,0,500);
  auto hist_ht_n = new TH1F("ht_n","",10,0,10);
  auto hist_ht_pt = new TH1F("ht_pt","",500,0,500);
  auto hist_met_n = new TH1F("met_n","",10,0,10);
  auto hist_met_pt = new TH1F("met_pt","",500,0,500);
  auto hist_jetak4_n = new TH1F("jetak4_n","",10,0,10);
  auto hist_jetak4_pt = new TH1F("jetak4_pt","",500,0,500);
  auto hist_jetak8_n = new TH1F("jetak8_n","",10,0,10);
  auto hist_jetak8_pt = new TH1F("jetak8_pt","",500,0,500);
  
  for(int event=0; event<chain->GetEntries(); event++) {

    if(event%100000==1) std::cout<<event<<" entries processed"<<std::endl;
    chain->GetEntry(event);

    if(ht_n==0) continue;
    
    hist_mu_n->Fill(mu_n);
    for(int muctr=0; muctr<mu_n; muctr++) {
      hist_mu_pt->Fill(mu_pt->at(muctr));
    }
    hist_eg_n->Fill(eg_n);
    for(int egctr=0; egctr<eg_n; egctr++) {
      hist_eg_pt->Fill(eg_pt->at(egctr));
    }
    hist_ht_n->Fill(ht_n);
    for(int htctr=0; htctr<ht_n; htctr++) {
      hist_ht_pt->Fill(ht_pt->at(htctr));
    }
    hist_met_n->Fill(met_n);
    for(int metctr=0; metctr<met_n; metctr++) {
      hist_met_pt->Fill(met_pt->at(metctr));
    }
    hist_jetak4_n->Fill(jetak4_n);
    for(int jetak4ctr=0; jetak4ctr<jetak4_n; jetak4ctr++) {
      hist_jetak4_pt->Fill(jetak4_pt->at(jetak4ctr));
    }
    hist_jetak8_n->Fill(jetak8_n);
    for(int jetak8ctr=0; jetak8ctr<jetak8_n; jetak8ctr++) {
      hist_jetak8_pt->Fill(jetak8_pt->at(jetak8ctr));
    }

    mu_pt->clear();
    eg_pt->clear();
    ht_pt->clear();
    met_pt->clear();
    jetak4_pt->clear();
    jetak8_pt->clear();
  }

  outfile->Write();
  outfile->Close();
  return -1;
}
