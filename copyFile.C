int copyFile() {

    auto chain = new TChain("demo/TrigObjects");
  //chain->Add("../out_1.root");
  chain->Add("/pnfs/iihe/cms/store/user/asahasra/EphemeralHLTPhysics1/DisplacedLeptonHLT_EphemeralHLTSkim1_asahasra/201216_220732/0000/out_*.root");

  std::cout<<chain->GetEntries()<<std::endl;

  int n;
  chain->SetBranchAddress("jetak8_n", &n);

  auto newfile = new TFile("skimmedFile_jetak8.root","RECREATE");
  auto newtree = chain->CloneTree(0);

  for(int event=0; event<chain->GetEntries(); event++) {

    if(event%100000==1) std::cout<<event<<" entries processed"<<std::endl;
    chain->GetEntry(event);

    if(n>0) {
      newtree->Fill();
    }
  }

  newfile->Write();
  return -1;
}
