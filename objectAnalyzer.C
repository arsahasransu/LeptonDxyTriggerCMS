int objectAnalyzer() {

  auto chain = new TChain("demo/TrigObjects");
  chain->Add("/pnfs/iihe/cms/store/user/asahasra/EphemeralHLTPhysics1/DisplacedLeptonHLT_EphemeralHLTSkim1_asahasra/201216_220732/0000/out_*.root");

  std::cout<<chain->GetEntries()<<std::endl;

  return -1;
}
