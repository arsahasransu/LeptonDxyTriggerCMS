int HQplotsAll() {


  // GENERATOR LEPTON PROPERTIES
  vector<double> gen_binslog10d0{-4.0,-3.6,-3.2,-2.8,-2.4,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0};
  int gen_nbinslog10d0 = gen_binslog10d0.size()-1;
  vector<double> gen_binspt{10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,34,38,42,46,50};
  int gen_nbinspt = gen_binspt.size()-1;
  seltext[0] = " ";
  seltext[1] = " ";
  
  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();
  
  file.push_back(sig3cmfile);
  name.push_back("gennoselgeneg");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("gennoselgeneg");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig1mfile);
  name.push_back("gennoselgeneg");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig3mfile);
  name.push_back("gennoselgeneg");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "log10d0", gen_nbinslog10d0, &gen_binslog10d0[0], true, false, false, (float []){6e-5,6}, (float []){0.525,0.7,0.775,0.99}, (float []){0.015,0.015}, true, "gen. lepton log_{10}d_{0} [log_{10}cm]", "normalized events / log_{10}cm", true, "genleptond0");
  //comparesamevariable(file, name, "subleadpt", gen_nbinspt, &gen_binspt[0], true, false, false, (float []){6e-5,6}, (float []){0.525,0.7,0.775,0.99}, (float []){0.015,0.015}, true, "gen. lepton_{2} p_{T} [GeV]", "normalized events / GeV", true, "genleptoneg2pt");
    
  // ANGULAR MATCHING - BEFORE AND AFTER PROMPT CORRECTION
 vector<double> bar_deta{-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.45,-0.4,-0.35,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1,1.4};
  int bar_deta_nbins = 49;
  vector<double> bar_qdphi{-0.8,-0.7,-0.6,-0.5,-0.45,-0.4,-0.35,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.4,0.45,0.5};
  int bar_qdphi_nbins = bar_qdphi.size()-1;
  seltext[0] = "gen e: p_{T}>10 GeV, |#eta_{prompt}|<1.479";
  seltext[1] = "#geq1 unseeded e/#gamma";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kRed-5);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selbarAnoselusgeneltrigebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kRed-2);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "dEta", bar_deta_nbins, &bar_deta[0], true, false, false, (float []){8e-5,0.9}, (float []){0.525,0.72,0.775,0.97}, (float []){-0.9,0.1}, true, "#Delta#eta(gen,SC)", "normalized events / unit", true, "genbasicptgt10selbarAnoselusgeneltrigebus_dEta");
  //comparesamevariable(file, name, "qdPhi", bar_qdphi_nbins, &bar_qdphi[0], true, false, false, (float []){2e-4,0.9}, (float []){0.15,0.72,0.4,0.97}, (float []){-0.75,0.02}, true, "charge#times#Delta#phi(gen,SC)", "normalized events / unit", true, "genbasicptgt10selbarAnoselusgeneltrigebus_qdPhi");
  //comparesamevariable(file, name, "dPromptEta", bar_deta_nbins, &bar_deta[0], true, false, false, (float []){8e-5,0.9}, (float []){0.525,0.72,0.775,0.97}, (float []){0.2,0.01}, true, "#Delta#eta_{prompt}(gen,SC)", "normalized events / unit", true, "genbasicptgt10selbarAnoselusgeneltrigebus_dPromptEta");
  //comparesamevariable(file, name, "qdPromptPhi", bar_qdphi_nbins, &bar_qdphi[0], true, false, false, (float []){2e-4,0.9}, (float []){0.15,0.72,0.4,0.97}, (float []){-0.75,0.02}, true, "charge#times#Delta#phi_{prompt}(gen,SC)", "normalized events / unit", true, "genbasicptgt10selbarAnoselusgeneltrigebus_qdPromptPhi");
  
  seltext[0] = "gen e: p_{T}>10 GeV, 1.479<|#eta_{prompt}|<2.4";
  seltext[1] = "#geq1 unseeded e/#gamma";  

  file.clear();
  name.clear();
  legend.clear();
  coloropt.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();

  file.push_back(sig3cmfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig30cmfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig1mfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kRed-5);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);

  file.push_back(sig3mfile);
  name.push_back("genbasicptgt10selecAnoselusgeneltrigeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kRed-2);
  histtype.push_back("hist same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("l");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "dEta", bar_deta_nbins, &bar_deta[0], true, false, false, (float []){8e-5,0.9}, (float []){0.525,0.72,0.775,0.97}, (float []){-0.9,0.1}, true, "#Delta#eta(gen,SC)", "normalized events / unit", true, "genbasicptgt10selecAnoselusgeneltrigeeus_dEta");
  //comparesamevariable(file, name, "qdPhi", bar_qdphi_nbins, &bar_qdphi[0], true, false, false, (float []){2e-4,0.9}, (float []){0.15,0.72,0.4,0.97}, (float []){-0.75,0.02}, true, "charge#times#Delta#phi(gen,SC)", "normalized events / unit", true, "genbasicptgt10selecAnoselusgeneltrigeeus_qdPhi");  
  //comparesamevariable(file, name, "dPromptEta", bar_deta_nbins, &bar_deta[0], true, false, false, (float []){8e-5,0.9}, (float []){0.525,0.72,0.775,0.97}, (float []){0.2,0.01}, true, "#Delta#eta_{prompt}(gen,SC)", "normalized events / unit", true, "genbasicptgt10selecAnoselusgeneltrigeeus_dPromptEta");
  //comparesamevariable(file, name, "qdPromptPhi", bar_qdphi_nbins, &bar_qdphi[0], true, false, false, (float []){2e-4,0.9}, (float []){0.15,0.72,0.4,0.97}, (float []){-0.75,0.02}, true, "charge#times#Delta#phi_{prompt}(gen,SC)", "normalized events / unit", true, "genbasicptgt10selecAnoselusgeneltrigeeus_qdPromptPhi");

  // KINEMATIC PROPERTIES OF THE SIGNAL LEPTON BEFORE AND AFTER GENMCH COMPARED TO THE DATA AND PROMPT MC
  vector<double> genbasicselptgt10_leadebsinin_binsleadebsinin{0.0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03};
  int genbasicselptgt10_leadebsinin_nbinsleadebsinin = genbasicselptgt10_leadebsinin_binsleadebsinin.size()-1;

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
  name.push_back("basicselusrecoebus");
  legend.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);

  file.push_back(dyfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  file.push_back(sig3cmfile);
  name.push_back("basicselusrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig3cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("c#tau = 3 cm, gen. matched");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  legendEntries = legend;    

  seltext[0] = "#geq2 unseeded e/#gamma";
  seltext[1] = "p_{T}>10 GeV, |#eta|<1.48, No ID";  
  //comparesamevariable(file, name, "subleadegin5x5noiseclnd", genbasicselptgt10_leadebsinin_nbinsleadebsinin, &genbasicselptgt10_leadebsinin_binsleadebsinin[0], true, false, false, (float []){3e-3,1.4}, (float []){0.52,0.65,0.77,0.99}, (float []){0.015,0.015}, true, "e/#gamma_{2} #sigmai#etai#eta", "normalized events / unit", true, "compare_genmchToData_EB_subleadegin5x5noiseclnd");
  seltext[0] = "#geq1 unseeded e/#gamma";
  seltext[1] = "p_{T}>10 GeV, |#eta|<1.48, No ID";  
  //comparesamevariable(file, name, "leadegin5x5noiseclnd", genbasicselptgt10_leadebsinin_nbinsleadebsinin, &genbasicselptgt10_leadebsinin_binsleadebsinin[0], true, false, false, (float []){3e-3,1.4}, (float []){0.52,0.65,0.77,0.99}, (float []){0.015,0.015}, true, "e/#gamma_{1} #sigmai#etai#eta", "normalized events / unit", true, "compare_genmchToData_EB_leadegin5x5noiseclnd");
  
  vector<double> genbasicselptgt10_leadeesinin_binsleadeesinin{0.0,0.004,0.008,0.012,0.016,0.02,0.024,0.028,0.032,0.036,0.04,0.044,0.048,0.052,0.056,0.06,0.064,0.068,0.072,0.076,0.08};
  int genbasicselptgt10_leadeesinin_nbinsleadeesinin = genbasicselptgt10_leadeesinin_binsleadeesinin.size()-1;
  seltext[0] = "#geq1 unseeded e/#gamma";
  seltext[1] = "p_{T}>10 GeV, 1.479<|#eta|<2.5, No ID";  

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
  name.push_back("basicselusrecoeeus");
  legend.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
    
  file.push_back(dyfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  file.push_back(sig3cmfile);
  name.push_back("basicselusrecoeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed-7);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig3cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("c#tau = 3 cm, gen. matched");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "subleadegin5x5noiseclnd", genbasicselptgt10_leadeesinin_nbinsleadeesinin, &genbasicselptgt10_leadeesinin_binsleadeesinin[0], true, false, false, (float []){3e-3,9}, (float []){0.52,0.65,0.77,0.99}, (float []){0.05,0.015}, true, "e/#gamma_{2} #sigmai#etai#eta", "normalized events / unit", true, "compare_genmchToData_EE_subleadegin5x5noiseclnd");
  //comparesamevariable(file, name, "leadegin5x5noiseclnd", genbasicselptgt10_leadeesinin_nbinsleadeesinin, &genbasicselptgt10_leadeesinin_binsleadeesinin[0], true, false, false, (float []){3e-3,9}, (float []){0.52,0.65,0.77,0.99}, (float []){0.05,0.015}, true, "e/#gamma_{1} #sigmai#etai#eta", "normalized events / unit", true, "compare_genmchToData_EE_leadegin5x5noiseclnd");

  // CHECK THE Z WINDOW COMPOSITION - PLOTS NOT FINAL, CONCLUSION - Z WINDOW CANNOT BE OBTAINED BY JUST EVENT SELECTIONS ON ECAL OBJECTS
  vector<double> genbasicselptgt10_leadeb_binsmee{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120};
  int genbasicselptgt10_leadeb_nbinsmee = genbasicselptgt10_leadeb_binsmee.size()-1;
  vector<double> leadsubleadeb_binsdR{-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0};
  int leadsubleadeb_nbinsdR = leadsubleadeb_binsdR.size()-1;

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
  name.push_back("basicselusrecoebus");
  legend.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dyfile);
  name.push_back("basicselusrecoebus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  file.push_back(datafile);
  name.push_back("selelevetoidusrecoebus");
  legend.push_back("2018 data");
  coloropt.push_back(kGray);
  histtype.push_back("p e same");
  markerstyle.push_back(24);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dyfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue-2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoebus");
  legend.push_back("2018 data");
  coloropt.push_back(kGray);
  histtype.push_back("p e same");
  markerstyle.push_back(22);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(datafile);
  name.push_back("ptgt25neq2noidusrecoebus");
  legend.push_back("2018 data");
  coloropt.push_back(kGray);
  histtype.push_back("p e same");
  markerstyle.push_back(26);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dyfile);
  name.push_back("ptgt25neq2noidusrecoebus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue-4);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  legendEntries = legend;    
  //comparesamevariable(file, name, "leadsubleadM", genbasicselptgt10_leadeb_nbinsmee, &genbasicselptgt10_leadeb_binsmee[0], true, false, false, (float []){1, 3e5}, (float []){-1,0.65,0.77,0.99}, (float []){-1,0.015}, false, "e/#gamma_{1} M(e_{1},e_{2}) [GeV]", "events / 1 GeV", true, "basicTocut1_EB_mee");
  //comparesamevariable(file, name, "leadsubleaddR", leadsubleadeb_nbinsdR, &leadsubleadeb_binsdR[0], true, false, false, (float []){1, 3e7}, (float []){-1,0.65,0.77,0.99}, (float []){-1,0.015}, false, "e/#gamma_{1} #DeltaR(e_{1},e_{2})", "events / 0.1 units", true, "basicTocut1_EB_dR");
  
  vector<double> promptdata_binspt{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120};
  int promptdata_nbinspt = promptdata_binspt.size()-1;
  vector<double> promptdata_binseta{-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
  int promptdata_nbinseta = promptdata_binseta.size()-1;

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
  name.push_back("basicselusrecoegus");
  legend.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dyfile);
  name.push_back("basicselusrecoegus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  file.push_back(datafile);
  name.push_back("selelevetoidusrecoegus");
  legend.push_back("2018 data");
  coloropt.push_back(kGray);
  histtype.push_back("p e same");
  markerstyle.push_back(24);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dyfile);
  name.push_back("genptgt10Abasicselusgenmchrecoegus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue-2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  file.push_back(datafile);
  name.push_back("selelevetozwindidusrecoegus");
  legend.push_back("2018 data");
  coloropt.push_back(kGray);
  histtype.push_back("p e same");
  markerstyle.push_back(22);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(datafile);
  name.push_back("ptgt25neq2noidusrecoegus");
  legend.push_back("2018 data");
  coloropt.push_back(kGray);
  histtype.push_back("p e same");
  markerstyle.push_back(26);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dyfile);
  name.push_back("ptgt25neq2noidusrecoegus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue-4);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  legendEntries = legend;    
  //comparesamevariable(file, name, "leadegpt", promptdata_nbinspt, &promptdata_binspt[0], true, false, false, (float []){1,3e5}, (float []){-1,0.65,0.77,0.99}, (float []){-1,0.015}, false, "e/#gamma_{1} p_{T} [GeV]", "events / 1 GeV", true, "basicTocut1_EB_eg1pt");
  //comparesamevariable(file, name, "subleadegpt", promptdata_nbinspt, &promptdata_binspt[0], true, false, false, (float []){1,3e5}, (float []){-1,0.65,0.77,0.99}, (float []){-1,0.015}, false, "e/#gamma_{2} p_{T} [GeV]", "events / 1 GeV", true, "basicTocut1_EB_eg2pt");
  //comparesamevariable(file, name, "leadegeta", promptdata_nbinseta, &promptdata_binseta[0], true, false, false, (float []){1,3e5}, (float []){-1,0.65,0.77,0.99}, (float []){-1,0.015}, false, "e/#gamma_{1} #eta", "events / 0.1", true, "basicTocut1_EB_eg1eta");
  //comparesamevariable(file, name, "subleadegeta", promptdata_nbinseta, &promptdata_binseta[0], true, false, false, (float []){1,3e5}, (float []){-1,0.65,0.77,0.99}, (float []){-1,0.015}, false, "e/#gamma_{2} #eta", "events / 0.1", true, "basicTocut1_EB_eg2eta");
  
  // COMPARE VARIABLES BETWEEN THE SIGNAL BENCHMARKS, PROMPT MC, AND DATA
  vector<double> genbasicselptgt10_leadeb_binshoe{0.0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.06,0.08,0.1,0.14,0.18,0.22,0.26,0.3,0.35,0.4};
  int genbasicselptgt10_leadeb_nbinshoe = genbasicselptgt10_leadeb_binshoe.size()-1;
  vector<double> genbasicselptgt10_leadeb_binspxlmch{-2,-1,0};
  int genbasicselptgt10_leadeb_nbinspxlmch = genbasicselptgt10_leadeb_binspxlmch.size()-1;
  vector<double> genbasicselptgt10_leadeb_binsmhits{-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
  int genbasicselptgt10_leadeb_nbinsmhits = genbasicselptgt10_leadeb_binsmhits.size()-1;
  seltext[0] = "#geq1 unseeded e/#gamma";
  seltext[1] = "p_{T}>10 GeV, |#eta|<1.48, No ID";  

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
  name.push_back("basicselusrecoebus");
  legend.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dyfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  file.push_back(sig3cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig1mfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig3mfile);
  name.push_back("genptgt10Abasicselusgenmchrecoebus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "leadegin5x5noiseclnd", genbasicselptgt10_leadebsinin_nbinsleadebsinin, &genbasicselptgt10_leadebsinin_binsleadebsinin[0], true, false, false, (float []){3e-3,9}, (float []){0.52,0.6,0.77,0.99}, (float []){0.001,2}, true, "e/#gamma_{1} #sigma_{i#etai#eta}", "normalized events / 0.001 unit", true, "genptgt10Abasicselusgenmchrecoebus_EB_leadegin5x5noiseclnd");
  //comparesamevariable(file, name, "leadeghovereoversupcluse", genbasicselptgt10_leadeb_nbinshoe, &genbasicselptgt10_leadeb_binshoe[0], true, false, false, (float []){2e-4,2}, (float []){0.52,0.6,0.77,0.99}, (float []){0.05,0.3}, true, "e/#gamma_{1} H/E", "normalized events / unit", true, "genptgt10Abasicselusgenmchrecoebus_EB_leadeghovere");
  //comparesamevariable(file, name, "leadegpixelmchvar_s2", genbasicselptgt10_leadeb_nbinspxlmch, &genbasicselptgt10_leadeb_binspxlmch[0], true, false, true, (float []){5e-5,5e5}, (float []){0.52,0.6,0.77,0.99}, (float []){0.05,0.3}, true, "e/#gamma_{1} pixel match s2", "normalized events", true, "genptgt10Abasicselusgenmchrecoebus_EB_leadegpxlmchvar");
  //comparesamevariable(file, name, "leadegmhits", genbasicselptgt10_leadeb_nbinsmhits, &genbasicselptgt10_leadeb_binsmhits[0], true, false, false, (float []){2e-5,2}, (float []){0.52,0.6,0.77,0.99}, (float []){0.05,0.3}, true, "e/#gamma_{1} missing hits", "normalized events / 1 unit", true, "genptgt10Abasicselusgenmchrecoebus_EB_leadegmhits");
  
  vector<double> genbasicselptgt10_leadee_binshoe{0.0,0.01,0.02,0.03,0.04,0.06,0.08,0.1,0.14,0.18,0.22,0.26,0.3,0.35,0.4};
  int genbasicselptgt10_leadee_nbinshoe = genbasicselptgt10_leadee_binshoe.size()-1;
  vector<double> genbasicselptgt10_leadee_binspxlmch{-2,-1,0};
  int genbasicselptgt10_leadee_nbinspxlmch = genbasicselptgt10_leadee_binspxlmch.size()-1;
  vector<double> genbasicselptgt10_leadee_binsmhits{-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
  int genbasicselptgt10_leadee_nbinsmhits = genbasicselptgt10_leadee_binsmhits.size()-1;
  seltext[0] = "#geq1 unseeded e/#gamma";
  seltext[1] = "p_{T}>10 GeV, 1.48<|#eta|<2.5, No ID";  

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
  name.push_back("basicselusrecoeeus");
  legend.push_back("2018 data");
  coloropt.push_back(kBlack);
  histtype.push_back("p e");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
    
  file.push_back(dyfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("DY #rightarrow ee MC");
  coloropt.push_back(kBlue);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(24);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);

  file.push_back(sig3cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 cm");
  coloropt.push_back(kRed+3);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig30cmfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 30 cm");
  coloropt.push_back(kRed);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig1mfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 1 m");
  coloropt.push_back(kOrange+2);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  file.push_back(sig3mfile);
  name.push_back("genptgt10Abasicselusgenmchrecoeeus");
  legend.push_back("#chi^{#pm} #rightarrow #chi^{0}l#nu, c#tau = 3 m");
  coloropt.push_back(kOrange+1);
  histtype.push_back("hist e1 same");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("le");
  scale.push_back(1);
  
  legendEntries = legend;    
  //comparesamevariable(file, name, "leadegin5x5noiseclnd", genbasicselptgt10_leadeesinin_nbinsleadeesinin, &genbasicselptgt10_leadeesinin_binsleadeesinin[0], true, false, false, (float []){3e-3,9}, (float []){0.525,0.6,0.775,0.99}, (float []){0.004, 2}, true, "e/#gamma_{1} #sigma_{i#etai#eta}", "normalized events / 0.004 unit", true, "genptgt10Abasicselusgenmchrecoeeus_EE_leadegin5x5noiseclnd");
  //comparesamevariable(file, name, "leadeghovereoversupcluse", genbasicselptgt10_leadee_nbinshoe, &genbasicselptgt10_leadee_binshoe[0], true, false, false, (float []){7e-4,2}, (float []){0.52,0.6,0.77,0.99}, (float []){0.02,0.5}, true, "e/#gamma_{1} H/E", "normalized events / 1 unit", true, "genptgt10Abasicselusgenmchrecoeeus_EE_leadeghovere");
  //comparesamevariable(file, name, "leadegpixelmchvar_s2", genbasicselptgt10_leadee_nbinspxlmch, &genbasicselptgt10_leadee_binspxlmch[0], true, false, true, (float []){2e-5,2}, (float []){0.52,0.6,0.77,0.99}, (float []){0.05,0.3}, true, "e/#gamma_{1} pixel match s2", "normalized events", true, "genptgt10Abasicselusgenmchrecoeeus_EE_leadegpxlmchvar");
  //comparesamevariable(file, name, "leadegmhits", genbasicselptgt10_leadee_nbinsmhits, &genbasicselptgt10_leadee_binsmhits[0], true, false, false, (float []){2e-5,2}, (float []){0.52,0.6,0.77,0.99}, (float []){0.05,0.3}, true, "e/#gamma_{1} missing hits", "normalized events / 1 unit", true, "genptgt10Abasicselusgenmchrecoeeus_EE_leadegmhits");


}
