from CRABClient.UserUtilities import config
config = config()

config.General.requestName = "JetMET2023C_crabRun240606_Ntuplizer_"
config.General.workArea = "JetMET2023C_crabRunNtuplizer"
config.General.transferLogs = False
config.General.transferOutputs = True

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "TriggerAnalyzerConfig_cfg.py"
config.JobType.inputFiles = ['Cert_Collisions2023_eraC_367095_368823_Golden.json']

config.Data.inputDataset= "/JetMET0/Run2023C-19Dec2023-v1/MINIAOD"
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.lumiMask = 'Cert_Collisions2023_eraC_367095_368823_Golden.json'
config.Data.publication = False
config.Data.outputDatasetTag = "JetMET2023C_crabRun240606_Ntuplizer"

config.Site.storageSite = "T2_UK_SGrid_RALPP"
