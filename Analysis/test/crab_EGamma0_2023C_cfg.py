from CRABClient.UserUtilities import config
config = config()

config.General.requestName = "EGamma0_2023C_crabRun250204_Ntuplizer_"
config.General.workArea = "EGamma0_2023C_crabRunNtuplizer"
config.General.transferLogs = False
config.General.transferOutputs = True

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "TriggerAnalyzerConfig_cfg.py"

config.Data.inputDataset= "/EGamma0/Run2023C-22Sep2023_v4-v1/MINIAOD"
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.totalUnits = 100
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_eraC_367095_368823_Golden.json'
config.Data.publication = False
config.Data.outputDatasetTag = "EGamma0_2023C_crabRun250204_Ntuplizer"

config.Site.storageSite = "T2_UK_SGrid_RALPP"
