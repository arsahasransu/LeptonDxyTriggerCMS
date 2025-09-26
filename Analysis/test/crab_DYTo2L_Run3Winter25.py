from CRABClient.UserUtilities import config
config = config()

config.General.requestName = "DYTo2L_crabRun250923_Ntuplizer_"
config.General.workArea = "DYTo2L_2025C_crabRunNtuplizer"
config.General.transferLogs = False
config.General.transferOutputs = True

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "TriggerAnalyzerConfig_cfg.py"

config.Data.inputDataset= "/DYto2L-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter25MiniAOD-142X_mcRun3_2025_realistic_v9-v3/MINIAODSIM"
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.totalUnits = 1000
config.Data.publication = False
config.Data.outputDatasetTag = "DYTo2L_2025C_crabRun250923_Ntuplizer"

config.Site.storageSite = "T2_UK_SGrid_RALPP"
