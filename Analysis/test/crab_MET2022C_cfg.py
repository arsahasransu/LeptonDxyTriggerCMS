from CRABClient.UserUtilities import config, getUsernameFromCRIC
config = config()

config.General.requestName = "MET2022C_crabRun230210Ntuplizer_"
config.General.workArea = "MET2022C_crabRunNtuplizer"
config.General.transferLogs = False
config.General.transferOutputs = True

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "TriggerAnalyzerConfig_cfg.py"
config.JobType.inputFiles = ['Cert_Collisions2022_355100_362760_Golden.json']

config.Data.inputDataset= "/MET/Run2022C-10Dec2022-v2/MINIAOD"
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
config.Data.LumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'
config.Data.publication = False
config.Data.outputDatasetTag = "MET2022C_crabRun230210Ntuplizer"
config.Data.ignoreLocality = True

config.Site.whitelist = ['T2_US*','T2_CH*']
config.Site.storageSite = "T2_BE_IIHE"
