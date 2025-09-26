import FWCore.ParameterSet.Config as cms

process = cms.Process("TriggerAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch///store/mc/Run3Winter25MiniAOD/DYto2L-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/142X_mcRun3_2025_realistic_v9-v3/110000/b1c66854-0a3a-4e60-9454-07b44f63e2cc.root'
    )
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("DYTo2L_Run3Winter25_EXOLLPTRG_Nano.root")
                               )

process.demo = cms.EDAnalyzer('TriggerAnalyzerMiniAOD2024',
                              isMC = cms.bool(False)
)


process.p = cms.Path(process.demo)
