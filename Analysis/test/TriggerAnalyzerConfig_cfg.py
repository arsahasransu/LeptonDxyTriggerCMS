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
        'file:/eos/user/a/asahasra/TestData/JetMET0_2023C_19Dec2023_v1_ReReco_MINIAOD.root'
    )
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("DiPhoton10_trigNtuples.root")
                               )

process.demo = cms.EDAnalyzer('TriggerAnalyzerMiniAOD2024'
)


process.p = cms.Path(process.demo)
