import FWCore.ParameterSet.Config as cms

process = cms.Process("TriggerAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:125a2199-8dfe-499f-b616-a28fe4d018c3.root'
    )
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("DiPhoton10_trigNtuples.root")
                               )

process.demo = cms.EDAnalyzer('TriggerAnalyzerRAWMiniAOD'
)


process.p = cms.Path(process.demo)
