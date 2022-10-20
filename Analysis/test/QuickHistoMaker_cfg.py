import FWCore.ParameterSet.Config as cms

process = cms.Process("TriggerAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:19856a6e-7ac9-4b16-bbf4-12dd69ee6589.root'
        'file:output_HLT2_after.root'
    )
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("DiPhoton10_trigNtuples_after.root")
                               )

process.demo = cms.EDAnalyzer('QuickHistoMaker'
)


process.p = cms.Path(process.demo)
