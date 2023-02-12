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
        'file:/afs/cern.ch/work/a/asahasra/private/Run3DataWithNewTriggers/CMSSW_12_6_3/src/7b8f72a2-5390-40e2-b9b8-db9d726de30e.root'
    )
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("DiPhoton10_trigNtuples.root")
                               )

process.demo = cms.EDAnalyzer('TriggerAnalyzerMiniAOD'
)


process.p = cms.Path(process.demo)
