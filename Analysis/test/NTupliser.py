import FWCore.ParameterSet.Config as cms

process = cms.Process("DEMO")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/pnfs/pp.rl.ac.uk/data/cms/store/user/asahasra/SingletTripletHDMToDisplacedL_M200deltaM20ctau1m_TuneCP5_14TeV-madgraph-pythia8/STHDM_M200deltaM20ctau1m_Run3Winter24MiniAOD_250418/250417_165811/0000/STHDM_Run3Winter24Miniaod_1.root'
    )
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("NTuples.root")
                               )

process.demo = cms.EDAnalyzer('TriggerAnalyzerMiniAOD2024'
)


process.p = cms.Path(process.demo)
