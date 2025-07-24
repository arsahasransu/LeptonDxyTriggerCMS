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
        'file:/pnfs/pp.rl.ac.uk/data/cms/store/user/asahasra/EGamma0/EGamma0_2023C_250601ReRECO/250601_092049/0000/EGamma0_2023_rereco_miniaod_100.root',
        # 'root://cms-xrd-global.cern.ch///store/data/Run2023C/EGamma0/MINIAOD/22Sep2023_v4-v1/410000/c607dcf9-2c85-4d9c-92fa-994427f93c35.root'
    )
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("EGamma0_EXOLLPTRG_Nano.root")
                               )

process.demo = cms.EDAnalyzer('TriggerAnalyzerMiniAOD2024',
                                isMC = cms.bool(False)
)


process.p = cms.Path(process.demo)
