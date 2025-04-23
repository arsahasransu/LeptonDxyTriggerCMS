import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('inpfilelist',
                    'input.txt',
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.string,
                    'Input file list')
options.register('outfile',
                    'output.root',
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.string,
                    'Output file name')
options.parseArguments()

process = cms.Process("DEMO")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Read input file names from a text file
with open(options.inpfilelist, "r") as f:
    file_names = [line.strip() for line in f if line.strip()]

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring( *file_names )
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(options.outfile),
                               )

process.demo = cms.EDAnalyzer('TriggerAnalyzerMiniAOD2024'
)


process.p = cms.Path(process.demo)
