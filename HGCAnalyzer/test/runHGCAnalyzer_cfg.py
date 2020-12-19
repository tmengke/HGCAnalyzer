import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Phase2C11_cff import Phase2C11
process = cms.Process('HGCGeomAnalysis',Phase2C11)
process.load('Configuration.Geometry.GeometryExtended2026D62_cff')
process.load('Configuration.Geometry.GeometryExtended2026D62Reco_cff')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')    
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoLocalCalo.HGCalRecProducers.hgcalRecHitMapProducer_cfi')
process.load('SimCalorimetry.HGCalAssociatorProducers.layerClusterAssociatorByEnergyScore_cfi')

process.task = cms.Task(process.hgcalRecHitMapProducer)
process.lcscore =  cms.Task(process.layerClusterAssociatorByEnergyScore)
#process.sq = cms.Sequence(process.task)

from Configuration.AlCa.GlobalTag import GlobalTag
#Global Tag used for production in
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

options = VarParsing.VarParsing ('analysis')

options.inputFiles= 'file:/uscms_data/d3/tmengke/HGCDPG/CMSSW_11_2_0_pre6/src/29090.0_SingleGammaPt25Eta1p7_2p7+2026D62+SingleGammaPt25Eta1p7_2p7_GenSimHLBeamSpot+DigiTrigger+RecoGlobal+HARVESTGlobal/step3.root'
options.outputFile = 'test.root'
options.maxEvents = -1
 
options.parseArguments()

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile)
                                   )

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
		options.inputFiles
	)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.hgcAnalyzer = cms.EDAnalyzer('HGCAnalyzer',
				    setZside = cms.untracked.int32(1)  # 1 zplus, -1 zminus, 0 both
				    )
process.p = cms.Path(process.hgcAnalyzer,process.task,process.lcscore)
