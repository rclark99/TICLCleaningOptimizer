import FWCore.ParameterSet.Config as cms

from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
from Configuration.AlCa.GlobalTag import GlobalTag

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
from FWCore.ParameterSet.VarParsing import VarParsing
from utils import read_csv

# HLT menu
from HLTrigger.Configuration.HLT_75e33_cff import *

# VarParsing instance
options = VarParsing("analysis")

options.register(
    "parametersFile",
    "default/default_params.csv",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of parameters file (CSV: one row per candidate point).",
)

options.register(
    "nEvents",
    100,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Number of events",
)

options.register(
    "outputFile",
    "jet_validation_metrics.root",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output ROOT file",
)

options.parseArguments()

process = cms.Process("RECO3", Phase2C17I13M9)

process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T33", "")

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(options.nEvents),
    output=cms.optional.untracked.allowed(cms.int32, cms.PSet),
)

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
        ["root://cmseos.fnal.gov//store/user/rclark1/TICL/Optimize/step2_QCD_82516435_99.root"]
    ),
    secondaryFileNames=cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(0),
)

process.options = cms.untracked.PSet(
    IgnoreCompletely=cms.untracked.vstring(),
    Rethrow=cms.untracked.vstring(),
    accelerators=cms.untracked.vstring("*"),
    allowUnscheduled=cms.obsolete.untracked.bool,
    canDeleteEarly=cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules=cms.untracked.bool(True),
    dumpOptions=cms.untracked.bool(False),
    emptyRunLumiMode=cms.obsolete.untracked.string,
    eventSetup=cms.untracked.PSet(
        forceNumberOfConcurrentIOVs=cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs=cms.untracked.uint32(0),
    ),
    fileMode=cms.untracked.string("FULLMERGE"),
    forceEventSetupCacheClearOnNewRun=cms.untracked.bool(False),
    holdsReferencesToDeleteEarly=cms.untracked.VPSet(),
    makeTriggerResults=cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly=cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(0),
    numberOfConcurrentRuns=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0),
    numberOfThreads=cms.untracked.uint32(4),
    printDependencies=cms.untracked.bool(False),
    sizeOfStackForThreadsInKB=cms.optional.untracked.uint32,
    throwIfIllegalParameter=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(True),
)

# parameter scan
params = read_csv(options.parametersFile)
print(params)

cleaningModules   = []
candidateModules  = []
pfModules         = []
puppiModules      = []
ak4JetModules     = []
ak4CorrJetModules = []
jetTesterModules  = []
metricsModules    = []

for i, p in enumerate(params):

    betaContamMin = float(p[0])
    R0            = float(p[1])
    sigmaT        = float(p[2])
    sigmaDR       = float(p[3])
    sigmaZ        = float(p[4])
    tAbsCut       = float(p[5])
    zAbsCut       = float(p[6])
    tPower        = float(p[7])
    drPower       = float(p[8])
    zPower        = float(p[9])
    wmin          = float(p[10])

    # trackster cleaning 
    setattr(process, f"hltTiclTracksterCleaning{i}",
        process.hltTiclTracksterCleaning.clone(
            clue3DTracksters        = cms.InputTag("hltTiclTrackstersCLUE3DHigh"),
            linkedTracksters        = cms.InputTag("hltTiclTracksterLinks"),
            clue3DInLinkedIndices   = cms.InputTag("hltTiclTracksterLinks", "linkedTracksterIdToInputTracksterId"),
            cleaner = process.hltTiclTracksterCleaning.cleaner.clone(
                betaContamMin = cms.double(betaContamMin),
                R0            = cms.double(R0),
                sigmaT        = cms.double(sigmaT),
                sigmaDR       = cms.double(sigmaDR),
                sigmaZ        = cms.double(sigmaZ),
                tAbsCut       = cms.double(tAbsCut),
                zAbsCut       = cms.double(zAbsCut),
                tPower        = cms.double(tPower),
                drPower       = cms.double(drPower),
                zPower        = cms.double(zPower),
                wmin          = cms.double(wmin),
            ),
            labelLinkedOut = cms.string("cleanedLinkedTracksters"),
            labelMapOut    = cms.string("cleanedLinkedTracksterIdToInputTracksterId"),
        )
    )
    cleaningModules.append(getattr(process, f"hltTiclTracksterCleaning{i}"))

    # TICL candidate 
    setattr(process, f"hltTiclCandidate{i}",
        process.hltTiclCandidate.clone(
            general_tracksterlinks_collections = cms.VInputTag(
                cms.InputTag(f"hltTiclTracksterCleaning{i}", "cleanedLinkedTracksterIdToInputTracksterId")
            ),
            general_tracksters_collections = cms.VInputTag(
                cms.InputTag(f"hltTiclTracksterCleaning{i}", "cleanedLinkedTracksters")
            ),
        )
    )
    candidateModules.append(getattr(process, f"hltTiclCandidate{i}"))

    # PF from TICL
    setattr(process, f"hltPfTICL{i}",
        process.hltPfTICL.clone(
            ticlCandidateSrc = cms.InputTag(f"hltTiclCandidate{i}")
        )
    )
    pfModules.append(getattr(process, f"hltPfTICL{i}"))

    # PUPPI 
    setattr(process, f"hltPFPuppi{i}",
        process.hltPFPuppi.clone(
            candName = cms.InputTag(f"hltPfTICL{i}")
        )
    )
    puppiModules.append(getattr(process, f"hltPFPuppi{i}"))

    # AK4 jets from PF candidates with puppi weights
    setattr(process, f"hltAK4PFPuppiJets{i}",
        process.hltAK4PFPuppiJets.clone(
            src       = cms.InputTag(f"hltPfTICL{i}"),
            srcWeights= cms.InputTag(f"hltPFPuppi{i}")
        )
    )
    ak4JetModules.append(getattr(process, f"hltAK4PFPuppiJets{i}"))

    # corrected jets 
    setattr(process, f"hltAK4PFPuppiJetsCorrected{i}",
        process.hltAK4PFPuppiJetsCorrected.clone(
            src        = cms.InputTag(f"hltAK4PFPuppiJets{i}"),
            correctors = cms.VInputTag("hltAK4PFPuppiJetCorrector"),
        )
    )
    ak4CorrJetModules.append(getattr(process, f"hltAK4PFPuppiJetsCorrected{i}"))

    # JetTester
    setattr(process, f"hltJetAnalyzerAK4PFPuppi{i}",
        process.hltJetAnalyzerAK4PFPuppi.clone(
            src            = cms.InputTag(f"hltAK4PFPuppiJets{i}"),
            JetCorrections = cms.InputTag("hltAK4PFPuppiJetCorrector"),
            srcGen         = cms.InputTag("ak4GenJetsNoNu"),
        )
    )
    jetTesterModules.append(getattr(process, f"hltJetAnalyzerAK4PFPuppi{i}"))

    # validation metrics
    setattr(process, f"jetValidationMetrics{i}",
        cms.EDAnalyzer(
            "JetValidationMetrics",
            jets               = cms.InputTag(f"hltAK4PFPuppiJetsCorrected{i}"),
            genjets            = cms.InputTag("ak4GenJetsNoNu"),
            recoJetPtThreshold = cms.double(30.0),
            matchGenPtThreshold= cms.double(20.0),
            RThreshold         = cms.double(0.4),
            mightGet           = cms.optional.untracked.vstring,
        )
    )
    metricsModules.append(getattr(process, f"jetValidationMetrics{i}"))

# output
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))

process.optimizeTask = cms.Task(
    *cleaningModules,
    *candidateModules,
    *pfModules,
    *puppiModules,
    *ak4JetModules,
    *ak4CorrJetModules,
    *jetTesterModules,
    *metricsModules,
)

process.TICLCleaningOptimize = cms.Path(cms.Sequence())
process.TICLCleaningOptimize.associate(process.optimizeTask)

process.schedule = cms.Schedule(process.TICLCleaningOptimize)

# customizations
process = setCrossingFrameOn(process)
process = customiseLogErrorHarvesterUsingOutputCommands(process)
process = customiseEarlyDelete(process)
