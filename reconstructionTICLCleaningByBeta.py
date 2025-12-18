import FWCore.ParameterSet.Config as cms

from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
from Configuration.AlCa.GlobalTag import GlobalTag

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
from FWCore.ParameterSet.VarParsing import VarParsing
from utils import read_csv

# HLT menu
# from HLTrigger.Configuration.HLT_75e33_cff import *

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

# options.register(
#     "outputFile",
#     "jet_validation_metrics.root",
#     VarParsing.multiplicity.singleton,
#     VarParsing.varType.string,
#     "Output ROOT file",
# )

options.parseArguments()

process = cms.Process("HLT3", Phase2C17I13M9)

process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TrackAssociator.default_cfi")
process.load("TrackingTools.MaterialEffects.Propagators_cff")
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.SimPhase2L1GlobalTriggerEmulator_cff')
process.load('L1Trigger.Configuration.Phase2GTMenus.SeedDefinitions.step1_2024.l1tGTMenu_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_75e33_cff')

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
    numberOfThreads=cms.untracked.uint32(1),
    printDependencies=cms.untracked.bool(False),
    sizeOfStackForThreadsInKB=cms.optional.untracked.uint32,
    throwIfIllegalParameter=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(True),
)

process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


# parameter scan
params = read_csv(options.parametersFile)
print(f"Loaded {len(params)} parameter points")

cleaningModules   = []
candidateModules  = []
pfModules         = []
puppiModules      = []
ak4JetModules     = []
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
        cms.EDProducer('TracksterCleaningProducer',
                clue3DTracksters        = cms.InputTag("hltTiclTrackstersCLUE3DHigh"),
                linkedTracksters        = cms.InputTag("hltTiclTracksterLinks"),
                clue3DInLinkedIndices   = cms.InputTag("hltTiclTracksterLinks", "linkedTracksterIdToInputTracksterId"),
                algo_verbosity = cms.int32(0),
                labelLinkedOut = cms.string('cleanedLinkedTracksters'),
                labelMapOut = cms.string('cleanedLinkedTracksterIdToInputTracksterId'),
                cleaner = cms.PSet(
                type = cms.string('Beta'),
                algo_verbosity = cms.int32(0),
                betaContamMin = cms.double(betaContamMin),
                R0 = cms.double(R0),
                useRawEnergy = cms.bool(True),
                epsE = cms.double(1e-06),
                epsDR = cms.double(1e-06),
                weightMode = cms.bool(True),
                emitDroppedAsStandalone = cms.bool(False),
                zAbsCut = cms.double(zAbsCut),
                tAbsCut = cms.double(tAbsCut),
                sigmaZ = cms.double(sigmaZ),
                sigmaT = cms.double(sigmaT),
                sigmaDR = cms.double(sigmaDR),
                zPower = cms.double(zPower),
                tPower = cms.double(tPower),
                drPower = cms.double(drPower),
                wmin = cms.double(wmin),
                doPruning = cms.bool(False),
                pruneWmin = cms.double(0.01),
                pruneUseSeparateKernels = cms.bool(False),
                sigmaZ_prune = cms.double(12.5),
                sigmaT_prune = cms.double(0.08),
                sigmaDR_prune = cms.double(0.08),
                zPower_prune = cms.double(1.5),
                tPower_prune = cms.double(0.5),
                drPower_prune = cms.double(0.5)
            )
        )
    )
    cleaningModules.append(getattr(process, f"hltTiclTracksterCleaning{i}"))

    # TICL candidate 
    setattr(process, f"hltTiclCandidate{i}",
        cms.EDProducer("TICLCandidateProducer",
            inferenceAlgo = cms.string('TracksterInferenceByPFN'),
            regressionAndPid = cms.bool(True),
            pluginInferenceAlgoTracksterInferenceByPFN = cms.PSet(
            algo_verbosity = cms.int32(0),
            onnxPIDModelPath = cms.FileInPath('RecoHGCal/TICL/data/ticlv5/onnx_models/PFN/linking/id_v0.onnx'),
            onnxEnergyModelPath = cms.FileInPath('RecoHGCal/TICL/data/ticlv5/onnx_models/PFN/linking/energy_v1.onnx'),
            inputNames = cms.vstring(
                'input',
                'input_tr_features'
            ),
            output_en = cms.vstring('enreg_output'),
            output_id = cms.vstring('pid_output'),
            eid_min_cluster_energy = cms.double(2.5),
            eid_n_layers = cms.int32(50),
            eid_n_clusters = cms.int32(10),
            doPID = cms.int32(1),
            doRegression = cms.int32(1),
            type = cms.string('TracksterInferenceByPFN')
            ),
            cutTk = cms.string('1.48 < abs(eta) < 3.0 && pt > 1. && quality("highPurity") && hitPattern().numberOfLostHits("MISSING_OUTER_HITS") < 5'),
            detector = cms.string('HGCAL'),
            egamma_tracksterlinks_collections = cms.VInputTag("hltTiclTracksterLinks"),
            egamma_tracksters_collections = cms.VInputTag("hltTiclTracksterCLUE3DHigh"),
            general_tracksterlinks_collections = cms.VInputTag(cms.InputTag(f"hltTiclTracksterCleaning{i}", "cleanedLinkedTracksterIdToInputTracksterId")),
            general_tracksters_collections = cms.VInputTag("hltTiclTrackstersCLUE3DHigh"),
            interpretationDescPSet = cms.PSet(
                algo_verbosity = cms.int32(0),
                cutTk = cms.string('1.48 < abs(eta) < 3.0 && pt > 1. && quality("highPurity") && hitPattern().numberOfLostHits("MISSING_OUTER_HITS") < 5'),
                delta_tk_ts_interface = cms.double(0.03),
                delta_tk_ts_layer1 = cms.double(0.02),
                timing_quality_threshold = cms.double(0.5),
                type = cms.string('General')
            ),
            layer_clusters = cms.InputTag("hltMergeLayerClusters"),
            layer_clustersTime = cms.InputTag("hltMergeLayerClusters","timeLayerCluster"),
            mightGet = cms.optional.untracked.vstring,
            muons = cms.InputTag("hltPhase2L3Muons"),
            original_masks = cms.VInputTag("hltMergeLayerClusters:InitialLayerClustersMask"),
            propagator = cms.string('PropagatorWithMaterial'),
            timingQualityThreshold = cms.double(0.5),
            timingSoA = cms.InputTag("mtdSoA"),
            tracks = cms.InputTag("hltGeneralTracks"),
            useMTDTiming = cms.bool(False),
            useTimingAverage = cms.bool(False)
        )
    )
    candidateModules.append(getattr(process, f"hltTiclCandidate{i}"))

    # PF from TICL
    setattr(process, f"hltPfTICL{i}",
        cms.EDProducer("PFTICLProducer",
            mightGet = cms.optional.untracked.vstring,
            muonSrc = cms.InputTag("hltPhase2L3Muons"),
            pfMuonAlgoParameters = cms.PSet(
                cosmicRejectionDistance = cms.double(1),
                eventFactorForCosmics = cms.double(10),
                eventFractionForCleaning = cms.double(0.5),
                eventFractionForRejection = cms.double(0.8),
                maxDPtOPt = cms.double(1),
                metFactorForCleaning = cms.double(4),
                metFactorForFakes = cms.double(4),
                metFactorForHighEta = cms.double(25),
                metFactorForRejection = cms.double(4),
                metSignificanceForCleaning = cms.double(3),
                metSignificanceForRejection = cms.double(4),
                minEnergyForPunchThrough = cms.double(100),
                minMomentumForPunchThrough = cms.double(100),
                minPtForPostCleaning = cms.double(20),
                ptErrorScale = cms.double(8),
                ptFactorForHighEta = cms.double(2),
                punchThroughFactor = cms.double(3),
                punchThroughMETFactor = cms.double(4),
                trackQuality = cms.string('highPurity')
            ),
            ticlCandidateSrc = cms.InputTag(f"hltTiclCandidate{i}"),
            timingQualityThreshold = cms.double(0.5),
            trackTimeErrorMap = cms.InputTag("tofPID","sigmat0"),
            trackTimeQualityMap = cms.InputTag("mtdTrackQualityMVA","mtdQualMVA"),
            trackTimeValueMap = cms.InputTag("tofPID","t0"),
            useMTDTiming = cms.bool(False),
            useTimingAverage = cms.bool(False)
        )
    )
    pfModules.append(getattr(process, f"hltPfTICL{i}"))

    # PUPPI 
    setattr(process, f"hltPFPuppi{i}",
        cms.EDProducer("PuppiProducer",
            DeltaZCut = cms.double(0.1),
            DeltaZCutForChargedFromPUVtxs = cms.double(0.2),
            EtaMaxCharged = cms.double(99999.0),
            EtaMaxPhotons = cms.double(2.5),
            EtaMinUseDeltaZ = cms.double(-1.0),
            MinPuppiWeight = cms.double(0.01),
            NumOfPUVtxsForCharged = cms.uint32(0),
            PUProxyValue = cms.InputTag("hltFixedGridRhoFastjetAll", "", "HLTX"),
            PtMaxCharged = cms.double(-1.0),
            PtMaxNeutrals = cms.double(200.0),
            PtMaxNeutralsStartSlope = cms.double(0.0),
            PtMaxPhotons = cms.double(-1.0),
            UseDeltaZCut = cms.bool(True),
            UseFromPVLooseTight = cms.bool(False),
            algos = cms.VPSet(
                cms.PSet(
                    EtaMaxExtrap = cms.double(2.0),
                    MedEtaSF = cms.vdouble(1.0, 1.0),
                    MinNeutralPt = cms.vdouble(0.5105, 0.821),
                    MinNeutralPtSlope = cms.vdouble(9.51e-06, 1.902e-05),
                    RMSEtaSF = cms.vdouble(1.0, 1.0),
                    etaMax = cms.vdouble(2.5, 3.5),
                    etaMin = cms.vdouble(0.0, 2.5),
                    ptMin = cms.vdouble(0.0, 0.0),
                    puppiAlgos = cms.VPSet(cms.PSet(
                        algoId = cms.int32(5),
                        applyLowPUCorr = cms.bool(True),
                        combOpt = cms.int32(0),
                        cone = cms.double(0.4),
                        rmsPtMin = cms.double(0.1),
                        rmsScaleFactor = cms.double(1.0),
                        useCharged = cms.bool(True)
                    ))
                ),
                cms.PSet(
                    EtaMaxExtrap = cms.double(2.0),
                    MedEtaSF = cms.vdouble(0.75),
                    MinNeutralPt = cms.vdouble(3.656),
                    MinNeutralPtSlope = cms.vdouble(5.072e-05),
                    RMSEtaSF = cms.vdouble(1.0),
                    etaMax = cms.vdouble(10.0),
                    etaMin = cms.vdouble(3.5),
                    ptMin = cms.vdouble(0.0),
                    puppiAlgos = cms.VPSet(cms.PSet(
                        algoId = cms.int32(5),
                        applyLowPUCorr = cms.bool(True),
                        combOpt = cms.int32(0),
                        cone = cms.double(0.4),
                        rmsPtMin = cms.double(0.5),
                        rmsScaleFactor = cms.double(1.0),
                        useCharged = cms.bool(False)
                    ))
                )
            ),
            applyCHS = cms.bool(True),
            candName = cms.InputTag(f"hltPfTICL{i}"),
            clonePackedCands = cms.bool(False),
            invertPuppi = cms.bool(False),
            puppiDiagnostics = cms.bool(False),
            puppiNoLep = cms.bool(False),
            useExistingWeights = cms.bool(False),
            useExp = cms.bool(False),
            usePUProxyValue = cms.bool(True),
            vertexName = cms.InputTag("hltOfflinePrimaryVertices", "", "HLTX"),
            vtxNdofCut = cms.int32(4),
            vtxZCut = cms.double(24)
        )
    )
    puppiModules.append(getattr(process, f"hltPFPuppi{i}"))

    # AK4 jets from PF candidates with puppi weights
    setattr(process, f"hltAK4PFPuppiJets{i}",
        cms.EDProducer("FastjetJetProducer",
            Active_Area_Repeats = cms.int32(1),
            GhostArea = cms.double(0.01),
            Ghost_EtaMax = cms.double(5.0),
            Rho_EtaMax = cms.double(4.4),
            applyWeight = cms.bool(True),
            doAreaDiskApprox = cms.bool(False),
            doAreaFastjet = cms.bool(True),
            doPUOffsetCorr = cms.bool(False),
            doPVCorrection = cms.bool(False),
            doRhoFastjet = cms.bool(False),
            inputEMin = cms.double(0.0),
            inputEtMin = cms.double(0.0),
            jetAlgorithm = cms.string('AntiKt'),
            jetPtMin = cms.double(5.0),
            jetType = cms.string('PFJet'),
            maxBadEcalCells = cms.uint32(9999999),
            maxBadHcalCells = cms.uint32(9999999),
            maxProblematicEcalCells = cms.uint32(9999999),
            maxProblematicHcalCells = cms.uint32(9999999),
            maxRecoveredEcalCells = cms.uint32(9999999),
            maxRecoveredHcalCells = cms.uint32(9999999),
            minSeed = cms.uint32(14327),
            rParam = cms.double(0.4),
            src = cms.InputTag(f"hltPfTICL{i}"),
            srcPVs = cms.InputTag(""),
            srcWeights = cms.InputTag(f"hltPFPuppi{i}"),
            useDeterministicSeed = cms.bool(True),
            voronoiRfact = cms.double(-0.9)
        )
    )
    ak4JetModules.append(getattr(process, f"hltAK4PFPuppiJets{i}"))

    # validation metrics
    setattr(process, f"jetValidationMetrics{i}",
        cms.EDAnalyzer("JetValidationMetrics",
            jets               = cms.InputTag(f"hltAK4PFPuppiJets{i}"),
            genjets            = cms.InputTag("ak4GenJetsNoNu"),
            recoJetPtThreshold = cms.double(30.0),
            matchGenPtThreshold= cms.double(20.0),
            RThreshold         = cms.double(0.4),
            mightGet           = cms.optional.untracked.vstring,
        )
    )
    metricsModules.append(getattr(process, f"jetValidationMetrics{i}"))

# output
outFile = getattr(options, "outputFile", "jet_validation_metrics.root")
process.TFileService = cms.Service("TFileService", fileName=cms.string(outFile))

process.optimizeTask = cms.Task(
    *cleaningModules,
    *candidateModules,
    *pfModules,
    *puppiModules,
    *ak4JetModules
)

process.metricsSeq = cms.Sequence()
for m in metricsModules:
    process.metricsSeq += m

process.TICLCleaningOptimize = cms.Path(process.metricsSeq)
process.TICLCleaningOptimize.associate(process.optimizeTask)

process.schedule = cms.Schedule(process.TICLCleaningOptimize)


# customizations
process = setCrossingFrameOn(process)
process = customiseLogErrorHarvesterUsingOutputCommands(process)
process = customiseEarlyDelete(process)
