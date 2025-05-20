#include "ApprovalTests.hpp"
#include "Parameters.h"
#include "catch2/catch_amalgamated.hpp"
#include <iomanip>
#include <iostream>
#include <string>

std::ostream &operator<<(std::ostream &os, const Parameters &params) {
  os << std::setprecision(18);
  os << "Current Parameters are: " << std::endl;

  os << "\n[Main]" << std::endl;
  os << "kB: " << params.kB << std::endl;
  os << "timeUnit: " << params.timeUnit << std::endl;
  os << "job: " << magic_enum::enum_name(params.job) << std::endl;
  os << "randomSeed: " << params.randomSeed << std::endl;
  os << "temperature: " << params.temperature << std::endl;
  os << "quiet: " << std::boolalpha << params.quiet << std::endl;
  os << "writeLog: " << std::boolalpha << params.writeLog << std::endl;
  os << "checkpoint: " << std::boolalpha << params.checkpoint << std::endl;
  os << "iniFilename: " << params.iniFilename << std::endl;
  os << "conFilename: " << params.conFilename << std::endl;
  os << "finiteDifference: " << params.finiteDifference << std::endl;
  os << "maxForceCalls: " << params.maxForceCalls << std::endl;
  os << "removeNetForce: " << std::boolalpha << params.removeNetForce
     << std::endl;
  os << "\n[Potential]" << std::endl;
  os << "potential: " << magic_enum::enum_name(params.potential) << std::endl;
  os << "MPIPollPeriod: " << params.MPIPollPeriod << std::endl;
  os << "MPIPotentialRank: " << params.MPIPotentialRank << std::endl;
  os << "LAMMPSLogging: " << std::boolalpha << params.LAMMPSLogging
     << std::endl;
  os << "LAMMPSThreads: " << params.LAMMPSThreads << std::endl;
  os << "EMTRasmussen: " << std::boolalpha << params.EMTRasmussen << std::endl;
  os << "LogPotential: " << std::boolalpha << params.LogPotential << std::endl;
  os << "extPotPath: " << params.extPotPath << std::endl;

  os << "\n[AMS]" << std::endl;
  os << "engine: " << params.engine << std::endl;
  os << "forcefield: " << params.forcefield << std::endl;
  os << "model: " << params.model << std::endl;
  os << "resources: " << params.resources << std::endl;
  os << "xc: " << params.xc << std::endl;
  os << "basis: " << params.basis << std::endl;

  os << "\n[AMS_ENV]" << std::endl;
  os << "amshome: " << params.amshome << std::endl;
  os << "scm_tmpdir: " << params.scm_tmpdir << std::endl;
  os << "scmlicense: " << params.scmlicense << std::endl;
  os << "scm_pythondir: " << params.scm_pythondir << std::endl;
  os << "amsbin: " << params.amsbin << std::endl;
  os << "amsresources: " << params.amsresources << std::endl;

  os << "\n[XTBPot]" << std::endl;
  os << "xtb_paramset: " << params.xtb_paramset << std::endl;
  os << "xtb_elec_temperature: " << params.xtb_elec_temperature << std::endl;
  os << "xtb_maxiter: " << params.xtb_maxiter << std::endl;
  os << "xtb_acc: " << params.xtb_acc << std::endl;

  os << "\n[Structure Comparison]" << std::endl;
  os << "distanceDifference: " << params.distanceDifference << std::endl;
  os << "neighborCutoff: " << params.neighborCutoff << std::endl;
  os << "checkRotation: " << std::boolalpha << params.checkRotation
     << std::endl;
  os << "indistinguishableAtoms: " << std::boolalpha
     << params.indistinguishableAtoms << std::endl;
  os << "energyDifference: " << params.energyDifference << std::endl;
  os << "removeTranslation: " << std::boolalpha << params.removeTranslation
     << std::endl;

  os << "\n[Process Search]" << std::endl;
  os << "processSearchMinimizeFirst: " << std::boolalpha
     << params.processSearchMinimizeFirst << std::endl;
  os << "processSearchMinimizationOffset: "
     << params.processSearchMinimizationOffset << std::endl;

  os << "\n[Saddle Search]" << std::endl;
  os << "saddleMaxJumpAttempts: " << params.saddleMaxJumpAttempts << std::endl;
  os << "saddleMaxIterations: " << params.saddleMaxIterations << std::endl;
  os << "saddleMethod: " << params.saddleMethod << std::endl;
  os << "saddleMinmodeMethod: " << params.saddleMinmodeMethod << std::endl;
  os << "saddleDisplaceType: " << params.saddleDisplaceType << std::endl;
  os << "saddleMaxEnergy: " << params.saddleMaxEnergy << std::endl;
  os << "saddleDisplaceMagnitude: " << params.saddleDisplaceMagnitude
     << std::endl;
  os << "saddleMaxSingleDisplace: " << params.saddleMaxSingleDisplace
     << std::endl;
  os << "saddleDisplaceRadius: " << params.saddleDisplaceRadius << std::endl;
  os << "saddleConvergedForce: " << params.saddleConvergedForce << std::endl;
  os << "saddlePerpForceRatio: " << params.saddlePerpForceRatio << std::endl;
  os << "saddleNonnegativeDisplacementAbort: " << std::boolalpha
     << params.saddleNonnegativeDisplacementAbort << std::endl;
  os << "saddleNonlocalCountAbort: " << params.saddleNonlocalCountAbort
     << std::endl;
  os << "saddleNonlocalDistanceAbort: " << params.saddleNonlocalDistanceAbort
     << std::endl;
  os << "saddleRemoveRotation: " << std::boolalpha
     << params.saddleRemoveRotation << std::endl;
  os << "saddleDynamicsTemperature: " << params.saddleDynamicsTemperature
     << std::endl;
  os << "saddleDynamicsStateCheckIntervalInput: "
     << params.saddleDynamicsStateCheckIntervalInput << std::endl;
  os << "saddleDynamicsStateCheckInterval: "
     << params.saddleDynamicsStateCheckInterval << std::endl;
  os << "saddleDynamicsRecordIntervalInput: "
     << params.saddleDynamicsRecordIntervalInput << std::endl;
  os << "saddleDynamicsRecordInterval: " << params.saddleDynamicsRecordInterval
     << std::endl;
  os << "saddleDynamicsLinearInterpolation: " << std::boolalpha
     << params.saddleDynamicsLinearInterpolation << std::endl;
  os << "saddleDynamicsMaxInitCurvature: "
     << params.saddleDynamicsMaxInitCurvature << std::endl;
  os << "saddleConfinePositive: " << std::boolalpha
     << params.saddleConfinePositive << std::endl;
  os << "saddleBowlBreakout: " << std::boolalpha << params.saddleBowlBreakout
     << std::endl;
  os << "saddleBowlActive: " << params.saddleBowlActive << std::endl;
  os << "saddleConfinePositiveMinForce: "
     << params.saddleConfinePositiveMinForce << std::endl;
  os << "saddleConfinePositiveScaleRatio: "
     << params.saddleConfinePositiveScaleRatio << std::endl;
  os << "saddleConfinePositiveBoost: " << params.saddleConfinePositiveBoost
     << std::endl;
  os << "saddleConfinePositiveMinActive: "
     << params.saddleConfinePositiveMinActive << std::endl;
  os << "saddleZeroModeAbortCurvature: " << params.saddleZeroModeAbortCurvature
     << std::endl;

  os << "\n[Optimizers]" << std::endl;
  os << "optMethod: " << magic_enum::enum_name(params.optMethod) << std::endl;
  os << "optConvergenceMetric: " << params.optConvergenceMetric << std::endl;
  os << "optConvergenceMetricLabel: " << params.optConvergenceMetricLabel
     << std::endl;
  os << "optMaxIterations: " << params.optMaxIterations << std::endl;
  os << "optMaxMove: " << params.optMaxMove << std::endl;
  os << "optConvergedForce: " << params.optConvergedForce << std::endl;
  os << "optTimeStepInput: " << params.optTimeStepInput << std::endl;
  os << "optTimeStep: " << params.optTimeStep << std::endl;
  os << "optMaxTimeStepInput: " << params.optMaxTimeStepInput << std::endl;
  os << "optMaxTimeStep: " << params.optMaxTimeStep << std::endl;
  os << "optLBFGSMemory: " << params.optLBFGSMemory << std::endl;
  os << "optLBFGSInverseCurvature: " << params.optLBFGSInverseCurvature
     << std::endl;
  os << "optLBFGSMaxInverseCurvature: " << params.optLBFGSMaxInverseCurvature
     << std::endl;
  os << "optLBFGSAutoScale: " << std::boolalpha << params.optLBFGSAutoScale
     << std::endl;
  os << "optLBFGSAngleReset: " << std::boolalpha << params.optLBFGSAngleReset
     << std::endl;
  os << "optLBFGSDistanceReset: " << std::boolalpha
     << params.optLBFGSDistanceReset << std::endl;
  os << "optQMSteepestDecent: " << std::boolalpha << params.optQMSteepestDecent
     << std::endl;
  os << "optCGNoOvershooting: " << std::boolalpha << params.optCGNoOvershooting
     << std::endl;
  os << "optCGKnockOutMaxMove: " << std::boolalpha
     << params.optCGKnockOutMaxMove << std::endl;
  os << "optCGLineSearch: " << std::boolalpha << params.optCGLineSearch
     << std::endl;
  os << "optCGLineConverged: " << params.optCGLineConverged << std::endl;
  os << "optCGLineSearchMaxIter: " << params.optCGLineSearchMaxIter
     << std::endl;
  os << "optCGMaxIterBeforeReset: " << params.optCGMaxIterBeforeReset
     << std::endl;
  os << "optSDAlpha: " << params.optSDAlpha << std::endl;
  os << "optSDTwoPoint: " << std::boolalpha << params.optSDTwoPoint
     << std::endl;

  os << "\n[Refine]" << std::endl;
  os << "refineOptMethod: " << magic_enum::enum_name(params.refineOptMethod)
     << std::endl;
  os << "refineThreshold: " << params.refineThreshold << std::endl;

  os << "\n[Dimer]" << std::endl;
  os << "dimerRotationAngle: " << params.dimerRotationAngle << std::endl;
  os << "dimerImproved: " << std::boolalpha << params.dimerImproved
     << std::endl;
  os << "dimerConvergedAngle: " << params.dimerConvergedAngle << std::endl;
  os << "dimerMaxIterations: " << params.dimerMaxIterations << std::endl;
  os << "dimerOptMethod: " << params.dimerOptMethod << std::endl;
  os << "dimerRotationsMax: " << params.dimerRotationsMax << std::endl;
  os << "dimerRotationsMin: " << params.dimerRotationsMin << std::endl;
  os << "dimerTorqueMax: " << params.dimerTorqueMax << std::endl;
  os << "dimerTorqueMin: " << params.dimerTorqueMin << std::endl;
  os << "dimerRemoveRotation: " << std::boolalpha << params.dimerRemoveRotation
     << std::endl;

  os << "\n[GPR Dimer]" << std::endl;
  os << "gprDimerRotationAngle: " << params.gprDimerRotationAngle << std::endl;
  os << "gprDimerConvergedAngle: " << params.gprDimerConvergedAngle
     << std::endl;
  os << "gprDimerRelaxConvAngle: " << params.gprDimerRelaxConvAngle
     << std::endl;
  os << "gprDimerInitRotationsMax: " << params.gprDimerInitRotationsMax
     << std::endl;
  os << "gprDimerRelaxRotationsMax: " << params.gprDimerRelaxRotationsMax
     << std::endl;
  os << "gprDimerDivisorTdimerGP: " << params.gprDimerDivisorTdimerGP
     << std::endl;
  os << "gprDimerMaxOuterIterations: " << params.gprDimerMaxOuterIterations
     << std::endl;
  os << "gprDimerMaxInnerIterations: " << params.gprDimerMaxInnerIterations
     << std::endl;
  os << "gprDimerMidpointMaxDisp: " << params.gprDimerMidpointMaxDisp
     << std::endl;
  os << "gprDimerRotOptMethod: " << params.gprDimerRotOptMethod << std::endl;
  os << "gprDimerTransOptMethod: " << params.gprDimerTransOptMethod
     << std::endl;
  os << "gprActiveRadius: " << params.gprActiveRadius << std::endl;
  os << "gprDimerSep: " << params.gprDimerSep << std::endl;
  os << "gprDimerConvStep: " << params.gprDimerConvStep << std::endl;
  os << "gprDimerMaxStep: " << params.gprDimerMaxStep << std::endl;
  os << "gprForceThreshold: " << params.saddleConvergedForce << std::endl;
  os << "gprDimerRatioAtLimit: " << params.gprDimerRatioAtLimit << std::endl;
  os << "gprDimerInitRotGP: " << std::boolalpha << params.gprDimerInitRotGP
     << std::endl;
  os << "gprDimerInitTransGP: " << std::boolalpha << params.gprDimerInitTransGP
     << std::endl;
  os << "gprDimerManyIterations: " << std::boolalpha
     << params.gprDimerManyIterations << std::endl;
  os << "gprDimerHyperOptMethod: " << params.gprDimerHyperOptMethod
     << std::endl;
  os << "gprDimerSigma2: " << params.gprDimerSigma2 << std::endl;
  os << "gprDimerJitterSigma2: " << params.gprDimerJitterSigma2 << std::endl;
  os << "gprDimerNoiseSigma2: " << params.gprDimerNoiseSigma2 << std::endl;
  os << "gprDimerPriorMu: " << params.gprDimerPriorMu << std::endl;
  os << "gprDimerPriorSigma2: " << params.gprDimerPriorSigma2 << std::endl;
  os << "gprDimerPriorNu: " << params.gprDimerPriorNu << std::endl;
  os << "gprOptCheckDerivatives: " << std::boolalpha
     << params.gprOptCheckDerivatives << std::endl;
  os << "gprOptMaxIterations: " << params.gprOptMaxIterations << std::endl;
  os << "gprOptTolFunc: " << params.gprOptTolFunc << std::endl;
  os << "gprOptTolSol: " << params.gprOptTolSol << std::endl;
  os << "gprOptLambdaLimit: " << params.gprOptLambdaLimit << std::endl;
  os << "gprOptLambdaInit: " << params.gprOptLambdaInit << std::endl;
  os << "gprUsePrune: " << std::boolalpha << params.gprUsePrune << std::endl;
  os << "gprPruneBegin: " << params.gprPruneBegin << std::endl;
  os << "gprPruneNVals: " << params.gprPruneNVals << std::endl;
  os << "gprPruneThreshold: " << params.gprPruneThreshold << std::endl;
  os << "gprReportLevel: " << params.gprReportLevel << std::endl;
  os << "gprDebugLevel: " << params.gprDebugLevel << std::endl;
  os << "gprDebugOutDir: " << params.gprDebugOutDir << std::endl;
  os << "gprDebugPosFile: " << params.gprDebugPosFile << std::endl;
  os << "gprDebugEnergyFile: " << params.gprDebugEnergyFile << std::endl;
  os << "gprDebugGradFile: " << params.gprDebugGradFile << std::endl;
  os << "gprDebugOutExt: " << params.gprDebugOutExt << std::endl;
  os << "gprDebugOffsetMidPoint: " << params.gprDebugOffsetMidPoint
     << std::endl;
  os << "gprDebugDy: " << params.gprDebugDy << std::endl;
  os << "gprDebugDz: " << params.gprDebugDz << std::endl;

  os << "\n[Surrogate]" << std::endl;
  os << "use_surrogate: " << params.use_surrogate << std::endl;
  os << "sub_job: " << magic_enum::enum_name(params.sub_job) << std::endl;
  os << "gp_uncertainity: " << params.gp_uncertainity << std::endl;
  os << "gp_linear_path_always: " << std::boolalpha
     << params.gp_linear_path_always << std::endl;
  os << "surrogatePotential: "
     << magic_enum::enum_name(params.surrogatePotential) << std::endl;

  os << "\n[CatLearn]" << std::endl;
  os << "catl_path: " << params.catl_path << std::endl;
  os << "catl_model: " << params.catl_model << std::endl;
  os << "catl_prior: " << params.catl_prior << std::endl;
  os << "catl_use_deriv: " << std::boolalpha << params.catl_use_deriv
     << std::endl;
  os << "catl_use_fingerprint: " << std::boolalpha
     << params.catl_use_fingerprint << std::endl;
  os << "catl_parallel: " << std::boolalpha << params.catl_parallel
     << std::endl;

  os << "\n[ASE ORCA]" << std::endl;
  os << "orca_path: " << params.orca_path << std::endl;
  os << "orca_nproc: " << params.orca_nproc << std::endl;
  os << "orca_sline: " << params.orca_sline << std::endl;

  os << "\n[Lanczos]" << std::endl;
  os << "lanczosTolerance: " << params.lanczosTolerance << std::endl;
  os << "lanczosMaxIterations: " << params.lanczosMaxIterations << std::endl;
  os << "lanczosQuitEarly: " << std::boolalpha << params.lanczosQuitEarly
     << std::endl;

  os << "\n[Prefactor]" << std::endl;
  os << "prefactorDefaultValue: " << params.prefactorDefaultValue << std::endl;
  os << "prefactorMaxValue: " << params.prefactorMaxValue << std::endl;
  os << "prefactorMinValue: " << params.prefactorMinValue << std::endl;
  os << "prefactorWithinRadius: " << params.prefactorWithinRadius << std::endl;
  os << "prefactorMinDisplacement: " << params.prefactorMinDisplacement
     << std::endl;
  os << "prefactorRate: " << params.prefactorRate << std::endl;
  os << "prefactorConfiguration: " << params.prefactorConfiguration
     << std::endl;
  os << "prefactorAllFreeAtoms: " << std::boolalpha
     << params.prefactorAllFreeAtoms << std::endl;
  os << "prefactorFilterScheme: " << params.prefactorFilterScheme << std::endl;
  os << "prefactorFilterFraction: " << params.prefactorFilterFraction
     << std::endl;

  os << "\n[Hessian]" << std::endl;
  os << "hessianAtomList: " << params.hessianAtomList << std::endl;
  os << "hessianZeroFreqValue: " << params.hessianZeroFreqValue << std::endl;

  os << "\n[Nudged Elastic Band]" << std::endl;
  os << "nebImages: " << params.nebImages << std::endl;
  os << "nebMaxIterations: " << params.nebMaxIterations << std::endl;
  os << "nebSpring: " << params.nebSpring << std::endl;
  os << "nebClimbingImageMethod: " << std::boolalpha
     << params.nebClimbingImageMethod << std::endl;
  os << "nebClimbingImageConvergedOnly: " << std::boolalpha
     << params.nebClimbingImageConvergedOnly << std::endl;
  os << "nebOldTangent: " << std::boolalpha << params.nebOldTangent
     << std::endl;
  os << "nebDoublyNudged: " << std::boolalpha << params.nebDoublyNudged
     << std::endl;
  os << "nebDoublyNudgedSwitching: " << std::boolalpha
     << params.nebDoublyNudgedSwitching << std::endl;
  os << "nebOptMethod: " << params.nebOptMethod << std::endl;
  os << "nebElasticBand: " << std::boolalpha << params.nebElasticBand
     << std::endl;
  os << "nebConvergedForce: " << params.nebConvergedForce << std::endl;
  os << "nebKSPMin: " << params.nebKSPMin << std::endl;
  os << "nebKSPMax: " << params.nebKSPMax << std::endl;
  os << "nebEnergyWeighted: " << std::boolalpha << params.nebEnergyWeighted
     << std::endl;

  os << "\n[Dynamics]" << std::endl;
  os << "mdTimeStepInput: " << params.mdTimeStepInput << std::endl;
  os << "mdTimeStep: " << params.mdTimeStep << std::endl;
  os << "mdTimeInput: " << params.mdTimeInput << std::endl;
  os << "mdTime: " << params.mdTime << std::endl;
  os << "mdSteps: " << params.mdSteps << std::endl;

  os << "\n[Parallel Replica]" << std::endl;
  os << "parrepRefineTransition: " << std::boolalpha
     << params.parrepRefineTransition << std::endl;
  os << "parrepAutoStop: " << std::boolalpha << params.parrepAutoStop
     << std::endl;
  os << "parrepDephaseLoopStop: " << std::boolalpha
     << params.parrepDephaseLoopStop << std::endl;
  os << "parrepDephaseTimeInput: " << params.parrepDephaseTimeInput
     << std::endl;
  os << "parrepDephaseTime: " << params.parrepDephaseTime << std::endl;
  os << "parrepDephaseLoopMax: " << params.parrepDephaseLoopMax << std::endl;
  os << "parrepStateCheckIntervalInput: "
     << params.parrepStateCheckIntervalInput << std::endl;
  os << "parrepStateCheckInterval: " << params.parrepStateCheckInterval
     << std::endl;
  os << "parrepRecordIntervalInput: " << params.parrepRecordIntervalInput
     << std::endl;
  os << "parrepRecordInterval: " << params.parrepRecordInterval << std::endl;
  os << "parrepCorrTimeInput: " << params.parrepCorrTimeInput << std::endl;
  os << "parrepCorrTime: " << params.parrepCorrTime << std::endl;

  os << "\n[TAD]" << std::endl;
  os << "tadLowT: " << params.tadLowT << std::endl;
  os << "tadMinPrefactor: " << params.tadMinPrefactor << std::endl;
  os << "tadConfidence: " << params.tadConfidence << std::endl;

  os << "\n[Thermostat]" << std::endl;
  os << "thermostat: " << params.thermostat << std::endl;
  os << "thermoAndersenAlpha: " << params.thermoAndersenAlpha << std::endl;
  os << "thermoAndersenTcolInput: " << params.thermoAndersenTcolInput
     << std::endl;
  os << "thermoAndersenTcol: " << params.thermoAndersenTcol << std::endl;
  os << "thermoNoseMass: " << params.thermoNoseMass << std::endl;
  os << "thermoLangevinFrictionInput: " << params.thermoLangevinFrictionInput
     << std::endl;
  os << "thermoLangevinFriction: " << params.thermoLangevinFriction
     << std::endl;

  os << "\n[Parallel Replica]" << std::endl;
  os << "repexcTemperatureDistribution: "
     << params.repexcTemperatureDistribution << std::endl;
  os << "repexcReplicas: " << params.repexcReplicas << std::endl;
  os << "repexcExchangeTrials: " << params.repexcExchangeTrials << std::endl;
  os << "repexcSamplingTimeInput: " << params.repexcSamplingTimeInput
     << std::endl;
  os << "repexcSamplingTime: " << params.repexcSamplingTime << std::endl;
  os << "repexcTemperatureHigh: " << params.repexcTemperatureHigh << std::endl;
  os << "repexcTemperatureLow: " << params.repexcTemperatureLow << std::endl;
  os << "repexcExchangePeriodInput: " << params.repexcExchangePeriodInput
     << std::endl;
  os << "repexcExchangePeriod: " << params.repexcExchangePeriod << std::endl;

  os << "\n[Hyperdynamics]" << std::endl;
  os << "biasPotential: " << params.biasPotential << std::endl;
  os << "bondBoostBALS: " << params.bondBoostBALS << std::endl;
  os << "bondBoostRMDTimeInput: " << params.bondBoostRMDTimeInput << std::endl;
  os << "bondBoostRMDTime: " << params.bondBoostRMDTime << std::endl;
  os << "bondBoostDVMAX: " << params.bondBoostDVMAX << std::endl;
  os << "bondBoostQRR: " << params.bondBoostQRR << std::endl;
  os << "bondBoostPRR: " << params.bondBoostPRR << std::endl;
  os << "bondBoostQcut: " << params.bondBoostQcut << std::endl;
  os << "basinHoppingDisplacement: " << params.basinHoppingDisplacement
     << std::endl;
  os << "basinHoppingInitialRandomStructureProbability: "
     << params.basinHoppingInitialRandomStructureProbability << std::endl;
  os << "basinHoppingPushApartDistance: "
     << params.basinHoppingPushApartDistance << std::endl;
  os << "basinHoppingSteps: " << params.basinHoppingSteps << std::endl;
  os << "basinHoppingQuenchingSteps: " << params.basinHoppingQuenchingSteps
     << std::endl;
  os << "basinHoppingSignificantStructure: " << std::boolalpha
     << params.basinHoppingSignificantStructure << std::endl;
  os << "basinHoppingSingleAtomDisplace: " << std::boolalpha
     << params.basinHoppingSingleAtomDisplace << std::endl;
  os << "basinHoppingDisplacementAlgorithm: "
     << params.basinHoppingDisplacementAlgorithm << std::endl;
  os << "basinHoppingDisplacementDistribution: "
     << params.basinHoppingDisplacementDistribution << std::endl;
  os << "basinHoppingSwapProbability: " << params.basinHoppingSwapProbability
     << std::endl;
  os << "basinHoppingJumpMax: " << params.basinHoppingJumpMax << std::endl;
  os << "basinHoppingJumpSteps: " << params.basinHoppingJumpSteps << std::endl;
  os << "basinHoppingAdjustDisplacement: " << std::boolalpha
     << params.basinHoppingAdjustDisplacement << std::endl;
  os << "basinHoppingAdjustPeriod: " << params.basinHoppingAdjustPeriod
     << std::endl;
  os << "basinHoppingAdjustFraction: " << params.basinHoppingAdjustFraction
     << std::endl;
  os << "basinHoppingTargetRatio: " << params.basinHoppingTargetRatio
     << std::endl;
  os << "basinHoppingWriteUnique: " << std::boolalpha
     << params.basinHoppingWriteUnique << std::endl;
  os << "basinHoppingStopEnergy (-DBL_MAX): " << params.basinHoppingStopEnergy
     << std::endl;

  os << "\n[Global Optimization]" << std::endl;
  os << "globalOptimizationMoveMethod: " << params.globalOptimizationMoveMethod
     << std::endl;
  os << "globalOptimizationDecisionMethod: "
     << params.globalOptimizationDecisionMethod << std::endl;
  os << "globalOptimizationSteps: " << params.globalOptimizationSteps
     << std::endl;
  os << "globalOptimizationBeta: " << params.globalOptimizationBeta
     << std::endl;
  os << "globalOptimizationAlpha: " << params.globalOptimizationAlpha
     << std::endl;
  os << "globalOptimizationMdmin: " << params.globalOptimizationMdmin
     << std::endl;
  os << "globalOptimizationTargetEnergy: "
     << params.globalOptimizationTargetEnergy << std::endl;

  os << "\n[Monte Carlo]" << std::endl;
  os << "monteCarloStepSize: " << params.monteCarloStepSize << std::endl;
  os << "monteCarloSteps: " << params.monteCarloSteps << std::endl;

  os << "\n[BGSD]" << std::endl;
  os << "alpha: " << params.alpha << std::endl;
  os << "beta: " << params.beta << std::endl;
  os << "gradientfinitedifference: " << params.gradientfinitedifference
     << std::endl;
  os << "Hforceconvergence: " << params.Hforceconvergence << std::endl;
  os << "grad2energyconvergence: " << params.grad2energyconvergence
     << std::endl;
  os << "grad2forceconvergence: " << params.grad2forceconvergence << std::endl;

  os << "\n[Debug]" << std::endl;
  os << "writeMovies: " << std::boolalpha << params.writeMovies << std::endl;
  os << "writeMoviesInterval: " << params.writeMoviesInterval << std::endl;
  return os;
}

std::vector<Parameters> getParameters() {
  // Return test data for Parameters
  return {Parameters()};
}

TEST_CASE("VerifyParameters") {
  ApprovalTests::Approvals::verifyAll("parameters", getParameters());
}
