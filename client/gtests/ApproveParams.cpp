/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
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
  os << "job: " << magic_enum::enum_name(params.main.job) << std::endl;
  os << "randomSeed: " << params.main.randomSeed << std::endl;
  os << "temperature: " << params.main.temperature << std::endl;
  os << "quiet: " << std::boolalpha << params.main.quiet << std::endl;
  os << "writeLog: " << std::boolalpha << params.main.writeLog << std::endl;
  os << "checkpoint: " << std::boolalpha << params.main.checkpoint << std::endl;
  os << "inpFilename: " << params.main.inpFilename << std::endl;
  os << "conFilename: " << params.main.conFilename << std::endl;
  os << "finiteDifference: " << params.main.finiteDifference << std::endl;
  os << "maxForceCalls: " << params.main.maxForceCalls << std::endl;
  os << "removeNetForce: " << std::boolalpha << params.main.removeNetForce
     << std::endl;

  os << "\n[Potential]" << std::endl;
  os << "potential: " << magic_enum::enum_name(params.pot.potential)
     << std::endl;
  os << "MPIPollPeriod: " << params.pot.MPIPollPeriod << std::endl;
  // TODO(rg): Should this be here? PotentialRank
  os << "MPIPotentialRank: " << params.MPIPotentialRank << std::endl;
  os << "LAMMPSLogging: " << std::boolalpha << params.pot.LAMMPSLogging
     << std::endl;
  os << "LAMMPSThreads: " << params.pot.LAMMPSThreads << std::endl;
  os << "EMTRasmussen: " << std::boolalpha << params.pot.EMTRasmussen
     << std::endl;
  os << "LogPotential: " << std::boolalpha << params.pot.LogPotential
     << std::endl;
  os << "extPotPath: " << params.pot.extPotPath << std::endl;

  os << "\n[AMS]" << std::endl;
  os << "engine: " << params.ams.engine << std::endl;
  os << "forcefield: " << params.ams.forcefield << std::endl;
  os << "model: " << params.ams.model << std::endl;
  os << "resources: " << params.ams.resources << std::endl;
  os << "xc: " << params.ams.xc << std::endl;
  os << "basis: " << params.ams.basis << std::endl;

  os << "\n[AMS_ENV]" << std::endl;
  os << "amshome: " << params.ams.amsenv.amshome << std::endl;
  os << "scm_tmpdir: " << params.ams.amsenv.scm_tmpdir << std::endl;
  os << "scmlicense: " << params.ams.amsenv.scmlicense << std::endl;
  os << "scm_pythondir: " << params.ams.amsenv.scm_pythondir << std::endl;
  os << "amsbin: " << params.ams.amsenv.amsbin << std::endl;
  os << "amsresources: " << params.ams.amsenv.amsresources << std::endl;

  os << "\n[XTBPot]" << std::endl;
  os << "xtb_paramset: " << params.pot.xtbp.paramset << std::endl;
  os << "xtb_elec_temperature: " << params.pot.xtbp.elec_temperature << std::endl;
  os << "xtb_maxiter: " << params.pot.xtbp.maxiter << std::endl;
  os << "xtb_acc: " << params.pot.xtbp.acc << std::endl;

  os << "\n[Structure Comparison]" << std::endl;
  os << "distanceDifference: " << params.structcomp.distanceDifference
     << std::endl;
  os << "neighborCutoff: " << params.structcomp.neighborCutoff << std::endl;
  os << "checkRotation: " << std::boolalpha << params.structcomp.checkRotation
     << std::endl;
  os << "indistinguishableAtoms: " << std::boolalpha
     << params.structcomp.indistinguishableAtoms << std::endl;
  os << "energyDifference: " << params.structcomp.energyDifference << std::endl;
  os << "removeTranslation: " << std::boolalpha
     << params.structcomp.removeTranslation << std::endl;

  os << "\n[Process Search]" << std::endl;
  os << "processSearchMinimizeFirst: " << std::boolalpha
     << params.procsearch.minimizeFirst << std::endl;
  os << "processSearchMinimizationOffset: "
     << params.procsearch.minimizationOffset << std::endl;

  os << "\n[Saddle Search]" << std::endl;
  os << "saddleMaxJumpAttempts: " << params.saddle.maxJumpAttempts << std::endl;
  os << "saddleMaxIterations: " << params.saddle.maxIterations << std::endl;
  os << "saddleMethod: " << params.saddle.method << std::endl;
  os << "saddleMinmodeMethod: " << params.saddle.minmodeMethod << std::endl;
  os << "saddleDisplaceType: " << params.saddle.displaceType << std::endl;
  os << "saddleMaxEnergy: " << params.saddle.maxEnergy << std::endl;
  os << "saddleDisplaceMagnitude: " << params.saddle.displaceMagnitude
     << std::endl;
  os << "saddleMaxSingleDisplace: " << params.saddle.maxSingleDisplace
     << std::endl;
  os << "saddleDisplaceRadius: " << params.saddle.displaceRadius << std::endl;
  os << "saddleConvergedForce: " << params.saddle.convergedForce << std::endl;
  os << "saddlePerpForceRatio: " << params.saddle.perpForceRatio << std::endl;
  os << "saddleNonnegativeDisplacementAbort: " << std::boolalpha
     << params.saddle.nonnegativeDisplacementAbort << std::endl;
  os << "saddleNonlocalCountAbort: " << params.saddle.nonlocalCountAbort
     << std::endl;
  os << "saddleNonlocalDistanceAbort: " << params.saddle.nonlocalDistanceAbort
     << std::endl;
  os << "saddleRemoveRotation: " << std::boolalpha
     << params.saddle.removeRotation << std::endl;
  os << "saddleDynamicsTemperature: " << params.saddle.dynamicsTemperature
     << std::endl;
  os << "saddleDynamicsStateCheckIntervalInput: "
     << params.saddle.dynamicsStateCheckIntervalInput << std::endl;
  os << "saddleDynamicsStateCheckInterval: "
     << params.saddle.dynamicsStateCheckInterval << std::endl;
  os << "saddleDynamicsRecordIntervalInput: "
     << params.saddle.dynamicsRecordIntervalInput << std::endl;
  os << "saddleDynamicsRecordInterval: " << params.saddle.dynamicsRecordInterval
     << std::endl;
  os << "saddleDynamicsLinearInterpolation: " << std::boolalpha
     << params.saddle.dynamicsLinearInterpolation << std::endl;
  os << "saddleDynamicsMaxInitCurvature: "
     << params.saddle.dynamicsMaxInitCurvature << std::endl;
  os << "saddleConfinePositive: " << std::boolalpha
     << params.saddle.confinePositive << std::endl;
  os << "saddleBowlBreakout: " << std::boolalpha << params.saddle.bowlBreakout
     << std::endl;
  os << "saddleBowlActive: " << params.saddle.bowlActive << std::endl;
  os << "saddleConfinePositiveMinForce: "
     << params.saddle.confinePositiveMinForce << std::endl;
  os << "saddleConfinePositiveScaleRatio: "
     << params.saddle.confinePositiveScaleRatio << std::endl;
  os << "saddleConfinePositiveBoost: " << params.saddle.confinePositiveBoost
     << std::endl;
  os << "saddleConfinePositiveMinActive: "
     << params.saddle.confinePositiveMinActive << std::endl;
  os << "saddleZeroModeAbortCurvature: " << params.saddle.zeroModeAbortCurvature
     << std::endl;

  os << "\n[Optimizers]" << std::endl;
  os << "optMethod: " << magic_enum::enum_name(params.optim.method)
     << std::endl;
  os << "optConvergenceMetric: " << params.optim.convergenceMetric << std::endl;
  os << "optConvergenceMetricLabel: " << params.optim.convergenceMetricLabel
     << std::endl;
  os << "optMaxIterations: " << params.optim.maxIterations << std::endl;
  os << "optMaxMove: " << params.optim.maxMove << std::endl;
  os << "optConvergedForce: " << params.optim.convergedForce << std::endl;
  os << "optTimeStepInput: " << params.optim.timeStepInput << std::endl;
  os << "optTimeStep: " << params.optim.timeStep << std::endl;
  os << "optMaxTimeStepInput: " << params.optim.maxTimeStepInput << std::endl;
  os << "optMaxTimeStep: " << params.optim.maxTimeStep << std::endl;
  os << "optLBFGSMemory: " << params.optim.LBFGSMemory << std::endl;
  os << "optLBFGSInverseCurvature: " << params.optim.LBFGSInverseCurvature
     << std::endl;
  os << "optLBFGSMaxInverseCurvature: " << params.optim.LBFGSMaxInverseCurvature
     << std::endl;
  os << "optLBFGSAutoScale: " << std::boolalpha << params.optim.LBFGSAutoScale
     << std::endl;
  os << "optLBFGSAngleReset: " << std::boolalpha << params.optim.LBFGSAngleReset
     << std::endl;
  os << "optLBFGSDistanceReset: " << std::boolalpha
     << params.optim.LBFGSDistanceReset << std::endl;
  os << "optQMSteepestDecent: " << std::boolalpha
     << params.optim.QMSteepestDecent << std::endl;
  os << "optCGNoOvershooting: " << std::boolalpha
     << params.optim.CGNoOvershooting << std::endl;
  os << "optCGKnockOutMaxMove: " << std::boolalpha
     << params.optim.CGKnockOutMaxMove << std::endl;
  os << "optCGLineSearch: " << std::boolalpha << params.optim.CGLineSearch
     << std::endl;
  os << "optCGLineConverged: " << params.optim.CGLineConverged << std::endl;
  os << "optCGLineSearchMaxIter: " << params.optim.CGLineSearchMaxIter
     << std::endl;
  os << "optCGMaxIterBeforeReset: " << params.optim.CGMaxIterBeforeReset
     << std::endl;
  os << "optSDAlpha: " << params.optim.SDAlpha << std::endl;
  os << "optSDTwoPoint: " << std::boolalpha << params.optim.SDTwoPoint
     << std::endl;

  os << "\n[Refine]" << std::endl;
  os << "refineOptMethod: "
     << magic_enum::enum_name(params.optim.refineOptMethod) << std::endl;
  os << "refineThreshold: " << params.optim.refineThreshold << std::endl;

  os << "\n[Dimer]" << std::endl;
  os << "dimerRotationAngle: " << params.dimer.rotationAngle << std::endl;
  os << "dimerImproved: " << std::boolalpha << params.dimer.improved
     << std::endl;
  os << "dimerConvergedAngle: " << params.dimer.convergedAngle << std::endl;
  os << "dimerMaxIterations: " << params.dimer.maxIterations << std::endl;
  os << "dimerOptMethod: " << params.dimer.optMethod << std::endl;
  os << "dimerRotationsMax: " << params.dimer.rotationsMax << std::endl;
  os << "dimerRotationsMin: " << params.dimer.rotationsMin << std::endl;
  os << "dimerTorqueMax: " << params.dimer.torqueMax << std::endl;
  os << "dimerTorqueMin: " << params.dimer.torqueMin << std::endl;
  os << "dimerRemoveRotation: " << std::boolalpha << params.dimer.removeRotation
     << std::endl;

  os << "\n[GPR Dimer]" << std::endl;
  os << "gprDimerRotationAngle: " << params.gprd.rotationAngle << std::endl;
  os << "gprDimerConvergedAngle: " << params.gprd.convergedAngle << std::endl;
  os << "gprDimerRelaxConvAngle: " << params.gprd.relaxConvAngle << std::endl;
  os << "gprDimerInitRotationsMax: " << params.gprd.initRotationsMax
     << std::endl;
  os << "gprDimerRelaxRotationsMax: " << params.gprd.relaxRotationsMax
     << std::endl;
  os << "gprDimerDivisorTdimerGP: " << params.gprd.divisorTdimerGP << std::endl;
  os << "gprDimerMaxOuterIterations: " << params.gprd.maxOuterIterations
     << std::endl;
  os << "gprDimerMaxInnerIterations: " << params.gprd.maxInnerIterations
     << std::endl;
  os << "gprDimerMidpointMaxDisp: " << params.gprd.midpointMaxDisp << std::endl;
  os << "gprDimerRotOptMethod: " << params.gprd.rotOptMethod << std::endl;
  os << "gprDimerTransOptMethod: " << params.gprd.transOptMethod << std::endl;
  os << "gprActiveRadius: " << params.gprd.activeRadius << std::endl;
  os << "gprDimerSep: " << params.gprd.dimerSep << std::endl;
  os << "gprDimerConvStep: " << params.gprd.convStep << std::endl;
  os << "gprDimerMaxStep: " << params.gprd.maxStep << std::endl;
  os << "gprForceThreshold: " << params.gprd.forceThreshold << std::endl;
  os << "gprDimerRatioAtLimit: " << params.gprd.ratioAtLimit << std::endl;
  os << "gprDimerInitRotGP: " << std::boolalpha << params.gprd.initRotGP
     << std::endl;
  os << "gprDimerInitTransGP: " << std::boolalpha << params.gprd.initTransGP
     << std::endl;
  os << "gprDimerManyIterations: " << std::boolalpha
     << params.gprd.manyIterations << std::endl;
  os << "gprDimerHyperOptMethod: " << params.gprd.hyperOptMethod << std::endl;
  os << "gprDimerSigma2: " << params.gprd.sigma2 << std::endl;
  os << "gprDimerJitterSigma2: " << params.gprd.jitterSigma2 << std::endl;
  os << "gprDimerNoiseSigma2: " << params.gprd.noiseSigma2 << std::endl;
  os << "gprDimerPriorMu: " << params.gprd.priorMu << std::endl;
  os << "gprDimerPriorSigma2: " << params.gprd.priorSigma2 << std::endl;
  os << "gprDimerPriorNu: " << params.gprd.priorNu << std::endl;
  os << "gprOptCheckDerivatives: " << std::boolalpha
     << params.gprd.optCheckDerivatives << std::endl;
  os << "gprOptMaxIterations: " << params.gprd.optMaxIterations << std::endl;
  os << "gprOptTolFunc: " << params.gprd.optTolFunc << std::endl;
  os << "gprOptTolSol: " << params.gprd.optTolSol << std::endl;
  os << "gprOptLambdaLimit: " << params.gprd.optLambdaLimit << std::endl;
  os << "gprOptLambdaInit: " << params.gprd.optLambdaInit << std::endl;
  os << "gprUsePrune: " << std::boolalpha << params.gprd.usePrune << std::endl;
  os << "gprPruneBegin: " << params.gprd.pruneBegin << std::endl;
  os << "gprPruneNVals: " << params.gprd.pruneNVals << std::endl;
  os << "gprPruneThreshold: " << params.gprd.pruneThreshold << std::endl;
  os << "gprReportLevel: " << params.gprd.reportLevel << std::endl;
  os << "gprDebugLevel: " << params.gprd.debugLevel << std::endl;
  os << "gprDebugOutDir: " << params.gprd.debugOutDir << std::endl;
  os << "gprDebugPosFile: " << params.gprd.debugPosFile << std::endl;
  os << "gprDebugEnergyFile: " << params.gprd.debugEnergyFile << std::endl;
  os << "gprDebugGradFile: " << params.gprd.debugGradFile << std::endl;
  os << "gprDebugOutExt: " << params.gprd.debugOutExt << std::endl;
  os << "gprDebugOffsetMidPoint: " << params.gprd.debugOffsetMidPoint
     << std::endl;
  os << "gprDebugDy: " << params.gprd.debugDy << std::endl;
  os << "gprDebugDz: " << params.gprd.debugDz << std::endl;

  os << "\n[Surrogate]" << std::endl;
  os << "use_surrogate: " << params.surrogate.use << std::endl;
  os << "sub_job: " << magic_enum::enum_name(params.surrogate.sub_job)
     << std::endl;
  os << "gp_uncertainity: " << params.surrogate.gp_uncertainity << std::endl;
  os << "gp_linear_path_always: " << std::boolalpha
     << params.surrogate.gp_linear_path_always << std::endl;
  os << "surrogatePotential: "
     << magic_enum::enum_name(params.surrogate.potential) << std::endl;

  os << "\n[CatLearn]" << std::endl;
  os << "catl_path: " << params.catl.path << std::endl;
  os << "catl_model: " << params.catl.model << std::endl;
  os << "catl_prior: " << params.catl.prior << std::endl;
  os << "catl_use_deriv: " << std::boolalpha << params.catl.use_deriv
     << std::endl;
  os << "catl_use_fingerprint: " << std::boolalpha
     << params.catl.use_fingerprint << std::endl;
  os << "catl_parallel: " << std::boolalpha << params.catl.parallel
     << std::endl;

  os << "\n[ASE ORCA]" << std::endl;
  os << "orca_path: " << params.aseorca.orca_path << std::endl;
  os << "orca_nproc: " << params.aseorca.orca_nproc << std::endl;
  os << "orca_sline: " << params.aseorca.simpleinput << std::endl;

  os << "\n[Lanczos]" << std::endl;
  os << "lanczosTolerance: " << params.lanczos.tolerance << std::endl;
  os << "lanczosMaxIterations: " << params.lanczos.maxIterations << std::endl;
  os << "lanczosQuitEarly: " << std::boolalpha << params.lanczos.quitEarly
     << std::endl;

  os << "\n[Prefactor]" << std::endl;
  os << "prefactorDefaultValue: " << params.prefactor.defaultValue << std::endl;
  os << "prefactorMaxValue: " << params.prefactor.maxValue << std::endl;
  os << "prefactorMinValue: " << params.prefactor.minValue << std::endl;
  os << "prefactorWithinRadius: " << params.prefactor.withinRadius << std::endl;
  os << "prefactorMinDisplacement: " << params.prefactor.minDisplacement
     << std::endl;
  os << "prefactorRate: " << magic_enum::enum_name(params.prefactor.rate)
     << std::endl;
  os << "prefactorConfiguration: "
     << magic_enum::enum_name(params.prefactor.configuration) << std::endl;
  os << "prefactorAllFreeAtoms: " << std::boolalpha
     << params.prefactor.allFreeAtoms << std::endl;
  os << "prefactorFilterScheme: "
     << magic_enum::enum_name(params.prefactor.filterScheme) << std::endl;
  os << "prefactorFilterFraction: " << params.prefactor.filterFraction
     << std::endl;

  os << "\n[Hessian]" << std::endl;
  os << "hessianAtomList: " << params.hessian.atomList << std::endl;
  os << "hessianZeroFreqValue: " << params.hessian.zeroFreqValue << std::endl;

  os << "\n[Nudged Elastic Band]" << std::endl;
  os << "nebImages: " << params.neb.images << std::endl;
  os << "nebMaxIterations: " << params.neb.maxIterations << std::endl;
  os << "nebSpring: " << params.neb.spring << std::endl;
  os << "nebClimbingImageMethod: " << std::boolalpha
     << params.neb.climbingImageMethod << std::endl;
  os << "nebClimbingImageConvergedOnly: " << std::boolalpha
     << params.neb.climbingImageConvergedOnly << std::endl;
  os << "nebOldTangent: " << std::boolalpha << params.neb.oldTangent
     << std::endl;
  os << "nebDoublyNudged: " << std::boolalpha << params.neb.doublyNudged
     << std::endl;
  os << "nebDoublyNudgedSwitching: " << std::boolalpha
     << params.neb.doublyNudgedSwitching << std::endl;
  os << "nebOptMethod: " << params.neb.optMethod << std::endl;
  os << "nebElasticBand: " << std::boolalpha << params.neb.elasticBand
     << std::endl;
  os << "nebConvergedForce: " << params.neb.convergedForce << std::endl;
  os << "nebKSPMin: " << params.neb.KSPMin << std::endl;
  os << "nebKSPMax: " << params.neb.KSPMax << std::endl;
  os << "nebEnergyWeighted: " << std::boolalpha << params.neb.energyWeighted
     << std::endl;

  os << "\n[Dynamics]" << std::endl;
  os << "mdTimeStepInput: " << params.md.timeStepInput << std::endl;
  os << "mdTimeStep: " << params.md.timeStep << std::endl;
  os << "mdTimeInput: " << params.md.timeInput << std::endl;
  os << "mdTime: " << params.md.time << std::endl;
  os << "mdSteps: " << params.md.steps << std::endl;

  os << "\n[Parallel Replica]" << std::endl;
  os << "parrepRefineTransition: " << std::boolalpha
     << params.parrep.refineTransition << std::endl;
  os << "parrepAutoStop: " << std::boolalpha << params.parrep.autoStop
     << std::endl;
  os << "parrepDephaseLoopStop: " << std::boolalpha
     << params.parrep.dephaseLoopStop << std::endl;
  os << "parrepDephaseTimeInput: " << params.parrep.dephaseTimeInput
     << std::endl;
  os << "parrepDephaseTime: " << params.parrep.dephaseTime << std::endl;
  os << "parrepDephaseLoopMax: " << params.parrep.dephaseLoopMax << std::endl;
  os << "parrepStateCheckIntervalInput: "
     << params.parrep.stateCheckIntervalInput << std::endl;
  os << "parrepStateCheckInterval: " << params.parrep.stateCheckInterval
     << std::endl;
  os << "parrepRecordIntervalInput: " << params.parrep.recordIntervalInput
     << std::endl;
  os << "parrepRecordInterval: " << params.parrep.recordInterval << std::endl;
  os << "parrepCorrTimeInput: " << params.parrep.corrTimeInput << std::endl;
  os << "parrepCorrTime: " << params.parrep.corrTime << std::endl;

  os << "\n[TAD]" << std::endl;
  os << "tadLowT: " << params.tad.lowT << std::endl;
  os << "tadMinPrefactor: " << params.tad.minPrefactor << std::endl;
  os << "tadConfidence: " << params.tad.confidence << std::endl;

  os << "\n[Thermostat]" << std::endl;
  os << "thermostat: " << params.thermostat.kind << std::endl;
  os << "thermoAndersenAlpha: " << params.thermostat.andersenAlpha << std::endl;
  os << "thermoAndersenTcolInput: " << params.thermostat.andersenTcolInput
     << std::endl;
  os << "thermoAndersenTcol: " << params.thermostat.andersenTcol << std::endl;
  os << "thermoNoseMass: " << params.thermostat.noseMass << std::endl;
  os << "thermoLangevinFrictionInput: "
     << params.thermostat.langevinFrictionInput << std::endl;
  os << "thermoLangevinFriction: " << params.thermostat.langevinFriction
     << std::endl;

  os << "\n[Parallel Replica]" << std::endl;
  os << "repexcTemperatureDistribution: "
     << params.repexc.temperatureDistribution << std::endl;
  os << "repexcReplicas: " << params.repexc.replicas << std::endl;
  os << "repexcExchangeTrials: " << params.repexc.exchangeTrials << std::endl;
  os << "repexcSamplingTimeInput: " << params.repexc.samplingTimeInput
     << std::endl;
  os << "repexcSamplingTime: " << params.repexc.samplingTime << std::endl;
  os << "repexcTemperatureHigh: " << params.repexc.temperatureHigh << std::endl;
  os << "repexcTemperatureLow: " << params.repexc.temperatureLow << std::endl;
  os << "repexcExchangePeriodInput: " << params.repexc.exchangePeriodInput
     << std::endl;
  os << "repexcExchangePeriod: " << params.repexc.exchangePeriod << std::endl;

  os << "\n[Hyperdynamics]" << std::endl;
  os << "biasPotential: " << params.bondBoost.biasPotential << std::endl;
  os << "bondBoostBALS: " << params.bondBoost.BALS << std::endl;
  os << "bondBoostRMDTimeInput: " << params.bondBoost.RMDTimeInput << std::endl;
  os << "bondBoostRMDTime: " << params.bondBoost.RMDTime << std::endl;
  os << "bondBoostDVMAX: " << params.bondBoost.DVMAX << std::endl;
  os << "bondBoostQRR: " << params.bondBoost.QRR << std::endl;
  os << "bondBoostPRR: " << params.bondBoost.PRR << std::endl;
  os << "bondBoostQcut: " << params.bondBoost.Qcut << std::endl;

  os << "basinHoppingDisplacement: " << params.bhop.displacement << std::endl;
  os << "basinHoppingInitialRandomStructureProbability: "
     << params.bhop.initialRandomStructureProbability << std::endl;
  os << "basinHoppingPushApartDistance: " << params.bhop.pushApartDistance
     << std::endl;
  os << "basinHoppingSteps: " << params.bhop.steps << std::endl;
  os << "basinHoppingQuenchingSteps: " << params.bhop.quenchingSteps
     << std::endl;
  os << "basinHoppingSignificantStructure: " << std::boolalpha
     << params.bhop.significantStructure << std::endl;
  os << "basinHoppingSingleAtomDisplace: " << std::boolalpha
     << params.bhop.singleAtomDisplace << std::endl;
  os << "basinHoppingDisplacementAlgorithm: "
     << params.bhop.displacementAlgorithm << std::endl;
  os << "basinHoppingDisplacementDistribution: "
     << params.bhop.displacementDistribution << std::endl;
  os << "basinHoppingSwapProbability: " << params.bhop.swapProbability
     << std::endl;
  os << "basinHoppingJumpMax: " << params.bhop.jumpMax << std::endl;
  os << "basinHoppingJumpSteps: " << params.bhop.jumpSteps << std::endl;
  os << "basinHoppingAdjustDisplacement: " << std::boolalpha
     << params.bhop.adjustDisplacement << std::endl;
  os << "basinHoppingAdjustPeriod: " << params.bhop.adjustPeriod << std::endl;
  os << "basinHoppingAdjustFraction: " << params.bhop.adjustFraction
     << std::endl;
  os << "basinHoppingTargetRatio: " << params.bhop.targetRatio << std::endl;
  os << "basinHoppingWriteUnique: " << std::boolalpha << params.bhop.writeUnique
     << std::endl;
  os << "basinHoppingStopEnergy (-DBL_MAX): " << params.bhop.stopEnergy
     << std::endl;

  os << "\n[Global Optimization]" << std::endl;
  os << "globalOptimizationMoveMethod: " << params.globopt.moveMethod
     << std::endl;
  os << "globalOptimizationDecisionMethod: " << params.globopt.decisionMethod
     << std::endl;
  os << "globalOptimizationSteps: " << params.globopt.steps << std::endl;
  os << "globalOptimizationBeta: " << params.globopt.beta << std::endl;
  os << "globalOptimizationAlpha: " << params.globopt.alpha << std::endl;
  os << "globalOptimizationMdmin: " << params.globopt.mdmin << std::endl;
  os << "globalOptimizationTargetEnergy: " << params.globopt.targetEnergy
     << std::endl;

  os << "\n[Monte Carlo]" << std::endl;
  os << "monteCarloStepSize: " << params.monte_carlo.stepSize << std::endl;
  os << "monteCarloSteps: " << params.monte_carlo.steps << std::endl;

  os << "\n[BGSD]" << std::endl;
  os << "alpha: " << params.bgsd.alpha << std::endl;
  os << "beta: " << params.bgsd.beta << std::endl;
  os << "gradientfinitedifference: " << params.bgsd.gradientfinitedifference
     << std::endl;
  os << "Hforceconvergence: " << params.bgsd.Hforceconvergence << std::endl;
  os << "grad2energyconvergence: " << params.bgsd.grad2energyconvergence
     << std::endl;
  os << "grad2forceconvergence: " << params.bgsd.grad2forceconvergence
     << std::endl;

  os << "\n[Debug]" << std::endl;
  os << "writeMovies: " << std::boolalpha << params.debug.writeMovies
     << std::endl;
  os << "writeMoviesInterval: " << params.debug.writeMoviesInterval
     << std::endl;
  return os;
}

std::vector<Parameters> getParameters() {
  // Return test data for Parameters
  return {Parameters()};
}

TEST_CASE("VerifyParameters") {
  ApprovalTests::Approvals::verifyAll("parameters", getParameters());
}
