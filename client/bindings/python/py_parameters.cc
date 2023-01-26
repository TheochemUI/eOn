#include "py_wrapper.hpp"

void py_parameters(py::module_ &m) {
    py::class_<Parameters>(m, "Parameters", py::dynamic_attr()) // dynamic incurs a penalty
        .def(py::init())
        /*
        ** Functions
        */

        .def("load", py::overload_cast<std::string>(&Parameters::load))
        // .def("load", py::overload_cast<FILE>(&Parameters::load))

        /*
        ** Parameters
        */

        // Constants
        .def_readwrite("kB", &Parameters::kB)
        .def_readwrite("timeUnit", &Parameters::timeUnit)

        // [Main] //
        .def_readwrite("job", &Parameters::job)
        .def_readwrite("randomSeed", &Parameters::randomSeed)
        .def_readwrite("temperature", &Parameters::temperature)
        .def_readwrite("quiet", &Parameters::quiet)
        .def_readwrite("writeLog", &Parameters::writeLog)
        .def_readwrite("checkpoint", &Parameters::checkpoint)
        .def_readwrite("iniFilename", &Parameters::iniFilename)
        .def_readwrite("conFilename", &Parameters::conFilename)
        .def_readwrite("finiteDifference", &Parameters::finiteDifference)
        .def_readwrite("maxForceCalls", &Parameters::maxForceCalls)
        .def_readwrite("removeNetForce", &Parameters::removeNetForce)

        // [Potential] //
        .def_readwrite("potential", &Parameters::potential)
        .def_readwrite("MPIPollPeriod", &Parameters::MPIPollPeriod)
        .def_readwrite("LAMMPSLogging", &Parameters::LAMMPSLogging)
        .def_readwrite("LAMMPSThreads", &Parameters::LAMMPSThreads)
        .def_readwrite("EMTRasmussen", &Parameters::EMTRasmussen)
        .def_readwrite("LogPotential", &Parameters::LogPotential)
        .def_readwrite("extPotPath", &Parameters::extPotPath)

        // [AMS] and [AMS_IO] //
        .def_readwrite("engine", &Parameters::engine) // MOPAC, ADF, BAND, REAXFF, FORCEFIELD
        .def_readwrite("forcefield", &Parameters::forcefield) // OPt.ff etc. (REAXFF)
        .def_readwrite("model", &Parameters::model)           // Model hamiltonian (MOPAC)
        .def_readwrite("resources", &Parameters::resources)   // DFTB
        .def_readwrite("xc", &Parameters::xc)                 // Exchange (BAND, ADF)
        .def_readwrite("basis", &Parameters::basis)           // Basis (BAND, ADF)

        // [AMS_ENV] //
        .def_readwrite("amshome", &Parameters::amshome)
        .def_readwrite("scm_tmpdir", &Parameters::scm_tmpdir)
        .def_readwrite("scmlicense", &Parameters::scmlicense)
        .def_readwrite("scm_pythondir", &Parameters::scm_pythondir)
        .def_readwrite("amsbin", &Parameters::amsbin)
        .def_readwrite("amsresources", &Parameters::amsresources)

        // [Structure Comparison] //
        .def_readwrite(
            "distanceDifference",
            &Parameters::distanceDifference) // The distance criterion for comparing geometries
        .def_readwrite(
            "neighborCutoff",
            &Parameters::neighborCutoff) // radius used in the local atomic structure analysis
        .def_readwrite("checkRotation", &Parameters::checkRotation)
        .def_readwrite("indistinguishableAtoms", &Parameters::indistinguishableAtoms)
        .def_readwrite("energyDifference", &Parameters::energyDifference)
        .def_readwrite("removeTranslation", &Parameters::removeTranslation)

        // [Process Search] //
        .def_readwrite("processSearchMinimizeFirst", &Parameters::processSearchMinimizeFirst)
        .def_readwrite(
            "processSearchMinimizationOffset",
            &Parameters::processSearchMinimizationOffset) // how far from the saddle to displace
                                                          // the minimization images

        // [Saddle Search] //
        .def_readwrite("saddleMaxJumpAttempts",
                       &Parameters::saddleMaxJumpAttempts) // number of displacements to reach a
                                                           // convex region;  if 0, a search is
                                                           // started after the displacement
        .def_readwrite("saddleMaxIterations",
                       &Parameters::saddleMaxIterations) // max iterations for saddle point
                                                         // searches and minimization
        .def_readwrite("saddleMethod", &Parameters::saddleMethod)
        .def_readwrite("saddleMinmodeMethod",
                       &Parameters::saddleMinmodeMethod) // algorithm to be used for lowest
                                                         // eigenmode determination
        .def_readwrite("saddleDisplaceType",
                       &Parameters::saddleDisplaceType) // displacement type to use
        .def_readwrite("saddleMaxEnergy",
                       &Parameters::saddleMaxEnergy) // energy above product state that will cause
                                                     // termination of the saddle point search
        .def_readwrite("saddleDisplaceMagnitude",
                       &Parameters::saddleDisplaceMagnitude) // norm of the displacement vector
        .def_readwrite(
            "saddleMaxSingleDisplace",
            &Parameters::saddleMaxSingleDisplace) // maximum value of displacement in x, y and z
                                                  // direction for atoms being displaced
        .def_readwrite("saddleDisplaceRadius",
                       &Parameters::saddleDisplaceRadius) // atoms within this radius of the
                                                          // displacement atoms are also displaced
        .def_readwrite("saddleConvergedForce",
                       &Parameters::saddleConvergedForce) // force convergence criterion required
                                                          // for a saddle point search
        .def_readwrite(
            "saddlePerpForceRatio",
            &Parameters::saddlePerpForceRatio) // proportion to keep of the perpendicular force
                                               // when the lowest eigenvalue is positive
        .def_readwrite("saddleNonnegativeDisplacementAbort",
                       &Parameters::saddleNonnegativeDisplacementAbort) // abort the saddle search
                                                                        // if the displacement does
                                                                        // not have a negative mode
        .def_readwrite(
            "saddleNonlocalCountAbort",
            &Parameters::saddleNonlocalCountAbort) // abort the search if this many atoms move more
                                                   // than NonlocalDistanceAbort
        .def_readwrite(
            "saddleNonlocalDistanceAbort",
            &Parameters::saddleNonlocalDistanceAbort) // abort the search if NonlocalCountAbort
                                                      // atoms move more than this distance
        .def_readwrite("saddleRemoveRotation",
                       &Parameters::saddleRemoveRotation) // remove dominant rotational component
                                                          // when system is translated
        .def_readwrite("saddleDynamicsTemperature",
                       &Parameters::saddleDynamicsTemperature) // temperature for dynamics saddle
                                                               // search method
        .def_readwrite("saddleDynamicsStateCheckIntervalInput",
                       &Parameters::saddleDynamicsStateCheckIntervalInput)
        .def_readwrite("saddleDynamicsStateCheckInterval",
                       &Parameters::saddleDynamicsStateCheckInterval) // how often to minimize
        .def_readwrite("saddleDynamicsRecordIntervalInput",
                       &Parameters::saddleDynamicsRecordIntervalInput)
        .def_readwrite("saddleDynamicsRecordInterval", &Parameters::saddleDynamicsRecordInterval)
        .def_readwrite("saddleDynamicsLinearInterpolation",
                       &Parameters::saddleDynamicsLinearInterpolation)
        .def_readwrite("saddleDynamicsMaxInitCurvature",
                       &Parameters::saddleDynamicsMaxInitCurvature)

        .def_readwrite("saddleConfinePositive", &Parameters::saddleConfinePositive) // undocumented
        .def_readwrite("saddleBowlBreakout", &Parameters::saddleBowlBreakout)       // undocumented
        .def_readwrite("saddleBowlActive", &Parameters::saddleBowlActive)           // undocumented
        .def_readwrite("saddleConfinePositiveMinForce",
                       &Parameters::saddleConfinePositiveMinForce) // undocumented
        .def_readwrite("saddleConfinePositiveScaleRatio",
                       &Parameters::saddleConfinePositiveScaleRatio) // undocumented
        .def_readwrite("saddleConfinePositiveBoost",
                       &Parameters::saddleConfinePositiveBoost) // undocumented
        .def_readwrite("saddleConfinePositiveMinActive",
                       &Parameters::saddleConfinePositiveMinActive) // undocumented
        .def_readwrite("saddleZeroModeAbortCurvature", &Parameters::saddleZeroModeAbortCurvature)

        // [Optimizer] //
        .def_readwrite("optMethod", &Parameters::optMethod)
        .def_readwrite("optConvergenceMetric",
                       &Parameters::optConvergenceMetric) // norm, max_atom, max_component
        .def_readwrite("optConvergenceMetricLabel", &Parameters::optConvergenceMetricLabel)
        .def_readwrite("optMaxIterations",
                       &Parameters::optMaxIterations) // maximum iterations for saddle point
                                                      // searches and minimization
        .def_readwrite(
            "optMaxMove",
            &Parameters::optMaxMove) // maximum displacement vector for a step during optimization
        .def_readwrite("optConvergedForce",
                       &Parameters::optConvergedForce) // force convergence criterion required for
                                                       // an optimization
        .def_readwrite("optTimeStepInput", &Parameters::optTimeStepInput)
        .def_readwrite("optTimeStep", &Parameters::optTimeStep) // time step size used in quickmin
        .def_readwrite("optMaxTimeStepInput",
                       &Parameters::optMaxTimeStepInput) // maximum time step for FIRE.
        .def_readwrite("optMaxTimeStep",
                       &Parameters::optMaxTimeStep) // maximum time step for FIRE.
        .def_readwrite(
            "optLBFGSMemory",
            &Parameters::optLBFGSMemory) // number of previous forces to keep in the bfgs memory
        .def_readwrite("optLBFGSInverseCurvature", &Parameters::optLBFGSInverseCurvature)
        .def_readwrite("optLBFGSMaxInverseCurvature", &Parameters::optLBFGSMaxInverseCurvature)
        .def_readwrite("optLBFGSAutoScale", &Parameters::optLBFGSAutoScale)
        .def_readwrite("optLBFGSAngleReset", &Parameters::optLBFGSAngleReset)
        .def_readwrite("optLBFGSDistanceReset", &Parameters::optLBFGSDistanceReset)
        .def_readwrite("optQMSteepestDecent",
                       &Parameters::optQMSteepestDecent) // if set the velocity will always be set
                                                         // to zero in quickmin
        .def_readwrite(
            "optCGNoOvershooting",
            &Parameters::optCGNoOvershooting) // if set it is ensured that the approximate line
                                              // search in conjugate gradients never overshoot the
                                              // minimum along the search line
        .def_readwrite(
            "optCGKnockOutMaxMove",
            &Parameters::optCGKnockOutMaxMove) // if set the old search direction is nullified when
                                               // steps larger than the optMaxMove are conducted
        .def_readwrite("optCGLineSearch",
                       &Parameters::optCGLineSearch) // if set full line search is conducted
        .def_readwrite("optCGLineConverged",
                       &Parameters::optCGLineConverged) // convergence criteria for line search,
                                                        // ratio between force component along
                                                        // search line and the orthogonal part
        .def_readwrite(
            "optCGLineSearchMaxIter",
            &Parameters::optCGLineSearchMaxIter) // maximal nr of iterations during line search
        .def_readwrite("optCGMaxIterBeforeReset",
                       &Parameters::optCGMaxIterBeforeReset) // max nr of cg steps before reset, if
                                                             // 0 no resetting is done
        .def_readwrite("optSDAlpha", &Parameters::optSDAlpha)
        .def_readwrite("optSDTwoPoint", &Parameters::optSDTwoPoint)

        // [Dimer] //
        .def_readwrite("dimerRotationAngle",
                       &Parameters::dimerRotationAngle) // finite difference rotation angle
        .def_readwrite("dimerImproved",
                       &Parameters::dimerImproved) // turn on the improved dimer method
        .def_readwrite(
            "dimerConvergedAngle",
            &Parameters::dimerConvergedAngle) // stop rotating when angle drops below this value
        .def_readwrite("dimerMaxIterations",
                       &Parameters::dimerMaxIterations) // maximum number of rotation iterations
        .def_readwrite(
            "dimerOptMethod",
            &Parameters::dimerOptMethod) // method to determine the next rotation direction
        .def_readwrite("dimerRotationsMax", &Parameters::dimerRotationsMax) // old
        .def_readwrite("dimerRotationsMin", &Parameters::dimerRotationsMin) // old
        .def_readwrite("dimerTorqueMax", &Parameters::dimerTorqueMax)       // old
        .def_readwrite("dimerTorqueMin", &Parameters::dimerTorqueMin)       // old
        .def_readwrite("dimerRemoveRotation",
                       &Parameters::dimerRemoveRotation) // remove dominant rotational component
                                                         // when estimating the eigenmode

        // [GPR Dimer] //
        .def_readwrite("gprDimerRotationAngle",
                       &Parameters::gprDimerRotationAngle) // finite difference rotation angle
        .def_readwrite("gprDimerConvergedAngle",
                       &Parameters::gprDimerConvergedAngle) // stop rotating when angle drops below
                                                            // this value {T_anglerot_init}
        .def_readwrite(
            "gprDimerRelaxConvAngle",
            &Parameters::gprDimerRelaxConvAngle) // stop rotating when angle drops below this value
                                                 // during relaxation {T_anglerot_gp}
        .def_readwrite("gprDimerInitRotationsMax",
                       &Parameters::gprDimerInitRotationsMax) // {num_iter_initrot}
        .def_readwrite("gprDimerRelaxRotationsMax",
                       &Parameters::gprDimerRelaxRotationsMax) // {num_iter_rot_gp}
        .def_readwrite("gprDimerDivisorTdimerGP",
                       &Parameters::gprDimerDivisorTdimerGP) // {divisor_T_dimer_gp}
        .def_readwrite(
            "gprDimerMaxOuterIterations",
            &Parameters::gprDimerMaxOuterIterations) // maximum number of outer iterations or new
                                                     // sets of observations {num_bigiter}
        .def_readwrite("gprDimerMaxInnerIterations",
                       &Parameters::gprDimerMaxInnerIterations) // maximum number of steps during
                                                                // the relaxation phase {num_iter}
        .def_readwrite("gprDimerMidpointMaxDisp",
                       &Parameters::gprDimerMidpointMaxDisp) // {disp_max}
        .def_readwrite("gprDimerRotOptMethod",
                       &Parameters::gprDimerRotOptMethod) // method to determine the next rotation
                                                          // direction {method_rot}
        .def_readwrite("gprDimerTransOptMethod",
                       &Parameters::gprDimerTransOptMethod) // method to determine the next
                                                            // rotation direction {method_trans}
        .def_readwrite("gprActiveRadius",
                       &Parameters::gprActiveRadius) // activation radius for inclusion in
                                                     // covariance matrix {actdist_fro}
        .def_readwrite("gprDimerSep",
                       &Parameters::gprDimerSep) // distance from the middle point of the dimer to
                                                 // the two images {dimer_sep}
        .def_readwrite(
            "gprDimerConvStep",
            &Parameters::gprDimerConvStep) // step length for convex regions {param_trans[0]}
        .def_readwrite("gprDimerMaxStep",
                       &Parameters::gprDimerMaxStep) // maximum step length {param_trans[1]}
        .def_readwrite(
            "gprForceThreshold",
            &Parameters::gprForceThreshold) // maximum component of the force acting on the middle
                                            // point of the dimer (stops when accurate force
                                            // components are below this value) {T_dimer}
        .def_readwrite(
            "gprDimerRatioAtLimit",
            &Parameters::gprDimerRatioAtLimit) // {ratio_at_limit} defines the limit for the ratio
                                               // of inter-atomic distances between the image and
                                               // its "nearest observed data point"
        .def_readwrite("gprDimerInitRotGP",
                       &Parameters::gprDimerInitRotGP) // initial rotations without GP (1) or with
                                                       // GP (0) {initrot_nogp}
        .def_readwrite("gprDimerInitTransGP",
                       &Parameters::gprDimerInitTransGP) // initial translations without GP (1) or
                                                         // with GP (0) {inittrans_nogp}
        .def_readwrite(
            "gprDimerManyIterations",
            &Parameters::gprDimerManyIterations) // indicates of the number of iterations is larger
                                                 // than required for dimer convergence on the
                                                 // accurate energy surface (1), otherwise (0) the
                                                 // relaxation phase is continued from the current
                                                 // dimer if the maximum iterations in the
                                                 // relaxation phase is reached {islarge_num_iter}
                                                 // GPR Params
        .def_readwrite("gprDimerHyperOptMethod",
                       &Parameters::gprDimerHyperOptMethod) // method to optimize hyperparameters
                                                            // {optimization_alg}
        .def_readwrite("gprDimerSigma2", &Parameters::gprDimerSigma2) // GPR variance {gp_sigma2}
        .def_readwrite("gprDimerJitterSigma2",
                       &Parameters::gprDimerJitterSigma2) // GPR jitter variance {jitter_sigma2}
        .def_readwrite("gprDimerNoiseSigma2",
                       &Parameters::gprDimerNoiseSigma2)                // noise Variance {sigma2}
        .def_readwrite("gprDimerPriorMu", &Parameters::gprDimerPriorMu) // prior mean {prior_mu}
        .def_readwrite("gprDimerPriorSigma2",
                       &Parameters::gprDimerPriorSigma2) // prior variance {prior_s2}
        .def_readwrite("gprDimerPriorNu",
                       &Parameters::gprDimerPriorNu) // prior degrees of freedom {prior_nu}
                                                     // GPR Optimization Parameters
        .def_readwrite("gprOptCheckDerivatives",
                       &Parameters::gprOptCheckDerivatives) // {check_derivative}
        .def_readwrite("gprOptMaxIterations", &Parameters::gprOptMaxIterations) // {max_iter}
        .def_readwrite("gprOptTolFunc", &Parameters::gprOptTolFunc)             // {tolerance_func}
        .def_readwrite("gprOptTolSol", &Parameters::gprOptTolSol)               // {tolerance_sol}
        .def_readwrite("gprOptLambdaLimit", &Parameters::gprOptLambdaLimit)     // {lambda_limit}
        .def_readwrite("gprOptLambdaInit", &Parameters::gprOptLambdaInit)       // {lambda}
        .def_readwrite("gprUsePrune", &Parameters::gprUsePrune)                 // {use_prune}
        .def_readwrite("gprPruneBegin", &Parameters::gprPruneBegin)             // {start_prune_at}
        .def_readwrite("gprPruneNVals", &Parameters::gprPruneNVals)             // {nprune_vals}
        .def_readwrite("gprPruneThreshold", &Parameters::gprPruneThreshold) // {prune_threshold}

        // GPR Debugging Parameters
        .def_readwrite("gprReportLevel", &Parameters::gprReportLevel)   // {report_level}
        .def_readwrite("gprDebugLevel", &Parameters::gprDebugLevel)     // {debug_level}
        .def_readwrite("gprDebugOutDir", &Parameters::gprDebugOutDir)   // {debug_output_dir}
        .def_readwrite("gprDebugPosFile", &Parameters::gprDebugPosFile) // {debug_output_file_R}
        .def_readwrite("gprDebugEnergyFile",
                       &Parameters::gprDebugEnergyFile)                   // {debug_output_file_E}
        .def_readwrite("gprDebugGradFile", &Parameters::gprDebugGradFile) // {debug_output_file_G}
        .def_readwrite("gprDebugOutExt",
                       &Parameters::gprDebugOutExt) // {debug_output_file_extension}
        .def_readwrite("gprDebugOffsetMidPoint",
                       &Parameters::gprDebugOffsetMidPoint)   // {debug_offset_from_mid_point}
        .def_readwrite("gprDebugDy", &Parameters::gprDebugDy) // {debug_dy}
        .def_readwrite("gprDebugDz", &Parameters::gprDebugDz) // {debug_dz}

        // [Lanczos] //
        .def_readwrite("lanczosTolerance",
                       &Parameters::lanczosTolerance) // difference between the lowest eignevalues
                                                      // of two successive iterations
        .def_readwrite("lanczosMaxIterations",
                       &Parameters::lanczosMaxIterations) // maximum number of iterations
        .def_readwrite("lanczosQuitEarly", &Parameters::lanczosQuitEarly)

        // [Prefactor] //
        .def_readwrite(
            "prefactorDefaultValue",
            &Parameters::prefactorDefaultValue) // default prefactor; calculate explicitly if zero
        .def_readwrite("prefactorMaxValue",
                       &Parameters::prefactorMaxValue) // max prefactor allowed
        .def_readwrite("prefactorMinValue",
                       &Parameters::prefactorMinValue) // min prefactor allowed
        .def_readwrite(
            "prefactorWithinRadius",
            &Parameters::prefactorWithinRadius) // atoms within this radius of the displaced atoms
                                                // are put in the Hessian, unless filterMode is
                                                // fraction
        .def_readwrite(
            "prefactorMinDisplacement",
            &Parameters::prefactorMinDisplacement) // atoms with displacement between min1 or min2
                                                   // and the saddle point are put in the Hessian
        .def_readwrite("prefactorRate", &Parameters::prefactorRate) // method to estimate prefactor
        .def_readwrite("prefactorConfiguration",
                       &Parameters::prefactorConfiguration) // configuration for which the
                                                            // frequencies should be determined
        .def_readwrite(
            "prefactorAllFreeAtoms",
            &Parameters::prefactorAllFreeAtoms) // use all free atom when determining the prefactor
        .def_readwrite(
            "prefactorFilterScheme",
            &Parameters::prefactorFilterScheme) // "cutoff" or "fraction", which use
                                                // prefactorMinDisplacement or
                                                // prefactorFilterFraction, respectively.
        .def_readwrite(
            "prefactorFilterFraction",
            &Parameters::prefactorFilterFraction) // Include atoms whose summed motion comprise
                                                  // more than prefactorFilterFraction in the
                                                  // prefactor calculation. Prioritizes atoms that
                                                  // move more.

        // [Hessian] //
        .def_readwrite("hessianAtomList", &Parameters::hessianAtomList)
        .def_readwrite("hessianZeroFreqValue", &Parameters::hessianZeroFreqValue)

        // [Nudged Elastic Band] //
        .def_readwrite("nebImages", &Parameters::nebImages)
        .def_readwrite("nebMaxIterations", &Parameters::nebMaxIterations)
        .def_readwrite("nebSpring", &Parameters::nebSpring)
        .def_readwrite("nebClimbingImageMethod", &Parameters::nebClimbingImageMethod)
        .def_readwrite("nebClimbingImageConvergedOnly", &Parameters::nebClimbingImageConvergedOnly)
        .def_readwrite("nebOldTangent", &Parameters::nebOldTangent)
        .def_readwrite("nebDoublyNudged", &Parameters::nebDoublyNudged)
        .def_readwrite("nebDoublyNudgedSwitching", &Parameters::nebDoublyNudgedSwitching)
        .def_readwrite("nebOptMethod", &Parameters::nebOptMethod)
        .def_readwrite("nebElasticBand", &Parameters::nebElasticBand)
        .def_readwrite("nebConvergedForce",
                       &Parameters::nebConvergedForce) // force convergence criterion required for
                                                       // an optimization

        // [Molecular Dynamics] //
        .def_readwrite("mdTimeStepInput", &Parameters::mdTimeStepInput)
        .def_readwrite("mdTimeStep", &Parameters::mdTimeStep)
        .def_readwrite("mdTimeInput", &Parameters::mdTimeInput)
        .def_readwrite("mdTime", &Parameters::mdTime)
        .def_readwrite("mdSteps", &Parameters::mdSteps)

        // [Parallel Replica] //
        .def_readwrite("parrepRefineTransition", &Parameters::parrepRefineTransition)
        .def_readwrite("parrepAutoStop", &Parameters::parrepAutoStop)
        .def_readwrite("parrepDephaseLoopStop", &Parameters::parrepDephaseLoopStop)
        .def_readwrite("parrepDephaseTimeInput", &Parameters::parrepDephaseTimeInput)
        .def_readwrite("parrepDephaseTime", &Parameters::parrepDephaseTime)
        .def_readwrite("parrepDephaseLoopMax", &Parameters::parrepDephaseLoopMax)
        .def_readwrite("parrepStateCheckIntervalInput", &Parameters::parrepStateCheckIntervalInput)
        .def_readwrite("parrepStateCheckInterval", &Parameters::parrepStateCheckInterval)
        .def_readwrite("parrepRecordIntervalInput", &Parameters::parrepRecordIntervalInput)
        .def_readwrite("parrepRecordInterval", &Parameters::parrepRecordInterval)
        .def_readwrite("parrepCorrTimeInput", &Parameters::parrepCorrTimeInput)
        .def_readwrite("parrepCorrTime", &Parameters::parrepCorrTime)

        // [Temperature Accelerated Dynamics] //
        .def_readwrite("tadLowT", &Parameters::tadLowT)
        .def_readwrite("tadMinPrefactor", &Parameters::tadMinPrefactor)
        .def_readwrite("tadConfidence", &Parameters::tadConfidence)

        // [Thermostat] //
        .def_readwrite("thermostat", &Parameters::thermostat)
        .def_readwrite("thermoAndersenAlpha", &Parameters::thermoAndersenAlpha)
        .def_readwrite("thermoAndersenTcolInput", &Parameters::thermoAndersenTcolInput)
        .def_readwrite("thermoAndersenTcol", &Parameters::thermoAndersenTcol)
        .def_readwrite("thermoNoseMass", &Parameters::thermoNoseMass)
        .def_readwrite("thermoLangevinFrictionInput", &Parameters::thermoLangevinFrictionInput)
        .def_readwrite("thermoLangevinFriction", &Parameters::thermoLangevinFriction)
        // std::vector<int> thermoAtoms;

        // [Replica Exchange] //
        .def_readwrite("repexcTemperatureDistribution", &Parameters::repexcTemperatureDistribution)
        .def_readwrite("repexcReplicas", &Parameters::repexcReplicas)
        .def_readwrite("repexcExchangeTrials", &Parameters::repexcExchangeTrials)
        .def_readwrite("repexcSamplingTimeInput", &Parameters::repexcSamplingTimeInput)
        .def_readwrite("repexcSamplingTime", &Parameters::repexcSamplingTime)
        .def_readwrite("repexcTemperatureHigh", &Parameters::repexcTemperatureHigh)
        .def_readwrite("repexcTemperatureLow", &Parameters::repexcTemperatureLow)
        .def_readwrite("repexcExchangePeriodInput", &Parameters::repexcExchangePeriodInput)
        .def_readwrite("repexcExchangePeriod", &Parameters::repexcExchangePeriod)

        // [Bond Boost] //
        .def_readwrite("biasPotential", &Parameters::biasPotential)
        .def_readwrite("bondBoostBALS", &Parameters::bondBoostBALS)
        .def_readwrite("bondBoostRMDTimeInput", &Parameters::bondBoostRMDTimeInput)
        .def_readwrite("bondBoostRMDTime", &Parameters::bondBoostRMDTime)
        .def_readwrite("bondBoostDVMAX", &Parameters::bondBoostDVMAX)
        .def_readwrite("bondBoostQRR", &Parameters::bondBoostQRR)
        .def_readwrite("bondBoostPRR", &Parameters::bondBoostPRR)
        .def_readwrite("bondBoostQcut", &Parameters::bondBoostQcut)

        // [Basin Hopping] //
        .def_readwrite("basinHoppingDisplacement", &Parameters::basinHoppingDisplacement)
        .def_readwrite("basinHoppingInitialRandomStructureProbability",
                       &Parameters::basinHoppingInitialRandomStructureProbability)
        .def_readwrite("basinHoppingPushApartDistance", &Parameters::basinHoppingPushApartDistance)
        .def_readwrite("basinHoppingSteps", &Parameters::basinHoppingSteps)
        .def_readwrite("basinHoppingQuenchingSteps", &Parameters::basinHoppingQuenchingSteps)
        .def_readwrite("basinHoppingSignificantStructure",
                       &Parameters::basinHoppingSignificantStructure)
        .def_readwrite("basinHoppingSingleAtomDisplace",
                       &Parameters::basinHoppingSingleAtomDisplace)
        .def_readwrite("basinHoppingDisplacementAlgorithm",
                       &Parameters::basinHoppingDisplacementAlgorithm)
        .def_readwrite("basinHoppingDisplacementDistribution",
                       &Parameters::basinHoppingDisplacementDistribution)
        .def_readwrite("basinHoppingSwapProbability", &Parameters::basinHoppingSwapProbability)
        .def_readwrite("basinHoppingJumpMax", &Parameters::basinHoppingJumpMax)
        .def_readwrite("basinHoppingJumpSteps", &Parameters::basinHoppingJumpSteps)
        .def_readwrite("basinHoppingAdjustDisplacement",
                       &Parameters::basinHoppingAdjustDisplacement)
        .def_readwrite("basinHoppingAdjustPeriod", &Parameters::basinHoppingAdjustPeriod)
        .def_readwrite("basinHoppingAdjustFraction", &Parameters::basinHoppingAdjustFraction)
        .def_readwrite("basinHoppingTargetRatio", &Parameters::basinHoppingTargetRatio)
        .def_readwrite("basinHoppingWriteUnique", &Parameters::basinHoppingWriteUnique)
        .def_readwrite("basinHoppingStopEnergy", &Parameters::basinHoppingStopEnergy)

        // [Global Optimization] //
        .def_readwrite("globalOptimizationMoveMethod", &Parameters::globalOptimizationMoveMethod)
        .def_readwrite("globalOptimizationDecisionMethod",
                       &Parameters::globalOptimizationDecisionMethod)
        .def_readwrite("globalOptimizationSteps", &Parameters::globalOptimizationSteps)
        .def_readwrite("globalOptimizationBeta", &Parameters::globalOptimizationBeta)
        .def_readwrite("globalOptimizationAlpha", &Parameters::globalOptimizationAlpha)
        .def_readwrite("globalOptimizationMdmin", &Parameters::globalOptimizationMdmin)
        .def_readwrite("globalOptimizationTargetEnergy",
                       &Parameters::globalOptimizationTargetEnergy)

        // [Monte Carlo] //
        .def_readwrite("monteCarloStepSize", &Parameters::monteCarloStepSize)
        .def_readwrite("monteCarloSteps", &Parameters::monteCarloSteps)

        // [BGSD] //

        .def_readwrite("alpha", &Parameters::alpha)
        .def_readwrite("beta", &Parameters::beta)
        .def_readwrite("gradientfinitedifference", &Parameters::gradientfinitedifference)
        .def_readwrite("Hforceconvergence", &Parameters::Hforceconvergence)
        .def_readwrite("grad2energyconvergence", &Parameters::grad2energyconvergence)
        .def_readwrite("grad2forceconvergence", &Parameters::grad2forceconvergence)

        // MPI stuff, not actually specified in config file
        // it is used to pass information to the GPAW MPI potential.
        .def_readwrite("MPIPotentialRank", &Parameters::MPIPotentialRank)
#ifdef EONMPI
        .def_readwrite("_Comm", &Parameters::_Comm) MPIClientComm
#endif
        // [Debug] //
        .def_readwrite("writeMovies", &Parameters::writeMovies)
        .def_readwrite("writeMoviesInterval", &Parameters::writeMoviesInterval)
        /*
        ** Python helpers
        */

        .def("__repr__", [](const Parameters &a) { return "<Parameter object>"; });
}
