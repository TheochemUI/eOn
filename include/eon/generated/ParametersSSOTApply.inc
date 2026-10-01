// AUTO-GENERATED from schema/eon_params.capnp — do not edit.
// Regenerate: python tools/params_ssot/codegen.py
// Included at the end of apply_ssot_defaults. GD is in scope.
  {
    auto &o = ParametersLoadAccess::rgpot_options(p);
    o.backend = std::string(GD::RGPOTPOT_BACKEND);
    o.basis = std::string(GD::RGPOTPOT_BASIS);
    o.theory = std::string(GD::RGPOTPOT_THEORY);
    o.scf_type = std::string(GD::RGPOTPOT_SCF_TYPE);
    o.functional = std::string(GD::RGPOTPOT_FUNCTIONAL);
    o.cutoff_ry = GD::RGPOTPOT_CUTOFF_RY;
    o.charge = GD::RGPOTPOT_CHARGE;
    o.multiplicity = GD::RGPOTPOT_MULTIPLICITY;
    o.engine_path = std::string(GD::RGPOTPOT_ENGINE_PATH);
    o.engine_library = std::string(GD::RGPOTPOT_ENGINE_LIBRARY);
    o.engine_root = std::string(GD::RGPOTPOT_ENGINE_ROOT);
    o.title = std::string(GD::RGPOTPOT_TITLE);
    o.memory_mb = GD::RGPOTPOT_MEMORY_MB;
    o.scratch_dir = std::string(GD::RGPOTPOT_SCRATCH_DIR);
    o.input_block = std::string(GD::RGPOTPOT_INPUT_BLOCK);
    o.permanent_dir = std::string(GD::RGPOTPOT_PERMANENT_DIR);
    o.params_path = std::string(GD::RGPOTPOT_PARAMS_PATH);
    o.ranks_per_image = GD::RGPOTPOT_RANKS_PER_IMAGE;
    o.model_path = std::string(GD::RGPOTPOT_MODEL_PATH);
    o.device = std::string(GD::RGPOTPOT_DEVICE);
    o.length_unit = std::string(GD::RGPOTPOT_LENGTH_UNIT);
    o.extensions_directory = std::string(GD::RGPOTPOT_EXTENSIONS_DIRECTORY);
    o.check_consistency = GD::RGPOTPOT_CHECK_CONSISTENCY;
    o.uncertainty_threshold = GD::RGPOTPOT_UNCERTAINTY_THRESHOLD;
    o.torch_determinism_strict = GD::RGPOTPOT_TORCH_DETERMINISM_STRICT;
    o.xtb_paramset = std::string(GD::RGPOTPOT_XTB_PARAMSET);
    o.xtb_accuracy = GD::RGPOTPOT_XTB_ACCURACY;
    o.xtb_electronic_temperature = GD::RGPOTPOT_XTB_ELECTRONIC_TEMPERATURE;
    o.xtb_max_iterations = GD::RGPOTPOT_XTB_MAX_ITERATIONS;
    o.xtb_charge = GD::RGPOTPOT_XTB_CHARGE;
    o.xtb_uhf = GD::RGPOTPOT_XTB_UHF;
  }
