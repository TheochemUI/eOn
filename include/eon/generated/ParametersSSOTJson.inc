// AUTO-GENERATED from schema/eon_params.capnp — do not edit.
// Regenerate: python tools/params_ssot/codegen.py
// Included inside namespace eonc::config, after JSON_OPT.
inline void project_ssot_json_write(json &j, const Parameters &p) {
  {
    const auto &o = ParametersLoadAccess::rgpot_options(p);
    j["RgpotPot"] = {
        {"backend", o.backend},
        {"basis", o.basis},
        {"theory", o.theory},
        {"scf_type", o.scf_type},
        {"functional", o.functional},
        {"cutoff_ry", o.cutoff_ry},
        {"charge", o.charge},
        {"multiplicity", o.multiplicity},
        {"engine_path", o.engine_path},
        {"engine_library", o.engine_library},
        {"engine_root", o.engine_root},
        {"title", o.title},
        {"memory_mb", o.memory_mb},
        {"scratch_dir", o.scratch_dir},
        {"input_block", o.input_block},
        {"permanent_dir", o.permanent_dir},
        {"params_path", o.params_path},
        {"ranks_per_image", o.ranks_per_image},
        {"model_path", o.model_path},
        {"device", o.device},
        {"length_unit", o.length_unit},
        {"extensions_directory", o.extensions_directory},
        {"check_consistency", o.check_consistency},
        {"uncertainty_threshold", o.uncertainty_threshold},
        {"torch_determinism_strict", o.torch_determinism_strict},
        {"xtb_paramset", o.xtb_paramset},
        {"xtb_accuracy", o.xtb_accuracy},
        {"xtb_electronic_temperature", o.xtb_electronic_temperature},
        {"xtb_max_iterations", o.xtb_max_iterations},
        {"xtb_charge", o.xtb_charge},
        {"xtb_uhf", o.xtb_uhf},
        {"task_name", o.task_name}
    };
  }
}

inline void project_ssot_json_read(const json &j, Parameters &p) {
  if (j.contains("RgpotPot")) {
    auto &s = j.at("RgpotPot");
    auto &o = ParametersLoadAccess::rgpot_options(p);
    if (s.contains("backend")) o.backend = s.at("backend").get<decltype(o.backend)>();
    if (s.contains("basis")) o.basis = s.at("basis").get<decltype(o.basis)>();
    else if (s.contains("nwchem_basis")) o.basis = s.at("nwchem_basis").get<decltype(o.basis)>();
    if (s.contains("theory")) o.theory = s.at("theory").get<decltype(o.theory)>();
    else if (s.contains("nwchem_theory")) o.theory = s.at("nwchem_theory").get<decltype(o.theory)>();
    if (s.contains("scf_type")) o.scf_type = s.at("scf_type").get<decltype(o.scf_type)>();
    else if (s.contains("nwchem_scf_type")) o.scf_type = s.at("nwchem_scf_type").get<decltype(o.scf_type)>();
    if (s.contains("functional")) o.functional = s.at("functional").get<decltype(o.functional)>();
    else if (s.contains("cpmd_functional")) o.functional = s.at("cpmd_functional").get<decltype(o.functional)>();
    if (s.contains("cutOffRy")) o.cutoff_ry = s.at("cutOffRy").get<decltype(o.cutoff_ry)>();
    else if (s.contains("cutoff_ry")) o.cutoff_ry = s.at("cutoff_ry").get<decltype(o.cutoff_ry)>();
    else if (s.contains("cpmd_cut_off_ry")) o.cutoff_ry = s.at("cpmd_cut_off_ry").get<decltype(o.cutoff_ry)>();
    if (s.contains("charge")) o.charge = s.at("charge").get<decltype(o.charge)>();
    else if (s.contains("nwchem_charge")) o.charge = s.at("nwchem_charge").get<decltype(o.charge)>();
    if (s.contains("multiplicity")) o.multiplicity = s.at("multiplicity").get<decltype(o.multiplicity)>();
    else if (s.contains("nwchem_multiplicity")) o.multiplicity = s.at("nwchem_multiplicity").get<decltype(o.multiplicity)>();
    if (s.contains("engine_path")) o.engine_path = s.at("engine_path").get<decltype(o.engine_path)>();
    if (s.contains("engine_library")) o.engine_library = s.at("engine_library").get<decltype(o.engine_library)>();
    if (s.contains("engine_root")) o.engine_root = s.at("engine_root").get<decltype(o.engine_root)>();
    if (s.contains("title")) o.title = s.at("title").get<decltype(o.title)>();
    if (s.contains("memory_mb")) o.memory_mb = s.at("memory_mb").get<decltype(o.memory_mb)>();
    if (s.contains("scratch_dir")) o.scratch_dir = s.at("scratch_dir").get<decltype(o.scratch_dir)>();
    if (s.contains("input_block")) o.input_block = s.at("input_block").get<decltype(o.input_block)>();
    if (s.contains("permanent_dir")) o.permanent_dir = s.at("permanent_dir").get<decltype(o.permanent_dir)>();
    if (s.contains("params_path")) o.params_path = s.at("params_path").get<decltype(o.params_path)>();
    if (s.contains("ranks_per_image")) o.ranks_per_image = s.at("ranks_per_image").get<decltype(o.ranks_per_image)>();
    if (s.contains("model_path")) o.model_path = s.at("model_path").get<decltype(o.model_path)>();
    if (s.contains("device")) o.device = s.at("device").get<decltype(o.device)>();
    if (s.contains("length_unit")) o.length_unit = s.at("length_unit").get<decltype(o.length_unit)>();
    if (s.contains("extensions_directory")) o.extensions_directory = s.at("extensions_directory").get<decltype(o.extensions_directory)>();
    if (s.contains("check_consistency")) o.check_consistency = s.at("check_consistency").get<decltype(o.check_consistency)>();
    if (s.contains("uncertainty_threshold")) o.uncertainty_threshold = s.at("uncertainty_threshold").get<decltype(o.uncertainty_threshold)>();
    if (s.contains("torch_determinism_strict")) o.torch_determinism_strict = s.at("torch_determinism_strict").get<decltype(o.torch_determinism_strict)>();
    if (s.contains("paramset")) o.xtb_paramset = s.at("paramset").get<decltype(o.xtb_paramset)>();
    else if (s.contains("xtb_paramset")) o.xtb_paramset = s.at("xtb_paramset").get<decltype(o.xtb_paramset)>();
    if (s.contains("accuracy")) o.xtb_accuracy = s.at("accuracy").get<decltype(o.xtb_accuracy)>();
    else if (s.contains("xtb_accuracy")) o.xtb_accuracy = s.at("xtb_accuracy").get<decltype(o.xtb_accuracy)>();
    if (s.contains("electronic_temperature")) o.xtb_electronic_temperature = s.at("electronic_temperature").get<decltype(o.xtb_electronic_temperature)>();
    else if (s.contains("xtb_electronic_temperature")) o.xtb_electronic_temperature = s.at("xtb_electronic_temperature").get<decltype(o.xtb_electronic_temperature)>();
    if (s.contains("max_iterations")) o.xtb_max_iterations = s.at("max_iterations").get<decltype(o.xtb_max_iterations)>();
    else if (s.contains("xtb_max_iterations")) o.xtb_max_iterations = s.at("xtb_max_iterations").get<decltype(o.xtb_max_iterations)>();
    if (s.contains("uhf")) o.xtb_uhf = s.at("uhf").get<decltype(o.xtb_uhf)>();
    else if (s.contains("xtb_uhf")) o.xtb_uhf = s.at("xtb_uhf").get<decltype(o.xtb_uhf)>();
    if (s.contains("task_name")) o.task_name = s.at("task_name").get<decltype(o.task_name)>();
    if (s.contains("xtb_charge")) o.xtb_charge = s.at("xtb_charge").get<decltype(o.xtb_charge)>();
    else o.xtb_charge = static_cast<double>(o.charge);
  }
}
