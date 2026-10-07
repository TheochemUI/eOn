// AUTO-GENERATED from schema/eon_params.capnp — do not edit.
// Regenerate: python tools/params_ssot/codegen.py
// Included inside namespace eonc::config, after toLowerCase.
inline void project_ssot_ini(const INIReader &ini, Parameters &params) {
  if (params.potential_options().potential == PotType::RGPOT) {
    auto &o = ParametersLoadAccess::rgpot_options(params);
    const char *sec = "RgpotPot";
    o.backend = ini.Get(sec, "backend", o.backend);
    o.basis = ini.Get(sec, "basis", ini.Get(sec, "nwchem_basis", o.basis));
    o.theory = ini.Get(sec, "theory", ini.Get(sec, "nwchem_theory", o.theory));
    o.scf_type = ini.Get(sec, "scf_type", ini.Get(sec, "nwchem_scf_type", o.scf_type));
    o.functional = ini.Get(sec, "functional", ini.Get(sec, "cpmd_functional", o.functional));
    o.cutoff_ry = ini.GetReal(sec, "cutOffRy", ini.GetReal(sec, "cutoff_ry", ini.GetReal(sec, "cpmd_cut_off_ry", o.cutoff_ry)));
    o.charge = ini.GetInteger(sec, "charge", ini.GetInteger(sec, "nwchem_charge", o.charge));
    o.multiplicity = ini.GetInteger(sec, "multiplicity", ini.GetInteger(sec, "nwchem_multiplicity", o.multiplicity));
    o.engine_path = ini.Get(sec, "engine_path", o.engine_path);
    o.engine_library = ini.Get(sec, "engine_library", o.engine_library);
    o.engine_root = ini.Get(sec, "engine_root", o.engine_root);
    o.title = ini.Get(sec, "title", o.title);
    o.memory_mb = ini.GetInteger(sec, "memory_mb", o.memory_mb);
    o.scratch_dir = ini.Get(sec, "scratch_dir", o.scratch_dir);
    o.input_block = ini.Get(sec, "input_block", o.input_block);
    o.permanent_dir = ini.Get(sec, "permanent_dir", o.permanent_dir);
    o.params_path = ini.Get(sec, "params_path", o.params_path);
    o.ranks_per_image = ini.GetInteger(sec, "ranks_per_image", o.ranks_per_image);
    o.model_path = ini.Get(sec, "model_path", o.model_path);
    o.device = ini.Get(sec, "device", o.device);
    o.length_unit = ini.Get(sec, "length_unit", o.length_unit);
    o.extensions_directory = ini.Get(sec, "extensions_directory", o.extensions_directory);
    o.check_consistency = ini.GetBoolean(sec, "check_consistency", o.check_consistency);
    o.uncertainty_threshold = ini.GetReal(sec, "uncertainty_threshold", o.uncertainty_threshold);
    o.torch_determinism_strict = ini.GetBoolean(sec, "torch_determinism_strict", o.torch_determinism_strict);
    o.xtb_paramset = ini.Get(sec, "paramset", ini.Get(sec, "xtb_paramset", o.xtb_paramset));
    o.xtb_accuracy = ini.GetReal(sec, "accuracy", ini.GetReal(sec, "xtb_accuracy", o.xtb_accuracy));
    o.xtb_electronic_temperature = ini.GetReal(sec, "electronic_temperature", ini.GetReal(sec, "xtb_electronic_temperature", o.xtb_electronic_temperature));
    o.xtb_max_iterations = static_cast<int>(ini.GetInteger(sec, "max_iterations", ini.GetInteger(sec, "xtb_max_iterations", o.xtb_max_iterations)));
    o.xtb_uhf = static_cast<int>(ini.GetInteger(sec, "uhf", ini.GetInteger(sec, "xtb_uhf", o.xtb_uhf)));
    o.task_name = ini.Get(sec, "task_name", o.task_name);
    o.xtb_charge = ini.GetReal(sec, "xtb_charge", static_cast<double>(o.charge));
    const std::string be = toLowerCase(o.backend);
    if (be == "xtb" || be == "xtbpot" || be == "gfn" || be == "gfnxtb") {
      o.xtb_paramset = ini.Get("XTBPot", "paramset", o.xtb_paramset);
      o.xtb_accuracy = ini.GetReal("XTBPot", "accuracy", o.xtb_accuracy);
      o.xtb_electronic_temperature = ini.GetReal("XTBPot", "electronic_temperature", o.xtb_electronic_temperature);
      o.xtb_max_iterations = static_cast<int>(ini.GetInteger("XTBPot", "max_iterations", o.xtb_max_iterations));
      o.xtb_uhf = static_cast<int>(ini.GetInteger("XTBPot", "uhf", o.xtb_uhf));
      o.xtb_charge = ini.GetReal("XTBPot", "charge", o.xtb_charge);
    }
    if (be == "cpmd" || be == "cpmdc" || be == "cpmdpot") {
      o.functional = ini.Get("cpmd", "functional", ini.Get("cpmd", "cpmd_functional", o.functional));
      o.cutoff_ry = ini.GetReal("cpmd", "cutOffRy", ini.GetReal("cpmd", "cutoff_ry", ini.GetReal("cpmd", "cpmd_cut_off_ry", o.cutoff_ry)));
      o.charge = ini.GetInteger("cpmd", "charge", o.charge);
      o.multiplicity = ini.GetInteger("cpmd", "multiplicity", o.multiplicity);
      o.title = ini.Get("cpmd", "title", o.title);
      o.memory_mb = ini.GetInteger("cpmd", "memory_mb", o.memory_mb);
      o.input_block = ini.Get("cpmd", "input_block", o.input_block);
    }
  }
}
