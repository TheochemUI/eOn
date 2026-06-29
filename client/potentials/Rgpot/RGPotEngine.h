#pragma once

#include <memory>
#include <string>

struct RGPotEngineOptions {
  std::string backend;
  std::string basis{"sto-3g"};
  std::string theory{"scf"};
  std::string scf_type{"rhf"};
  std::string functional{"BLYP"};
  double cutoff_ry{70.0};
  int charge{0};
  int multiplicity{1};
  std::string engine_path;
  std::string engine_library;
  std::string engine_root;
  std::string title;
  int memory_mb{0};
  std::string scratch_dir;
  std::string input_block; // optional NWChem inputBlocks text
};

/** Opaque rgpot-backed engine (nwchemc / cpmdc). No eOn Potential.h here. */
class RGPotEngine {
public:
  explicit RGPotEngine(const RGPotEngineOptions &opt);
  ~RGPotEngine();
  RGPotEngine(const RGPotEngine &) = delete;
  RGPotEngine &operator=(const RGPotEngine &) = delete;

  [[nodiscard]] const std::string &backend() const noexcept { return backend_; }
  [[nodiscard]] bool available() const;
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, const double *box) const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  std::string backend_;
};
