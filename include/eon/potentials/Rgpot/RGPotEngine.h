#pragma once

#include <cstddef>
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
  std::string input_block;   // raw input text for NWChem or CPMD inputBlocks
  std::string permanent_dir; // CPMD FILEPATH for RESTART files (cpmdc)
  std::string params_path;   // optional CPMDParams message file (cpmdc)
  int ranks_per_image{0};    // cpmdc: ranks per calculator group, 0 = off
  // Metatomic (backend=metatomic): dlopen libmetatomic_engine.so
  std::string model_path;
  std::string device{"cpu"};
  std::string length_unit{"angstrom"};
  std::string extensions_directory;
  bool check_consistency{false};
  double uncertainty_threshold{-1.0};
  bool torch_determinism_strict{false};
  // XTB (backend=xtb): dlopen libxtb_engine.so — not native -Dwith_xtb link
  std::string xtb_paramset{"GFN2xTB"}; // GFNFF / GFN0xTB / GFN1xTB / GFN2xTB
  double xtb_accuracy{1.0};
  double xtb_electronic_temperature{300.0};
  int xtb_max_iterations{250};
  double xtb_charge{0.0};
  int xtb_uhf{0};
};

/** Opaque rgpot-backed engine (nwchemc / cpmdc / metatomic / xtb). */
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

  /// Number of calculator groups the MPI world is split into (1 when
  /// ranks_per_image is off) and the group this rank belongs to.
  [[nodiscard]] int calculatorGroups() const noexcept;
  [[nodiscard]] int calculatorIndex() const noexcept;
  /// Ranks in the MPI world the groups were bound on (1 without MPI).
  [[nodiscard]] int calculatorWorld() const noexcept;
  /// This process's rank in that world (0 without MPI).
  [[nodiscard]] int worldRank() const noexcept;
  /// Collective on MPI_COMM_WORLD: every rank leaves with world rank 0's
  /// bytes. A no-op on one rank.
  void broadcastFromDriver(void *data, std::size_t bytes) const;
  /// Registers MPI_Finalize at exit (once per process) when the world has
  /// more than one rank.
  void finalizeMpiAtExit() const;
  /// Collective on MPI_COMM_WORLD: every rank leaves with the energy and
  /// the 3N forces computed by the first rank of group `owner`. Only
  /// CPMD's parent rank holds the true forces, so this also keeps the
  /// ranks of one group on the same step.
  void shareResult(int owner, long N, double *F, double *U) const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  std::string backend_;
};
