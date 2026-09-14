#pragma once
/*
** Thin RGPOT-metatomic backend: dlopen libmetatomic_engine.so (no torch).
*/
#include "metatomic_c_abi.h"
#include <memory>
#include <string>

struct MetatomicEngineOptions {
  std::string model_path;
  std::string device{"cpu"};
  std::string length_unit{"angstrom"};
  std::string extensions_directory;
  bool check_consistency{false};
  double uncertainty_threshold{-1.0};
  std::string engine_path; // optional explicit .so path
  bool torch_determinism_strict{false};
};

class MetatomicEngineLoader {
public:
  explicit MetatomicEngineLoader(const MetatomicEngineOptions &opt);
  ~MetatomicEngineLoader();
  MetatomicEngineLoader(const MetatomicEngineLoader &) = delete;
  MetatomicEngineLoader &operator=(const MetatomicEngineLoader &) = delete;

  [[nodiscard]] bool available() const noexcept { return m_pot != nullptr; }
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) const;

private:
  void *m_lib{nullptr};
  RgpotMtaPot *m_pot{nullptr};
  using create_fn = RgpotMtaPot *(*)(const RgpotMtaConfig *, char *, size_t);
  using destroy_fn = void (*)(RgpotMtaPot *);
  using force_fn = int (*)(RgpotMtaPot *, long, const double *, const int *,
                           double *, double *, double *, const double *);
  create_fn m_create{nullptr};
  destroy_fn m_destroy{nullptr};
  force_fn m_force{nullptr};
};
