#pragma once
/*
** Thin RGPOT-xtb backend: dlopen libxtb_engine.so (ISO_C_BINDING engine).
** Host does not link libxtb — same pattern as MetatomicEngineLoader.
*/
#include "xtb_c_abi.h"
#include <memory>
#include <string>

struct XTBEngineOptions {
  int method{RGPOT_XTB_METHOD_GFN2}; // RGPOT_XTB_METHOD_*
  double accuracy{1.0};
  double electronic_temperature{300.0};
  int max_iterations{250};
  double charge{0.0};
  int uhf{0};
  std::string engine_path;
};

class XTBEngineLoader {
public:
  explicit XTBEngineLoader(const XTBEngineOptions &opt);
  ~XTBEngineLoader();
  XTBEngineLoader(const XTBEngineLoader &) = delete;
  XTBEngineLoader &operator=(const XTBEngineLoader &) = delete;

  [[nodiscard]] bool available() const noexcept { return m_pot != nullptr; }
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) const;

private:
  void *m_lib{nullptr};
  RgpotXtbPot *m_pot{nullptr};
  using create_fn = RgpotXtbPot *(*)(const RgpotXtbConfig *, char *, size_t);
  using destroy_fn = void (*)(RgpotXtbPot *);
  using force_fn = int (*)(RgpotXtbPot *, long, const double *, const int *,
                           double *, double *, double *, const double *);
  create_fn m_create{nullptr};
  destroy_fn m_destroy{nullptr};
  force_fn m_force{nullptr};
};
