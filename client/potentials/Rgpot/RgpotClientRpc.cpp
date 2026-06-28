/*
** Cap'n Proto client to rgpot potserv. Isolated TU: defines Potential from schema.
*/
#include "RgpotClientRpc.h"

#include <capnp/ez-rpc.h>
#include <kj/debug.h>

#include <cstring>
#include <memory>
#include <stdexcept>

#include "Potentials.capnp.h"

namespace {

struct Holder {
  std::unique_ptr<capnp::EzRpcClient> client;
  Potential::Client pot{nullptr};
};

} // namespace

void *rgpotClientCreate(const std::string &host, int port) {
  auto *h = new Holder;
  h->client = std::make_unique<capnp::EzRpcClient>(host, port);
  h->pot = h->client->getMain<Potential>();
  return h;
}

void rgpotClientDestroy(void *holder) {
  delete static_cast<Holder *>(holder);
}

bool rgpotClientConfigure(void *holder, const RgpotConfigureSpec &spec,
                          std::string *message_out) {
  auto *h = static_cast<Holder *>(holder);
  auto &waitScope = h->client->getWaitScope();

  auto req = h->pot.configureRequest();
  auto cfg = req.initConfig();
  if (spec.backend == "CPMD" || spec.backend == "cpmd") {
    auto cp = cfg.initCpmd();
    cp.setFunctional(spec.cpmd_functional);
    cp.setTask(spec.cpmd_task);
    cp.setCutOffRy(spec.cpmd_cut_off_ry);
    cp.setCharge(spec.cpmd_charge);
    cp.setMultiplicity(spec.cpmd_multiplicity);
  } else {
    auto nw = cfg.initNwchem();
    nw.setBasis(spec.nwchem_basis);
    nw.setTheory(spec.nwchem_theory);
    nw.setScfType(spec.nwchem_scf_type);
    nw.setCharge(spec.nwchem_charge);
    nw.setMultiplicity(spec.nwchem_multiplicity);
    if (!spec.nwchem_input_blocks.empty()) {
      auto blocks = nw.initInputBlocks(
          static_cast<unsigned>(spec.nwchem_input_blocks.size()));
      for (unsigned i = 0; i < blocks.size(); ++i)
        blocks.set(i, spec.nwchem_input_blocks[i]);
    }
  }

  auto result = req.send().wait(waitScope);
  if (message_out) {
    *message_out = result.getMessage().cStr();
  }
  return result.getOk();
}

void rgpotClientCalculate(void *holder, long n_atoms, const double *positions,
                          const int *atomic_nrs, const double *box_9,
                          double *energy_out, double *forces_out) {
  auto *h = static_cast<Holder *>(holder);
  auto &waitScope = h->client->getWaitScope();

  auto req = h->pot.calculateRequest();
  auto fip = req.initFip();
  const unsigned npos = static_cast<unsigned>(n_atoms * 3);
  auto pos = fip.initPos(npos);
  for (unsigned i = 0; i < npos; ++i) {
    pos.set(i, positions[i]);
  }
  auto atmnrs = fip.initAtmnrs(static_cast<unsigned>(n_atoms));
  for (long i = 0; i < n_atoms; ++i) {
    atmnrs.set(static_cast<unsigned>(i), atomic_nrs[i]);
  }
  auto box = fip.initBox(9);
  for (unsigned i = 0; i < 9; ++i) {
    box.set(i, box_9[i]);
  }
  fip.setLengthUnit("angstrom");
  fip.setEnergyUnit("eV");

  auto result = req.send().wait(waitScope);
  auto pres = result.getResult();
  *energy_out = pres.getEnergy();
  auto forces = pres.getForces();
  if (forces.size() != npos) {
    throw std::runtime_error("rgpot calculate: force vector size mismatch");
  }
  for (unsigned i = 0; i < npos; ++i) {
    forces_out[i] = forces[i];
  }
}
