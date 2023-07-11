#include "../SurrogatePotential.h"

#ifdef WITH_CATLEARN
#include "../potentials/CatLearnPot/CatLearnPot.h"
#endif

#ifdef WITH_GPR_OPTIM
#include "potentials/GPRPotential/GPRPotential.h"
#endif

namespace helpers::create {
template <typename... Args>
std::shared_ptr<SurrogatePotential>
makeSurrogatePotential(PotType a_ptype, std::shared_ptr<Parameters> a_params,
                       Args &&...a_args) {
  switch (a_ptype) {
    // TODO: Every potential must know their own type
#ifdef WITH_CATLEARN
  case PotType::CatLearn: {
    return (std::make_shared<CatLearnPot>(a_params));
    break;
  }
#endif
#ifdef WITH_GPR_OPTIM
  case PotType::GPR_Optim: {
    return (std::make_shared<GPRPotential>(a_params,
                                           std::forward<Args>(a_args)...));
    break;
  }
#endif
  default:
    throw std::runtime_error(
        "No known surrogate potential could be constructed");
    break;
  }
}
} // namespace helpers::create
