#include "../SurrogatePotential.h"

#ifdef WITH_CATLEARN
#include "../potentials/CatLearnPot/CatLearnPot.h"
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
  default:
    throw std::runtime_error(
        "No known surrogate potential could be constructed");
    break;
  }
}
} // namespace helpers::create
