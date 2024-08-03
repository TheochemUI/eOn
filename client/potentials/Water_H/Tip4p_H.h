/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/

#pragma once
#include "../../Potential.h"
#include "../Water/Water.hpp"

extern "C" {
/** External function implemented in Fortran. Calculate interactions between
atoms using the potH-H2O force field.
@param[in] R           Array of length of the three vectors connecting the
H atom to the individual atoms in the water molecule in Angstrom. R(0)=H-O,
R(1)=H-H1 and R(2)=H-H2
@param[out] U           Pointer to energy in kJ/mol.
@param[out] F           Array of the norms for the three forces in
kJ/mol/Angstrom. F(0)=H-O, F(1)=H-H1 and F(2)=H-H2. Forces are . Direction of
the vectors for the distances and forces are equvalent
*/
void poth2oh_(double *R, double *U, double *F);
void setup_(void);
}
namespace eonc {
/** H-Water potential.*/
class Tip4p_H : public Potential<Tip4p_H> {

public:
  Tip4p_H() {
    setup_();
    tip4p_pot = std::make_shared<Tip4p>();
  }
  void forceImpl(const ForceInput &, ForceOut *) override final;

private:
  std::shared_ptr<Tip4p> tip4p_pot;
};
} // namespace eonc
