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
#include <string>
#include <string_view>
#include <unordered_map>

namespace eonc {
// clang-format off
enum class Element {
  Unknown = 0, H = 1, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl,
  Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb,
  Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La,
  Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir,
  Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, COUNTER};
// clang-format on

struct ElementProperties {
  const size_t atomicNumber;
  const double atomicMass;
};

const std::unordered_map<Element, ElementProperties> elementData = {
    {Element::Unknown, {0, 0.0}}, {Element::H, {1, 1.008}},
    {Element::He, {2, 4.0026}},   {Element::Li, {3, 6.94}},
    {Element::Be, {4, 9.0122}},   {Element::B, {5, 10.81}},
    {Element::C, {6, 12.011}},    {Element::N, {7, 14.007}},
    {Element::O, {8, 15.999}},    {Element::F, {9, 18.998}},
    {Element::Ne, {10, 20.180}},  {Element::Na, {11, 22.990}},
    {Element::Mg, {12, 24.305}},  {Element::Al, {13, 26.982}},
    {Element::Si, {14, 28.085}},  {Element::P, {15, 30.974}},
    {Element::S, {16, 32.06}},    {Element::Cl, {17, 35.45}},
    {Element::Ar, {18, 39.948}},  {Element::K, {19, 39.098}},
    {Element::Ca, {20, 40.078}},  {Element::Sc, {21, 44.956}},
    {Element::Ti, {22, 47.867}},  {Element::V, {23, 50.942}},
    {Element::Cr, {24, 51.996}},  {Element::Mn, {25, 54.938}},
    {Element::Fe, {26, 55.845}},  {Element::Co, {27, 58.933}},
    {Element::Ni, {28, 58.693}},  {Element::Cu, {29, 63.546}},
    {Element::Zn, {30, 65.38}},   {Element::Ga, {31, 69.723}},
    {Element::Ge, {32, 72.630}},  {Element::As, {33, 74.922}},
    {Element::Se, {34, 78.971}},  {Element::Br, {35, 79.904}},
    {Element::Kr, {36, 83.798}},  {Element::Rb, {37, 85.468}},
    {Element::Sr, {38, 87.62}},   {Element::Y, {39, 88.906}},
    {Element::Zr, {40, 91.224}},  {Element::Nb, {41, 92.906}},
    {Element::Mo, {42, 95.95}},   {Element::Tc, {43, 98}},
    {Element::Ru, {44, 101.07}},  {Element::Rh, {45, 102.91}},
    {Element::Pd, {46, 106.42}},  {Element::Ag, {47, 107.87}},
    {Element::Cd, {48, 112.41}},  {Element::In, {49, 114.82}},
    {Element::Sn, {50, 118.71}},  {Element::Sb, {51, 121.76}},
    {Element::Te, {52, 127.60}},  {Element::I, {53, 126.90}},
    {Element::Xe, {54, 131.29}},  {Element::Cs, {55, 132.91}},
    {Element::Ba, {56, 137.33}},  {Element::La, {57, 138.91}},
    {Element::Ce, {58, 140.12}},  {Element::Pr, {59, 140.91}},
    {Element::Nd, {60, 144.24}},  {Element::Pm, {61, 145}},
    {Element::Sm, {62, 150.36}},  {Element::Eu, {63, 151.96}},
    {Element::Gd, {64, 157.25}},  {Element::Tb, {65, 158.93}},
    {Element::Dy, {66, 162.50}},  {Element::Ho, {67, 164.93}},
    {Element::Er, {68, 167.26}},  {Element::Tm, {69, 168.93}},
    {Element::Yb, {70, 173.05}},  {Element::Lu, {71, 174.97}},
    {Element::Hf, {72, 178.49}},  {Element::Ta, {73, 180.95}},
    {Element::W, {74, 183.84}},   {Element::Re, {75, 186.21}},
    {Element::Os, {76, 190.23}},  {Element::Ir, {77, 192.22}},
    {Element::Pt, {78, 195.08}},  {Element::Au, {79, 196.97}},
    {Element::Hg, {80, 200.59}},  {Element::Tl, {81, 204.38}},
    {Element::Pb, {82, 207.2}},   {Element::Bi, {83, 208.98}},
    {Element::Po, {84, 209}},     {Element::At, {85, 210}},
    {Element::Rn, {86, 222}},     {Element::Fr, {87, 223}},
    {Element::Ra, {88, 226}},     {Element::Ac, {89, 227}},
    {Element::Th, {90, 232.04}},  {Element::Pa, {91, 231.04}},
    {Element::U, {92, 238.03}}};

std::string mass2atom(double atomicmass);
size_t symbol2atomicNumber(const std::string_view &symbol);
std::string atomicNumber2symbol(size_t n);

} // namespace eonc
