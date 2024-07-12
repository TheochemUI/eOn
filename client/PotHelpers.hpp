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

#include "C_Structs.h"

namespace eonc::pot {
// Typically this is done by the caller, however it is here as a sanity check
void zeroForceOut(const size_t &nAtoms, ForceOut *efvd);

template <typename T> class registry {
public:
  static size_t count;
  static size_t forceCalls;
  static T *head;
  T *prev;
  T *next;

protected:
  registry() {
    ++count;
    prev = nullptr;
    next = head;
    head = static_cast<T *>(this);
    if (next) {
      next->prev = head;
    }
  }

  registry(const registry &) {
    ++count;
    prev = nullptr;
    next = head;
    head = static_cast<T *>(this);
    if (next) {
      next->prev = head;
    }
  }

  ~registry() {
    --count;
    if (prev) {
      prev->next = next;
    }
    if (next) {
      next->prev = prev;
    }
    if (head == this) {
      head = next;
    }
  }

public:
  static void incrementForceCalls() { ++forceCalls; }
};

template <typename T> size_t registry<T>::count = 0;
template <typename T> size_t registry<T>::forceCalls = 0;
template <typename T> T *registry<T>::head = nullptr;

} // namespace eonc::pot
