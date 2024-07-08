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
// Floating Point Trapping. It is platform specific!
// This causes the program to crash on divison by zero,
// invalid operations, and overflows.
void enableFPE(void);
