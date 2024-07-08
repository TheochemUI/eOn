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
#include <vector>

int getBundleSize(void);
std::vector<std::string> unbundle(int number);
void bundle(int number, const std::vector<std::string> &filenames,
            std::vector<std::string> *bundledFilenames);
void deleteUnbundledFiles(const std::vector<std::string> &unbundledFilenames);
