// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2016 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Code within namespace sha256 are originally from Stephan Brumme.
// see http://create.stephan-brumme.com/disclaimer.html

#pragma once

#include <cstdint>
#include <cstddef>
#include <string>

namespace qlat
{  //

namespace sha256
{  //

const size_t HashBytes = 32;

const size_t HashValues = HashBytes / 4;

void processBlock(uint32_t newHash[8], const uint32_t oldHash[8],
                  const uint8_t data[64]);

void processInput(uint32_t hash[8], const uint32_t oldHash[8],
                  const uint64_t numBytes, const uint8_t* input,
                  const size_t inputSize);

void setInitialHash(uint32_t hash[8]);

void computeHash(uint32_t hash[8], const void* data, const size_t size);

void rawHashFromHash(uint8_t rawHash[HashBytes],
                     const uint32_t hash[HashValues]);

std::string showRawHash(const uint8_t rawHash[HashBytes]);

std::string showHash(const uint32_t hash[8]);

}  // namespace sha256

}  // namespace qlat
