// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#ifndef O2_FRAMEWORK_BASICOPS_H_
#define O2_FRAMEWORK_BASICOPS_H_
#include <array>
#include <string_view>
#include "CommonConstants/MathConstants.h"

namespace o2::framework
{
enum BasicOp : unsigned int {
  LogicalAnd, // 2-ar operations
  LogicalOr,
  Addition,
  Subtraction,
  Division,
  Multiplication,
  BitwiseAnd,
  BitwiseOr,
  BitwiseXor,
  LessThan,
  LessThanOrEqual,
  GreaterThan,
  GreaterThanOrEqual,
  Equal,
  NotEqual,
  Atan2, // 2-ar functions
  Power,
  Sqrt, // 1-ar functions
  Exp,
  Log,
  Log10,
  Sin,
  Cos,
  Tan,
  Asin,
  Acos,
  Atan,
  Abs,
  Round,
  BitwiseNot,
  Conditional // 3-ar functions
};

static constexpr std::array<std::string_view, BasicOp::Conditional + 1> mapping{
  "&&",
  "||",
  "+",
  "-",
  "/",
  "*",
  "&",
  "|",
  "^",
  "<",
  "<=",
  ">",
  ">=",
  "==",
  "!=",
  "natan2",
  "npow",
  "nsqrt",
  "nexp",
  "nlog",
  "nlog10",
  "nsin",
  "ncos",
  "ntan",
  "nasin",
  "nacos",
  "natan",
  "nabs",
  "nround",
  "nbitwise_not",
  "ifnode"};

static constexpr std::array<std::string_view, 9> mathConstants{
  "Almost0",
  "Epsilon",
  "Almost1",
  "VeryBig",
  "PI",
  "TwoPI",
  "PIHalf",
  "PIThird",
  "PIQuarter"};

static constexpr std::array<float, 9> mathConstantsValues{
  o2::constants::math::Almost0,
  o2::constants::math::Epsilon,
  o2::constants::math::Almost1,
  o2::constants::math::VeryBig,
  o2::constants::math::PI,
  o2::constants::math::TwoPI,
  o2::constants::math::PIHalf,
  o2::constants::math::PIThird,
  o2::constants::math::PIQuarter};
} // namespace o2::framework

#endif // O2_FRAMEWORK_BASICOPS_H_
