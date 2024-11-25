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

#define BOOST_TEST_MODULE Test GeneratorPythia8Param class
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <CommonUtils/ConfigurableParam.h>
#include <Generators/GeneratorPythia8Param.h>
#include <boost/property_tree/ptree.hpp>
#include "CCDB/BasicCCDBManager.h"

// Tests various aspects of the
// ConfigurableParamPromoter class, which is used to promote
// Pythia8GenConfig to a configurable param
BOOST_AUTO_TEST_CASE(pythia8_Pythia8GenConfig)
{
  o2::conf::ConfigurableParam::updateFromString(
    "GeneratorPythia8.config=Foo;GeneratorPythia8.includePartonEvent=true");

  using o2::eventgen::GeneratorPythia8Param;

  BOOST_CHECK_EQUAL(GeneratorPythia8Param::Instance().config, std::string("Foo"));
  BOOST_CHECK_EQUAL(GeneratorPythia8Param::Instance().includePartonEvent, true);

  BOOST_CHECK_EQUAL(GeneratorPythia8Param::Instance().includePartonEvent, o2::conf::ConfigurableParam::getValueAs<bool>("GeneratorPythia8.includePartonEvent"));
  // setValue - getValue
  o2::conf::ConfigurableParam::setValue("GeneratorPythia8.config", "Baz");
  BOOST_CHECK_EQUAL(o2::conf::ConfigurableParam::getValueAs<std::string>("GeneratorPythia8.config"), std::string("Baz"));
  BOOST_CHECK_EQUAL(GeneratorPythia8Param::Instance().config, std::string("Baz"));

  // member provenance
  BOOST_CHECK_EQUAL(GeneratorPythia8Param::Instance().getMemberProvenance("config"), o2::conf::ConfigurableParam::EParamProvenance::kRT);
  BOOST_CHECK_EQUAL(GeneratorPythia8Param::Instance().getMemberProvenance("verbose"), o2::conf::ConfigurableParam::EParamProvenance::kCODE);

  // config detach
  auto config_copy = GeneratorPythia8Param::Instance().detach();
  BOOST_CHECK_EQUAL(config_copy.config, std::string("Baz"));
  BOOST_CHECK_EQUAL(config_copy.includePartonEvent, true);

  // file IO
  TFile tmp_file("GeneratorParamConfig_tmp.root", "RECREATE");

  GeneratorPythia8Param::Instance().serializeTo(&tmp_file);
  // modify the instance to some intermediate fluent value
  o2::conf::ConfigurableParam::setValue("GeneratorPythia8.includePartonEvent", "0");
  BOOST_CHECK_EQUAL(config_copy.includePartonEvent, true);
  BOOST_CHECK_EQUAL(GeneratorPythia8Param::Instance().includePartonEvent, false);
  tmp_file.Close();

  // read back
  TFile tmp_file2("GeneratorParamConfig_tmp.root", "READ");
  const_cast<GeneratorPythia8Param&>(GeneratorPythia8Param::Instance()).initFrom(&tmp_file2);
  BOOST_CHECK_EQUAL(GeneratorPythia8Param::Instance().includePartonEvent, true);
  tmp_file2.Close();

  // CCDB IO
  std::string ccdbUrl = "http://ccdb-test.cern.ch:8080";
  bool hostReachable = false;
  o2::ccdb::CcdbApi api;
  api.init(ccdbUrl);
  std::string pathA = "/Generators/UnitTest/Pythia8/GeneratorPythia8Param";
  std::map<std::string, std::string> md;
  long start = 1000, stop = 2000;
  api.storeAsTFileAny(&GeneratorPythia8Param::Instance(), pathA, md, start, stop);

  // modify the instance to some intermediate fluent value
  o2::conf::ConfigurableParam::setValue("GeneratorPythia8.includePartonEvent", "0");

  auto returnedobj = api.retrieveFromTFileAny<o2::eventgen::GeneratorPythia8Param>(pathA, md, (start + stop) / 2);
  GeneratorPythia8Param::Instance().printKeyValues();
};
