// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   TPCDigitRootFileWriterSpec.cxx
/// @author Matthias Richter, Sandro Wenzel
/// @since  2018-04-19
/// @brief  Processor spec for a ROOT file writer for TPC digits

#include "TPCDigitRootWriterSpec.h"
#include "Framework/CallbackService.h"
#include "TPCBase/Sector.h"
#include "TPCBase/Digit.h"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <memory> // for make_shared, make_unique, unique_ptr
#include <stdexcept>

using namespace o2::framework;

namespace o2
{
namespace TPC
{

template <typename T>
TBranch* getOrMakeBranch(TTree& tree, const char* basename, int sector, T* ptr)
{
  std::stringstream stream;
  stream << basename << "_" << sector;
  const auto brname = stream.str().c_str();
  if (auto br = tree.GetBranch(brname)) {
    br->SetAddress((void*)&ptr);
	return br;
  }
  // otherwise make it
  return tree.Branch(brname, ptr);
}

/// create the processor spec
/// describing a processor aggregating digits for various TPC sectors and writing them to file
/// MC truth information is also aggregated and written out
DataProcessorSpec getTPCDigitRootFileWriterSpec()
{
  auto initFunction = [](InitContext& ic) {
    // get the option from the init context
    auto filename = ic.options().get<std::string>("tpc-digit-outfile");
    auto treename = ic.options().get<std::string>("treename");

    auto outputfile = std::make_shared<TFile>(filename.c_str(), "RECREATE");
    auto outputtree = std::make_shared<TTree>(treename.c_str(), treename.c_str());

    // container for incoming digits
    auto digits = std::make_shared<std::vector<o2::TPC::Digit>>();

    // we make the branch name on the fly since we don't know which sectors arrive

    // the callback to be set as hook at stop of processing for the framework
    auto finishWriting = [outputfile, outputtree]() {
      // correctly verify and set the entry numbers of the tree
      // to be consistent with the entries in the branches

      outputtree->Write();
      outputfile->Close();
    };
    ic.services().get<CallbackService>().set(CallbackService::Id::Stop, finishWriting);

    // set up the processing function
    // using by-copy capture of the worker instance shared pointer
    // the shared pointer makes sure to clean up the instance when the processing
    // function gets out of scope
    auto processingFct = [outputfile, outputtree, digits](ProcessingContext& pc) {
      auto indata = pc.inputs().get<std::vector<o2::TPC::Digit>>("input");
      auto sectorptr = pc.inputs().get<int>("sector");
      const int sector = *sectorptr.get();

      *digits.get() = std::move(indata);
      // connect this to a particular branch
      auto br = getOrMakeBranch(*outputtree.get(), "TPCDigit", sector, digits.get());
      br->Fill();
      br->ResetAddress();

      // outputtree->Fill();
    };

    // return the actual processing function as a lambda function using variables
    // of the init function
    return processingFct;
  };

  return DataProcessorSpec{ "TPCDigitWriter",
                            { InputSpec{ "input", "TPC", "DIGITS", 0, Lifetime::Timeframe }, // digit input
                              InputSpec{ "sector", "TPC", "SECTOR", 0, Lifetime::Condition } },
                            {}, // no output
                            AlgorithmSpec(initFunction),
                            Options{
                              { "outfile", VariantType::String, "tpctracks.root", { "Name of the input file" } },
                              { "treename", VariantType::String, "Tracks", { "Name of tree for tracks" } },
                            } };
}
} // end namespace TPC
} // end namespace o2




