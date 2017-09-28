// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// @file   runSim.cxx
/// @author Sebastian Klewin
///

#include <boost/program_options.hpp>
#include <iostream>
#include <sys/wait.h>

#include "../../../../macro/run_sim_tpc.C"
#include "../../../../macro/run_digi_tpc.C"
#include "../../../../macro/run_clus_tpc.C"

namespace bpo = boost::program_options;

int main(int argc, char *argv[])
{
  // Arguments parsing
  bpo::variables_map vm; 
  bpo::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Produce help message.")
    ("mode,m",      bpo::value<std::string>()->default_value("sim"),    R"(mode of processing, "sim", "digi", "clus" or "all".)")
    ("nEvents,n",   bpo::value<int>()->default_value(2),                "number of events to simulate.")
    ("mcEngine,e",  bpo::value<std::string>()->default_value("TGeant3"), "MC generator to be used.")
    ("continuous,c", bpo::value<int>()->default_value(1),                "Running in continuous mode 1 - Triggered mode 0");
  bpo::store(parse_command_line(argc, argv, desc), vm);
  bpo::notify(vm);

  // help
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return EXIT_SUCCESS;
  }

  // Actual "work"
  const int events = vm["nEvents"].as<int>();
  const std::string engine = vm["mcEngine"].as<std::string>();
  const std::string mode = vm["mode"].as<std::string>();
  const int isContinuous = vm["continuous"].as<int>();

  std::cout << "####" << std::endl;
  std::cout << "#### Starting TPC simulation tool for" << std::endl;
  std::cout << "#### " << events << " events and " << engine << " as MC engine" << std::endl;
  std::cout << "####" << std::endl;
  std::cout << std::endl << std::endl;


  if (mode == "sim") {
    runSimTpc(events,engine);
  } else if (mode == "digi") {
    runDigiTpc(events,engine, isContinuous);
  } else if (mode == "clus") {
    runClusTpc(events,engine, isContinuous);
  } else if (mode == "all") {
    int status;
    pid_t PID = fork();
    if (PID == -1) { std::cout << "ERROR" << std::endl; return EXIT_FAILURE;}
    if (PID == 0)  { runSimTpc(events,engine); return EXIT_SUCCESS;}
    else waitpid(PID,&status,0);
    
    PID = fork();
    if (PID == -1) { std::cout << "ERROR" << std::endl; return EXIT_FAILURE;}
    if (PID == 0)  { runDigiTpc(events,engine,isContinuous); return EXIT_SUCCESS;}
    else waitpid(PID,&status,0);

    PID = fork();
    if (PID == -1) { std::cout << "ERROR" << std::endl; return EXIT_FAILURE;}
    if (PID == 0)  { runClusTpc(events,engine,isContinuous); return EXIT_SUCCESS;}
    else waitpid(PID,&status,0);
  } else {
      std::cout << "Mode was not recognised" << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
