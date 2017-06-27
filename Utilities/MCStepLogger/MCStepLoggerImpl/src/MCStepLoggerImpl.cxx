#include <TVirtualMCApplication.h>
#include <TVirtualMC.h>
#include <dlfcn.h>
#include <iostream>
#include <map>
#include <set>

typedef void (TVirtualMCApplication::*stepptr)(void);

// function pointers, pointing to original symbols
// In this case, we are retrieving the Symbol from libBase which
// contains FairMCApplication
static void* libHandle = dlopen("libBase.so", RTLD_NOW);

class StepLogger
{
  int stepcounter = 0;

  std::set<int> trackset;
  std::set<int> pdgset;
  std::map<int, int> volumetosteps;
  std::map<int, char const*> idtovolname;

 public:
  void AddStep(TVirtualMC* mc)
  {
    stepcounter++;
    trackset.insert(mc->GetStack()->GetCurrentTrackNumber());
    pdgset.insert(mc->TrackPid());
    int copyNo;
    auto id = mc->CurrentVolID(copyNo);
    volumetosteps[id]++;
    if (idtovolname.find(id) == idtovolname.end()) {
      idtovolname.insert(std::pair<int, char const*>(id, mc->CurrentVolName()));
    }
  }

  ~StepLogger()
  {
    std::cerr << "did " << stepcounter << " steps \n";
    std::cerr << "transported " << trackset.size() << " different tracks \n";
    std::cerr << "transported " << pdgset.size() << " different types \n";
    // summarize steps per volume
    for (auto& p : volumetosteps) {
      std::cout << " VolName " << idtovolname[p.first] << " COUNT " << p.second << "\n";
    }
  }
};

// will be deconstructed at the end
static StepLogger logger;

extern "C" void dispatchOriginalStepping(TVirtualMCApplication* app)
{
  typedef void (TVirtualMCApplication::*StepMethodType)(void);
  static StepMethodType origMethod = nullptr;
  if (origMethod == nullptr) {
    void* symbolAddress = dlsym(libHandle, "_ZN17FairMCApplication8SteppingEv");

    // hack since C++ does not allow casting to C++ member function pointers
    // thanks to gist.github.com/mooware/1174572
    memcpy(&origMethod, &symbolAddress, sizeof(&symbolAddress));
  }
  (app->*origMethod)();
}

extern "C" void performLogging(TVirtualMCApplication* app)
{
  static TVirtualMC* mc = TVirtualMC::GetMC();
  logger.AddStep(mc);
}
