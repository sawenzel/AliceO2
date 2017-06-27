// to be able to hook into FairMCApplication
class FairMCApplication
{
 public:
  void Stepping();
};

// to be able to hook into AliMC
class AliMC
{
 public:
  void Stepping();
};

class TVirtualMCApplication;
extern "C" void performLogging(TVirtualMCApplication*);
extern "C" void dispatchOriginalStepping(TVirtualMCApplication*);

void FairMCApplication::Stepping()
{
  performLogging(reinterpret_cast<TVirtualMCApplication*>(this));
  dispatchOriginalStepping(reinterpret_cast<TVirtualMCApplication*>(this));
}

void AliMC::Stepping()
{
  performLogging(reinterpret_cast<TVirtualMCApplication*>(this));
  dispatchOriginalStepping(reinterpret_cast<TVirtualMCApplication*>(this));
}
