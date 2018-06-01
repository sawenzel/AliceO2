#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <fcntl.h>
#include <SimConfig/SimConfig.h>
#include <sys/wait.h>
#include <vector>

const char* serverlogname = "serverlog";
const char* workerlogname = "workerlog";
const char* mergerlogname = "mergerlog";

// helper exec to launch the devices
int main(int argc, char* argv[]) {
  std::string binpath(getenv("O2_ROOT"));
  binpath += "/bin";
  auto home = getenv("HOME");

  std::stringstream configss;
  configss << home << "/alisw_new/O2/run/primary-server.json";

  auto& conf = o2::conf::SimConfig::Instance();
  if (!conf.resetFromArguments(argc, argv)) {
    return 1;
  }

  std::vector<int> childpids;

  // the server
  int pid = fork();
  if (pid == 0) {
    int fd = open(serverlogname, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
    dup2(fd, 1); // make stdout go to file
    dup2(fd, 2); // make stderr go to file - you may choose to not do this
                 // or perhaps send stderr to another file
    close(fd); // fd no longer needed - the dup'ed handles are sufficient

    const std::string name("O2PrimaryServerDeviceRunner");
    const std::string path = binpath + "/" + name;

    // copy all arguments into a common vector
    const int Nargs = argc + 6;
    const char* arguments[Nargs];
    arguments[0] = name.c_str();
    arguments[1] = "--control"; arguments[2] = "static";
    arguments[3] = "--id"; arguments[4] = "primary-server";
    arguments[5] = "--mq-config";
    arguments[6] = configss.str().c_str();
    for (int i = 1; i < argc; ++i) {
      arguments[6 + i] = argv[i];
    }
    arguments[argc + 6] = nullptr;
    for (int i = 0; i < Nargs; ++i) {
      std::cerr << arguments[i] << "\n";
    }
    std::cerr << "$$$$";

    //execl(path.c_str(), name.c_str(), "--control", "static", "--id", "primary-server", "--mq-config",
    //      configss.str().c_str(), argv, (char*)0);
    execv(path.c_str(), (char *const *) arguments);
    return 0;
  }
  else {
	childpids.push_back(pid);
    std::cout << "Spawning particle server on PID " << pid << "; Redirect output to " << serverlogname << "\n";
  }

  auto f = getenv("ALICE_NSIMWORKERS");
  int nworkers = 1;
  if (f) {
    nworkers = atoi(f);
  }
  for (int id = 0; id < nworkers; ++id) {
    // the workers
    std::stringstream workerlogss;
    workerlogss << workerlogname << id;

	pid = fork();
    if (pid == 0) {
      int fd = open(workerlogss.str().c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
      dup2(fd, 1); // make stdout go to file
      dup2(fd, 2); // make stderr go to file - you may choose to not do this
                   // or perhaps send stderr to another file
      close(fd);   // fd no longer needed - the dup'ed handles are sufficient

      const std::string name("O2SimDeviceRunner");
      const std::string path = binpath + "/" + name;
      execl(path.c_str(), name.c_str(), "--control", "static", "--id", "worker", "--mq-config", configss.str().c_str(),
            (char*)0);
      return 0;
    } else {
      childpids.push_back(pid);
      std::cout << "Spawning sim worker " << id << " on PID " << pid << "; Redirect output to " << workerlogss.str() << "\n";
    }
  }

  // the hit merger
  int status, cpid;
  pid = fork();
  if (pid == 0) {
    int fd = open(mergerlogname, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
    dup2(fd, 1); // make stdout go to file
    dup2(fd, 2); // make stderr go to file - you may choose to not do this
                 // or perhaps send stderr to another file
    close(fd);   // fd no longer needed - the dup'ed handles are sufficient

    const std::string name("O2HitMergerRunner");
    const std::string path = binpath + "/" + name;
    execl(path.c_str(), name.c_str(), "--control", "static", "--id", "hitmerger", "--mq-config", configss.str().c_str(),
          (char*)0);
    return 0;
  } else {
    std::cout << "Spawning hit merger on PID " << pid << "; Redirect output to " << mergerlogname << "\n";
	childpids.push_back(pid);
  }

  // wait on all children
  for (auto p : childpids) {
    if ((cpid = wait(&status)) == pid) {
      printf("Process %d returned\n", pid);
    }
  }

  return 0;
}
