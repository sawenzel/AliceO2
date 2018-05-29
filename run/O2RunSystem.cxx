#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <fcntl.h>
#include <sys/wait.h>

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
    execl(path.c_str(), name.c_str(), "--control", "static", "--id", "primary-server", "--mq-config",
          configss.str().c_str(), (char*)0);
    return 0;
  }
  else {
    std::cout << "Spawning particle server on PID " << pid << "; Redirect output to " << serverlogname << "\n";
  }

  int nworkers = 1;
  if (argc > 1) {
    nworkers = atoi(argv[1]);
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
    if ((cpid = wait(&status)) == pid) {
      printf("Child %d returned\n", pid);
    }
  }
  return 0;
}
