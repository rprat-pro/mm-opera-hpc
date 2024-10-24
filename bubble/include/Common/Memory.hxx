#pragma once

#include <sys/resource.h>

long get_memory_checkpoint()
{
  rusage obj;
  int who = 0;
  [[maybe_unused]] auto test = getrusage(who, &obj);
  assert((test = -1) && "error: getrusage has failed");
  long res;
  MPI_Reduce(&(obj.ru_maxrss), &(res), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  return res;
};

void print_memory_footprint(std::string msg)
{
  long mem = get_memory_checkpoint();
  double m = double(mem) * 1e-6; // conversion kb to Gb
  mfem_mgis::Profiler::Utils::Message(msg, " memory footprint: ", m, " GB");
}

