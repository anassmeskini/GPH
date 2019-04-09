#include "core/AvaiLPSolver.h"
#include "core/Common.h"
#include "core/LPSolver.h"
#include "core/MIP.h"
#include "core/MPSReader.h"
#include "core/Timer.h"
#include "fmt/format.h"
#include <memory>

int
main(int argc, char** argv)
{
  MIP<double> mip;
  std::string filename("mip.mps");

  if (argc == 2)
    filename = std::string(argv[1]);

  try {
    auto t0 = Timer::now();
    mip = mpsreader::parse(filename);
    auto t1 = Timer::now();

    fmt::print("Reading the problem took: {}\n", Timer::seconds(t1, t0));
  } catch (const std::exception& ex) {
    std::cout << ex.what();
    return 1;
  }

  std::unique_ptr<LPSolver<double>> solver(new AvaiLPSolver(mip));
  try {
    auto t0 = Timer::now();
    LPResult result = solver->solve();
    auto t1 = Timer::now();

    fmt::print("Solving the LP took: {}\n", Timer::seconds(t1, t0));

    fmt::print("LP solver return status: {}\n", to_str(result.status));

    if (result.status == LPResult::OPTIMAL) {
      fmt::print("obj: {}\n", result.obj);
    }

  } catch (...) {
    fmt::print("Solver raised an exception\n");
  }
  return 0;
}
