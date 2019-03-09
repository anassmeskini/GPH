#ifndef _AVAILABLE_SOLVER_HPP
#define _AVAILABLE_SOLVER_HPP

#include "CPXSolver.h"
#include "DummySolver.h"

#ifdef CONCERT_CPLEX_FOUND
using AvaiLPSolver = CPXSolver;
#else
using AvaiLPSolver = DSolver;
#endif

#endif
