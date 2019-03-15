#ifndef _AVAILABLE_SOLVER_HPP
#define _AVAILABLE_SOLVER_HPP

#include "CPXSolver.h"
#include "DummySolver.h"
#include "SPXSolver.h"

#ifdef CONCERT_CPLEX_FOUND
using AvaiLPSolver = CPXSolver;
#elif SOPLEX_FOUND
using AvaiLPSolver = SPXSolver;
#endif

#endif
