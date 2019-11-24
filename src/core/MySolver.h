#ifndef _AVAILABLE_SOLVER_HPP
#define _AVAILABLE_SOLVER_HPP

#include "interfaces/CPXSolver.h"
#include "interfaces/GLPKSolver.h"
#include "interfaces/SPXSolver.h"

#ifdef CONCERT_CPLEX_FOUND
using MySolver = CPXSolver;
#elif SOPLEX_FOUND
using MySolver = SPXSolver;
#elif GLPK_FOUND
using MySolver = GLPKSolver;
#else
static_assert(0, "No LP solver available");
#endif

#endif
