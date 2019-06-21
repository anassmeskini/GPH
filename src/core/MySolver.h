#ifndef _AVAILABLE_SOLVER_HPP
#define _AVAILABLE_SOLVER_HPP

#include "interface/CPXSolver.h"
#include "interface/GLPKSolver.h"
#include "interface/SPXSolver.h"

#ifdef CONCERT_CPLEX_FOUND
using MySolver = CPXSolver;
#elif SOPLEX_FOUND
using MySolver = SPXSolver;
#elif GLPK_FOUND
using MySolver = GLPKSolver;
#endif

#endif
