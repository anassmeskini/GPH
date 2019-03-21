#ifndef _AVAILABLE_SOLVER_HPP
#define _AVAILABLE_SOLVER_HPP

#include "interface/CPXSolver.h"
#include "interface/GLPKSolver.h"
#include "interface/SPXSolver.h"

#ifdef CONCERT_CPLEX_FOUND
using AvaiLPSolver = CPXSolver;
#elif SOPLEX_FOUND
using AvaiLPSolver = SPXSolver;
#elif GLPK_FOUND
using AvaiLPSolver = GLPKSolver;
#endif

#endif
