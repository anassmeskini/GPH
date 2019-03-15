#include "LPSolver.h"

std::string to_str(LPResult::Status st)
{
	switch (st)
	{
	case LPResult::OPTIMAL:
		return "optimal";
	case LPResult::UNBOUNDED:
		return "unbounded";
	case LPResult::INFEASIBLE:
		return "infeasible";
	case LPResult::OTHER:
		return "unknown status";
	}

	assert(0);
	return "";
}
