#include "horizon_tracking.h"

namespace LP_MP
{
    int main()
    {
        std::string path = "Chain-Small-5x1.uai";
        //bool success = LP_MP::UAIMaxPotInput::ParseProblem(path);
        using solver_type = Solver<FMC_HORIZON_TRACKING, LP, StandardVisitor>;
        solver_type solver;
        const bool success = UAIMaxPotInput::ParseProblemGridAndDecomposeToChains(solver, path);
    }
}