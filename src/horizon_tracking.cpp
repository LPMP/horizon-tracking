#include "horizon_tracking.h"

namespace LP_MP
{
    int main()
    {
        std::string path = "Grid-Small-5x5.uai";
        using solver_type = Solver<FMC_HORIZON_TRACKING, LP, StandardVisitor>;
        solver_type solver;
        const bool success = UAIMaxPotInput::ParseProblemGridAndDecomposeToChains(solver, path);
        solver.Solve();
    }
}