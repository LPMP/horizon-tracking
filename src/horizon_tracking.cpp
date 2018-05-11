#include "horizon_tracking.h"

using namespace LP_MP;

int main()
{
    std::string path = "Grid-Small-5x5.uai";
    using solver_type = Solver<LP<FMC_HORIZON_TRACKING>, StandardVisitor>;
    solver_type solver;
    const bool success = UAIMaxPotInput::ParseProblemGridAndDecomposeToChains(solver, path);
    solver.Solve();
}
