#include "horizon_tracking.h"
#include "tree_decomposition.hxx"

using namespace LP_MP;

int main(int argc, char** argv)
{
    using solver_type = Solver<LP_subgradient_ascent<FMC_HORIZON_TRACKING>, StandardVisitor>;
    solver_type solver(argc, argv);
    solver.ReadProblem(UAIMaxPotInput::ParseProblemGridAndDecomposeToChains<solver_type>);
    solver.Solve();
}