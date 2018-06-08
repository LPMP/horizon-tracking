#include "horizon_tracking.h"
#include "LP_conic_bundle.hxx"

using namespace LP_MP;

int main(int argc, char** argv)
{
    using solver_type = Solver<LP_conic_bundle<FMC_HORIZON_TRACKING>, StandardVisitor>;
    solver_type solver(argc, argv);
    solver.ReadProblem(UAIMaxPotInput::ParseProblemGridAndDecomposeToChains<solver_type>);
    solver.Solve();
}