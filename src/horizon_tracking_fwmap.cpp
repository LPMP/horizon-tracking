#include "horizon_tracking.h"
#include "LP_FWMAP.hxx"

using namespace LP_MP;

int main(int argc, char** argv)
{
    using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
    solver_type solver(argc, argv);
    solver.ReadProblem(UAIMaxPotInput::ParseProblemGridAndDecomposeToChains<solver_type>);
    solver.Solve();
    solver.GetLP().write_back_reparametrization();
    // mark mrf factors
    solver.GetLP();
}
