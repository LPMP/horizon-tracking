#include "test.h"
#include "horizon_tracking.h"
#include "horizon_tracking_constructor.hxx"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include "LP_FWMAP.hxx"

using namespace LP_MP;

int main()
{
    {
        char *argv[] = {"horizon_tracking_test", "-i", "../../test/Chain-3x1.uai", "--maxIter", "5", "--timeout", "60", "-v", "0", NULL};
        int argc = sizeof(argv) / sizeof(char*) - 1;
        using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
        solver_type solver(argc, argv);
        solver.ReadProblem(UAIMaxPotInput::ParseProblemGridAndDecomposeToChains<solver_type>);
        solver.Solve();
        test(std::abs(solver.lower_bound() - 17) <= eps);
    }
    {
        char *argv[] = {"horizon_tracking_test", "-i", "../../test/Grid-5x5.uai", "--maxIter", "50", "--timeout", "60", "-v", "1", NULL};
        int argc = sizeof(argv) / sizeof(char*) - 1;
        using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
        solver_type solver(argc, argv);
        solver.ReadProblem(UAIMaxPotInput::ParseProblemGridAndDecomposeToChains<solver_type>);
        solver.Solve();
        test(std::abs(solver.lower_bound() - 4.307381) <= eps);
    }
}
