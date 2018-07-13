#include "test.h"
#include "horizon_tracking.h"
#include "solver.hxx"
#include "LP_FWMAP.hxx"
#include "test_max_potential.hxx"

using namespace LP_MP;

int main()
{
    {
        using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
        solver_type solver(solver_options);
        UAIMaxPotInput::ParseProblemStringGridAndDecomposeToChains<solver_type>(chain_uai_input, solver);
        solver.Solve();
        test(std::abs(solver.lower_bound() - 17) <= eps);
        solver.GetLP().write_back_reparametrization();
        test(std::abs(solver.GetLP().original_factors_lower_bound() - 17) <= eps);
    }
    {
        using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
        solver_type solver(solver_options);
        UAIMaxPotInput::ParseProblemStringGridAndDecomposeToChains<solver_type>(grid_uai_input, solver);
        solver.Solve();
        test(std::abs(solver.lower_bound() - 4.307381) <= eps);
        solver.GetLP().write_back_reparametrization();
        test(std::abs(solver.GetLP().original_factors_lower_bound() - 4.307381) <= eps);
    }
}
