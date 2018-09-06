#include "test.h"
#include "horizon_tracking.h"
#include "horizon_tracking_uai_input.h"
#include "horizon_tracking_chain_constructor.hxx"
#include "solver.hxx"
#include "LP_FWMAP.hxx"
#include "data_max_potential_chains_test.hxx"
#include <type_traits>
#include "horizon_tracking_primal_rounding.hxx"

using namespace LP_MP;

int main()
{
    // Test on 3 Node Chain:
    {
        using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING_CHAINS>, StandardVisitor>;
        solver_type solver(solver_options_small);
        auto input = horizon_tracking_uai_input::parse_string(chain_uai_input_small);
        construct_horizon_tracking_problem_on_grid_to_chains(input, solver, solver.template GetProblemConstructor<0>());
        solver.Solve();
        test(solver.lower_bound(), 49, 0);
        solver.GetLP().write_back_reparametrization();
        test(solver.GetLP().original_factors_lower_bound(), 49, 0);
        round_primal_solution(solver);

        test(solver.primal_cost(), 49, 0);
    }

    // Test on 5 Node Chain:
    {
        using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING_CHAINS>, StandardVisitor>;
        solver_type solver(solver_options_medium);
        auto input = horizon_tracking_uai_input::parse_string(chain_uai_input_medium);
        construct_horizon_tracking_problem_on_grid_to_chains(input, solver, solver.template GetProblemConstructor<0>());
        order_nodes_by_label_space_cadinality(solver.template GetProblemConstructor<0>());
        solver.Solve();
        test(solver.lower_bound(), 0.085579, 0);
        solver.GetLP().write_back_reparametrization();
        test(solver.GetLP().original_factors_lower_bound(), 0.085579, 0);
        round_primal_solution(solver);
        test(solver.GetLP().original_factors_lower_bound(), 0.085579, 0);
        test(solver.primal_cost(), 0.085579, 0);
    }
}
