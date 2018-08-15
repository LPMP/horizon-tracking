#include "horizon_tracking.h"
#include "tree_decomposition.hxx"

using namespace LP_MP;

int main(int argc, char** argv)
{
    using solver_type = Solver<LP_subgradient_ascent<FMC_HORIZON_TRACKING>, StandardVisitor>;
    solver_type solver(argc, argv);
    auto input = horizon_tracking_uai_input::parse_file(solver.get_input_file());
    construct_horizon_tracking_problem_on_grid(input, solver, solver.template GetProblemConstructor<0>());
    order_nodes_by_label_space_cadinality(solver.template GetProblemConstructor<0>());
    solver.Solve();
}
