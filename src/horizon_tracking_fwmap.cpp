#include "horizon_tracking.h"
#include "LP_FWMAP.hxx"

using namespace LP_MP;

int main(int argc, char** argv)
{
    using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
    solver_type solver(argc, argv);

    auto input = horizon_tracking_uai_input::parse_file(solver.get_input_file());
    construct_horizon_tracking_problem_on_grid(input, solver, solver.template GetProblemConstructor<0>());
    order_nodes_by_label_space_cadinality(solver.template GetProblemConstructor<0>());

    solver.Solve();
    solver.GetLP().write_back_reparametrization();
    // mark mrf factors
    auto numF = solver.GetLP().GetNumberOfFactors();
    std::vector<FactorTypeAdapter*> factors;
    for (auto i = 0; i < numF; i++)
    {
        auto currentF = solver.GetLP().GetFactor(i);
        if (dynamic_cast<FMC_HORIZON_TRACKING::UnaryFactor*>(currentF) || 
            dynamic_cast<FMC_HORIZON_TRACKING::PairwiseFactor*>(currentF))
            factors.push_back(currentF);
    }
    solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::backward>(factors.begin(), factors.end(),std::numeric_limits<INDEX>::max()-2);
    solver.RegisterPrimal();
    solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::forward>(factors.begin(), factors.end(),std::numeric_limits<INDEX>::max()-1);
    solver.RegisterPrimal();
    solver.WritePrimal();
    std::cout<<"\n\n Primal Cost: "<<solver.primal_cost();
    std::cout<<"\n Percentage duality gap: "<<100.0 * (solver.primal_cost() - solver.lower_bound()) / solver.lower_bound() <<"%\n\n";
}
 
