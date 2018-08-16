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
    
    // send messages from botleneck potentials down to mrf
    auto chain_constructor = solver.template GetProblemConstructor<0>();
    auto max_potential_factors = chain_constructor.max_graph_factors();
    for(auto* f : max_potential_factors) {
        f->send_messages_to_left<FMC_HORIZON_TRACKING::max_chain_to_max_potential_message_container>(1.0);
    }

    auto max_potential_on_chain_factors = chain_constructor.max_chain_factors();
    for(auto* f : max_potential_on_chain_factors) {
        f->send_messages_to_left<FMC_HORIZON_TRACKING::pairwise_max_message_container>(1.0);
    }

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
 
