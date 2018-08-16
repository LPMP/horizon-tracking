#include "test.h"
#include "horizon_tracking.h"
#include "solver.hxx"
#include "LP_FWMAP.hxx"
#include "data_max_potential_grids_test.hxx"
#include <type_traits>

using namespace LP_MP;

int main()
{
    // Test on 2x2 Grid:
    {
        using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
        solver_type solver(solver_options_small);
        auto input = horizon_tracking_uai_input::parse_string(grid_uai_input_small);
        construct_horizon_tracking_problem_on_grid(input, solver, solver.template GetProblemConstructor<0>());
        order_nodes_by_label_space_cadinality(solver.template GetProblemConstructor<0>());
        solver.Solve();
        test(std::abs(solver.lower_bound() -  0.242203) <= eps);
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
        test(std::abs(solver.GetLP().original_factors_lower_bound() - 0.242203) <= eps);
        auto numF = solver.GetLP().GetNumberOfFactors();
        std::vector<FactorTypeAdapter*> factors;
        for (auto i = 0; i < numF; i++)
        {
            auto currentF = solver.GetLP().GetFactor(i);
            if (dynamic_cast<FMC_HORIZON_TRACKING::UnaryFactor*>(currentF) || 
                dynamic_cast<FMC_HORIZON_TRACKING::PairwiseFactor*>(currentF))
                factors.push_back(currentF);
        }
        test(solver.primal_cost() == std::numeric_limits<REAL>::infinity());
        solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::forward>(factors.begin(), factors.end(), 1);
        solver.RegisterPrimal();
        solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::backward>(factors.begin(), factors.end(), 2);
        solver.RegisterPrimal();
        test(solver.primal_cost() < std::numeric_limits<REAL>::max());
        test((solver.primal_cost() - 0.242203) <= eps);
    }

    // Test on 3x3 Grid:
    {
        using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
        solver_type solver(solver_options_small);
        auto input = horizon_tracking_uai_input::parse_string(grid_uai_input_3x3);
        construct_horizon_tracking_problem_on_grid(input, solver, solver.template GetProblemConstructor<0>());
        order_nodes_by_label_space_cadinality(solver.template GetProblemConstructor<0>());
        solver.Solve();
        test(std::abs(solver.lower_bound() -  0.590735) <= eps);
        solver.GetLP().write_back_reparametrization();
        auto chain_constructor = solver.template GetProblemConstructor<0>();
        auto max_potential_factors = chain_constructor.max_graph_factors();
        for(auto* f : max_potential_factors) {
            f->send_messages_to_left<FMC_HORIZON_TRACKING::max_chain_to_max_potential_message_container>(1.0);
        }

        auto max_potential_on_chain_factors = chain_constructor.max_chain_factors();
        for(auto* f : max_potential_on_chain_factors) {
            f->send_messages_to_left<FMC_HORIZON_TRACKING::pairwise_max_message_container>(1.0);
        }
        test(std::abs(solver.GetLP().original_factors_lower_bound() - 0.590735) <= eps);
        auto numF = solver.GetLP().GetNumberOfFactors();
        std::vector<FactorTypeAdapter*> factors;
        for (auto i = 0; i < numF; i++)
        {
            auto currentF = solver.GetLP().GetFactor(i);
            if (dynamic_cast<FMC_HORIZON_TRACKING::UnaryFactor*>(currentF) || 
                dynamic_cast<FMC_HORIZON_TRACKING::PairwiseFactor*>(currentF))
                factors.push_back(currentF);
        }
        test(solver.primal_cost() == std::numeric_limits<REAL>::infinity());
        solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::forward>(factors.begin(), factors.end(), 1);
        solver.RegisterPrimal();
        solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::backward>(factors.begin(), factors.end(), 2);
        solver.RegisterPrimal();
        test(solver.primal_cost() < std::numeric_limits<REAL>::max());
        test((solver.primal_cost() -  0.590735) <= eps);
    }

    // Test on 5x5 Grid:
    {
        using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
        solver_type solver(solver_options_medium);
        auto input = horizon_tracking_uai_input::parse_string(grid_uai_input_medium);
        construct_horizon_tracking_problem_on_grid(input, solver, solver.template GetProblemConstructor<0>());
        order_nodes_by_label_space_cadinality(solver.template GetProblemConstructor<0>());
        solver.Solve();
        test(std::abs(solver.lower_bound() - 4.307381) <= eps);
        solver.GetLP().write_back_reparametrization();
                auto chain_constructor = solver.template GetProblemConstructor<0>();
        auto max_potential_factors = chain_constructor.max_graph_factors();
        for(auto* f : max_potential_factors) {
            f->send_messages_to_left<FMC_HORIZON_TRACKING::max_chain_to_max_potential_message_container>(1.0);
        }

        auto max_potential_on_chain_factors = chain_constructor.max_chain_factors();
        for(auto* f : max_potential_on_chain_factors) {
            f->send_messages_to_left<FMC_HORIZON_TRACKING::pairwise_max_message_container>(1.0);
        }
        test(std::abs(solver.GetLP().original_factors_lower_bound() - 4.307381) <= eps);
        auto numF = solver.GetLP().GetNumberOfFactors();
        std::vector<FactorTypeAdapter*> factors;
        for (auto i = 0; i < numF; i++)
        {
            auto currentF = solver.GetLP().GetFactor(i);
            if (dynamic_cast<FMC_HORIZON_TRACKING::UnaryFactor*>(currentF) || 
                dynamic_cast<FMC_HORIZON_TRACKING::PairwiseFactor*>(currentF))
                factors.push_back(currentF);
        }
        test(solver.primal_cost() == std::numeric_limits<REAL>::infinity());
        solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::backward>(factors.begin(), factors.end(), 1);
        solver.RegisterPrimal();
        solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::forward>(factors.begin(), factors.end(), 2);
        solver.RegisterPrimal();
        test(solver.primal_cost() < std::numeric_limits<REAL>::max());
        test((solver.primal_cost() - 4.307381) <= eps);
    }
}
