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
        UAIMaxPotInput::ParseProblemStringGridAndDecomposeToChains<solver_type>(grid_uai_input_small, solver);
        solver.Solve();
        test(std::abs(solver.lower_bound() -  0.242203) <= eps);
        solver.GetLP().write_back_reparametrization();
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
        UAIMaxPotInput::ParseProblemStringGridAndDecomposeToChains<solver_type>(grid_uai_input_3x3, solver);
        solver.Solve();
        test(std::abs(solver.lower_bound() -  0.590735) <= eps);
        solver.GetLP().write_back_reparametrization();
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
        UAIMaxPotInput::ParseProblemStringGridAndDecomposeToChains<solver_type>(grid_uai_input_medium, solver);
        solver.Solve();
        test(std::abs(solver.lower_bound() - 4.307381) <= eps);
        solver.GetLP().write_back_reparametrization();
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
