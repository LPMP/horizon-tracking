#include "test.h"
#include "horizon_tracking.h"
#include "solver.hxx"
#include "LP_FWMAP.hxx"
#include "test_max_potential.hxx"
#include <type_traits>

using namespace LP_MP;

int main()
{
    using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
    solver_type solver(solver_options);
    UAIMaxPotInput::ParseProblemStringGridAndDecomposeToChains<solver_type>(chain_uai_input, solver);
    solver.Solve();
    test(std::abs(solver.lower_bound() - 49) <= eps);
    std::cout<<"OriginalFactorLowerBound:"<<solver.GetLP().original_factors_lower_bound()<<std::endl;
    solver.GetLP().write_back_reparametrization();
    test(std::abs(solver.GetLP().original_factors_lower_bound() - 49) <= eps);
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
    solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::backward>(factors.begin(), factors.end());
    solver.RegisterPrimal();
    test(std::abs(solver.primal_cost() - 49) <= eps);
}