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
    auto numF = solver.GetLP().GetNumberOfFactors();
    std::vector<FactorTypeAdapter*> factors;
    for (auto i = 0; i < numF; i++)
    {
        auto currentF = solver.GetLP().GetFactor(i);
        if (dynamic_cast<FMC_HORIZON_TRACKING::UnaryFactor*>(currentF) || 
            dynamic_cast<FMC_HORIZON_TRACKING::PairwiseFactor*>(currentF))
            factors.push_back(currentF);
    }
    solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::forward>(factors.begin(), factors.end(), 1);
    solver.RegisterPrimal();
    solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::backward>(factors.begin(), factors.end(), 2);
    solver.RegisterPrimal();
    solver.WritePrimal();
    std::cout<<"\n Primal Cost: "<<solver.primal_cost();
}
