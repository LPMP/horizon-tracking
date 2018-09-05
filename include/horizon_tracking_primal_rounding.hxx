#ifndef LP_MP_HORIZON_TRACKING_PRIMAL_ROUNDING_HXX
#define LP_MP_HORIZON_TRACKING_PRIMAL_ROUNDING_HXX

#include "horizon_tracking.h"
#include "solver.hxx"
#include "LP_MP.h"

using namespace LP_MP;

template<typename SOLVER>
std::vector<FactorTypeAdapter*> get_mrf_factors(SOLVER& solver)
{
    std::vector<FactorTypeAdapter*> mrf_factors;
    for (auto i = 0; i < solver.GetLP().GetNumberOfFactors(); i++) {
        auto f = solver.GetLP().GetFactor(i);      
        if (dynamic_cast<FMC_HORIZON_TRACKING_CHAINS::UnaryFactor*>(f) || 
            dynamic_cast<FMC_HORIZON_TRACKING_CHAINS::PairwiseFactor*>(f))
            mrf_factors.push_back(f);
    }    
    return mrf_factors;
}

template<typename SOLVER>
void round_primal_solution(SOLVER& solver)
{
    //TO ADDRESS: Changing this primal pass make primal solution on 5x5 grid test very weak.
    solver.GetLP().ComputePassAndPrimal(std::numeric_limits<INDEX>::max()-4);
    auto mrf_factors = get_mrf_factors(solver);
    solver.GetLP().template ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::forward>(mrf_factors.begin(), mrf_factors.end(), std::numeric_limits<INDEX>::max()-2);
    solver.RegisterPrimal();
    solver.GetLP().template ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::backward>(mrf_factors.begin(), mrf_factors.end(), std::numeric_limits<INDEX>::max()-1);
    solver.RegisterPrimal();
}

#endif //LP_MP_HORIZON_TRACKING_PRIMAL_ROUNDING_HXX