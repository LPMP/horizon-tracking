#include "test.h"
#include "horizon_tracking.h" // rename to horizon_tracking.h
#include "horizon_tracking_constructor.hxx"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"

using namespace LP_MP;

int main()
{
    std::string path = "Chain-Small-5x1.uai";
    //bool success = LP_MP::UAIMaxPotInput::ParseProblem(path);
    using solver_type = Solver<FMC_HORIZON_TRACKING, LP, StandardVisitor>;
    solver_type solver;
    const bool success = UAIMaxPotInput::ParseProblemGridAndDecomposeToChains(solver, path);
    test(success == 1);
}
