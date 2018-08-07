#include "test.h"
#include "horizon_tracking_factors.hxx"

using namespace LP_MP;

int main()
{
    // problem with zero linear potentials
    const int numNodes = 3;
    std::vector<INDEX> numLabels = {3, 2, 3};
    std::vector<std::array<INDEX,2>> potential_size{{3,2}, {2,3}};

    three_dimensional_variable_array<REAL> LinearPairwisePotentials(potential_size.begin(), potential_size.end());
    three_dimensional_variable_array<REAL> MaxPairwisePotentials(potential_size.begin(), potential_size.end());
    std::vector<std::array<INDEX, 2>> messagePassingSchedule;
    std::vector<INDEX> numEdgesForNode({1, 2, 1});

    messagePassingSchedule.push_back({0, 1});
    messagePassingSchedule.push_back({2, 1});

    MaxPairwisePotentials(0,0,0) = 1;
    MaxPairwisePotentials(0,1,0) = 2;
    MaxPairwisePotentials(0,2,0) = 1;
    MaxPairwisePotentials(0,0,1) = 3;
    MaxPairwisePotentials(0,1,1) = 4;
    MaxPairwisePotentials(0,2,1) = 2;

    MaxPairwisePotentials(1,0,0) = 1.5;
    MaxPairwisePotentials(1,0,1) = 5;
    MaxPairwisePotentials(1,0,2) = 1;
    MaxPairwisePotentials(1,1,0) = 3;
    MaxPairwisePotentials(1,1,1) = 1;
    MaxPairwisePotentials(1,1,2) = 6;

    {
        max_potential_on_tree tree = max_potential_on_tree(MaxPairwisePotentials, LinearPairwisePotentials, numLabels, messagePassingSchedule, numEdgesForNode);
        REAL lb = tree.GetMaxPotLowerBound();
        REAL ub = tree.GetMaxPotUpperBound();
        tree.MaximizePotentialAndComputePrimal();
        std::vector<std::array<REAL, 3>> rootMarginals = tree.max_potential_marginals();
        REAL bestSumMarginal = INFINITY;
        for (const auto& currentM : rootMarginals)
        {
            bestSumMarginal = fmin(currentM[0] + currentM[1], bestSumMarginal);
        }
        test(bestSumMarginal == 1);
    }

    // now add linear potentials
    LinearPairwisePotentials(0,0,0) = 100;
    LinearPairwisePotentials(0,1,0) = 100;
    LinearPairwisePotentials(0,2,0) = 100;
    LinearPairwisePotentials(0,0,1) = 100;
    LinearPairwisePotentials(0,1,1) = 100;
    LinearPairwisePotentials(0,2,1) = 0;

    LinearPairwisePotentials(1,0,0) = 100;
    LinearPairwisePotentials(1,0,1) = 100;
    LinearPairwisePotentials(1,0,2) = 100;
    LinearPairwisePotentials(1,1,0) = 100;
    LinearPairwisePotentials(1,1,1) = 100;
    LinearPairwisePotentials(1,1,2) = 0;

    {
        max_potential_on_tree tree = max_potential_on_tree(MaxPairwisePotentials, LinearPairwisePotentials, numLabels, messagePassingSchedule, numEdgesForNode);
        REAL lb = tree.GetMaxPotLowerBound();
        REAL ub = tree.GetMaxPotUpperBound();
        tree.MaximizePotentialAndComputePrimal();
        std::vector<std::array<REAL, 3>> rootMarginals = tree.max_potential_marginals();
        REAL bestSumMarginal = INFINITY;
        for (const auto& currentM : rootMarginals)
        {
            bestSumMarginal = fmin(currentM[0] + currentM[1], bestSumMarginal);
        }
        test(bestSumMarginal == 6);
    }

    LinearPairwisePotentials(0,2,0) = 1.5;
    LinearPairwisePotentials(1,0,2) = 1.5;

    {
        max_potential_on_tree tree = max_potential_on_tree(MaxPairwisePotentials, LinearPairwisePotentials, numLabels, messagePassingSchedule, numEdgesForNode);
        REAL lb = tree.GetMaxPotLowerBound();
        REAL ub = tree.GetMaxPotUpperBound();
        tree.MaximizePotentialAndComputePrimal();
        std::vector<std::array<REAL, 3>> rootMarginals = tree.max_potential_marginals();
        REAL bestSumMarginal = INFINITY;
        for (const auto& currentM : rootMarginals)
        {
            bestSumMarginal = fmin(currentM[0] + currentM[1], bestSumMarginal);
        }
        test(bestSumMarginal == 4);
    }
}
