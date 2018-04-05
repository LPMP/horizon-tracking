#include "test.h"
#include "horizon_tracking_factors.hxx"

using namespace LP_MP;

int main()
{
    // problem with zero linear potentials
    const int numNodes = 3;
    std::vector<INDEX> numLabels = {3, 2, 3};

    std::vector<matrix<REAL>> LinearPairwisePotentials(numNodes - 1);
    for (int i = 0; i < numNodes - 1; i++)
    {
        LinearPairwisePotentials[i] = matrix<REAL>(numLabels[i], numLabels[i + 1]);
    }

    std::vector<matrix<REAL>> MaxPairwisePotentials(numNodes - 1);

    MaxPairwisePotentials[0] = matrix<REAL>(3,2);
    MaxPairwisePotentials[0](0,0) = 1;
    MaxPairwisePotentials[0](1,0) = 2;
    MaxPairwisePotentials[0](2,0) = 1;
    MaxPairwisePotentials[0](0,1) = 3;
    MaxPairwisePotentials[0](1,1) = 4;
    MaxPairwisePotentials[0](2,1) = 2;

    MaxPairwisePotentials[1] = matrix<REAL>(2,3);
    MaxPairwisePotentials[1](0,0) = 1.5;
    MaxPairwisePotentials[1](0,1) = 5;
    MaxPairwisePotentials[1](0,2) = 1;
    MaxPairwisePotentials[1](1,0) = 3;
    MaxPairwisePotentials[1](1,1) = 1;
    MaxPairwisePotentials[1](1,2) = 6;

    {
        max_potential_on_chain chain = max_potential_on_chain(MaxPairwisePotentials, LinearPairwisePotentials, numLabels);
        REAL objective = chain.LowerBound();
        test(objective == 1);
    }

    // now add linear potentials
    LinearPairwisePotentials[0](0,0) = 100;
    LinearPairwisePotentials[0](1,0) = 100;
    LinearPairwisePotentials[0](2,0) = 100;
    LinearPairwisePotentials[0](0,1) = 100;
    LinearPairwisePotentials[0](1,1) = 100;
    LinearPairwisePotentials[0](2,1) = 0;

    LinearPairwisePotentials[1](0,0) = 100;
    LinearPairwisePotentials[1](0,1) = 100;
    LinearPairwisePotentials[1](0,2) = 100;
    LinearPairwisePotentials[1](1,0) = 100;
    LinearPairwisePotentials[1](1,1) = 100;
    LinearPairwisePotentials[1](1,2) = 0;

    {
        max_potential_on_chain chain = max_potential_on_chain(MaxPairwisePotentials, LinearPairwisePotentials, numLabels);
        REAL objective = chain.LowerBound();
        test(objective == 6);
    }

    LinearPairwisePotentials[0](2,0) = 1.5;
    LinearPairwisePotentials[1](0,2) = 1.5;

    {
        max_potential_on_chain chain = max_potential_on_chain(MaxPairwisePotentials, LinearPairwisePotentials, numLabels);
        REAL objective = chain.LowerBound();
        test(objective == 4);
    }

}
