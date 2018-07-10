#include "test.h"
#include "horizon_tracking_factors.hxx"

using namespace LP_MP;

int main()
{
    // problem with zero linear potentials
    const int numNodes = 3;
    std::vector<INDEX> numLabels = {3, 2, 3};
    std::vector<std::array<INDEX,2>> potential_size{{3,2}, {2,3}};

    tensor3_variable<REAL> LinearPairwisePotentials(potential_size.begin(), potential_size.end());
    tensor3_variable<REAL> MaxPairwisePotentials(potential_size.begin(), potential_size.end());

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
        std::vector<std::vector<std::array<REAL,2>>> all_marginals;
        max_potential_on_chain chain = max_potential_on_chain(MaxPairwisePotentials, LinearPairwisePotentials, numLabels, 0);
        chain.MaximizePotentialAndComputePrimal();
        std::vector<std::array<REAL,3>> current_chain_marginals = chain.max_potential_marginals();
        std::vector<std::array<REAL,2>> current_chain_marginals_max;
        for (auto current_marginal_item : current_chain_marginals)
        {
            current_chain_marginals_max.push_back({current_marginal_item[0], current_marginal_item[1]});   // Ignoring the third column in the first iteration. 
        }
        all_marginals.push_back(current_chain_marginals_max);
        max_potential_on_graph graph = max_potential_on_graph(all_marginals);
        graph.MaximizePotentialAndComputePrimal();
        chain.set_max_potential_index(graph.max_potential_index(0));
        chain.MaximizePotentialAndComputePrimal();
        test(graph.EvaluatePrimal() + chain.EvaluatePrimal() == 1);
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
        std::vector<std::vector<std::array<REAL,2>>> all_marginals;
        max_potential_on_chain chain = max_potential_on_chain(MaxPairwisePotentials, LinearPairwisePotentials, numLabels, 0);
        chain.MaximizePotentialAndComputePrimal();
        std::vector<std::array<REAL,3>> current_chain_marginals = chain.max_potential_marginals();
        std::vector<std::array<REAL,2>> current_chain_marginals_max;
        for (auto current_marginal_item : current_chain_marginals)
        {
            current_chain_marginals_max.push_back({current_marginal_item[0], current_marginal_item[1]});
        }
        all_marginals.push_back(current_chain_marginals_max);
        max_potential_on_graph graph = max_potential_on_graph(all_marginals);
        graph.MaximizePotentialAndComputePrimal();
        chain.set_max_potential_index(graph.max_potential_index(0));
        chain.MaximizePotentialAndComputePrimal();
        test(graph.EvaluatePrimal() + chain.EvaluatePrimal() == 6);
    }

    LinearPairwisePotentials(0,2,0) = 1.5;
    LinearPairwisePotentials(1,0,2) = 1.5;

    {
        std::vector<std::vector<std::array<REAL,2>>> all_marginals;
        max_potential_on_chain chain = max_potential_on_chain(MaxPairwisePotentials, LinearPairwisePotentials, numLabels, 0);
        chain.MaximizePotentialAndComputePrimal();
        std::vector<std::array<REAL,3>> current_chain_marginals = chain.max_potential_marginals();
        std::vector<std::array<REAL,2>> current_chain_marginals_max;
        for (auto current_marginal_item : current_chain_marginals)
        {
            current_chain_marginals_max.push_back({current_marginal_item[0], current_marginal_item[1]});   
        }
        all_marginals.push_back(current_chain_marginals_max);
        max_potential_on_graph graph = max_potential_on_graph(all_marginals);
        graph.MaximizePotentialAndComputePrimal();
        chain.set_max_potential_index(graph.max_potential_index(0));
        chain.MaximizePotentialAndComputePrimal();
        test(graph.EvaluatePrimal() + chain.EvaluatePrimal() == 7);
    }
}
