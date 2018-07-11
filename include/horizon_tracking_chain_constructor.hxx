#ifndef LP_MP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX
#define LP_MP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX

#include "horizon_tracking_uai_constructor.hxx"

namespace LP_MP {

template<typename MRF_CONSTRUCTOR, typename MAX_CHAIN_FACTOR, typename MAX_POTENTIAL_FACTOR, typename PAIRWISE_MAX_CHAIN_MESSAGE, typename MAX_CHAIN_MAX_POTENTIAL_MESSAGE>
class max_chain_constructor : public MRF_CONSTRUCTOR {
public:
using FMC = typename MRF_CONSTRUCTOR::FMC;
using mrf_constructor = MRF_CONSTRUCTOR;
using max_chain_factor_container = MAX_CHAIN_FACTOR;
using max_potential_factor_container = MAX_POTENTIAL_FACTOR;
using pairwise_max_factor_message_container = PAIRWISE_MAX_CHAIN_MESSAGE;
using max_chain_max_potential_message_container = MAX_CHAIN_MAX_POTENTIAL_MESSAGE;
using mrf_constructor::mrf_constructor;

template <typename ITERATOR, typename TENSORSIZE>
max_chain_factor_container* add_max_chain(ITERATOR node_var_begin, ITERATOR node_var_end,
                                    TENSORSIZE potSizeBegin, TENSORSIZE potSizeEnd, 
                                    const tensor3_variable<REAL>& maxPairwisePotentials,
                                    INDEX chainIndex,
                                    factor_tree<FMC>* t = nullptr)
{
    std::vector<INDEX> numLabels;
    for(auto it = node_var_begin; it!=node_var_end; ++it) {
        const INDEX i = (*it);
        numLabels.push_back( this->GetNumberOfLabels(i) );
    }
    tensor3_variable<REAL> linearPairwisePotentials(potSizeBegin, potSizeEnd);
    for(std::size_t n1=0; n1<maxPairwisePotentials.dim1(); ++n1) {
        for(std::size_t l1=0; l1<maxPairwisePotentials.dim2(n1); ++l1) {
            for(std::size_t l2=0; l2<maxPairwisePotentials.dim3(n1); ++l2) {
                linearPairwisePotentials(n1, l1, l2) = 0.0;
            }
        }
    }
    auto* chain_factor = this->lp_->template add_factor<max_chain_factor_container>(maxPairwisePotentials, linearPairwisePotentials, numLabels, chainIndex);

    INDEX pairwise_index=0;
    INDEX node1_index = 0;
    INDEX node2_index = 1;
    for(auto it = node_var_begin; std::next(it, 1)!=node_var_end; ++it, ++pairwise_index, ++node1_index, ++node2_index) {
        const INDEX i = (*it);
        const INDEX j = *std::next(it, 1);
        auto* pairwise_factor = this->GetPairwiseFactor(i,j);
        auto* msg = this->lp_->template add_message<pairwise_max_factor_message_container>(pairwise_factor, chain_factor, pairwise_index, node1_index, node2_index);

        if(t != nullptr) {
            t->add_message(msg, Chirality::right); 
        }
    }

    return chain_factor;
}

template<typename ITERATOR>
max_potential_factor_container* add_max_potential(ITERATOR max_chain_begin, ITERATOR max_chain_end, factor_tree<FMC>* t = nullptr)
{
    std::vector<std::vector<std::array<REAL,2>>> all_marginals;
    for(auto max_chain_it = max_chain_begin; max_chain_it!=max_chain_end; ++max_chain_it) {
        auto* f = (*max_chain_it)->GetFactor();
        f->MaximizePotentialAndComputePrimal();
        std::vector<std::array<REAL,3>> current_chain_marginals = f->max_potential_marginals();
        std::vector<std::array<REAL,2>> current_chain_marginals_max;
        for (auto current_marginal_item : current_chain_marginals) {
            current_chain_marginals_max.push_back({current_marginal_item[0], current_marginal_item[1]});   // Ignoring the third column in the first iteration. 
        }
        all_marginals.push_back(current_chain_marginals_max);
    }

    auto* max_factor = this->lp_->template add_factor<max_potential_factor_container>(all_marginals);
    for(auto max_chain_it = max_chain_begin; max_chain_it!=max_chain_end; ++max_chain_it) {
        const auto chain_index = std::distance(max_chain_begin, max_chain_it);
        auto* current_chain = *max_chain_it;

        auto* msg = this->lp_->template add_message<max_chain_max_potential_message_container>(current_chain, max_factor, chain_index);

        if(t != nullptr) {
            t->add_message(msg, Chirality::right);
        } 
    }
    return max_factor;
}
};

template<typename FMC, typename MAX_CHAIN_FACTOR, typename MAX_POTENTIAL_FACTOR, typename MAX_CHAIN_MAX_POTENTIAL_MESSAGE>
class max_disjoint_chain_constructor {
public:
using max_chain_factor_container = MAX_CHAIN_FACTOR;
using max_potential_factor_container = MAX_POTENTIAL_FACTOR;
using max_chain_max_potential_message_container = MAX_CHAIN_MAX_POTENTIAL_MESSAGE;

max_chain_factor_container* add_max_chain(const tensor3_variable<REAL>& linearPairwisePotentials, 
                                    const tensor3_variable<REAL>& maxPairwisePotentials,
                                    const std::vector<INDEX>& numLabels,
                                    INDEX chainIndex,
                                    factor_tree<FMC>* t = nullptr)
{
    auto* f = this->lp_->template add_factor<max_chain_factor_container>(maxPairwisePotentials, linearPairwisePotentials, numLabels, chainIndex, false);
    return f;
}

template<typename ITERATOR>
max_potential_factor_container* add_max_potential(ITERATOR max_chain_begin, ITERATOR max_chain_end, factor_tree<FMC>* t = nullptr)
{
    std::vector<std::vector<std::array<REAL,2>>> all_marginals;
    for(auto max_chain_it = max_chain_begin; max_chain_it!=max_chain_end; ++max_chain_it) {
        auto* f = (*max_chain_it)->GetFactor();
        f->MaximizePotentialAndComputePrimal();
        std::vector<std::array<REAL,3>> current_chain_marginals = f->max_potential_marginals();
        std::vector<std::array<REAL,2>> current_chain_marginals_max;
        for (auto current_marginal_item : current_chain_marginals)
        {
            current_chain_marginals_max.push_back({current_marginal_item[0], current_marginal_item[1]});   // Ignoring the third column in the first iteration. 
        }
        all_marginals.push_back(current_chain_marginals_max);
    }

    auto* m = this->lp_->template add_factor<max_potential_factor_container>(all_marginals);
    for(auto max_chain_it = max_chain_begin; max_chain_it!=max_chain_end; ++max_chain_it) {
        const auto i = std::distance(max_chain_begin, max_chain_it);
        auto* c = *max_chain_it;

        auto* msg = this->lp_->template add_message<max_chain_max_potential_message_container>(c, m, i);

        if(t != nullptr) {
            t->add_message(msg, Chirality::right);
        } 
    }
    return m;
}
};

namespace UAIMaxPotInput {
template<typename SOLVER>
   bool ParseProblemChains(const std::string& filename, SOLVER& s)
   {
        using FMC = typename SOLVER::FMC;
        auto& chain_constructor = s.template GetProblemConstructor<0>();

        std::cout << "parsing " << filename << "\n";
        pegtl::file_parser problem(filename);
        std::vector<MaxPotInput> input;
        bool read_suc = problem.parse< grammar, action >(input);
        assert(read_suc);
        assert(input.size() == 2); // One max potential and one linear potential field
       
        INDEX numNodes = input[0].number_of_variables_;
        std::vector<INDEX> numLabels;
        std::vector<typename std::remove_reference_t<decltype(chain_constructor)>::max_chain_factor_container*> max_chain_potentials;
        factor_tree<FMC> tree;

        std::vector<std::vector<INDEX>> functionTableSizes;
        //TODO: Add functionality for multiple disjoint chains.
        for (INDEX currentNodeIndex = 0; currentNodeIndex != numNodes - 1; ++currentNodeIndex)
        {
            INDEX l1Size = input[0].cardinality_[currentNodeIndex];
            INDEX l2Size = input[0].cardinality_[currentNodeIndex + 1];
            //Assuming equivalent max potential and linear potential graphs:
            assert(input[0].cardinality_[currentNodeIndex] == input[1].cardinality_[currentNodeIndex]);
            assert(input[0].cardinality_[currentNodeIndex + 1] == input[1].cardinality_[currentNodeIndex + 1]);
            functionTableSizes.push_back(std::vector<INDEX>{l1Size, l2Size});

            numLabels.push_back(l1Size);
            if (currentNodeIndex == numNodes - 2)
                numLabels.push_back(l2Size);
        }

        tensor3_variable<REAL> linearPairwisePotentials(functionTableSizes.begin(), functionTableSizes.end());
        tensor3_variable<REAL> maxPairwisePotentials(functionTableSizes.begin(), functionTableSizes.end());
            
        // Populate potentials:
        for (INDEX currentPotentialsIndex = 0; currentPotentialsIndex < input.size(); currentPotentialsIndex++)
        { 
            INDEX cliqueIndexChain = 0;
            // process the pairwise potentials first:
            for (INDEX cliqueIndex = 0; cliqueIndex < input[currentPotentialsIndex].clique_scopes_.size(); cliqueIndex++)
            {
                const auto& currentCliqueScope = input[currentPotentialsIndex].clique_scopes_[cliqueIndex];
                if (currentCliqueScope.size() == 1)
                    continue; // Add unaries later.

                INDEX dim23InputIndex = 0;
                for (INDEX l1 = 0; l1 < input[currentPotentialsIndex].cardinality_[currentCliqueScope[0]]; l1++)
                {
                    for (INDEX l2 = 0; l2 < input[currentPotentialsIndex].cardinality_[currentCliqueScope[1]]; l2++)
                    {
                        if (currentPotentialsIndex == 0)
                            linearPairwisePotentials(cliqueIndexChain, l1, l2) = input[currentPotentialsIndex].function_tables_[cliqueIndex][dim23InputIndex];
                        
                        else
                            maxPairwisePotentials(cliqueIndexChain, l1, l2) = input[currentPotentialsIndex].function_tables_[cliqueIndex][dim23InputIndex];

                        dim23InputIndex++;
                    }
                }
                cliqueIndexChain++;
            }
        }

        for (INDEX currentPotentialsIndex = 1; currentPotentialsIndex < input.size(); currentPotentialsIndex++)
        { 
            for (INDEX cliqueIndex = 0; cliqueIndex < input[currentPotentialsIndex].clique_scopes_.size(); cliqueIndex++)
            {
                const auto& currentCliqueScope = input[currentPotentialsIndex].clique_scopes_[cliqueIndex];
                if (currentCliqueScope.size() != 1)
                    continue;
                
                INDEX edgeIndexToRight = currentCliqueScope[0];
                if (edgeIndexToRight < numNodes - 1)
                {
                    // Send to the edges on the right:
                    INDEX numStatesCurrentNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0]];
                    INDEX numStatesNextNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0] + 1];
                    for (INDEX l1 = 0; l1 < numStatesCurrentNode; l1++)
                    {
                        for (INDEX l2 = 0; l2 < numStatesNextNode; l2++)
                        {
                            if (currentPotentialsIndex == 0)
                                linearPairwisePotentials(edgeIndexToRight, l1, l2) += input[currentPotentialsIndex].function_tables_[cliqueIndex][l1];
                            else
                                maxPairwisePotentials(edgeIndexToRight, l1, l2) += input[currentPotentialsIndex].function_tables_[cliqueIndex][l1];
                        }
                    }
                }
                else
                {
                    // Send to the left edges for the corner nodes:
                    INDEX edgeIndexToLeft = edgeIndexToRight - 1;
                    INDEX numStatesCurrentNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0]];
                    INDEX numStatesPrevNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0] - 1];
                    for (INDEX l2 = 0; l2 < numStatesCurrentNode; l2++)
                    {
                        for (INDEX l1 = 0; l1 < numStatesPrevNode; l1++)
                        {
                            if (currentPotentialsIndex == 0)
                                linearPairwisePotentials(edgeIndexToLeft, l1, l2) += input[currentPotentialsIndex].function_tables_[cliqueIndex][l1];
                            else
                                maxPairwisePotentials(edgeIndexToLeft, l1, l2) += input[currentPotentialsIndex].function_tables_[cliqueIndex][l1];
                        }
                    }
                }
            }
        }

        auto* f = chain_constructor.add_max_chain(linearPairwisePotentials, maxPairwisePotentials, numLabels, 0, &tree);
        max_chain_potentials.push_back(f);


        chain_constructor.add_max_potential(max_chain_potentials.begin(), max_chain_potentials.end(), &tree);
        tree.init();
        s.GetLP().add_tree(tree);
        s.GetLP().construct_decomposition();

        return read_suc;
   }

   template<typename SOLVER>
   bool ParseProblemGridAndDecomposeToChains(const std::string& filename, SOLVER& s)
   {
        std::cout << "parsing " << filename << "\n";
        pegtl::file_parser problem(filename);
        std::vector<MaxPotInput> input;
        bool read_suc = problem.parse< grammar, action >(input);
        assert(read_suc);
        assert(input.size() == 2); // One max potential and one linear potential field
         build_bottleneck_labelling_problem_grid_chains(input, s);
        return read_suc;
   }

    template<typename SOLVER>
    bool ParseProblemStringGridAndDecomposeToChains(const std::string& instance, SOLVER& s)
    {
        std::cout << "parsing string\n";
        std::vector<MaxPotInput> input;
        bool read_suc = pegtl::parse<grammar, action>(instance,"",input);
        if(read_suc) {
            assert(input.size() == 2); // One max potential and one linear potential field
            build_bottleneck_labelling_problem_grid_chains(input, s);
        }
        return read_suc;
    }

    template<typename SOLVER>
    void build_bottleneck_labelling_problem_grid_chains(std::vector<MaxPotInput> input, SOLVER& s)
    {
        auto& chain_constructor = s.template GetProblemConstructor<0>();

        using FMC = typename SOLVER::FMC;
        build_mrf(chain_constructor, input[0]);

        auto trees = chain_constructor.compute_forest_cover();
        for(auto& tree : trees) {
            s.GetLP().add_tree(tree);
        }
        
        INDEX numNodes = input[0].number_of_variables_;
        INDEX xEdgeDistance = 1;
        INDEX yEdgeDistance = 0;
        INDEX horizLastNode = 0;
        INDEX horizLength = 0;
        INDEX vertLastNode = 0;
        INDEX vertLength = 0;
        for (const auto& currentCliqueScope : input[0].clique_scopes_) {
            if (currentCliqueScope.size() == 1)
                continue;

            assert(currentCliqueScope.size() <= 2);
            INDEX currentEdgeDistance = currentCliqueScope[1] - currentCliqueScope[0];
            if (yEdgeDistance == 0 && currentEdgeDistance != xEdgeDistance)
                yEdgeDistance = currentEdgeDistance;
                
            else if (currentEdgeDistance != xEdgeDistance)
                assert(yEdgeDistance == currentEdgeDistance);
            
            if (currentCliqueScope[0] == horizLastNode && currentEdgeDistance == xEdgeDistance) {
                horizLength++;
                horizLastNode = currentCliqueScope[1];
            }

            if (currentCliqueScope[0] == vertLastNode && currentEdgeDistance == yEdgeDistance) {
                vertLength++;
                vertLastNode = currentCliqueScope[1];
            }
        }
        INDEX gridSizeX = horizLength + 1;
        INDEX gridSizeY = vertLength + 1;
        INDEX numChains = gridSizeX + gridSizeY;
        if (gridSizeX == 1 || gridSizeY == 1)
            numChains = 1;

        std::vector<std::set<INDEX>> chains(numChains);
        for (INDEX currentChain = 0; currentChain < numChains; currentChain++) {
            if (currentChain < gridSizeY) {
                for (INDEX i = 0; i < gridSizeX; i++)
                    chains[currentChain].insert(i + currentChain * gridSizeX);
            }

            else {
                for (INDEX i = 0; i < gridSizeY; i++)
                    chains[currentChain].insert(i * gridSizeX + currentChain - gridSizeY);
            }
        }

        std::vector<typename std::remove_reference_t<decltype(chain_constructor)>::max_chain_factor_container*> max_chain_potentials;
        factor_tree<FMC> tree;
        INDEX chainIndex = 0;
        for (const auto& currentChain : chains) {
            std::vector<std::vector<INDEX>> functionTableSizes;

            for (auto currentNodeItr = currentChain.begin(); currentNodeItr != std::prev(currentChain.end()); ++currentNodeItr) {
                INDEX l1Size = input[0].cardinality_[*currentNodeItr];
                INDEX l2Size = input[0].cardinality_[*std::next(currentNodeItr, 1)];
                //Assuming equivalent max potential and linear potential graphs:
                assert(input[0].cardinality_[*currentNodeItr] == input[1].cardinality_[*currentNodeItr]);
                assert(input[0].cardinality_[*std::next(currentNodeItr, 1)] == input[1].cardinality_[*std::next(currentNodeItr, 1)]);
                functionTableSizes.push_back(std::vector<INDEX>{l1Size, l2Size});
            }

            tensor3_variable<REAL> maxPairwisePotentials(functionTableSizes.begin(), functionTableSizes.end());
            
            // Populate max potentials:
            for (INDEX currentPotentialsIndex = 1; currentPotentialsIndex < input.size(); currentPotentialsIndex++) { 
                INDEX cliqueIndexChain = 0;
                // process the pairwise potentials first:
                for (INDEX cliqueIndex = 0; cliqueIndex < input[currentPotentialsIndex].clique_scopes_.size(); cliqueIndex++) {
                    const auto& currentCliqueScope = input[currentPotentialsIndex].clique_scopes_[cliqueIndex];
                    if (currentCliqueScope.size() == 1)
                        continue; // Add unaries later.
                    
                    if (currentChain.count(currentCliqueScope[0]) == 0 || 
                        currentChain.count(currentCliqueScope[1]) == 0) 
                        continue;   // Current clique is not present in the current chain.

                    INDEX dim23InputIndex = 0;
                    for (INDEX l1 = 0; l1 < input[currentPotentialsIndex].cardinality_[currentCliqueScope[0]]; l1++) {
                        for (INDEX l2 = 0; l2 < input[currentPotentialsIndex].cardinality_[currentCliqueScope[1]]; l2++) {
                            maxPairwisePotentials(cliqueIndexChain, l1, l2) = input[currentPotentialsIndex].function_tables_[cliqueIndex][dim23InputIndex];
                            dim23InputIndex++;
                        }
                    }
                    cliqueIndexChain++;
                }
            }

            for (INDEX currentPotentialsIndex = 1; currentPotentialsIndex < input.size(); currentPotentialsIndex++) { 
                for (INDEX cliqueIndex = 0; cliqueIndex < input[currentPotentialsIndex].clique_scopes_.size(); cliqueIndex++) {
                    const auto& currentCliqueScope = input[currentPotentialsIndex].clique_scopes_[cliqueIndex];
                    if (currentCliqueScope.size() != 1)
                        continue;
                    
                    if (currentChain.count(currentCliqueScope[0]) == 0) 
                        continue;   // Current clique is not present in the current chain.

                    INDEX chainDelta = *std::next(currentChain.begin(), 1) - *currentChain.begin();
                    INDEX edgeIndexToRight = (currentCliqueScope[0] - *currentChain.begin()) / chainDelta;
                    if (edgeIndexToRight < currentChain.size() - 1) {
                        // Send to the edges on the right:
                        INDEX numStatesCurrentNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0]];
                        INDEX numStatesNextNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0] + chainDelta];
                        for (INDEX l1 = 0; l1 < numStatesCurrentNode; l1++) {
                            for (INDEX l2 = 0; l2 < numStatesNextNode; l2++) {
                                maxPairwisePotentials(edgeIndexToRight, l1, l2) += input[currentPotentialsIndex].function_tables_[cliqueIndex][l1];
                            }
                        }
                    } else {
                        // Send to the left edges for the corner nodes:
                        INDEX edgeIndexToLeft = edgeIndexToRight - 1;
                        INDEX numStatesCurrentNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0]];
                        INDEX numStatesPrevNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0] - chainDelta];
                        for (INDEX l2 = 0; l2 < numStatesCurrentNode; l2++) {
                            for (INDEX l1 = 0; l1 < numStatesPrevNode; l1++) {
                                maxPairwisePotentials(edgeIndexToLeft, l1, l2) += input[currentPotentialsIndex].function_tables_[cliqueIndex][l1];
                            }
                        }
                    }
                }
            }

            auto* f = chain_constructor.add_max_chain(currentChain.begin(), currentChain.end(), functionTableSizes.begin(), functionTableSizes.end(), maxPairwisePotentials, chainIndex++, &tree);
            max_chain_potentials.push_back(f);

        }

        auto* f = chain_constructor.add_max_potential(max_chain_potentials.begin(), max_chain_potentials.end(), &tree);
        s.GetLP().add_tree(tree);
        s.GetLP().construct_decomposition();
   }
}
}


#endif // LP_MP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX
