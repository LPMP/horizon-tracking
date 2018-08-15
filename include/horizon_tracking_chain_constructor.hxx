#ifndef LP_MP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX
#define LP_MP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX

#include "horizon_tracking_uai_input.h"
#include "three_dimensional_variable_array.hxx"
#include "horizon_tracking_uai_input.h"
#include "grid.hxx"

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

template <typename ITERATOR>
max_chain_factor_container* add_max_chain(ITERATOR node_var_begin, ITERATOR node_var_end,
                                    const three_dimensional_variable_array<REAL>& maxPairwisePotentials,
                                    const std::size_t chainIndex,
                                    factor_tree<FMC>* t = nullptr)
{
    std::vector<INDEX> num_labels;
    for(auto it = node_var_begin; it!=node_var_end; ++it) {
        const INDEX i = (*it);
        num_labels.push_back( this->get_number_of_labels(i) );
    }
    three_dimensional_variable_array<REAL> linearPairwisePotentials(maxPairwisePotentials);
    for(std::size_t n1=0; n1<maxPairwisePotentials.dim1(); ++n1) {
        for(std::size_t l1=0; l1<maxPairwisePotentials.dim2(n1); ++l1) {
            for(std::size_t l2=0; l2<maxPairwisePotentials.dim3(n1); ++l2) {
                linearPairwisePotentials(n1, l1, l2) = 0.0;
            }
        }
    }
    auto* chain_factor = this->lp_->template add_factor<max_chain_factor_container>(maxPairwisePotentials, linearPairwisePotentials, num_labels, chainIndex);

    INDEX pairwise_index=0;
    INDEX node1_index = 0;
    INDEX node2_index = 1;
    for(auto it = node_var_begin; std::next(it, 1)!=node_var_end; ++it, ++pairwise_index, ++node1_index, ++node2_index) {
        const INDEX i = (*it);
        const INDEX j = *std::next(it, 1);
        auto* pairwise_factor = this->get_pairwise_factor(i,j);
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

template<typename SOLVER, typename HORIZON_TRACKING_CONSTRUCTOR>
void construct_horizon_tracking_problem_on_grid(const horizon_tracking_input& input, SOLVER& solver, HORIZON_TRACKING_CONSTRUCTOR& chain_constructor)
{
    // construct mrf part
    chain_constructor.construct(input.mrf); 
    auto trees = chain_constructor.compute_forest_cover();
    for(auto& tree : trees) {
        solver.GetLP().add_tree(tree);
    }

    for(const auto& bottleneck_potential : input.bottleneck_potentials) {
        // we assume that bottleneck potential and mrf potential have same number of variables and in same order. TODO: Check for this!

        // check whether bottleneck potentials come from grid graph and if so, construct horizontal and vertical chains
        auto grid = recognize_grid(bottleneck_potential.pairwise_indices);

        // allocate space for max potentials on chains
        std::vector<three_dimensional_variable_array<REAL>> max_pairwise_potentials_on_chains(grid.number_of_chains());
        for(std::size_t i=0; i<grid.number_of_chains(); ++i) {
            const auto nodes = grid.chain(i);
            assert(nodes.size() > 0);
            std::vector<std::array<std::size_t,2>> function_table_size;
            function_table_size.reserve(nodes.size()-1);
            for(auto node_it = nodes.begin(); node_it!=std::prev(nodes.end()); ++node_it) {
                const std::size_t l1Size = bottleneck_potential.cardinality(*node_it);
                const std::size_t l2Size = bottleneck_potential.cardinality(*std::next(node_it, 1));
                function_table_size.push_back({l1Size, l2Size});
            }
            max_pairwise_potentials_on_chains[i].resize(function_table_size.begin(), function_table_size.end());
        }

        // Populate max potentials
        for(std::size_t i=0; i<bottleneck_potential.no_pairwise_factors(); ++i) {
            const auto pairwise_variables = bottleneck_potential.get_pairwise_variables(i);
            const auto [chain_number, chain_position] = grid.edge_to_chain(pairwise_variables[0], pairwise_variables[1]);
            assert(chain_number < max_pairwise_potentials_on_chains.size());

            auto pairwise_potential = bottleneck_potential.get_pairwise_potential(i);
            for(std::size_t l1=0; l1<bottleneck_potential.cardinality(pairwise_variables[0]); ++l1) {
                for(std::size_t l2=0; l2<bottleneck_potential.cardinality(pairwise_variables[1]); ++l2) {
                    max_pairwise_potentials_on_chains[chain_number](chain_position, l1, l2) = pairwise_potential(l1,l2);
                }
            }
        }

        // add chain potentials
        using FMC = typename SOLVER::FMC;
        factor_tree<FMC> tree;
        std::vector<typename std::remove_reference_t<decltype(chain_constructor)>::max_chain_factor_container*> max_chain_potentials;
        max_chain_potentials.reserve(grid.number_of_chains());

        for(std::size_t i=0; i<grid.number_of_chains(); ++i) {
            const auto nodes = grid.chain(i);
            auto* f = chain_constructor.add_max_chain(nodes.begin(), nodes.end(), max_pairwise_potentials_on_chains[i], i, &tree);
            max_chain_potentials.push_back(f);
        }

        auto* f = chain_constructor.add_max_potential(max_chain_potentials.begin(), max_chain_potentials.end(), &tree);
        solver.GetLP().add_tree(tree);
    }
}

class cardinality_traversal {
    public:
        cardinality_traversal(std::size_t num_nodes) {
            nodes.resize(num_nodes);
            nodes_covered.resize(num_nodes);
        }

        struct adjacency_list {
            std::size_t num_labels;
            std::vector<std::size_t> neighbors;
        };

        struct node {
            std::size_t node_index;
            std::size_t num_labels;
            bool operator<(const node& o) const { return node_index < o.node_index; } // or the other way around?
        };

        void add_node(std::size_t index, std::size_t num_labels) {
            nodes[index].num_labels = num_labels;
            nodes_covered[index] = false;
            if (num_labels < min_num_labels)
            {
                min_num_labels = num_labels;
                min_label_node_index = index;
            }
        }

        void add_edge(std::size_t i, std::size_t j) {
            nodes[i].neighbors.push_back(j);
            nodes[j].neighbors.push_back(i);
        }

        std::vector<std::size_t> get_traversal_order() const {
            std::vector<std::size_t> traversal_order;
            std::priority_queue<node, std::vector<node>> nodeIndexAndNodeLabels;

            nodeIndexAndNodeLabels.push({min_label_node_index, min_num_labels});
            nodes_covered[min_label_node_index] = true;

            while (!nodeIndexAndNodeLabels.empty()) {
                auto currentBestNode = nodeIndexAndNodeLabels.top();
                nodeIndexAndNodeLabels.pop();
                traversal_order.push_back(currentBestNode.node_index);

                for (const auto& currentNeighbourIndex : nodes[currentBestNode.node_index].neighbors) {
                    if (nodes_covered[currentNeighbourIndex])
                        continue;
                    
                    nodeIndexAndNodeLabels.push({currentNeighbourIndex, nodes[currentNeighbourIndex].num_labels});
                    nodes_covered[currentNeighbourIndex] = true;
                }
            }
            return traversal_order;
        }

    private:

    std::vector<adjacency_list> nodes;
    mutable std::vector<bool> nodes_covered;
    std::size_t min_label_node_index;
    std::size_t min_num_labels = std::numeric_limits<std::size_t>::max();
};

template<typename MRF_CONSTRUCTOR>
void order_nodes_by_label_space_cadinality(MRF_CONSTRUCTOR& mrf)
{
    cardinality_traversal traversal(mrf.get_number_of_variables());
    for(std::size_t i=0; i<mrf.get_number_of_variables(); ++i) {
        const auto no_labels = mrf.get_number_of_labels(i);
        traversal.add_node(i, no_labels); 
    }

    for(std::size_t p=0; p<mrf.get_number_of_pairwise_factors(); ++p) {
        const auto [i,j] = mrf.get_pairwise_variables(p);
        traversal.add_edge(i,j);
    }

    std::vector<std::size_t> order = traversal.get_traversal_order();
    std::vector<std::size_t> inverse_order(order.size());
    for(std::size_t i=0; i<order.size(); ++i) {
        inverse_order[order[i]] = i;
    }

    for (std::size_t i=0; i<order.size()-1; ++i ) {
        auto* f1 = mrf.get_unary_factor(order[i]);
        auto* f2 = mrf.get_unary_factor(order[i+1]);
        mrf.get_lp()->AddFactorRelation(f1,f2);
    }

    for(std::size_t p=0; p<mrf.get_number_of_pairwise_factors(); ++p) {
        const auto [i,j] = mrf.get_pairwise_variables(p);
        auto* f_p = mrf.get_pairwise_factor(p);
        auto f_i = mrf.get_unary_factor(i);
        auto f_j = mrf.get_unary_factor(j);
        if (inverse_order[i] < inverse_order[j]) {
            mrf.get_lp()->AddFactorRelation(f_i, f_p);
            mrf.get_lp()->AddFactorRelation(f_p, f_j);
        } else {
            mrf.get_lp()->AddFactorRelation(f_j, f_p);
            mrf.get_lp()->AddFactorRelation(f_p, f_i);
        }
    }
}



template<typename FMC, typename MAX_CHAIN_FACTOR, typename MAX_POTENTIAL_FACTOR, typename MAX_CHAIN_MAX_POTENTIAL_MESSAGE>
class max_disjoint_chain_constructor {
public:
using max_chain_factor_container = MAX_CHAIN_FACTOR;
using max_potential_factor_container = MAX_POTENTIAL_FACTOR;
using max_chain_max_potential_message_container = MAX_CHAIN_MAX_POTENTIAL_MESSAGE;

max_chain_factor_container* add_max_chain(const three_dimensional_variable_array<REAL>& linearPairwisePotentials, 
                                    const three_dimensional_variable_array<REAL>& maxPairwisePotentials,
                                    const std::vector<INDEX>& num_labels,
                                    INDEX chainIndex,
                                    factor_tree<FMC>* t = nullptr)
{
    auto* f = this->lp_->template add_factor<max_chain_factor_container>(maxPairwisePotentials, linearPairwisePotentials, num_labels, chainIndex, false);
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

/*
namespace UAIMaxPotInput {
template<typename SOLVER>
    // TODO: needed?
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
        std::vector<INDEX> num_labels;
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

            num_labels.push_back(l1Size);
            if (currentNodeIndex == numNodes - 2)
                num_labels.push_back(l2Size);
        }

        three_dimensional_variable_array<REAL> linearPairwisePotentials(functionTableSizes.begin(), functionTableSizes.end());
        three_dimensional_variable_array<REAL> maxPairwisePotentials(functionTableSizes.begin(), functionTableSizes.end());
            
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

        auto* f = chain_constructor.add_max_chain(linearPairwisePotentials, maxPairwisePotentials, num_labels, 0, &tree);
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

            three_dimensional_variable_array<REAL> maxPairwisePotentials(functionTableSizes.begin(), functionTableSizes.end());
            
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
*/
} // namespace LP_MP 

#endif // LP_MP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX
