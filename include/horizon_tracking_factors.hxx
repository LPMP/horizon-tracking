#ifndef LP_MP_HORIZON_TRACKING_FACTORS_HXX
#define LP_MP_HORIZON_TRACKING_FACTORS_HXX

#include "vector.hxx"
#include "mrf_problem_construction.hxx"

namespace LP_MP {

    class max_factor {
        public:
            template<typename ITERATOR>
            max_factor(ITERATOR size_begin, ITERATOR size_end)
            :
                max_potential(size_begin, size_end),
                linear_potential(size_begin, size_end),
                primal(std::distance(size_begin, size_end)),
                max_potential_sorted(max_potential.total_size())
            {}

            template<typename MATRIX>
            void assign_max_potential_matrix(const INDEX entry, const MATRIX& costs)
            {
                INDEX c = 0;
                for(INDEX i=0; i<costs.dim1(); ++i) {
                    for(INDEX j=0; j<costs.dim2(); ++j) {
                        max_potential[entry][c++] = costs(i,j);
                    }
                }
            }

            REAL EvaluatePrimal() const
            {
                REAL linear_cost = 0.0;
                REAL max_cost = -std::numeric_limits<REAL>::infinity();
                for(INDEX i=0; i<max_potential.size(); ++i) {
                    linear_cost += linear_potential[i][ primal[i] ];
                    max_cost = std::max(max_potential[i][ primal[i] ], max_cost);
                }
                return linear_cost + max_cost;
            }

            INDEX no_variables() const { return max_potential.size(); }
            INDEX no_labels(const INDEX i) const { return max_potential[i].size(); }

            template<typename LAMBDA>
            REAL combinatorial_search(LAMBDA f) const
            {
                std::vector<bool> variable_covered(max_potential.size(), false);
                INDEX no_variables_covered = 0; 
                std::vector<REAL> linear_cost(no_variables(), 0.0);

                REAL lb = std::numeric_limits<REAL>::infinity();
                REAL total_linear_cost = 0.0;

                for(const auto& e : max_potential_sorted) {
                    const REAL cost = e.cost;
                    const REAL var = e.var;
                    const REAL label = e.label;
                    if(linear_potential[var][label] < linear_cost[var]) {
                        linear_cost[var] = linear_potential[var][label];
                        total_linear_cost += -linear_cost[var] + linear_potential[var][label]; 
                    }
                    if(variable_covered[var] == false) {
                        variable_covered[var] = true;
                        f(lb, var, label);
                        no_variables_covered++;
                    }
                    if(cost + total_linear_cost < lb) {
                        if(no_variables_covered == no_variables()) {
                            lb = cost + total_linear_cost;
                        }
                        f(lb, var, label);
                    }
                }

                return lb; 
            }

            void MaximizePotentialAndComputePrimal()
            {
                auto record_primal = [&primal](const REAL lb, const INDEX var, const INDEX label) {
                    primal[var] = label;
                };

                combinatorial_search(record_primal);
            }

            REAL LowerBound() const
            {
                auto no_op = [](auto, auto, auto) {};
                return combinatorial_search(no_op);
            }

            void sort_max_potentials()
            {
                INDEX c=0;
                for(INDEX var=0; var<max_potential.size(); ++var) {
                    for(INDEX label=0; label<max_potential[var].size(); ++label) {
                        max_potential_sorted[c++] = {max_potential[var][label], var, label};
                    }
                }
                std::sort(max_potential_sorted.begin(), max_potential_sorted.end(), [](auto& a, auto& b) { return a.cost < b.cost; });
            }

            void init_primal()
            {
                std::fill(primal.begin(), primal.end(), std::numeric_limits<INDEX>::max());
            }

            auto export_variables() { std::tie(max_potential, linear_potential); }

            two_dim_variable_array<REAL> max_potential;
            two_dim_variable_array<REAL> linear_potential;

        private:
            mutable REAL max_potential;
            mutable bool max_potential_valid = false;
            vector<INDEX> primal; 

            struct max_potential_entry {
                REAL cost;
                INDEX variable;
                INDEX label;
            };
            vector<max_potential_entry> max_potential_sorted;

    };

    class max_potential_on_chain {
        public:
            template<typename MRF_CONSTRUCTOR, typename ITERATOR>
            void setup_mcf(MFC_CONSTRUCTOR& mrf, ITERATOR var_begin, ITERATOR var_end)
            {
                INDEX no_mini_nodes = 0;
                INDEX no_nodes = std::distance(var_begin, var_end);
                INDEX no_edges = 0;

                 // For source node:
                INDEX no_labels_first = mrf.GetNumberOfLabels(*var_begin); 
                no_edges += no_labels_first;
                no_mini_nodes += 1;

                // For nodes of the chain:
                for(auto it=var_begin; std::next(it)!=var_end; ++it) {
                    const INDEX i = *it;
                    const INDEX j = *std::next(it);
                    no_mini_nodes += mrf.GetNumberOfLabels(j);
                    no_edges += mrf.GetNumberOfLabels(i) * mrf.GetNumberOfLabels(j);
                }

                // For terminal node:
                INDEX no_labels_end = mrf.GetNumberOfLabels(*(var_end)); //TODO: Check var_end if end inclusive
                no_edges +=  no_labels_end;
                no_mini_nodes += 1;

                MCF::SSP mcf(no_mini_nodes, no_edges);

                INDEX mini_node_counter = 0;
                
                //Adding additional edges from source to first node:
                for (INDEX xj = 0; xj < no_labels_first; ++xj){
                    mcf.add_edge(0, xj, 0, 1, 0.0);
                }

                mini_node_counter += no_labels_first;

                for(auto it=var_begin; std::next(it)!=var_end; ++it) {
                    const INDEX i = *it;
                    const INDEX j = *std::next(it);
                    const INDEX no_labels_i = mrf.GetNumberOfLabels(i);

                    if(i<j) {
                        for(INDEX xi=0; xi<mrf.GetNumberOfLabels(i); ++xi) {
                            for(INDEX xj=0; xj<mrf.GetNumberOfLabels(j); ++xj) {
                                mcf.add_edge(mini_node_counter + xi, mini_node_counter + no_labels_i * xj, 0, 1, 0.0);
                            }
                        }
                    } else {
                        assert(j<i);
                        for(INDEX xj=0; xj<mrf.GetNumberOfLabels(j); ++xj) {
                            for(INDEX xi=0; xi<mrf.GetNumberOfLabels(i); ++xi) {
                                mcf.add_edge(mini_node_counter + xi, mini_node_counter + no_labels_i * xj, 0, 1, 0.0);
                            }
                        }
                    }
                    mini_node_counter += no_labels_i;
                }

                // Adding additional edges from last node to terminal:
                for (INDEX xi = 0; xi < no_labels_end; ++xi){
                    mcf.add_edge(mini_node_counter + xi, no_mini_nodes - 1, 0, 1, 0.0);
                }

                REAL best_solution_cost = INFINITY;
                vector<INDEX> solution;                 // Stores the optimal label for each node index
                std::map<INDEX, INDEX> predecessors;    // For each mini-node, stores the incoming arc which minimizes the distance to that mini-node
                vector<INDEX> distances_from_source(mini_node_counter, INFINITY);
                distances_from_source[0] = 0;           // Source node is 0 distance away from itself
                vector<INDEX> sorting_order = GetPairwisePotsSortingOrder(all_pairwise_pots);
                std::queue<std::pair<INDEX, INDEX>> nodes_to_update;
                for(const auto& current_edge_index : sorting_order)
                {
                    nodes_to_update.push(all_pairwise_pots[current_edge_index].index_mini_node1, current_edge_index);
                    all_pairwise_pots[current_edge_index].added = true;         // True allows the arc to be considered for shortest path
                    bool found_path = UpdateDistances(nodes_to_update, all_pairwise_pots, distances_from_source, predecessors);
                    if (found_path)
                    {
                        REAL current_solution_cost = all_pairwise_pots[current_edge_index].max_pot + distances_from_source[mini_node_counter - 1];
                        if (current_solution_cost < best_solution_cost)
                        {
                            INDEX current_mini_node_index = mini_node_counter - 1;  // Assuming that the terminal node is the ending node.

                            for (INDEX i = no_nodes - 1; i >= 0; i--)
                            {
                                solution[i] = all_pairwise_pots[predecessors[current_mini_node_index]].l1;
                                current_mini_node_index = all_pairwise_pots[predecessors[current_mini_node_index]].index_mini_node1;
                            }
                        }
                    }
                }
            }

            REAL LowerBound() const
            {

            }

            vector<INDEX> GetPairwisePotsSortingOrder(const vector<max_linear_pairwise_pot>& pots) 
            {
                vector<INDEX> idx(pots.size());
                std::iota(idx.begin(), idx.end(), 0);

                std::sort(idx.begin(), idx.end(),
                    [&pots](INDEX i1, INDEX i2) {return pots[i1].max_pot < pots[i2].max_pot;});

                return idx;
            }

            // TODO: Take unary cost into account as well
            bool UpdateDistances(std::queue<std::pair<INDEX, INDEX>>& nodes_to_update, const vector<max_linear_pairwise_pot>& all_pairwise_pots, vector<REAL>& distances_from_source, std::map<INDEX, INDEX>& predecessors)
            {
                bool reached_terminal = false;
                while(!nodes_to_update.empty())
                {
                    std::pair<INDEX, INDEX> currentPair = nodes_to_update.front();
                    nodes_to_update.pop();
                    INDEX node_to_update = currentPair.first;
                    INDEX edge_index_start = currentPair.second;
                    bool anyAdded = false;

                    for (INDEX edge_index = edge_index_start; edge_index < all_pairwise_pots.size(); edge_index++)
                    {
                        // Do not consider this potential as it has not been added through sorting yet.
                        if (!all_pairwise_pots[edge_index].is_added)
                            continue;

                        // Assuming that the pairwise pots are stored contiguously for each left node.
                        if (all_pairwise_pots[edge_index].index_mini_node1 != node_to_update)
                        {
                            if (!anyAdded)
                                continue;
                            break;
                        }
                        anyAdded = true;
                        max_linear_pairwise_pot current_pot = all_pairwise_pots[edge_index];
                        REAL new_potential_node_2 = distances_from_source[current_pot.index_mini_node1]
                                             + current_pot.linear_pot;
                        if (distances_from_source[current_pot.index_mini_node2] > new_potential_node_2)
                        {
                            distances_from_source[current_pot.index_mini_node2] = new_potential_node_2;
                            nodes_to_update.push(std::pair<INDEX, INDEX>(current_pot.index_mini_node2, edge_index));
                            predecessors[current_pot.index_mini_node2] = edge_index;

                            //Assuming that if we reached the last index of distances that must be the terminal node.
                            if (current_pot.index_mini_node2 == distances_from_source.size() - 1)
                                reached_terminal = true;
                        }
                    }
                }
                return reached_terminal;
            }
        private:
            struct max_linear_pairwise_pot {
                REAL max_pot;
                REAL linear_pot;
                const INDEX i,j; // pairwise potential variables
                const INDEX l1,l2; // pairwise potential labels 
                const INDEX index_mini_node1, index_mini_node2; 
                // Assuming each label is a mini-node.
                bool is_added = false;  //TODO: Set this to true, when a pairwise term is added to the chain based on sorting of max pot.
            };
            vector<max_linear_pairwise_pot> all_pairwise_pots;
    };

    class pairwise_max_factor_message {
        public:
            pairwise_max_factor_message(const INDEX _entry) : entry(_entry) {}

            template<typename FACTOR, typename MSG>
            void RepamRight(FACTOR& r, const MSG& msgs)
            {
                for(INDEX i=0; i<r[entry].size(); ++i) {
                    r.linear_potential[entry][i] += msgs[i];
                }
            }

            template<typename FACTOR, typename MSG>
            void RepamLeft(FACTOR& l, const MSG& msgs)
            {
                INDEX c=0;
                for(INDEX i=0; i<l.dim1(); ++i) {
                    for(INDEX j=0; j<l.dim2(); ++j) {
                        l.cost(i,j) += msgs[c++];
                    }
                } 
            }

            template<typename LEFT_FACTOR, typename MSG>
            void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
            {
                vector<REAL> m(l.size());
                INDEX c=0;
                for(INDEX i=0; i<l.dim1(); ++i) {
                    for(INDEX j=0; j<l.dim2(); ++j) {
                        m[c++] = l(i,j);
                    }
                }
                msg -= omega*m; 
            }

            template<typename RIGHT_FACTOR, typename MSG>
            void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0)
            {
                assert(false);
            }

            template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
            bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
            {
                const INDEX left_primal = l.primal()[0]*l.dim1() + l.primal()[1];
                if(left_primal < l.size()) {
                    const bool changed = (left_primal != r.primal[entry]);
                    l.primal()[0] = r.primal[entry] / l.dim1();
                    l.primal()[1] = r.primal[entry] % l.dim1();
                    return changed;
                } else {
                    return false;
                }
            }

            template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
            bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
            {
                const INDEX left_primal = l.primal()[0]*l.dim1() + l.primal()[1];
                if(r.primal() < r.no_labels(entry)) {
                    const bool changed = (left_primal != r.primal[entry]);
                    r.primal[entry] = left_primal;
                    return changed;
                } else {
                    return false;
                }
            }

            template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
            bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
            {
                const INDEX left_primal = l.primal()[0]*l.dim1() + l.primal()[1];
                return left_primal == r.primal[entry];
            } 


        private:
            const INDEX entry;
    };
}

#endif // LP_MP_HORIZON_TRACKING_FACTORS_HXX
