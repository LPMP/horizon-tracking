#ifndef LP_MP_HORIZON_TRACKING_FACTORS_HXX
#define LP_MP_HORIZON_TRACKING_FACTORS_HXX

#include "vector.hxx"
#include "mrf_problem_construction.hxx"

namespace LP_MP {
    /*
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

    }; */

    struct max_linear_pairwise_pot {
        REAL max_pot;
        REAL linear_pot;
        INDEX i,j; // pairwise potential variables
        INDEX l1,l2; // pairwise potential labels 
        INDEX index_mini_node1, index_mini_node2; 
        // Assuming each label is a mini-node.
        bool is_added = false;
    };

    class max_potential_on_chain {
        public:
            max_potential_on_chain(std::vector<max_linear_pairwise_pot> allPairwisePots, int numNodes, int numMiniNodes)
            {
                all_pairwise_pots = allPairwisePots;
                NumNodes = numNodes;
                NumMiniNodes = numMiniNodes;
            }

            void Solve()
            {
                solution.assign(NumNodes, 0);
                InsertTerminalNode();
                REAL best_solution_cost = INFINITY;
                std::map<INDEX, INDEX> predecessors;    // For each mini-node, stores the incoming arc which minimizes the distance to that mini-node
                std::vector<REAL> distances_from_source(NumMiniNodes, INFINITY);
                for (int i = 0; i < all_pairwise_pots.size(); i++)
                {
                    if (all_pairwise_pots[i].i > 0)
                        break;

                    distances_from_source[all_pairwise_pots[i].l1] = 0;    
                }

                std::vector<INDEX> sorting_order = GetPairwisePotsSortingOrder(all_pairwise_pots);
                std::queue<std::pair<INDEX, INDEX>> nodes_to_update;
                for(const auto& current_edge_index : sorting_order)
                {
                    nodes_to_update.push(std::pair<INDEX, INDEX>(all_pairwise_pots[current_edge_index].index_mini_node1, current_edge_index));
                    all_pairwise_pots[current_edge_index].is_added = true;         // True allows the arc to be considered for shortest path
                    bool found_path = UpdateDistances(nodes_to_update, all_pairwise_pots, distances_from_source, predecessors);
                    if (found_path)
                    {
                        REAL current_solution_cost = all_pairwise_pots[current_edge_index].max_pot + distances_from_source[NumMiniNodes - 1];
                        if (current_solution_cost < best_solution_cost)
                        {
                            best_solution_cost = current_solution_cost;
                            INDEX current_mini_node_index = NumMiniNodes - 1;  // Assuming that the terminal node is the ending node.
                            for (int i = NumNodes - 2; i >= 0; i--)
                            {
                                solution[i] = all_pairwise_pots[predecessors[current_mini_node_index]].l1;
                                current_mini_node_index = all_pairwise_pots[predecessors[current_mini_node_index]].index_mini_node1;
                            }
                        }
                    }
                }
                solutionObjective = best_solution_cost;
                RemoveTerminalNode();
            }

            std::vector<INDEX> GetSolution() const
            {
                return solution;
            }

            REAL GetSolutionObjective() const
            {
                return solutionObjective;
            }

            REAL LowerBound() const
            {

            }

        private:
            std::vector<max_linear_pairwise_pot> all_pairwise_pots;
            int NumNodes, NumMiniNodes;
            std::vector<INDEX> solution;
            REAL solutionObjective;

            void InsertTerminalNode()
            {
                INDEX lastNodeNumStates = all_pairwise_pots[all_pairwise_pots.size() - 1].l2 + 1;
                for (int l1 = 0; l1 < lastNodeNumStates; l1++)
                {
                    max_linear_pairwise_pot currentPot;
                    currentPot.i = NumNodes;
                    currentPot.j = NumNodes + 1;
                    currentPot.index_mini_node1 = NumMiniNodes - lastNodeNumStates + l1;
                    currentPot.index_mini_node2 = NumMiniNodes;
                    currentPot.l1 = l1;
                    currentPot.l2 = 0;
                    currentPot.max_pot = 0;
                    currentPot.linear_pot = 0;
                    all_pairwise_pots.push_back(currentPot);
                }
                NumNodes = NumNodes + 1;
                NumMiniNodes = NumMiniNodes + 1;
            }

            void RemoveTerminalNode()
            {
                for (int currentEdgeToRemove = all_pairwise_pots.size() - 1; currentEdgeToRemove >= 0; currentEdgeToRemove++)
                {
                    if(all_pairwise_pots[currentEdgeToRemove].j < NumNodes - 1)
                        break;

                    all_pairwise_pots.erase(all_pairwise_pots.begin() + currentEdgeToRemove);
                }

                NumNodes = NumNodes - 1;
                NumMiniNodes = NumMiniNodes - 1;
            }

            bool UpdateDistances(std::queue<std::pair<INDEX, INDEX>>& nodes_to_update, const std::vector<max_linear_pairwise_pot>& all_pairwise_pots, std::vector<REAL>& distances_from_source, std::map<INDEX, INDEX>& predecessors)
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
                            if (current_pot.index_mini_node2 == NumMiniNodes - 1)
                                reached_terminal = true;
                        }
                    }
                }
                return reached_terminal;
            }


            std::vector<INDEX> GetPairwisePotsSortingOrder(const std::vector<max_linear_pairwise_pot>& pots) 
            {
                std::vector<INDEX> idx(pots.size());
                std::iota(idx.begin(), idx.end(), 0);

                std::sort(idx.begin(), idx.end(),
                    [&pots](INDEX i1, INDEX i2) {return pots[i1].max_pot < pots[i2].max_pot;});

                return idx;
            }


    };

    /*
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
    */
}

#endif // LP_MP_HORIZON_TRACKING_FACTORS_HXX
