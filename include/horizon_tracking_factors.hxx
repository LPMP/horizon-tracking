#ifndef LP_MP_HORIZON_TRACKING_FACTORS_HXX
#define LP_MP_HORIZON_TRACKING_FACTORS_HXX

#include "vector.hxx"
#include <queue>
#include <limits>

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


    class max_potential_on_chain {
        struct MaxPairwisePotential {
            REAL value;
            INDEX n1, n2; 
            INDEX l1, l2; 
            bool isAdded = false;
        };

        public:       
            max_potential_on_chain(const tensor3_variable<REAL>& maxPairwisePotentials, const tensor3_variable<REAL>& linearPairwisePotentials, const std::vector<INDEX>& numLabels)
            :
                LinearPairwisePotentials(linearPairwisePotentials),
                NumNodes(numLabels.size()),
                NumLabels(numLabels)
            {
                assert(maxPairwisePotentials.dim1() + 1 == numLabels.size()); 
                assert(maxPairwisePotentials.dim1() == linearPairwisePotentials.dim1());

                for(std::size_t n1=0; n1<maxPairwisePotentials.dim1(); ++n1) {
                    assert(maxPairwisePotentials.dim2(n1) == linearPairwisePotentials.dim2(n1) && maxPairwisePotentials.dim3(n1) == linearPairwisePotentials.dim3(n1));
                    for(std::size_t i=0; i<maxPairwisePotentials.dim2(n1); ++i) {
                        for(std::size_t j=0; j<maxPairwisePotentials.dim3(n1); ++j) {
                            MaxPairwisePotentials.push_back( {maxPairwisePotentials(n1,i,j), n1, n1+1, i, j, false} );
                        }
                    }
                }
            }


            REAL LowerBound()
            {
                Solve();
                return solutionObjective;
                // compute optimal solution and return its cost
            }

            REAL EvaluatePrimal()
            {
                return solutionObjective;
                // return cost of current solution
            }

            void MaximizePotentialAndComputePrimal() 
            {
                Solve();
                // compute optimal solution and store it
            }

            template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
            template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( LinearPairwisePotentials.data() ); }

            auto export_variables() { return std::tie( LinearPairwisePotentials.data() ); }
            void init_primal() {}
            
            template<typename ARRAY>
            void apply(ARRAY& a) const 
            { 
                std::size_t offset = 0;
                for(INDEX n=0; n<solution.size()-1; ++n) {
                    a[offset + solution[n]* LinearPairwisePotentials.dim3(n) + solution[n+1]];
                    offset += LinearPairwisePotentials.dim2(n) * LinearPairwisePotentials.dim3(n); 
                }
            }

            template<typename EXTERNAL_SOLVER>
            void construct_constraints(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v) const
            {
                assert(false);
            }
            template<typename EXTERNAL_SOLVER>
            void convert_primal(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v)
            {
                assert(false);
            } 

            tensor3_variable<REAL> LinearPairwisePotentials;
            std::vector<INDEX> solution;

        private:
            std::vector<MaxPairwisePotential> MaxPairwisePotentials;
            std::vector<INDEX> NumLabels;

            int NumNodes;
            REAL solutionObjective;

            void InsertTerminalNode()
            {
                INDEX lastNodeNumStates = NumLabels[NumNodes - 1];
                for (int l1 = 0; l1 < lastNodeNumStates; l1++)
                {
                    MaxPairwisePotential currentPot;
                    currentPot.n1 = NumNodes - 1;
                    currentPot.n2 = NumNodes;
                    currentPot.l1 = l1;
                    currentPot.l2 = 0;
                    currentPot.value = 0;
                    MaxPairwisePotentials.push_back(currentPot);
                }
                NumLabels.push_back(1); // Terminal Node has one label
                NumNodes = NumNodes + 1;
            }

            void RemoveTerminalNode()
            {
                for (int currentEdgeToRemove = MaxPairwisePotentials.size() - 1; currentEdgeToRemove >= 0; currentEdgeToRemove--)
                {
                    if (MaxPairwisePotentials[currentEdgeToRemove].n2 < NumNodes - 1)
                        break;

                    MaxPairwisePotentials.erase(MaxPairwisePotentials.begin() + currentEdgeToRemove);
                }

                NumNodes = NumNodes - 1;
                NumLabels.pop_back();
            }

            void Solve()
            {
                solution.assign(NumNodes, 0);
                InsertTerminalNode();
                REAL bestSolutionCost = INFINITY;
                std::vector<INDEX> predecessors(NumNodes);    // For each node, store the label of the previous node
                std::vector<std::vector<REAL> > distanceFromSource(NumNodes);
                REAL initialValue = 0;
                for (INDEX i = 0 ; i < NumNodes ; i++ )
                {
                    if (i > 0)
                        initialValue = std::numeric_limits<REAL>::infinity();

                    distanceFromSource[i].resize(NumLabels[i], initialValue);
                }

                std::vector<INDEX> sortingOrder = GetPairwisePotsSortingOrder(MaxPairwisePotentials);
                std::queue<INDEX> edgesToUpdate;
                for(const auto& currentEdgeToInsert : sortingOrder)
                {
                    edgesToUpdate.push(currentEdgeToInsert);
                    MaxPairwisePotentials[currentEdgeToInsert].isAdded = true;         // True allows the arc to be considered for shortest path
                    bool foundPath = UpdateDistances(edgesToUpdate, distanceFromSource, predecessors);
                    if (foundPath)
                    {
                        REAL currentCost = MaxPairwisePotentials[currentEdgeToInsert].value + distanceFromSource[NumNodes - 1][0];
                        if (currentCost < bestSolutionCost)
                        {
                            bestSolutionCost = currentCost;
                            for (int currentNodeToBackTrack = NumNodes - 2; currentNodeToBackTrack >= 0; currentNodeToBackTrack--)
                            {
                                solution[currentNodeToBackTrack] = predecessors[currentNodeToBackTrack + 1];
                            }
                        }
                    }
                }
                solutionObjective = bestSolutionCost;
                RemoveTerminalNode();
            }

            bool UpdateDistances(std::queue<INDEX>& edgesToUdate, std::vector<std::vector<REAL> >& distanceFromSource, std::vector<INDEX>& predecessors)
            {
                bool reachedTerminal = false;
                while(!edgesToUdate.empty())
                {
                    INDEX currentEdge = edgesToUdate.front();
                    edgesToUdate.pop();
                    auto currentMaxPot = MaxPairwisePotentials[currentEdge];
                    
                    INDEX n1 = currentMaxPot.n1;
                    INDEX n2 = currentMaxPot.n2;
                    INDEX l1 = currentMaxPot.l1;
                    INDEX l2 = currentMaxPot.l2;

                    REAL currentDistanceTon2l2 = distanceFromSource[n2][l2];
                    REAL currentLinearPot = 0;
                    if (n2 < NumNodes - 1) // As LinearPairwisePotentials does not contain potentials from last node to terminal node
                        currentLinearPot = LinearPairwisePotentials(n1,l1, l2);

                    REAL offeredDistanceTon2l2 = distanceFromSource[n1][l1] + currentLinearPot;
                    if (offeredDistanceTon2l2 >= currentDistanceTon2l2)
                        continue;

                    distanceFromSource[n2][l2] = offeredDistanceTon2l2;
                    predecessors[n2] = l1;
                    
                    if (n2 == NumNodes - 1)
                    {
                        reachedTerminal = true;
                        continue;
                    }

                    INDEX n3 = n2 + 1;

                    // The distance of n2, l2 has been updated so add all of its immediate children to the queue to be inspected.
                    INDEX firstEdgeToConsider = currentEdge + (NumLabels[n2] - 1 - l2) + (NumLabels[n1] - 1 - l1) * NumLabels[n2] + l2 * NumLabels[n3] + 1;
                    for (INDEX l3 = 0, currentEdgeToConsider = firstEdgeToConsider; l3 < NumLabels[n2]; l3++, currentEdgeToConsider++)
                    {
                        // Do not consider this potential as it has not been added through sorting yet.
                        if (!MaxPairwisePotentials[currentEdgeToConsider].isAdded)
                            continue;
                        
                        edgesToUdate.push(currentEdgeToConsider);
                    }
                }
                return reachedTerminal;
            }


            std::vector<INDEX> GetPairwisePotsSortingOrder(const std::vector<MaxPairwisePotential>& pots) 
            {
                std::vector<INDEX> idx(pots.size());
                std::iota(idx.begin(), idx.end(), 0);

                std::sort(idx.begin(), idx.end(),
                    [&pots](INDEX i1, INDEX i2) {return pots[i1].value < pots[i2].value;});

                return idx;
            }
    };

class unary_max_potential_on_chain_message {
public:
    template<typename FACTOR, typename MSG>
    void RepamRight(FACTOR& r, const MSG& msgs)
    {
        if(variable < r.LinearPairwisePotentials.dim1()) {
            for(std::size_t i=0; i<r.LinearPairwisePotentials.dim2(variable); ++i) {
                for(std::size_t j=0; j<r.LinearPairwisePotentials.dim3(variable); ++j) {
                    r.LinearPairwisePotentials(variable,i,j) += msgs[i];
                }
            }
        } else {
            for(std::size_t i=0; i<r.LinearPairwisePotentials.dim2(variable-1); ++i) {
                for(std::size_t j=0; j<r.LinearPairwisePotentials.dim3(variable-1); ++j) {
                    r.LinearPairwisePotentials(variable,i,j) += msgs[j];
                }
            } 
        }
    }

    template<typename FACTOR, typename MSG>
    void RepamLeft(FACTOR& l, const MSG& msgs)
    {
        for(INDEX i=0; i<l.size(); ++i) {
            l[i] += msgs[i];
        } 
    }

    template<typename LEFT_FACTOR, typename MSG>
    void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
    {
        msg -= omega*l;
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0)
    {
        assert(false);
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
    {
        if(l.primal() < l.size()) {
            const bool changed = (l.primal() != r.solution[variable]);
            l.primal() = r.solution[variable];
            return changed;
        } else {
            return false;
        }
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
    {
        if(r.primal() < l.size()) {
            const bool changed = (l.primal() != r.solution[variable]);
            r.solution[variable] = l.primal();
            return changed;
        } else {
            return false;
        }
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
        return l.primal() == r.solution[variable];
    } 

private:
    const std::size_t variable;
};



// }


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
