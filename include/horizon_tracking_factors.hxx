#ifndef LP_MP_HORIZON_TRACKING_FACTORS_HXX
#define LP_MP_HORIZON_TRACKING_FACTORS_HXX

#include "vector.hxx"
#include <queue>
#include <limits>
#include "two_dimensional_variable_array.hxx"
#include <unordered_set>
#include <cmath>

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

    struct EdgePriority
    {
        REAL value;
        INDEX index;
    };
    
    bool operator< (const EdgePriority& lhs, const EdgePriority& rhs)
    {
        return lhs.value < rhs.value;
    }

    class max_potential_on_graph {

        struct MaxPotentialElement {
            REAL value;
            INDEX edgeIndex;
            INDEX labelIndex;
            bool isAdded = false;
        };

        public:       
            max_potential_on_graph(const two_dim_variable_array<REAL>& maxPotentials, const two_dim_variable_array<REAL>& linearPotentials)
            : LinearPotentials(linearPotentials)
            {
                assert(maxPotentials.size() == linearPotentials.size());

                for(std::size_t currentEdgeIndex = 0; currentEdgeIndex < maxPotentials.size(); ++currentEdgeIndex)
                 {
                    assert(maxPotentials[currentEdgeIndex].size() == linearPotentials[currentEdgeIndex].size()); 
                    
                    for(std::size_t currentLabelIndex = 0; currentLabelIndex < maxPotentials[currentEdgeIndex].size(); ++currentLabelIndex )
                    {
                        MaxPotentials.push_back( { maxPotentials[currentEdgeIndex][currentLabelIndex], currentEdgeIndex, currentLabelIndex, false } );
                    }
                }
            }


            REAL LowerBound() const
            {
                Solve();
                return solutionObjective;
                // compute optimal solution and return its cost
            }

            REAL EvaluatePrimal() const
            {
                return solutionObjective;
                // return cost of current solution
            }
 
            void MaximizePotentialAndComputePrimal() 
            {
                Solve();
                // compute optimal solution and store it
            }

        private:
            std::vector<MaxPotentialElement> MaxPotentials;
            two_dim_variable_array<REAL> LinearPotentials;
            mutable std::vector<INDEX> solution;
            mutable REAL solutionObjective;

            void Solve() const
            {
                std::vector<INDEX> sortingOrder = GetMaxPotsSortingOrder(MaxPotentials);
                std::unordered_set<INDEX> addedEdges;
                INDEX numEdges = LinearPotentials.size();
                REAL s = 0;
                std::vector<INDEX> l(numEdges, INFINITY);
                double bestObjective = INFINITY;
                std::vector<INDEX> bestSolutionEdgeLabels(numEdges);
                std::vector<INDEX> currentSolutionEdgeLabels(numEdges);

                for(const auto& currentElementToInsert : sortingOrder)
                {
                    INDEX currentEdgeIndex = MaxPotentials[currentElementToInsert].edgeIndex;
                    INDEX currentLabelIndex = MaxPotentials[currentElementToInsert].labelIndex;
                    REAL currentLinearCost = LinearPotentials[currentEdgeIndex][currentLabelIndex];
                    REAL currentMaxCost =  MaxPotentials[currentElementToInsert].value;

                    // If the edge is not yet covered:
                    if (addedEdges.count(currentEdgeIndex) == 0)  
                    {
                        addedEdges.insert(currentEdgeIndex);
                        s += currentLinearCost;
                        l[currentEdgeIndex] = currentLinearCost;
                        currentSolutionEdgeLabels[currentEdgeIndex] = currentLabelIndex;
                    }
                    
                    // If edge has been added, but current label has lower linear cost. We have two scenarios:
                    // 1. Graph has not been covered completely in which case we want to increase our max pot. threshold anyway. Thus, if we are gaining 
                    //      an improvement in linear cost take it.
                    // 2. Graph is covered completely, but we want to see adding current label can possibly give us any benefit, which will be 
                    //    checked by the 3rd if condition.
                    if (currentLinearCost < l[currentEdgeIndex])  
                    {
                        s = s - l[currentEdgeIndex] + currentLinearCost;
                        l[currentEdgeIndex] = currentLinearCost;
                        currentSolutionEdgeLabels[currentEdgeIndex] = currentLabelIndex;
                    }

                    // Found another solution which is better than the previous one, in which case mark current solution as the best so far.
                    if (addedEdges.size() == numEdges && bestObjective > s + currentMaxCost) 
                    {
                        bestObjective = s + currentMaxCost;
                        bestSolutionEdgeLabels = currentSolutionEdgeLabels; //TODO : what to do with the arg-min i.e. the labels?
                    }
                }

                solutionObjective  = bestObjective;
                solution = bestSolutionEdgeLabels;
            }

            std::vector<INDEX> GetMaxPotsSortingOrder(const std::vector<MaxPotentialElement>& pots) const
            {
                std::vector<INDEX> idx(pots.size());
                std::iota(idx.begin(), idx.end(), 0);

                std::sort(idx.begin(), idx.end(),
                    [&pots](INDEX i1, INDEX i2) {return pots[i1].value < pots[i2].value;});

                return idx;
            }

    };

class ShortestPathTreeInChain {
        struct Node
        {
            bool parentAssigned = false;
            INDEX parentLabel;
            std::unordered_set<INDEX> childLabels;
            std::unordered_set<INDEX> possibleChildLabels;
            REAL distance = std::numeric_limits<REAL>::infinity();
            bool status = true; // 1 -> Closed.
        };

        struct SpTreeEdge
        {
            INDEX n1, l1, n2, l2;
            REAL maxPotValue;
        };

        bool IsEqual(const SpTreeEdge& lhs, const SpTreeEdge& rhs)
        {
            return lhs.n1 == rhs.n1 && lhs.n2 == rhs.n2 &&
                    lhs.l1 == rhs.l1 && lhs.l2 == rhs.l2 &&
                    lhs.maxPotValue == rhs.maxPotValue;
        }

        struct SpTreeEdgeCompByn1
        {
            bool operator()(const SpTreeEdge& lhs, const SpTreeEdge& rhs)
            {
                return lhs.n1 < rhs.n1;
            }
        };

        struct SpTreeEdgeCompByMaxPot
        {
            bool operator()(const SpTreeEdge& lhs, const SpTreeEdge& rhs)
            {
                return lhs.maxPotValue > rhs.maxPotValue;
            }
        };

    private:
        std::vector<std::vector<Node>> ChainNodes;
        std::set<SpTreeEdge, SpTreeEdgeCompByMaxPot> treeEdgeIndicesAndMaxPots;
        std::priority_queue<SpTreeEdge, std::vector<SpTreeEdge>, SpTreeEdgeCompByn1> edgesToRemove; 

    public:
        ShortestPathTreeInChain(INDEX numNodes)
        {   
            ChainNodes.resize(numNodes);
        }

        void SetLabels(INDEX nodeIndex, INDEX numLabels)
        {
            ChainNodes[nodeIndex].resize(numLabels);
            if (nodeIndex == 0)
            {
                ChainNodes[0][0].distance = 0;
                assert(numLabels == 1);
            }
        }

        void SetPossibleChildLabels(INDEX n1, INDEX l2Max)
        {
            for (INDEX l1 = 0; l1 < ChainNodes[n1].size(); l1++)
            {
                for (INDEX l2 = 0; l2 < l2Max; l2++)
                {
                    ChainNodes[n1][l1].possibleChildLabels.insert(l2);
                }
            }
        }

        SpTreeEdge GetMaxPotValueInTree() const
        {
            return *treeEdgeIndicesAndMaxPots.begin();
        }

        bool CheckParentForShortestPath(INDEX currentNode, INDEX currentLabel, INDEX previousLabel, REAL currentEdgeDistance, REAL currentMaxPotValue)
        {
            assert(currentNode > 0);
            REAL currentOfferedDistance = ChainNodes[currentNode - 1][previousLabel].distance + currentEdgeDistance;

            if (!ChainNodes[currentNode][currentLabel].parentAssigned ||
                currentOfferedDistance < ChainNodes[currentNode][currentLabel].distance)
            {
                SetParent(currentNode, currentLabel, previousLabel, currentOfferedDistance, currentMaxPotValue); 
                return true;
            }
            return false;
        }

        void SetParent(INDEX currentNode, INDEX currentLabel, INDEX parentLabel, REAL newDistance, REAL maxPotValue)
        {
            if (ChainNodes[currentNode][currentLabel].parentAssigned)
            {
                INDEX previousParent = ChainNodes[currentNode][currentLabel].parentLabel;
                if (ChainNodes[currentNode - 1][previousParent].childLabels.count(currentLabel) > 0)
                    ChainNodes[currentNode - 1][previousParent].childLabels.erase(currentLabel);
            }

            ChainNodes[currentNode][currentLabel].parentLabel = parentLabel;
            ChainNodes[currentNode][currentLabel].distance = newDistance;
            ChainNodes[currentNode - 1][parentLabel].childLabels.insert(currentLabel);
            ChainNodes[currentNode][currentLabel].parentAssigned = true;
            treeEdgeIndicesAndMaxPots.insert({currentNode - 1, parentLabel, currentNode, currentLabel, maxPotValue});
        }

        void SetStatus(INDEX n, INDEX l, bool status)
        {
            ChainNodes[n][l].status = status;
        }

        bool GetStatusOfNode(INDEX n, INDEX l)
        {
            return ChainNodes[n][l].status;
        } 

        REAL GetDistance(INDEX n, INDEX l) const
        {
            return ChainNodes[n][l].distance;
        }

        REAL IncreaseDistance(INDEX n, INDEX l, REAL value)
        {
            ChainNodes[n][l].distance += value;
        }

        INDEX GetParentLabel(INDEX n, INDEX l) const
        {
            return ChainNodes[n][l].parentLabel;
        }

        std::unordered_set<std::array<INDEX, 2>> GetLocallyAffectedNodes()
        {
            std::unordered_set<std::array<INDEX, 2>> locallyAffectedVertices;
            while (edgesToRemove.size() > 0)
            {
                SpTreeEdge currentEdgeToRemove = edgesToRemove.top();
                edgesToRemove.pop();
                if (locallyAffectedVertices.count({currentEdgeToRemove.n2, currentEdgeToRemove.l2}) > 0)    // Already explored.
                    continue;

                GetDescendants(currentEdgeToRemove.n2, currentEdgeToRemove.l2, locallyAffectedVertices);
            }
            return locallyAffectedVertices;
        }

        void GetDescendants(INDEX n, INDEX l, std::unordered_set<std::array<INDEX, 2>>& descendants)
        {
            std::queue<std::array<INDEX, 2>> toExplore;
            toExplore.push({n, l});
            while (toExplore.size() > 0)
            {
                std::array<INDEX, 2> nodeAndLabel = toExplore.front();
                descendants.insert({nodeAndLabel[0], nodeAndLabel[1]});
                toExplore.pop();
                for (const auto& currentChildLabel: ChainNodes[nodeAndLabel[0]][nodeAndLabel[1]].childLabels)
                {
                    if (descendants.count({nodeAndLabel[0] + 1, currentChildLabel}) > 0)
                        continue; // TODO: Check if Already explored

                    toExplore.push({nodeAndLabel[0] + 1, currentChildLabel});    
                }
            }
        }

        bool CheckAndPopMaxPotInTree(INDEX n1, INDEX l1, INDEX n2, INDEX l2, REAL maxPotValue)
        {
            SpTreeEdge edge{n1, l1, n2, l2, maxPotValue};
            for (auto it = treeEdgeIndicesAndMaxPots.begin(); it != treeEdgeIndicesAndMaxPots.end(); ++it)
            {
                if (IsEqual(*it, edge))
                {
                    edgesToRemove.push(edge);
                    treeEdgeIndicesAndMaxPots.erase(it);
                    if (ChainNodes[n1][l1].childLabels.count(l2) > 0)
                        ChainNodes[n1][l1].childLabels.erase(l2);

                    ChainNodes[n2][l2].parentAssigned = false;    
                    return 1;
                }
            }

            return 0;
        }

        bool isEdgeDeleted(INDEX n1, INDEX l1, INDEX l2)
        {
            return ChainNodes[n1][l1].possibleChildLabels.count(l2) == 0;
        }

        void RemovePossibleEdge(INDEX n1, INDEX l1, INDEX l2)
        {
            ChainNodes[n1][l1].possibleChildLabels.erase(l2);
        }

        bool CheckPathToTerminal()
        {
            std::queue<std::array<REAL, 2>> q;
            q.push({ChainNodes.size() - 1, 0});
            while (q.size() > 0)
            {   
                std::array<REAL, 2> e = q.front();
                q.pop();

                if (e[0] == 0)
                    return 1;

                if (!ChainNodes[e[0]][e[1]].parentAssigned)
                    return 0;

                if (e[0] > 0)
                    q.push({e[0] - 1, ChainNodes[e[0]][e[1]].parentLabel});
            }
            return 1;
        }
    };



    class max_potential_on_chain {
        struct MaxPairwisePotential {
            REAL value;
            INDEX n1, n2; 
            INDEX l1, l2; 
            bool isAdded = false;
        };

        struct AffectedVertex
        {
            INDEX n;
            INDEX l;
            INDEX parentLabel;
            REAL delta;
            REAL newDistance;
            REAL maxPot;
            bool operator<(const AffectedVertex& rhs) const
            {
                if (delta > rhs.delta)
                    return 1;

                if (delta == rhs.delta)
                    return newDistance > rhs.newDistance;

                return 0;
            }           
        };
        
        public:       
            max_potential_on_chain(const tensor3_variable<REAL>& maxPairwisePotentials, const tensor3_variable<REAL>& linearPairwisePotentials, const std::vector<INDEX>& numLabels, bool useEdgeDeletion = true)
            :
                LinearPairwisePotentials(linearPairwisePotentials),
                MaxPairwisePotentials(maxPairwisePotentials),
                NumNodes(numLabels.size()),
                NumLabels(numLabels),
                UseEdgeDeletion(useEdgeDeletion)
            {
                assert(maxPairwisePotentials.dim1() + 1 == numLabels.size()); 
                assert(maxPairwisePotentials.dim1() == linearPairwisePotentials.dim1());

                for(std::size_t n1=0; n1<maxPairwisePotentials.dim1(); ++n1) {
                    assert(maxPairwisePotentials.dim2(n1) == linearPairwisePotentials.dim2(n1) && maxPairwisePotentials.dim3(n1) == linearPairwisePotentials.dim3(n1));
                    for(std::size_t i=0; i<maxPairwisePotentials.dim2(n1); ++i) {
                        for(std::size_t j=0; j<maxPairwisePotentials.dim3(n1); ++j) {
                            MaxPotentials1D.push_back( {maxPairwisePotentials(n1,i,j), n1, n1+1, i, j, false} );
                        }
                    }
                }

                MaxPotsSortingOrder = GetPairwisePotsSortingOrder(MaxPotentials1D);
                solution_.assign(NumNodes, 0);
                InsertTerminalNode();
            }


            REAL LowerBound() const
            {
                Solve();
                return solutionObjective;
                // compute optimal solution and return its cost
            }

            REAL EvaluatePrimal() const
            {
                return solutionObjective;
                // return cost of current solution
            }

            REAL objective() const  //TODO: Duplicate to EvaluatePrimal.
            {
                // TODO: I think we dont need this now:
                // const auto *std::min_element(max_potential_marginals_.begin(), max_potential_marginals_.end(), [](auto& m1, auto& m2) { return m1[1] + m1[2] < m2[1] + m2[2]; })
                return solutionObjective;
            }

            void MaximizePotentialAndComputePrimal() 
            {
                Solve();    
                // compute optimal solution and store it
            }

            template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
            template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( LinearPairwisePotentials.data() ); }

            auto export_variables() { return std::tie( ); }
            void init_primal() {}
            
            template<typename ARRAY>
            void apply(ARRAY& a) const 
            { 
                std::size_t offset = 0;
                for(INDEX n=0; n<solution_.size()-1; ++n) {
                    a[offset + solution_[n]* LinearPairwisePotentials.dim3(n) + solution_[n+1]];
                    offset += LinearPairwisePotentials.dim2(n) * LinearPairwisePotentials.dim3(n); 
                }
            }

            template<typename EXTERNAL_SOLVER>
            void construct_constraints(EXTERNAL_SOLVER& s) const
            {
                assert(false);
            }
            template<typename EXTERNAL_SOLVER>
            void convert_primal(EXTERNAL_SOLVER& s)
            {
                assert(false);
            } 

            tensor3_variable<REAL> LinearPairwisePotentials;

            // setter
            INDEX& solution(const std::size_t i) { assert(i < solution_.size()); return solution_[i]; }
            // getter
            INDEX solution(const std::size_t i) const { assert(i < solution_.size()); return solution_[i]; }
        
            std::array<REAL,3>& max_potential_marginal(const std::size_t i) { assert(i < max_potential_marginals_.size()); return max_potential_marginals_[i]; }
            std::array<REAL,3> max_potential_marginal(const std::size_t i) const { assert(i < max_potential_marginals_.size()); return max_potential_marginals_[i]; }

        protected:
                std::vector<INDEX> NumLabels;
                mutable std::vector<INDEX> solution_;

        private:
            mutable std::vector<MaxPairwisePotential> MaxPotentials1D;
            std::vector<INDEX> MaxPotsSortingOrder;
            tensor3_variable<REAL> MaxPairwisePotentials;
            mutable std::vector<std::array<REAL,3>> max_potential_marginals_; // (i) max potential, (ii) minimum linear potential, (iii) cost of configuration TODO: Sort these 
            bool MaxPotMarginalsInitialized = false;

            // cost = max_potential_marginals_[1] + max_potential_marginals_[2] TODO: remove max potential value from objective calculation
            bool UseEdgeDeletion;

            int NumNodes;
            mutable REAL solutionObjective = std::numeric_limits<REAL>::infinity();

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
                    currentPot.value = std::numeric_limits<REAL>::min();
                    MaxPotentials1D.push_back(currentPot);
                }
                NumLabels.push_back(1); // Terminal Node has one label
                NumNodes = NumNodes + 1;
            }

            void RemoveTerminalNode()
            {
                for (int currentEdgeToRemove = MaxPotentials1D.size() - 1; currentEdgeToRemove >= 0; currentEdgeToRemove--)
                {
                    if (MaxPotentials1D[currentEdgeToRemove].n2 < NumNodes - 1)
                        break;

                    MaxPotentials1D.erase(MaxPotentials1D.begin() + currentEdgeToRemove);
                }

                NumNodes = NumNodes - 1;
                NumLabels.pop_back();
            }

            void Solve() const
            {
                if (UseEdgeDeletion)
                    SolveByEdgeDeletion();
                else
                    SolveByEdgeAddition();

                MaxPotMarginalsInitialized = true;
            }

            void SolveByEdgeAddition() const
            {
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

                INDEX currentMaxPotIndex = 0;
                for(const auto& currentEdgeToInsert : MaxPotsSortingOrder)
                {
                    MaxPotentials1D[currentEdgeToInsert].isAdded = true;         // True allows the arc to be considered for shortest path
                    bool foundPath = UpdateDistances(currentEdgeToInsert, distanceFromSource, predecessors);
                    if (foundPath)
                    {
                        REAL currentCost =  distanceFromSource[NumNodes - 1][0]; 

                        if (!MaxPotMarginalsInitialized)
                            max_potential_marginals_.push_back({MaxPotentials1D[currentEdgeToInsert].value, currentCost, 0});
                        
                        else
                        {
                            assert(MaxPotentials1D[currentEdgeToInsert].value == max_potential_marginals_[currentMaxPotIndex][0]);
                            currentCost += max_potential_marginals_[currentMaxPotIndex][2];
                            max_potential_marginals_[currentMaxPotIndex][1] = distanceFromSource[NumNodes - 1][0];
                            currentMaxPotIndex++;
                        }

                        if (currentCost < bestSolutionCost)     //TODO: So the solution we will get will not consider the max potential values, even for first iteration.
                        {
                            bestSolutionCost = currentCost;
                            for (int currentNodeToBackTrack = NumNodes - 2; currentNodeToBackTrack >= 0; currentNodeToBackTrack--)
                            {
                                solution_[currentNodeToBackTrack] = predecessors[currentNodeToBackTrack + 1];
                            }
                        }
                    }
                }

                //TODO: I think we dont need this now:
                // *std::min_element(max_potential_marginals_.begin(), max_potential_marginals_.end(), [](auto& m1, auto& m2) { return m1[1] + m1[2] < m2[1] + m2[2]; })
                solutionObjective = bestSolutionCost;
            }


            bool UpdateDistances(INDEX edgeToUpdate, std::vector<std::vector<REAL> >& distanceFromSource, std::vector<INDEX>& predecessors) const
            {
                bool reachedTerminal = false;
                std::queue<EdgePriority> pQueue;  //TODO: Priority queue probably does not offer any benefit for topological sort shortest path.
                auto currentMaxPot = MaxPotentials1D[edgeToUpdate];
                
                INDEX n1 = currentMaxPot.n1;
                INDEX n2 = currentMaxPot.n2;
                INDEX l1 = currentMaxPot.l1;
                INDEX l2 = currentMaxPot.l2;
                REAL currentLinearPot = 0;
                if (n2 < NumNodes - 1) // As LinearPairwisePotentials does not contain potentials from last node to terminal node
                        currentLinearPot = LinearPairwisePotentials(n1, l1, l2);
                
                REAL offeredDistanceTon2l2 = distanceFromSource[n1][l1] + currentLinearPot;

                pQueue.push(EdgePriority{offeredDistanceTon2l2, edgeToUpdate});

                while(!pQueue.empty())
                {
                    EdgePriority currentEdgeStruct = pQueue.front();
                    pQueue.pop();
                    INDEX currentEdge = currentEdgeStruct.index;
                    REAL offeredDistanceTon2l2 = currentEdgeStruct.value;
                    auto currentMaxPot = MaxPotentials1D[currentEdge];
                    
                    INDEX n1 = currentMaxPot.n1;
                    INDEX n2 = currentMaxPot.n2;
                    INDEX l1 = currentMaxPot.l1;
                    INDEX l2 = currentMaxPot.l2;

                    REAL currentDistanceTon2l2 = distanceFromSource[n2][l2];
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
                    for (INDEX l3 = 0, currentEdgeToConsider = firstEdgeToConsider; l3 < NumLabels[n3]; l3++, currentEdgeToConsider++)
                    {
                        // Do not consider this potential as it has not been added through sorting yet.
                        if (!MaxPotentials1D[currentEdgeToConsider].isAdded)
                            continue;
                        
                        auto childMaxPot = MaxPotentials1D[currentEdgeToConsider];
                        assert(childMaxPot.l1 == l2);
                        assert(childMaxPot.l2 == l3); // Might fail if some of the pairwise potentials are not present thus causing jumps!
                        REAL currentLinearPot = 0;
                        if (n3 < NumNodes - 1) // As LinearPairwisePotentials does not contain potentials from last node to terminal node
                            currentLinearPot = LinearPairwisePotentials(n2, l2, l3);
                        
                        REAL offeredDistanceTon3l3 = offeredDistanceTon2l2 + currentLinearPot;

                        pQueue.push(EdgePriority{offeredDistanceTon3l3, currentEdgeToConsider});
                    }
                }
                return reachedTerminal;
            }

            std::vector<INDEX> GetPairwisePotsSortingOrder(const std::vector<MaxPairwisePotential>& pots) const
            {
                std::vector<INDEX> idx(pots.size());
                std::iota(idx.begin(), idx.end(), 0);
                std::sort(idx.begin(), idx.end(), [&pots](INDEX i1, INDEX i2) {return pots[i1].value < pots[i2].value;});
                return idx;
            }

            ShortestPathTreeInChain FindAndInitializeSPTree() const
            {
                ShortestPathTreeInChain spTree(1 + NumNodes);   // Also includes source node
                spTree.SetLabels(0, 1);
                spTree.SetPossibleChildLabels(0, NumLabels[0]);

                for (INDEX i = 0; i < NumNodes; i++)
                {
                    spTree.SetLabels(i + 1, NumLabels[i]);

                    if (i < NumNodes - 1)
                        spTree.SetPossibleChildLabels(i + 1, NumLabels[i + 1]);
                }
                
                for (INDEX n = 0; n < NumNodes; n++)
                {
                    for (INDEX l = 0; l < NumLabels[n]; l++)
                    {
                        INDEX prevNumLabels = 1;
                    if (n > 0)
                            prevNumLabels = NumLabels[n - 1];

                        for (INDEX prevLabel = 0; prevLabel < prevNumLabels; prevLabel++)
                        {
                            REAL currentLinearPot = 0;
                            REAL currentMaxPot = std::numeric_limits<REAL>::min();
                            if (n < NumNodes - 1 && n > 0)   // Source and Terminal node potentials are not present in it.
                            {
                                currentLinearPot = LinearPairwisePotentials(n - 1, prevLabel, l);
                                currentMaxPot = MaxPairwisePotentials(n - 1, prevLabel, l);
                            }

                            spTree.CheckParentForShortestPath(n + 1, l, prevLabel, currentLinearPot, currentMaxPot);
                        }
                    }
                }

                return spTree;
            }

            void SolveByEdgeDeletion() const
            {
                ShortestPathTreeInChain spTree = FindAndInitializeSPTree();
                REAL treeMaxPotValue = spTree.GetMaxPotValueInTree().maxPotValue;
                INDEX currentMaxPotIndex = 0;
                if (!MaxPotMarginalsInitialized)
                    marginals_.insert(marginals_.begin(), {treeMaxPotValue, spTree.GetDistance(NumNodes, 0), 0});
                else
                {
                    assert(marginals_[currentMaxPotIndex][0] == treeMaxPotValue);
                    marginals_[currentMaxPotIndex][1] = spTree.GetDistance(NumNodes, 0);
                }

                currentCost = marginals_[currentMaxPotIndex][2] + spTree.GetDistance(NumNodes, 0);
                if (currentCost < solutionObjective)
                {
                    solutionObjective = currentCost; 
                    INDEX currentLabel = 0; // for terminal node;
                    for (int currentNodeToBackTrack = NumNodes - 1; currentNodeToBackTrack >= 1; currentNodeToBackTrack--) // Exclude source node.
                    {
                        solution_[currentNodeToBackTrack - 1] = spTree.GetParentLabel(currentNodeToBackTrack + 1, currentLabel);
                        currentLabel = solution_[currentNodeToBackTrack];
                    }
                }

                currentMaxPotIndex++;

                for (int i = MaxPotsSortingOrder.size() - 1; i >= 0; i--)
                {   
                    auto currentMaxPotEdge = MaxPotentials1D[sortingOrder[i]];
                    bool wasTreeEdge = spTree.CheckAndPopMaxPotInTree(
                        currentMaxPotEdge.n1 + 1, currentMaxPotEdge.l1, currentMaxPotEdge.n2 + 1, currentMaxPotEdge.l2, currentMaxPotEdge.value); // +1 due to source node.
                    spTree.RemovePossibleEdge(currentMaxPotEdge.n1 + 1, currentMaxPotEdge.l1, currentMaxPotEdge.l2);

                    if (wasTreeEdge)
                    {
                        std::priority_queue<AffectedVertex> locallyAffectedPQueue;
                        std::unordered_set<std::array<INDEX, 2>> locallyAffectedVertices = spTree.GetLocallyAffectedNodes();
                        for (auto const& currentAffectedNode : locallyAffectedVertices)
                        {   
                            spTree.SetStatus(currentAffectedNode[0], currentAffectedNode[1], 0);    //mark as open
                            INDEX bestParentLabel;
                            REAL bestParentMaxPotValue; 
                            REAL bestParentDistance = std::numeric_limits<REAL>::infinity();
                            for (INDEX prevLabel = 0; prevLabel < NumLabels[currentAffectedNode[0] - 2]; prevLabel++)
                            {
                                REAL currentLinearPot = 0;
                                REAL currentMaxPot = std::numeric_limits<REAL>::min();
                                if (currentAffectedNode[0] - 2 < NumNodes - 2)
                                {
                                    currentLinearPot = LinearPairwisePotentials(currentAffectedNode[0] - 2, prevLabel, currentAffectedNode[1]);
                                    currentMaxPot = MaxPairwisePotentials(currentAffectedNode[0] - 2, prevLabel, currentAffectedNode[1]);
                                }

                                if (locallyAffectedVertices.count({currentAffectedNode[0] - 1, prevLabel}) > 0 ||
                                    spTree.isEdgeDeleted(currentAffectedNode[0] - 1, prevLabel, currentAffectedNode[1])
                                    || !spTree.GetStatusOfNode(currentAffectedNode[0] - 1, prevLabel))
                                    continue;
 
                                REAL db = spTree.GetDistance(currentAffectedNode[0] - 1, prevLabel); // parent distance

                                REAL newDistance = db + currentLinearPot;
                                if (newDistance < bestParentDistance)
                                {
                                    bestParentLabel = prevLabel;
                                    bestParentDistance = newDistance;
                                    bestParentMaxPotValue = currentMaxPot;
                                }
                            }

                            if (!std::isinf(bestParentDistance))
                            {
                                REAL currentDistance = spTree.GetDistance(currentAffectedNode[0], currentAffectedNode[1]);
                                locallyAffectedPQueue.push({currentAffectedNode[0], currentAffectedNode[1], bestParentLabel, bestParentDistance - currentDistance, bestParentDistance, bestParentMaxPotValue});
                            }
                        }

                        //STEP 3:
                        std::unordered_set<std::array<INDEX, 2>> toRemoveFromQueue;
                        while (locallyAffectedPQueue.size() > 0)
                        {
                            AffectedVertex bestVertex = locallyAffectedPQueue.top();
                            locallyAffectedPQueue.pop();

                            if (toRemoveFromQueue.count({bestVertex.n, bestVertex.l}) > 0)
                                continue;
                            
                            spTree.SetParent(bestVertex.n, bestVertex.l, bestVertex.parentLabel, bestVertex.newDistance, bestVertex.maxPot);

                            std::unordered_set<std::array<INDEX, 2>> descendants;
                            spTree.GetDescendants(bestVertex.n, bestVertex.l, descendants);
                            for (const auto& currentDesc : descendants)
                            {
                                if (currentDesc[0] != bestVertex.n || currentDesc[1] != bestVertex.l)
                                    spTree.IncreaseDistance(currentDesc[0], currentDesc[1], bestVertex.delta);

                                spTree.SetStatus(currentDesc[0], currentDesc[1], 1); 
                                toRemoveFromQueue.insert({currentDesc[0], currentDesc[1]});
                            }
                            // Relax outgoing edges of just consolidated vertices.
                            for (const auto& currentDesc : descendants)
                            {
                                INDEX nextNode = currentDesc[0] + 1;
                                if (nextNode > NumNodes)
                                    continue;  
                                
                                for (INDEX nextLabel = 0; nextLabel < NumLabels[nextNode - 1]; nextLabel++)
                                {
                                    REAL linearPotentialValue = 0;
                                    if (currentDesc[0] - 1 < NumNodes - 2)
                                    {
                                        REAL linearPotentialValue = LinearPairwisePotentials(currentDesc[0] - 1, currentDesc[1], nextLabel);
                                    }
                                    
                                    if (spTree.GetStatusOfNode(nextNode, nextLabel) ||
                                    spTree.isEdgeDeleted(currentDesc[0], currentDesc[1], nextLabel))
                                        continue;
                                    
                                    REAL newDistance = spTree.GetDistance(currentDesc[0], currentDesc[1]) + linearPotentialValue;
                                    REAL delta = newDistance - spTree.GetDistance(nextNode, nextLabel);
                                    locallyAffectedPQueue.push({nextNode, nextLabel, currentDesc[1], delta, newDistance});
                                }
                            }                            
                        }

                        if (!spTree.CheckPathToTerminal())
                            return;

                        REAL treeMaxPotValue = spTree.GetMaxPotValueInTree().maxPotValue;
                        
                        if (!MaxPotMarginalsInitialized)
                            marginals_.insert(marginals_.begin(), {treeMaxPotValue, spTree.GetDistance(NumNodes, 1, 0)}); // push front to maintain the increasing max pot order.
                        else
                        {
                            assert(marginals_[currentMaxPotIndex][0] == treeMaxPotValue);
                            marginals_[currentMaxPotIndex][1] = spTree.GetDistance(NumNodes, 1, 0);
                        }

                        currentCost = marginals_[currentMaxPotIndex][2] + spTree.GetDistance(NumNodes, 0);
                        if (currentCost < solutionObjective)
                        {
                            solutionObjective = currentCost; 
                            INDEX currentLabel = 0; // for terminal node;
                            for (int currentNodeToBackTrack = NumNodes - 1; currentNodeToBackTrack >= 1; currentNodeToBackTrack--) // Exclude source node.
                            {
                                solution_[currentNodeToBackTrack - 1] = spTree.GetParentLabel(currentNodeToBackTrack + 1, currentLabel); //TODO: Verify this!
                                currentLabel = solution_[currentNodeToBackTrack];
                            }
                        }

                        currentMaxPotIndex++;                        
                    }
                }
            }
    };

class LabelStateSpace {
    private:
        REAL MaxPotentialUpperBound;
        REAL MaxPotentialLowerBound;
        bool lowerBoundAdded = false;
        bool cleared = false;
        std::map<REAL, REAL> potentials; 

    public:
        LabelStateSpace()
        {}

        LabelStateSpace(REAL maxPotLowerBound, REAL maxPotUpperBound) :
        MaxPotentialLowerBound(maxPotLowerBound), MaxPotentialUpperBound(maxPotUpperBound)
        {
            assert(MaxPotentialLowerBound <= MaxPotentialUpperBound);
        }

        const std::map<REAL, REAL>& GetPotentials() const
        {
            assert(!cleared);
            return potentials;
        }

        void CheckAndInsert(REAL maxPot, REAL linearPot)
        {
            assert(!cleared);
            if (maxPot > MaxPotentialUpperBound)
                return; // This edge has value greater than ub, and we hope to find something <= ub, thus discarding this.

            else if (maxPot < MaxPotentialLowerBound)
            {
                if (!lowerBoundAdded)
                {
                    potentials.insert(std::pair<REAL, REAL>(MaxPotentialLowerBound, linearPot));
                    lowerBoundAdded = true;
                }
                else
                {
                    auto itr = potentials.find(MaxPotentialLowerBound);
                    assert(itr != potentials.end());
                    itr->second = fmin(linearPot, itr->second);
                }
            }

            else
            {
                auto lbKey = potentials.lower_bound(maxPot);
                if (lbKey != potentials.end() || potentials.size() > 0)
                {
                    if (lbKey->first == maxPot) // Current key already exists in the map
                    {
                        lbKey->second = fmin(linearPot, lbKey->second);
                        return;
                    }
                    else if (lbKey != potentials.begin()) // There is a key before the current element which needs to be checked with current max pot.
                    {
                        assert(maxPot > std::prev(lbKey)->first);
                        if (!(linearPot < std::prev(lbKey)->second)) // Inserting has no benefit as linearPot is not strictly less than previous linear pot.
                            return;
                    }
                }
                potentials.insert(lbKey, std::pair<REAL, REAL>(maxPot, linearPot)); // Insert with hint 
                if (maxPot == MaxPotentialLowerBound)
                    lowerBoundAdded = true;
            }
        }

        REAL GetMaxPotLowerBound() const
        {
            return MaxPotentialLowerBound;
        }

        REAL GetMaxPotUpperBound() const
        {
            return MaxPotentialUpperBound;
        }

        void TakeUnion(LabelStateSpace newStateSpace)
        {
            assert(!isCleared());
            assert(!newStateSpace.isCleared());
            const auto& newPotentials = newStateSpace.GetPotentials();
            for (auto const& currentPot : newPotentials)
            {
                CheckAndInsert(currentPot.first, currentPot.second);
            }
        }

        static LabelStateSpace MergeStateSpaceFromDifferentNeighbours(LabelStateSpace s1, LabelStateSpace s2)
        {
            assert(!s1.isCleared());
            assert(!s2.isCleared());
            const auto& s1Potentials = s1.GetPotentials();
            const auto& s2Potentials = s2.GetPotentials();

            if (s1Potentials.size() == 0)
                return s2;

            if (s2Potentials.size() == 0)
                return s1;

            LabelStateSpace merged = LabelStateSpace(s1.GetMaxPotLowerBound(), s1.GetMaxPotUpperBound());
            assert(s1.GetMaxPotLowerBound() == s2.GetMaxPotLowerBound());
            assert(s1.GetMaxPotUpperBound() == s2.GetMaxPotUpperBound());
            for (auto const& currentS1Pot : s1Potentials)
            {
                for (auto const& currentS2Pot : s2Potentials)
                {
                    merged.CheckAndInsert(fmax(currentS1Pot.first, currentS2Pot.first), currentS1Pot.second + currentS2Pot.second);
                    //TODO: Probably can be made more efficient, also validate the logic.
                }
            }
            return merged;
        }

        void ClearPotentials()
        {
            potentials.clear();
            cleared = true;
        }

        bool isCleared() const
        {
            return cleared;
        }

};

class max_potential_on_tree {
    public:       
        max_potential_on_tree(const tensor3_variable<REAL>& maxPairwisePotentials, const tensor3_variable<REAL>& linearPairwisePotentials,
         const std::vector<INDEX>& numLabels, const std::vector<std::array<INDEX, 2>>& messagePassingSchedule, const std::vector<INDEX>& numEdgesForNode)
        :
            LinearPairwisePotentials(linearPairwisePotentials),
            MaxPairwisePotentials(maxPairwisePotentials),
            NumLabels(numLabels),
            NumNodes(numLabels.size()),
            NumEdgesForNode(numEdgesForNode),
            MessagePassingSchedule(messagePassingSchedule)
        {
            assert(maxPairwisePotentials.dim1() + 1 == numLabels.size()); 
            assert(maxPairwisePotentials.dim1() == linearPairwisePotentials.dim1());
        }

        // Call this function whenever linear/max potentials get changed so the bounds need to be recomputed.
        void PotsChanged()
        {
            boundsDirty = true;
        }

        REAL GetMaxPotLowerBound()
        {
            if (boundsDirty)
            {
                RecomputeBounds();
                boundsDirty = false;
            }
            return MaxPotentialLowerBound;
        }

        REAL GetMaxPotUpperBound()
        {
            if (boundsDirty)
            {
                RecomputeBounds();
                boundsDirty = false;
            }
            return MaxPotentialUpperBound;
        }

        // Returns the marginals of the root node.
        std::vector<std::array<REAL, 2>> ComputeMarginalsForBothPotentials(REAL maxPotLowerBound, REAL maxPotUpperBound) const
        {
            std::vector<std::vector<LabelStateSpace>> messages(NumNodes);
            for (INDEX i = 0; i < NumNodes; i++)
            {
                messages[i].resize(NumLabels[i]);
                for (INDEX l = 0; l < NumLabels[i]; l++) 
                    messages[i][l] = LabelStateSpace(maxPotLowerBound, maxPotUpperBound);
            }

            std::vector<INDEX> totalSentAndReceivedMessages(NumNodes, 0);
            std::vector<bool> messageReceived(NumNodes, false);
            INDEX edgeIndex = 0;
            for (const auto & currentEdge : MessagePassingSchedule) // TODO: Assuming the edge index is the n1 of current edge.
            {
                INDEX tail = currentEdge[0];
                INDEX head = currentEdge[1];
                bool isTailLeaf = false;
                if (!messageReceived[tail]) // Only leaf nodes can send a message (only one) and only before receiving any.
                {
                    isTailLeaf = true;
                    assert(totalSentAndReceivedMessages[tail] == 0);
                }
                for (INDEX lh = 0; lh < NumLabels[head]; lh++)
                {
                    LabelStateSpace lhStateSpace = GetHeadLabelStateSpaceFromCurrentTail(messages[tail], edgeIndex, lh, maxPotLowerBound, maxPotUpperBound, head < tail, isTailLeaf);
                    LabelStateSpace lhPrevStateSpace = messages[head][lh];
                    messages[head][lh] = LabelStateSpace::MergeStateSpaceFromDifferentNeighbours(lhStateSpace, lhPrevStateSpace);
                }

                totalSentAndReceivedMessages[head]++;
                totalSentAndReceivedMessages[tail]++;

                // Clear the potentials of all nodes except root node to conserve memory as they have 'messaged' their beliefs already.
                if (currentEdge != MessagePassingSchedule.back())
                {
                    if (totalSentAndReceivedMessages[tail] == NumEdgesForNode[tail])
                    {   
                        for (INDEX lt = 0; lt < NumLabels[tail]; lt++)
                            messages[tail][lt].ClearPotentials();
                    }
                    if (totalSentAndReceivedMessages[head] == NumEdgesForNode[head])
                    {
                        for (INDEX lh = 0; lh < NumLabels[head]; lh++)
                            messages[head][lh].ClearPotentials();
                    }
                }
                edgeIndex++;
            }

            // Merge the marginals of all the labels of root node.
            const auto& lastEdge = MessagePassingSchedule.back();
            INDEX rootNode = lastEdge[1];
            LabelStateSpace mergedRootNodeStateSpace;
            for (INDEX l = 0; l < NumLabels[rootNode]; l++)
            {
                mergedRootNodeStateSpace.TakeUnion(messages[rootNode][l]);
            }
            auto rootPots = mergedRootNodeStateSpace.GetPotentials();
            std::vector<std::array<REAL, 2>> rootPotsArray;
            for (const auto &currentPair : rootPots)
            {
                rootPotsArray.push_back({currentPair.first, currentPair.second});
            }
            return rootPotsArray;
        }

        // setter
        INDEX& solution(const std::size_t i) { assert(i < solution_.size()); return solution_[i]; }
        // getter
        INDEX solution(const std::size_t i) const { assert(i < solution_.size()); return solution_[i]; }

    private:
        std::vector<INDEX> solution_;
        tensor3_variable<REAL> MaxPairwisePotentials;
        tensor3_variable<REAL> LinearPairwisePotentials;
        INDEX NumNodes;
        std::vector<INDEX> NumLabels;
        std::vector<INDEX> NumEdgesForNode;
        
        std::vector<std::array<INDEX, 2>> MessagePassingSchedule;
        mutable REAL MaxPotentialLowerBound;    // Computed by max potential message passing.
        mutable REAL MaxPotentialUpperBound;
        mutable REAL LinearPotentialLowerBound; // Computed by conventional message passing.
        mutable REAL LinearPotentialUpperBound; // TODO: Can be computed from max potential message passing and by used to prune paths longer than this bound.

        bool boundsDirty = true;

        void RecomputeBounds() const
        {
            std::array<REAL, 2> bounds = MessagePassingForOnePotential(LinearPairwisePotentials, MaxPairwisePotentials, true);
            LinearPotentialLowerBound = bounds[0];
            MaxPotentialUpperBound = bounds[1];

            bounds = MessagePassingForOnePotential(MaxPairwisePotentials, LinearPairwisePotentials, false);
            MaxPotentialLowerBound = bounds[0];
            LinearPotentialUpperBound = bounds[1]; //TODO: If this is not useful remove it and compute max potential lb only once and store it, as it will not change.
        }

        LabelStateSpace GetHeadLabelStateSpaceFromCurrentTail(const std::vector<LabelStateSpace>& tailMessages, INDEX edgeIndex, INDEX lh, REAL maxPotLB, REAL maxPotUB, bool isReverse, bool isTailLeaf = false) const 
        {
            LabelStateSpace lhStateSpace(maxPotLB, maxPotUB);
            for (INDEX lt = 0; lt < tailMessages.size(); lt++)
            {
                REAL edgeMaxPot = !isReverse ? MaxPairwisePotentials(edgeIndex, lt, lh) : MaxPairwisePotentials(edgeIndex, lh, lt);
                REAL edgeLinearPot = !isReverse ? LinearPairwisePotentials(edgeIndex, lt, lh) : LinearPairwisePotentials(edgeIndex, lh, lt);
                if (!isTailLeaf)
                {
                    for (const auto& currentMessageTolt: tailMessages[lt].GetPotentials()) // Iterator over all messages incoming to l1.
                    {
                        lhStateSpace.CheckAndInsert(fmax(currentMessageTolt.second, edgeMaxPot), currentMessageTolt.first + edgeLinearPot);
                    }
                }
                else
                {
                    lhStateSpace.CheckAndInsert(edgeMaxPot, edgeLinearPot);                    
                }
            }
            return lhStateSpace;
        }

        std::array<REAL, 2> MessagePassingForOnePotential(const tensor3_variable<REAL>& mainPairwisePots, const tensor3_variable<REAL>& otherPairwisePots, bool doConventional) const
        {
            std::vector<std::vector<REAL>> mainMessages(NumNodes);
            std::vector<std::vector<REAL>> otherMessages(NumNodes);
            std::vector<bool> messageSent(NumNodes, false);
            std::vector<bool> messageReceived(NumNodes, false);
            for (INDEX i = 0; i < NumNodes; i++)
            {
                mainMessages[i].resize(NumLabels[i], 0);
                otherMessages[i].resize(NumLabels[i], 0);
            }

            INDEX edgeIndex = 0;
            for (const auto & currentEdge : MessagePassingSchedule) // TODO: Assuming the edge index is the n1 of current edge.
            {
                INDEX tail = currentEdge[0];
                INDEX head = currentEdge[1];
                if (!messageReceived[tail]) // Only leaf nodes can send a message (only one) and only before receiving any.
                {
                    assert(!messageSent[tail]);
                }
                for (INDEX lh = 0; lh < mainMessages[head].size(); lh++)
                {
                    std::array<REAL, 2> minMessage = ComputeMessageValue(mainMessages[tail], otherMessages[tail], edgeIndex, lh, mainPairwisePots, otherPairwisePots, doConventional, head < tail);
                    if (doConventional)
                    {
                        mainMessages[head][lh] += minMessage[0];
                        otherMessages[head][lh] = fmax(otherMessages[head][lh], minMessage[1]);
                    }
                    else
                    {
                        mainMessages[head][lh] = fmax(mainMessages[head][lh], minMessage[0]);                                            
                        otherMessages[head][lh] += minMessage[1];
                    }
                }
                messageSent[tail] = true;
                messageReceived[head] = true;
                edgeIndex++;
            }

            const auto& lastEdge = MessagePassingSchedule.back();
            INDEX rootNode = lastEdge[1];
            assert(messageReceived[rootNode]);
            assert(!messageSent[rootNode]);
            REAL mainBound = INFINITY;
            REAL otherBound = INFINITY;

            for (int l = 0; l < NumLabels[rootNode]; l++)
            {
                if (mainMessages[rootNode][l] < mainBound)
                {
                    mainBound = mainMessages[rootNode][l];
                    otherBound = otherMessages[rootNode][l];
                }
            }

            return std::array<REAL, 2>({mainBound, otherBound});
        }

        std::array<REAL, 2> ComputeMessageValue(const std::vector<REAL>& tailMainMessages, const std::vector<REAL>& tailOtherMessages, INDEX edgeIndex, INDEX lh, 
        const tensor3_variable<REAL>& mainPairwisePots, const tensor3_variable<REAL>& otherPairwisePots, bool doConventional, bool isReverse) const
        {
            REAL bestMainMessage = INFINITY;
            REAL bestOtherMessage = INFINITY;
            INDEX bestlt;
            for (INDEX lt = 0; lt < tailMainMessages.size(); lt++)
            {
                REAL currentEdgePot = !isReverse ? mainPairwisePots(edgeIndex, lt, lh) :  mainPairwisePots(edgeIndex, lh, lt);
                REAL currentIncomingMsg = tailMainMessages[lt];

                if (doConventional)
                {
                    REAL msg = currentIncomingMsg + currentEdgePot;
                    if (msg < bestMainMessage)
                    {
                        bestMainMessage = msg;
                        bestlt = lt;
                    }
                }

                else
                {
                    REAL msg = fmax(currentIncomingMsg, currentEdgePot);
                    if (msg < bestMainMessage)
                    {
                        bestMainMessage = msg;
                        bestlt = lt;
                    }
                }
            }

            if (doConventional)
                bestOtherMessage = fmax(!isReverse ? otherPairwisePots(edgeIndex, bestlt, lh) :  otherPairwisePots(edgeIndex, lh, bestlt), tailOtherMessages[bestlt]);

            else
                bestOtherMessage = (!isReverse ? otherPairwisePots(edgeIndex, bestlt, lh) : otherPairwisePots(edgeIndex, lh, bestlt)) + tailOtherMessages[bestlt];

            return std::array<REAL, 2>({bestMainMessage, bestOtherMessage});
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
                for(INDEX i=0; i<r.LinearPairwisePotentials[entry].size(); ++i) {
                    r.LinearPairwisePotentials[entry][i] += msgs[i];
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
                    assert(false);
                    //const bool changed = (left_primal != r.solution[entry]);
                    //l.primal()[0] = r.solution[entry] / l.dim1();
                    //l.primal()[1] = r.solution[entry] % l.dim1();
                    //return changed;
                } else {
                    return false;
                }
            }

            template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
            bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
            {
                const INDEX left_primal = l.primal()[0]*l.dim1() + l.primal()[1];
                if(r.solution() < r.NumLabels[entry]) {
                    assert(false);
                    //const bool changed = (left_primal != r.solution[entry]);
                    //r.solution[entry] = left_primal;
                    //return changed;
                } else {
                    return false;
                }
            }

            template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
            bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
            {
                const INDEX left_primal = l.primal()[0]*l.dim1() + l.primal()[1];
                return left_primal == r.solution[entry];
            } 

            template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
            void construct_constraints(SOLVER& s, 
                    LEFT_FACTOR& l, typename SOLVER::vector l_left_msg_variables, typename SOLVER::vector l_right_msg_variables, typename SOLVER::matrix l_pairwise_variables, 
                    RIGHT_FACTOR& r)
            {
            }

        private:
            const INDEX entry;
    };
    
}

#endif // LP_MP_HORIZON_TRACKING_FACTORS_HXX
