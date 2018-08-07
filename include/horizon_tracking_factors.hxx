#ifndef LP_MP_HORIZON_TRACKING_FACTORS_HXX
#define LP_MP_HORIZON_TRACKING_FACTORS_HXX

#include "vector.hxx"
#include <queue>
#include <limits>
#include "two_dimensional_variable_array.hxx"
#include "three_dimensional_variable_array.hxx"
#include <unordered_set>
#include <cmath>

namespace LP_MP {

    class max_potential_on_graph {

        struct MaxPotentialElement {
            REAL value;
            INDEX tableIndex;
            INDEX labelIndex;
        };

        public:
            // indices correspond to: chain index, entry in chain, max pot value, linear cost.
            max_potential_on_graph(const std::vector<std::vector<std::array<REAL,2>>>& marginals_collection)
            : marginals_collection_(marginals_collection)
            {
                for(std::size_t currentTableIndex = 0; currentTableIndex < marginals_collection_.size(); ++currentTableIndex)
                 {
                    for(std::size_t currentLabel = 0; currentLabel < marginals_collection_[currentTableIndex].size(); ++currentLabel )
                    {
                        MaxPotentials.push_back( { marginals_collection_[currentTableIndex][currentLabel][0], currentTableIndex, currentLabel } );
                    }
                }
                SortingOrder = GetMaxPotsSortingOrder(MaxPotentials);
            }

            void init_primal() 
            {
                primal_computed = false;
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

            template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
            template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

            auto export_variables() { return std::tie( ); }
                        
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

            INDEX& max_potential_index(const INDEX tableIndex) 
            {
                assert(tableIndex < max_potential_index_.size()); 
                return max_potential_index_[tableIndex]; 
            }
            INDEX max_potential_index(const INDEX tableIndex) const
            {
                assert(tableIndex < max_potential_index_.size());
                if (!primal_computed)
                    Solve();

                return max_potential_index_[tableIndex]; 
            }

            std::vector<std::vector<std::array<REAL,2>>>& marginals_collection() 
            {
                return marginals_collection_; 
            }
            std::vector<std::vector<std::array<REAL,2>>> marginals_collection() const 
            {
                return marginals_collection_; 
            }

        private:
            std::vector<std::vector<std::array<REAL,2>>> marginals_collection_;

            // For primal recovery: Compute slack and send message to chains. Contains: max pot value, linear pot value, chain index which can receive message (TODO: Can it be more than one?)
            // TODO: Not needed during optimization only required for the last pass.
            mutable std::vector<std::vector<std::array<REAL,2>>> graph_marginals; // same layout as marginals_collection_ but contains the sum of linear pots overall all tables (chains).

            std::vector<MaxPotentialElement> MaxPotentials;
            mutable std::vector<INDEX> max_potential_index_;
            mutable REAL solutionObjective;
            mutable bool primal_computed = false;
            std::vector<INDEX> SortingOrder;

            void Solve() const
            {
                INDEX numCovered = 0;
                INDEX numTables = marginals_collection_.size();
                std::vector<bool> coveredTables(numTables, false);
                REAL s = 0;
                std::vector<REAL> l(numTables, INFINITY);
                double bestObjective = INFINITY;
                std::vector<INDEX> bestLabelsForTables(numTables);
                std::vector<INDEX> lablesForTables(numTables);
                graph_marginals.clear();
                graph_marginals.resize(numTables);

                for(const auto& currentElementToInsert : SortingOrder)
                {
                    INDEX currentTableIndex = MaxPotentials[currentElementToInsert].tableIndex;
                    INDEX currentLabelIndex = MaxPotentials[currentElementToInsert].labelIndex;
                    REAL currentLinearCost = marginals_collection_[currentTableIndex][currentLabelIndex][1];
                    REAL currentMaxCost =  MaxPotentials[currentElementToInsert].value;
                    assert(currentMaxCost == marginals_collection_[currentTableIndex][currentLabelIndex][0]);

                    // If the edge is not yet covered:
                    if (!coveredTables[currentTableIndex])  
                    {
                        coveredTables[currentTableIndex] = true;
                        numCovered++;
                        s += currentLinearCost;
                        l[currentTableIndex] = currentLinearCost;
                        lablesForTables[currentTableIndex] = currentLabelIndex;
                    }
                    
                    // If edge has been added, but current label has lower linear cost. We have two scenarios:
                    // 1. Graph has not been covered completely in which case we want to increase our max pot. threshold anyway. Thus, if we are gaining 
                    //      an improvement in linear cost take it.
                    // 2. Graph is covered completely, but we want to see adding current label can possibly give us any benefit, which will be 
                    //    checked by the 3rd if condition.
                    if (currentLinearCost < l[currentTableIndex])  
                    {
                        s = s - l[currentTableIndex] + currentLinearCost;
                        l[currentTableIndex] = currentLinearCost;
                        lablesForTables[currentTableIndex] = currentLabelIndex;
                    }

                    if (numCovered == numTables) 
                    {
                        // Build all the marginals even the ones which cannot be optimal, as they are used for downward messages:
                        graph_marginals[currentTableIndex].push_back({currentMaxCost, s});                  

                        // Found another solution which is better than the previous one, in which case mark current solution as the best so far.
                        if (bestObjective > s + currentMaxCost)
                        {
                            bestObjective = s + currentMaxCost;
                            solutionObjective = s + currentMaxCost;
                            bestLabelsForTables = lablesForTables;
                        }
                    }
                }

                max_potential_index_ = bestLabelsForTables; //TODO: The previous values for this solution will be lost! Do we need them?
                primal_computed = true;
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
            REAL distance = std::numeric_limits<REAL>::max();
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

        struct EdgePriority
        {
            REAL value;
            INDEX index;
        };
        
        struct EdgePriorityComparison
        {
            bool operator() (const EdgePriority& lhs, const EdgePriority& rhs)
            {
                return lhs.value < rhs.value;
            }
        };
        
        public:     
            three_dimensional_variable_array<REAL> LinearPairwisePotentials;  

            max_potential_on_chain(const three_dimensional_variable_array<REAL>& maxPairwisePotentials, const three_dimensional_variable_array<REAL>& linearPairwisePotentials, const std::vector<INDEX>& numLabels, INDEX chainIndex, bool useEdgeDeletion = false)
            :
                LinearPairwisePotentials(linearPairwisePotentials),
                MaxPairwisePotentials(maxPairwisePotentials),
                NumNodes(numLabels.size()),
                NumLabels(numLabels),
                ChainIndex(chainIndex),
                UseEdgeDeletion(useEdgeDeletion)
            {
                assert(maxPairwisePotentials.dim1() + 1 == numLabels.size()); 
                assert(maxPairwisePotentials.dim1() == linearPairwisePotentials.dim1());

                for(std::size_t n1=0; n1<maxPairwisePotentials.dim1(); ++n1) {
                    assert(maxPairwisePotentials.dim2(n1) == linearPairwisePotentials.dim2(n1) && maxPairwisePotentials.dim3(n1) == linearPairwisePotentials.dim3(n1));
                    for(std::size_t i=0; i<maxPairwisePotentials.dim2(n1); ++i) {
                        for(std::size_t j=0; j<maxPairwisePotentials.dim3(n1); ++j) {
                            MaxPotentials1D.push_back( {maxPairwisePotentials(n1,i,j), n1, n1+1, i, j} );
                        }
                    }
                }

                MaxPotsSortingOrder = GetPairwisePotsSortingOrder(MaxPotentials1D);
                solution_.assign(NumNodes, 0);
                init_primal();
                InsertTerminalNode();
            }

            REAL LowerBound() const
            {
                if (!max_potential_marginals_valid_)
                    Solve();
                    
                INDEX bestIndex = GetBestMarginal();
                return max_potential_marginals_[bestIndex][1] + max_potential_marginals_[bestIndex][2];
                // compute optimal solution and return its cost
            }

            REAL EvaluatePrimal() const
            {
                assert(max_potential_index_ != std::numeric_limits<std::size_t>::max());
                return max_potential_marginals_[max_potential_index_][1] + max_potential_marginals_[max_potential_index_][2];
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
                const bool labels_computed = *std::max_element(solution_.begin(), solution_.end()) < std::numeric_limits<INDEX>::max();
                const bool max_potential_index_computed = max_potential_index_ != std::numeric_limits<std::size_t>::max();
                
                if(max_potential_index_computed && labels_computed) return;
                if(max_potential_index_computed) { ComputeSolution(); return; }
                if(!max_potential_index_computed && !labels_computed) { 
                    if (!max_potential_marginals_valid_) {
                        Solve();
                        set_max_potential_index(GetBestMarginal());
                    }
                    ComputeSolution();
                }
                if(labels_computed && !max_potential_index_computed) {
                    assert(false);
                }
            }

            template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
            template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( LinearPairwisePotentials.data() ); }

            auto export_variables() { return std::tie( ); }
            void init_primal() 
            {
                max_potential_index_ = std::numeric_limits<std::size_t>::max();
                solutionObjective = std::numeric_limits<std::size_t>::max();
                std::fill(solution_.begin(), solution_.end(), std::numeric_limits<INDEX>::max());
            } 

            void invalidate_marginals()
            {
                max_potential_marginals_valid_ = false;
            }
            
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

            INDEX& solution(const std::size_t i) 
            {
                assert(i < solution_.size()); 
                return solution_[i]; 
            }

            INDEX solution(const std::size_t i) const 
            {
                assert(i < solution_.size());
                const bool labels_computed = *std::max_element(solution_.begin(), solution_.end()) < std::numeric_limits<std::size_t>::max();
                assert(labels_computed);
                return solution_[i]; 
            }

            void set_max_potential_index(const std::size_t index) const
            {
                assert(max_potential_marginals_valid_);
                if (max_potential_index_ == index)
                    return; // Can use this for optimization?

                max_potential_index_ = index;
            }

            std::size_t max_potential_index() const 
            {
                return max_potential_index_;
            }
        
            std::array<REAL,3>& max_potential_marginal(const std::size_t i) { assert(i < max_potential_marginals_.size()); return max_potential_marginals_[i]; }
            std::array<REAL,3> max_potential_marginal(const std::size_t i) const 
            {
                assert(i < max_potential_marginals_.size());
                return max_potential_marginals_[i];
            }

            std::vector<std::array<REAL,3>> max_potential_marginals() const 
            {
                return max_potential_marginals_;
            }

        protected:
                std::vector<INDEX> NumLabels;
                // primal solution
                mutable std::vector<INDEX> solution_;
                mutable std::size_t max_potential_index_;

        private:
            mutable std::vector<MaxPairwisePotential> MaxPotentials1D;
            std::vector<INDEX> MaxPotsSortingOrder;
            three_dimensional_variable_array<REAL> MaxPairwisePotentials;
            mutable std::vector<std::array<REAL,3>> max_potential_marginals_; // (i) max potential, (ii) minimum linear potential, (iii) cost of configuration 
            mutable bool max_potential_marginals_valid_ = false;
            mutable bool MaxPotMarginalsInitialized = false;
            INDEX ChainIndex;

            bool UseEdgeDeletion;

            int NumNodes;
            mutable REAL solutionObjective = std::numeric_limits<std::size_t>::max();

            // This function will only compute the marginals, not the solution objective nor the primal solution.
            void Solve() const
            {
                if (UseEdgeDeletion)
                    SolveByEdgeDeletion();
                else
                    SolveByEdgeAddition();

                MaxPotMarginalsInitialized = true;
                max_potential_marginals_valid_ = true;
            }

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
                    currentPot.value = std::numeric_limits<REAL>::lowest();
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

            INDEX GetBestMarginal() const
            {
                assert(max_potential_marginals_valid_);
                REAL bestCost = std::numeric_limits<std::size_t>::max();
                INDEX bestIndex;

                for (INDEX currentMaxPotIndex = 0; currentMaxPotIndex < max_potential_marginals_.size(); ++currentMaxPotIndex)
                {
                    auto currentMarginal = max_potential_marginals_[currentMaxPotIndex];
                    auto currentCost = currentMarginal[0] + currentMarginal[1] + currentMarginal[2];
                    if (currentCost < bestCost)
                    {
                        bestIndex = currentMaxPotIndex;
                        bestCost = currentCost;
                    }
                }
                return bestIndex;
            }
            
            // Computes the primal solution and its objective and stores it.
            // TODO: There should be a flag for knowing if we need to recompute or not
            void ComputeSolution() const
            {
                assert(max_potential_marginals_valid_);
                REAL maxPotThreshold = max_potential_marginals_[max_potential_index_][0];
                ShortestPathTreeInChain spTree = FindAndInitializeSPTree(maxPotThreshold);
                solutionObjective =  max_potential_marginals_[max_potential_index_][1] + max_potential_marginals_[max_potential_index_][2];
                // max_potential_marginals_[max_potential_index_][2] + spTree.GetDistance(NumNodes, 0);
                INDEX currentLabel = 0; // for terminal node;
                for (int currentNodeToBackTrack = NumNodes - 1; currentNodeToBackTrack >= 1; currentNodeToBackTrack--) // Exclude source node.
                {
                    solution_[currentNodeToBackTrack - 1] = spTree.GetParentLabel(currentNodeToBackTrack + 1, currentLabel);
                    currentLabel = solution_[currentNodeToBackTrack - 1]; //TODO: TEST!
                }
            }

            void SolveByEdgeAddition() const
            {
                REAL bestSolutionCost = INFINITY;
                std::vector<std::vector<REAL> > distanceFromSource(NumNodes);
                REAL initialValue = 0;
                for (INDEX i = 0 ; i < NumNodes ; i++ )
                {
                    if (i > 0)
                        initialValue = std::numeric_limits<REAL>::max();

                    distanceFromSource[i].resize(NumLabels[i], initialValue);
                }

                INDEX currentMaxPotIndex = 0;
                for(const auto& currentEdgeToInsert : MaxPotsSortingOrder)
                {
                    bool foundPath = UpdateDistances(currentEdgeToInsert, distanceFromSource, MaxPotentials1D[currentEdgeToInsert].value);

                    REAL currentLinearCost =  distanceFromSource[NumNodes - 1][0]; 

                    //TODO: Storing ALL possible max potentials, even the ones which are not feasible!                    
                    if (currentLinearCost == std::numeric_limits<REAL>::max())
                       continue;

                    // Insert the marginal, and do not increment the index if the max pot was already present
                    // at previous index in which case the marginal was not inserted and we only took min:
                    if (InsertMarginal<true>(MaxPotentials1D[currentEdgeToInsert].value, currentMaxPotIndex, currentLinearCost))
                        currentMaxPotIndex++;
                }
            }

            bool UpdateDistances(INDEX edgeToUpdate, std::vector<std::vector<REAL> >& distanceFromSource, REAL maxPotThresh) const
            {
                bool reachedTerminal = false;
                std::queue<INDEX> queue;
                queue.push(edgeToUpdate);

                while(!queue.empty())
                {
                    INDEX currentEdge = queue.front();
                    queue.pop();
                    const auto& currentMaxPot = MaxPotentials1D[currentEdge];
                    REAL currentLinearPot = 0;
                    const auto& n1 = currentMaxPot.n1;
                    const auto& n2 = currentMaxPot.n2;
                    const auto& l1 = currentMaxPot.l1;
                    const auto& l2 = currentMaxPot.l2;

                    if (n2 < NumNodes - 1) // As LinearPairwisePotentials does not contain potentials from last node to terminal node
                        currentLinearPot = LinearPairwisePotentials(n1, l1, l2);
                    
                    REAL offeredDistanceTon2l2 = distanceFromSource[n1][l1] + currentLinearPot;
                    auto currentDistanceTon2l2 = distanceFromSource[n2][l2];

                    if (offeredDistanceTon2l2 >= currentDistanceTon2l2)
                        continue;

                    distanceFromSource[n2][l2] = offeredDistanceTon2l2;
                    
                    if (n2 == NumNodes - 1)
                    {
                        reachedTerminal = true;
                        continue;
                    }

                    INDEX n3 = n2 + 1;

                    // The distance of n2, l2 has been updated so add all of its immediate children to the queue to be inspected.
                    INDEX firstEdgeToConsider = currentEdge + (NumLabels[n2] - 1 - l2) + (NumLabels[n1] - 1 - l1) * NumLabels[n2] + l2 * NumLabels[n3] + 1;
                    for (INDEX l3 = 0, currentEdgeToConsider = firstEdgeToConsider; l3 < NumLabels[n3]; ++l3, ++currentEdgeToConsider)
                    {
                        // Do not consider this potential as it has not been added through sorting yet.
                        if (MaxPotentials1D[currentEdgeToConsider].value > maxPotThresh)
                            continue;
                        
                        // auto childMaxPot = MaxPotentials1D[currentEdgeToConsider];
                        // assert(childMaxPot.l1 == l2);
                        // assert(childMaxPot.l2 == l3); // Might fail if some of the pairwise potentials are not present thus causing jumps!
                        queue.push(currentEdgeToConsider);
                    }
                }
                return reachedTerminal;
            }

            template <bool insertEnd>
            bool InsertMarginal(REAL maxPotValue, INDEX insertionIndex, REAL currentLinearCost) const
            {
                if (!MaxPotMarginalsInitialized)
                {
                    if(insertEnd)
                    {
                        if(max_potential_marginals_.size() > 0 && max_potential_marginals_.back()[0] == maxPotValue)
                            max_potential_marginals_.back()[1] = std::min(max_potential_marginals_.back()[1], currentLinearCost);
                        else
                            max_potential_marginals_.push_back({maxPotValue, currentLinearCost, 0});
                    }
                    else
                    {
                        if(max_potential_marginals_.size() > 0 && max_potential_marginals_.front()[0] == maxPotValue)
                            max_potential_marginals_.front()[1] = std::min(max_potential_marginals_.front()[1], currentLinearCost);
                        else
                            max_potential_marginals_.insert(max_potential_marginals_.begin(), {maxPotValue, currentLinearCost, 0});
                    }
                }
                else
                {                       
                    // Check if the current max potential value is also present in the marginals (can only be present at adjacent index as they are sorted),
                    // if yes just take the minimum of the linear costs.
                    INDEX adjacentIndex = insertionIndex + (insertEnd ? -1:1);
                    if (adjacentIndex >= 0 && adjacentIndex <= max_potential_marginals_.size() - 1 &&
                         max_potential_marginals_[adjacentIndex][0] == maxPotValue)
                    {
                        max_potential_marginals_[adjacentIndex][1] = std::min(max_potential_marginals_[adjacentIndex][1], currentLinearCost);
                        return false;
                    }
                    
                    assert(maxPotValue == max_potential_marginals_[insertionIndex][0]);
                    
                    max_potential_marginals_[insertionIndex][1] = currentLinearCost;
                    return true;
                }
            }

            std::vector<INDEX> GetPairwisePotsSortingOrder(const std::vector<MaxPairwisePotential>& pots) const
            {
                std::vector<INDEX> idx(pots.size());
                std::iota(idx.begin(), idx.end(), 0);
                std::sort(idx.begin(), idx.end(), [&pots](INDEX i1, INDEX i2) {return pots[i1].value < pots[i2].value;});
                return idx;
            }

            ShortestPathTreeInChain FindAndInitializeSPTree(REAL maxPotThreshold = INFINITY) const
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
                            REAL currentMaxPot = std::numeric_limits<REAL>::lowest();
                            if (n < NumNodes - 1 && n > 0)   // Source and Terminal node potentials are not present in it.
                            {
                                currentLinearPot = LinearPairwisePotentials(n - 1, prevLabel, l);
                                currentMaxPot = MaxPairwisePotentials(n - 1, prevLabel, l);
                            }
                            
                            if (currentMaxPot > maxPotThreshold)
                                continue;

                            spTree.CheckParentForShortestPath(n + 1, l, prevLabel, currentLinearPot, currentMaxPot);
                        }
                    }
                }

                return spTree;
            }

            void SolveByEdgeDeletion() const
            {
                // TODO: SolveByEdgeDeletion can potentially delete multiple edges at each iteration if they are not present in the shortest path
                // so by that way we won't get the marginals for all possible marginals, however we do need all marginals because in the 
                // subsequent iterations some unpopulated max potential values would be needed if now these values are the values of the 
                // shortest path due to reparameterization of the linear potentials. This function should not be used until this issue is addressed. So:
                assert(false);

                //TODO: Only seeing the tree Max pot value is not going to work, we need to see the marginals for ALL possible paths,
                ShortestPathTreeInChain spTree = FindAndInitializeSPTree();
                REAL treeMaxPotValue = spTree.GetMaxPotValueInTree().maxPotValue;
                INDEX currentMaxPotIndex = max_potential_marginals_.size() - 1;

                if (InsertMarginal<false>(treeMaxPotValue, currentMaxPotIndex, spTree.GetDistance(NumNodes, 0)))
                    currentMaxPotIndex--;

                for (int i = MaxPotsSortingOrder.size() - 1; i >= 0; i--)
                {   
                    auto currentMaxPotEdge = MaxPotentials1D[MaxPotsSortingOrder[i]];
                    bool wasTreeEdge = spTree.CheckAndPopMaxPotInTree(
                        currentMaxPotEdge.n1 + 1, currentMaxPotEdge.l1, currentMaxPotEdge.n2 + 1, currentMaxPotEdge.l2, currentMaxPotEdge.value); // +1 due to source node.
                    spTree.RemovePossibleEdge(currentMaxPotEdge.n1 + 1, currentMaxPotEdge.l1, currentMaxPotEdge.l2);

                    if (!wasTreeEdge)
                    {
                        if (InsertMarginal<false>(currentMaxPotEdge.value, currentMaxPotIndex, spTree.GetDistance(NumNodes, 0)))
                            currentMaxPotIndex--;
                    }
                    else
                    {
                        std::priority_queue<AffectedVertex> locallyAffectedPQueue;
                        std::unordered_set<std::array<INDEX, 2>> locallyAffectedVertices = spTree.GetLocallyAffectedNodes();
                        for (auto const& currentAffectedNode : locallyAffectedVertices)
                        {   
                            spTree.SetStatus(currentAffectedNode[0], currentAffectedNode[1], 0);    //mark as open
                            INDEX bestParentLabel;
                            REAL bestParentMaxPotValue; 
                            REAL bestParentDistance = std::numeric_limits<REAL>::max();
                            for (INDEX prevLabel = 0; prevLabel < NumLabels[currentAffectedNode[0] - 2]; prevLabel++)
                            {
                                REAL currentLinearPot = 0;
                                REAL currentMaxPot = std::numeric_limits<REAL>::lowest();
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

                            if (bestParentDistance < std::numeric_limits<REAL>::max()) //TODO Replace by max
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

                        if (InsertMarginal<false>(spTree.GetMaxPotValueInTree().maxPotValue, currentMaxPotIndex, spTree.GetDistance(NumNodes, 0)))
                            currentMaxPotIndex--;                      
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

                    //Need to insert ALL max pots between lb and ub even if for now they dont seem optimal, they CAN become optimal after reparametrization!
                    // else if (lbKey != potentials.begin()) // There is a key before the current element which needs to be checked with current max pot.
                    // {
                    //     assert(maxPot > std::prev(lbKey)->first);
                    //     if (!(linearPot < std::prev(lbKey)->second)) // Inserting has no benefit as linearPot is not strictly less than previous linear pot.
                    //         return;
                    // }
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

        //TODO: Optimize this:
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
        max_potential_on_tree(const three_dimensional_variable_array<REAL>& maxPairwisePotentials, const three_dimensional_variable_array<REAL>& linearPairwisePotentials,
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
            ComputeMaxPotLowerBound();
            solution_.assign(NumNodes, 0);
        }

        // Call this function whenever linear/max potentials get changed so the bounds need to be recomputed.
        void LinearPotsChanged()
        {
            boundsDirty = true;
        }

        REAL GetMaxPotLowerBound()
        {
            return MaxPotentialLowerBoundForAllTrees;
        }

        REAL GetMaxPotUpperBound()
        {
            if (boundsDirty)
            {
                ComputeMaxPotUpperBound();
                boundsDirty = false;
            }
            return MaxPotentialUpperBound;
        }

        REAL LowerBound() const
        {
            ComputeSolution();
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
            return solutionObjective;
        }

        void MaximizePotentialAndComputePrimal() 
        {
            ComputeSolution();
            // compute optimal solution and store it
        }

        INDEX& solution(const std::size_t i) { assert(i < solution_.size()); return solution_[i]; }
        INDEX solution(const std::size_t i) const { assert(i < solution_.size()); return solution_[i]; }

        void set_max_potential_index(const std::size_t index)
        {
            if (max_potential_index_ == index)
                return;

            max_potential_index_ = index;
            ComputeSolution();
        }

        std::size_t max_potential_index() const { return max_potential_index_; }
                        
        std::array<REAL,3>& max_potential_marginal(const std::size_t i) { assert(i < max_potential_marginals_.size()); return max_potential_marginals_[i]; }
        std::array<REAL,3> max_potential_marginal(const std::size_t i) const { assert(i < max_potential_marginals_.size()); return max_potential_marginals_[i]; }
        std::vector<std::array<REAL,3>> max_potential_marginals() const { return max_potential_marginals_; }

    private:
        mutable std::vector<INDEX> solution_;
        three_dimensional_variable_array<REAL> MaxPairwisePotentials;
        three_dimensional_variable_array<REAL> LinearPairwisePotentials;
        INDEX NumNodes;
        std::vector<INDEX> NumLabels;
        std::vector<INDEX> NumEdgesForNode;
        
        std::vector<std::array<INDEX, 2>> MessagePassingSchedule;
        static REAL MaxPotentialLowerBoundForAllTrees;    // Computed by max potential message passing.
        mutable REAL MaxPotentialUpperBound;
        mutable REAL LinearPotentialLowerBound; // Computed by conventional message passing.
        // mutable REAL LinearPotentialUpperBound; // TODO: Can be computed from max potential message passing and by used to prune paths longer than this bound.
        mutable std::size_t max_potential_index_;
        mutable std::vector<std::array<REAL,3>> max_potential_marginals_; // (i) max potential, (ii) minimum linear potential, (iii) cost of configuration 
        mutable bool MaxPotMarginalsInitialized = false;

        bool boundsDirty = true;
        mutable REAL solutionObjective = std::numeric_limits<REAL>::max();

        void ComputeMaxPotLowerBound() const
        {
            std::array<REAL, 2> bounds = MessagePassingForOnePotential(MaxPairwisePotentials, LinearPairwisePotentials, false);
            MaxPotentialLowerBoundForAllTrees = fmax(bounds[0], MaxPotentialLowerBoundForAllTrees); // To take max over all trees.
            // LinearPotentialUpperBound = bounds[1]; //TODO: If this is not useful remove it and compute max potential lb only once and store it, as it will not change.
            // Seems like LinearPotentialUpperBound wont help in anything, so we dont need to call this function again and again after updated linear pots.
        }

        void ComputeMaxPotUpperBound() const
        {
            std::array<REAL, 2> bounds = MessagePassingForOnePotential(LinearPairwisePotentials, MaxPairwisePotentials, true);
            LinearPotentialLowerBound = bounds[0];
            MaxPotentialUpperBound = bounds[1];
        }

        // Computes and stores the marginals of the root node.
        void Solve() const
        {
            std::vector<std::vector<LabelStateSpace>> messages(NumNodes);
            for (INDEX i = 0; i < NumNodes; i++)
            {
                messages[i].resize(NumLabels[i]);
                for (INDEX l = 0; l < NumLabels[i]; l++) 
                    messages[i][l] = LabelStateSpace(MaxPotentialLowerBoundForAllTrees, MaxPotentialUpperBound);
            }

            std::vector<INDEX> totalSentAndReceivedMessages(NumNodes, 0);
            std::vector<bool> messageReceived(NumNodes, false);
            INDEX edgeIndex = 0;
            for (const auto & currentEdge : MessagePassingSchedule)
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
                    LabelStateSpace lhStateSpace = GetHeadLabelStateSpaceFromCurrentTail(messages[tail], edgeIndex, lh, head < tail, isTailLeaf);
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
            INDEX maxPotIndex = 0;
            for (const auto &currentPair : rootPots)
            {
                if (!MaxPotMarginalsInitialized)
                    max_potential_marginals_.push_back({currentPair.first, currentPair.second, 0});
                else
                {
                    // Due to std::map the rootPots should be already sorted w.r.t increasing max pots:
                    assert(currentPair.first == max_potential_marginals_[maxPotIndex][0]); 
                    max_potential_marginals_[maxPotIndex][1] = currentPair.second;
                }
                maxPotIndex++;
            }
        }

        LabelStateSpace GetHeadLabelStateSpaceFromCurrentTail(const std::vector<LabelStateSpace>& tailMessages, INDEX edgeIndex, INDEX lh, bool isReverse, bool isTailLeaf = false) const 
        {
            LabelStateSpace lhStateSpace(MaxPotentialLowerBoundForAllTrees, MaxPotentialUpperBound);
            for (INDEX lt = 0; lt < tailMessages.size(); lt++)
            {
                // Assuming that the 2D matrix of potentials for each edge is stored always in the order in which tail always counts as l2 and head as l1.
                // TODO: Maybe do not need isReverse if the above condition is satisfied?
                assert(MaxPairwisePotentials.dim3(edgeIndex) == tailMessages.size());
                assert(LinearPairwisePotentials.dim3(edgeIndex) == tailMessages.size()); 

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

        // TODO: There should be a flag for knowing if we need to recompute or not
        void ComputeSolution() const
        {
            std::vector<std::vector<REAL>> messages(NumNodes);
            std::vector<bool> messageSent(NumNodes, false);
            std::vector<bool> messageReceived(NumNodes, false);
            for (INDEX i = 0; i < NumNodes; i++)
            {
                messages[i].resize(NumLabels[i], 0);
            }

            INDEX edgeIndex = 0;
            for (const auto & currentEdge : MessagePassingSchedule)
            {
                INDEX tail = currentEdge[0];
                INDEX head = currentEdge[1];
                if (!messageReceived[tail]) // Only leaf nodes can send a message (only one) and only before receiving any.
                {
                    assert(!messageSent[tail]);
                }
                for (INDEX lh = 0; lh < messages[head].size(); lh++)
                {
                    REAL minMessage = ComputeMessageValue(messages[tail], edgeIndex, lh, head < tail);
                    messages[head][lh] += minMessage;
                }
                messageSent[tail] = true;
                messageReceived[head] = true;
                edgeIndex++;
            }

            bool solutionObjectiveSet = false;

            for (auto it = MessagePassingSchedule.rbegin(); it != MessagePassingSchedule.rend(); ++it)           
            {
                edgeIndex--;
                const auto & currentEdge = *it;
                INDEX tail = currentEdge[1];
                INDEX head = currentEdge[0];
                REAL tailNodeBestMessage = INFINITY;
                for (INDEX lt = 0; lt < messages[tail].size(); lt++)
                {
                    if (messages[tail][lt] < tailNodeBestMessage)
                    {
                        solution_[tail] = lt;
                        tailNodeBestMessage = messages[tail][lt];
                    }
                }
                
                if (!solutionObjectiveSet) // Works for root node.
                {
                    solutionObjective = tailNodeBestMessage + max_potential_marginals_[max_potential_index_][2];
                    solutionObjectiveSet = true;
                }

                for (INDEX lh = 0; lh < messages[head].size(); lh++)
                {
                    REAL minMessage = ComputeMessageValue(messages[tail], edgeIndex, lh, head < tail, solution_[tail]);
                    messages[head][lh] += minMessage;
                }
            }
        }

        REAL ComputeMessageValue(const std::vector<REAL>& tailMessages, INDEX edgeIndex, 
        INDEX lh, bool isReverse, bool labelTail = -1) const
        {
            REAL bestMessageValue = INFINITY;
            INDEX bestlt;
            INDEX ltStart = 0;
            INDEX ltEnd = tailMessages.size() - 1;
            
            if (labelTail >= 0)
            {
                ltStart = labelTail;
                ltEnd = labelTail;
            }

            for (INDEX lt = ltStart; lt <= ltEnd; lt++)
            {
                REAL currentEdgeMaxPot = !isReverse ? MaxPairwisePotentials(edgeIndex, lt, lh) :  MaxPairwisePotentials(edgeIndex, lh, lt);
                if (currentEdgeMaxPot > max_potential_marginals_[max_potential_index_][2])
                    continue;

                REAL currentEdgePot = !isReverse ? LinearPairwisePotentials(edgeIndex, lt, lh) :  LinearPairwisePotentials(edgeIndex, lh, lt);
                REAL currentIncomingMsg = tailMessages[lt];

                REAL msg = currentIncomingMsg + currentEdgePot;
                if (msg < bestMessageValue)
                {
                    bestMessageValue = msg;
                    bestlt = lt;
                }
            }

            return bestMessageValue;
        }

        std::array<REAL, 2> MessagePassingForOnePotential(const three_dimensional_variable_array<REAL>& mainPairwisePots, const three_dimensional_variable_array<REAL>& otherPairwisePots, bool doConventional) const
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
            for (const auto & currentEdge : MessagePassingSchedule)
            {
                INDEX tail = currentEdge[0];
                INDEX head = currentEdge[1];
                if (!messageReceived[tail]) // Only leaf nodes can send a message (only one) and only before receiving any.
                {
                    assert(!messageSent[tail]);
                }
                for (INDEX lh = 0; lh < mainMessages[head].size(); lh++)
                {
                    std::array<REAL, 2> minMessage = ComputeMessageValuePair(mainMessages[tail], otherMessages[tail], edgeIndex, lh, mainPairwisePots, otherPairwisePots, doConventional, head < tail);
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

        std::array<REAL, 2> ComputeMessageValuePair(const std::vector<REAL>& tailMainMessages, const std::vector<REAL>& tailOtherMessages, INDEX edgeIndex, INDEX lh, 
        const three_dimensional_variable_array<REAL>& mainPairwisePots, const three_dimensional_variable_array<REAL>& otherPairwisePots, bool doConventional, bool isReverse) const
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
REAL max_potential_on_tree::MaxPotentialLowerBoundForAllTrees = std::numeric_limits<REAL>::lowest(); //TODO: Will need to be re-initialized for a new instance.


class unary_max_potential_on_chain_message {
    public:
        // UNUSED CLASS:
        unary_max_potential_on_chain_message(const std::size_t nodeIndex) : variable(nodeIndex) {}
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
            r.invalidate_marginals();
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
            if(r.primal()[variable] < l.size()) {
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

class pairwise_max_factor_tree_message {
    public:
        pairwise_max_factor_tree_message(const INDEX _pairwise_entry, INDEX _unary_1, INDEX _unary_2) :
            pairwise_entry(_pairwise_entry), 
            unary_1(_unary_1),
            unary_2(_unary_2)
            {}

        template<typename FACTOR, typename MSG>
        void RepamRight(FACTOR& r, const MSG& msgs)
        {
            INDEX i = 0;
            for(INDEX l1=0; l1<r.LinearPairwisePotentials.dim2(pairwise_entry); ++l1) {
                for(INDEX l2=0; l2<r.LinearPairwisePotentials.dim3(pairwise_entry); ++l2, ++i) {
                    r.LinearPairwisePotentials(pairwise_entry, l1, l2) += msgs[i];
                }
            }
            r.invalidate_marginals();
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
            bool changed_unary_1 = false;
            if(r.solution(unary_1) < l.dim1()) {
                changed_unary_1 = (r.solution(unary_1) != l.primal()[0]);
                l.primal()[0] = r.solution(unary_1);
            }

            bool changed_unary_2 = false;
            if(r.solution(unary_2) < l.dim2()) {
                changed_unary_2 = (r.solution(unary_2) != l.primal()[1]);
                l.primal()[1] = r.solution(unary_2);
            }

            return changed_unary_1 || changed_unary_2;
        }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
        {
            bool changed_unary_1 = false;
            if(l.primal()[0] < l.dim1()) {
                changed_unary_1 = (r.solution(unary_1) != l.primal()[0]);
                r.solution(unary_1) = l.primal()[0];
            }

            bool changed_unary_2 = false;
            if(l.primal()[1] < l.dim2()) {
                changed_unary_2 = (r.solution(unary_2) != l.primal()[1]);
                r.solution(unary_2) = l.primal()[1];
            }

            return changed_unary_1 || changed_unary_2;
        }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
        {
            return r.solution(unary_1) == l.primal()[0] && r.solution(unary_2) == l.primal()[1];
        } 

        template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
        void construct_constraints(SOLVER& s, 
                LEFT_FACTOR& l, typename SOLVER::vector l_left_msg_variables, typename SOLVER::vector l_right_msg_variables, typename SOLVER::matrix l_pairwise_variables, 
                RIGHT_FACTOR& r)
        {
        }

    private:
        const INDEX pairwise_entry;
        const INDEX unary_1, unary_2;
};

// left factor is chain/tree
// right factor is max_potential_on_graph
class max_factor_tree_graph_message {
    public:

        max_factor_tree_graph_message(const std::size_t _entry) : entry(_entry) {}

        template<typename FACTOR, typename MSG>
        void RepamRight(FACTOR& r, const MSG& msg)
        {
            assert(r.marginals_collection()[entry].size() == msg.size());
            for(std::size_t i=0; i<r.marginals_collection()[entry].size(); ++i) {
                r.marginals_collection()[entry][i][1] += msg[i];
            }
        }

        template<typename FACTOR, typename MSG>
        void RepamLeft(FACTOR& l, const MSG& msg)
        {
            assert(msg.size() == l.max_potential_marginals().size());
            for(std::size_t i=0; i<msg.size(); ++i) {
                l.max_potential_marginal(i)[2] += msg[i]; 
            }
        }

        template<typename LEFT_FACTOR, typename MSG>
        void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
        {
            vector<REAL> m(l.max_potential_marginals().size());
            for(std::size_t i=0; i<m.size(); ++i) {
                m[i] = l.max_potential_marginal(i)[1] + l.max_potential_marginal(i)[2];
            }
            const auto min = m.min();
            for(auto& x : m) { x-= min; }
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
            if(r.max_potential_index(entry) != std::numeric_limits<INDEX>::max()) {
                const bool changed = (l.max_potential_index() != r.max_potential_index(entry));
                l.set_max_potential_index(r.max_potential_index(entry));
                return changed;
            } else {
                return false;
            }
            // l.max_potential_index() = r.max_potential_index(entry);
        }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
        {
            if(l.max_potential_index() != std::numeric_limits<INDEX>::max()) {
                const bool changed = (r.max_potential_index(entry) != l.max_potential_index());
                r.max_potential_index(entry) = l.max_potential_index();
                return changed;
            } else {
                return false;
            }
        }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
        {
            return r.max_potential_index(entry) == l.max_potential_index();
        } 

        template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
        void construct_constraints(SOLVER& s, LEFT_FACTOR& l, RIGHT_FACTOR& r)
        {
        }

    private:
    const std::size_t entry; // TODO: change name?
};
    
}

#endif // LP_MP_HORIZON_TRACKING_FACTORS_HXX
