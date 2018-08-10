#ifndef  LP_MP_HORIZON_TRACKING_UAI_CONSTRUCTOR_HXX
#define  LP_MP_HORIZON_TRACKING_UAI_CONSTRUCTOR_HXX

#include "vector.hxx"
#include "parse_rules.h"
#include "pegtl/parse.hh"
#include "tree_decomposition.hxx"
#include "arboricity.h"
#include "mrf_problem_construction.hxx"
#include "tree_decomposition.hxx"
#include <string>

namespace LP_MP {
namespace UAIMaxPotInput {
   struct MaxPotInput {
      std::size_t number_of_variables_;
      std::size_t number_of_cliques_;
      std::vector<INDEX> cardinality_;
      std::vector<std::vector<INDEX>> clique_scopes_;
      std::vector<std::vector<REAL>> function_tables_;
   };

// For getting the traversal order of the nodes, in the order of increasing number of labels.
class UaiFileGraph {
    public:
        UaiFileGraph(std::size_t numNodes) {
            NumNodes = numNodes;
            Nodes.resize(numNodes);
            NodesCovered.resize(numNodes);
        }

        struct GraphNode {
            std::size_t numLabels;
            std::vector<std::size_t> neighbourIndices;
        };

        void AddNode(std::size_t index, std::size_t numLabels) {
            Nodes[index].numLabels = numLabels;
            NodesCovered[index] = false;
            if (numLabels < minimumNumLabels)
            {
                minimumNumLabels = numLabels;
                minLabelNodeIndex = index;
            }
        }

        void AddEdge(std::size_t n1Index, std::size_t n2Index) {
            Nodes[n1Index].neighbourIndices.push_back(n2Index);
            Nodes[n2Index].neighbourIndices.push_back(n1Index);
        }

        struct Node {
            std::size_t NodeIndex;
            std::size_t NumLabels;
        };

        struct NodeCompare {
            bool operator() (const Node& lhs, const Node& rhs) const
            {
                return lhs.NumLabels > rhs.NumLabels;
            }
        };

        std::vector<std::size_t> getTraversalOrder() const {
            std::vector<std::size_t> traversalOrder;
            std::priority_queue<Node, std::vector<Node>, NodeCompare> nodeIndexAndNodeLabels;

            nodeIndexAndNodeLabels.push({minLabelNodeIndex, minimumNumLabels});
            NodesCovered[minLabelNodeIndex] = true;

            while (!nodeIndexAndNodeLabels.empty()) {
                auto currentBestNode = nodeIndexAndNodeLabels.top();
                nodeIndexAndNodeLabels.pop();
                traversalOrder.push_back(currentBestNode.NodeIndex);

                for (const auto& currentNeighbourIndex : Nodes[currentBestNode.NodeIndex].neighbourIndices) {
                    if (NodesCovered[currentNeighbourIndex])
                        continue;
                    
                    nodeIndexAndNodeLabels.push({currentNeighbourIndex, Nodes[currentNeighbourIndex].numLabels});
                    NodesCovered[currentNeighbourIndex] = true;
                }
            }
            return traversalOrder;
        }

    private:

    std::size_t NumNodes;
    std::vector<GraphNode> Nodes;
    mutable std::vector<bool> NodesCovered;
    std::size_t minLabelNodeIndex;
    std::size_t minimumNumLabels = std::numeric_limits<std::size_t>::max();
};


   // import basic parsers
   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::opt_invisible;
   using Parsing::mand_invisible;
   using Parsing::positive_integer;
   using Parsing::real_number;
   
   struct init_line_markov : pegtl::seq< opt_whitespace, pegtl::string<'M','A','R','K','O','V'>, opt_whitespace > {};
   struct init_line_max_pot : pegtl::seq< opt_whitespace, pegtl::string<'M','A','X','-','P','O','T','E','N','T','I','A','L','S'>, opt_whitespace > {};
   struct number_of_variables : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace > {};
   // vector of integers denoting how many labels each variable has
   struct cardinality : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace > {};
   struct number_of_cliques : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace> {};
   // first is the number of variables in the clique, then the actual variables.
   // the clique_scopes should match number_of_clique_lines, each line consisting of a sequence of integers
   struct new_clique_scope : pegtl::seq< positive_integer > {};
   struct clique_scope : pegtl::seq< positive_integer > {};
   struct clique_scope_line : pegtl::seq< opt_whitespace, new_clique_scope, pegtl::plus< opt_whitespace, clique_scope >, opt_whitespace, pegtl::eol > {};
   struct clique_scopes_end
   {
      //template< pegtl::apply_mode A, pegtl::rewind_mode M, template< typename ... > class Action, template< typename ... > class Control, typename Input >
      template< pegtl::apply_mode A, template< typename ... > class Action, template< typename ... > class Control, typename Input >
         static bool match( Input &, std::vector<MaxPotInput>& input )
         {
            return input.back().number_of_cliques_ == input.back().clique_scopes_.size();
         }
   }; 
   struct clique_scopes : pegtl::until< clique_scopes_end, clique_scope_line > {};
   // a function table is begun by number of entries and then a list of real numbers. Here we record all the values in the real stack
   // do zrobienia: treat whitespace
   struct new_function_table : pegtl::seq< positive_integer > {};
   struct function_table_entry : pegtl::seq< real_number > {};
   struct function_tables_end
   {
      //template< pegtl::apply_mode A, pegtl::rewind_mode M, template< typename ... > class Action, template< typename ... > class Control, typename Input >
      template< pegtl::apply_mode A, template< typename ... > class Action, template< typename ... > class Control, typename Input >
         static bool match( Input &, std::vector<MaxPotInput>& input )
         {
            return input.back().number_of_cliques_ == input.back().function_tables_.size();
         }
   };
   struct function_table_end
   {
      //template< pegtl::apply_mode A, pegtl::rewind_mode M, template< typename ... > class Action, template< typename ... > class Control, typename Input >
      template< pegtl::apply_mode A, template< typename ... > class Action, template< typename ... > class Control, typename Input >
         static bool match( Input &, std::vector<MaxPotInput>& input )
         {
            auto& table = input.back().function_tables_.back();
            if(table.back() + 1 == table.size()) {
               table.resize(table.size()-1); // remove last entry which holds the end size of the table
               return true;
            } else {
               return false;
            }
         }
   };
   struct function_table : pegtl::seq< new_function_table, opt_invisible, pegtl::until< function_table_end, opt_invisible, function_table_entry >, opt_invisible > {};
   struct function_tables : pegtl::seq< opt_invisible, pegtl::until< function_tables_end, function_table >, opt_invisible > {};//,pegtl::seq<pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>, real_number, pegtl::plus<pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>, real_number>> {};
   //template<> struct control< function_table > : pegtl::change_state_and_action< function_table, ..., object_action > {};

   struct grammar : pegtl::seq<pegtl::seq<init_line_markov, pegtl::eol,
                    number_of_variables, pegtl::eol,
                    pegtl::plus< cardinality >, pegtl::eol,
                    number_of_cliques, pegtl::eol,
                    clique_scopes,
                    opt_invisible,
                    function_tables>, 
                    pegtl::plus<
                    init_line_max_pot, pegtl::eol,
                    number_of_variables, pegtl::eol,
                    pegtl::plus< cardinality >, pegtl::eol,
                    number_of_cliques, pegtl::eol,
                    clique_scopes,
                    opt_invisible,
                    function_tables
                   >> {};


   template< typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<> struct action<init_line_max_pot > {
      template<typename Input>
      static void apply(const Input& in, std::vector<MaxPotInput>& input) 
      {
         input.resize(input.size() + 1);
      }
   };

   template<> struct action<init_line_markov > {
      template<typename Input>
      static void apply(const Input& in, std::vector<MaxPotInput>& input) 
      {
         input.resize(input.size() + 1);
      }
   };

   template<> struct action< number_of_variables > {
      template<typename Input>
      static void apply(const Input& in, std::vector<MaxPotInput>& input) 
      {
         input.back().number_of_variables_ = std::stoul(in.string());
      }
   };

   template<> struct action< number_of_cliques > {
      template<typename Input>
      static void apply(const Input & in, std::vector<MaxPotInput>& input)
      {
         input.back().number_of_cliques_ = std::stoul(in.string()); 
      }
   };

   template<> struct action< cardinality > {
      template<typename Input>
      static void apply(const Input & in, std::vector<MaxPotInput>& input)
      {
         input.back().cardinality_.push_back(std::stoul(in.string()));
      }
   };

   template<> struct action< new_clique_scope > {
      template<typename Input>
      static void apply(const Input &, std::vector<MaxPotInput>& input)
      {
         input.back().clique_scopes_.push_back(std::vector<INDEX>(0));
      }
   };
   template<> struct action< clique_scope > {
      template<typename Input>
      static void apply(const Input & in, std::vector<MaxPotInput>& input)
      {
         input.back().clique_scopes_.back().push_back(std::stoul(in.string()));
         assert(input.back().clique_scopes_.back().back() < input.back().number_of_variables_);
      }
   };
   template<> struct action< new_function_table > {
      template<typename Input>
      static void apply(const Input & in, std::vector<MaxPotInput>& input)
      {
         const INDEX no_entries = std::stoul(in.string());
         std::vector<REAL> entries;
         entries.reserve(no_entries+1);
         entries.push_back(no_entries);
         input.back().function_tables_.push_back(std::move(entries));
      }
   };
   template<> struct action< function_table_entry > {
      template<typename Input>
      static void apply(const Input & in, std::vector<MaxPotInput>& input)
      {
         auto& table = input.back().function_tables_.back();
         table.push_back(std::stod(in.string()));
         std::swap(table.back(), *(table.rbegin()+1)); // exchange last element, which always denotes the final number of entries in the function table
      }
   };

   bool ParseProblem(const std::string& filename)
   {
      std::cout << "parsing " << filename << "\n";
      pegtl::file_parser problem(filename);
      std::vector<MaxPotInput> input;
      bool read_suc = problem.parse< grammar, action >(input);
      return read_suc;
   }

    template<typename MRF_CONSTRUCTOR>
    void build_mrf(MRF_CONSTRUCTOR& mrf, const MaxPotInput& input)
    {
        assert(input.number_of_cliques_ == input.clique_scopes_.size());
        assert(input.number_of_cliques_ == input.function_tables_.size());

        // only unary and pairwise potentials supported right now
        for(INDEX i=0; i<input.number_of_cliques_; ++i) {
            assert(input.clique_scopes_[i].size() < 3);
        }

        UaiFileGraph traversalOrder = UaiFileGraph(input.number_of_variables_);
        for (std::size_t i=0; i < input.number_of_variables_; ++i) {
            const auto noLabels = input.cardinality_[i];
            traversalOrder.AddNode(i, noLabels);
        }
        
        for (std::size_t i=0; i<input.number_of_cliques_; ++i) {
            if (input.clique_scopes_[i].size() == 2) 
                traversalOrder.AddEdge(input.clique_scopes_[i][0], input.clique_scopes_[i][1]);
        }

        std::vector<std::size_t> order = traversalOrder.getTraversalOrder();
        std::vector<std::size_t> origNodeTraversalPriority(input.number_of_variables_);
        // first input the unaries, as pairwise potentials need them to be able to link to them
        // add unary factors with cost zero for each variables. There are models where unaries are not explicitly added.
        // std::size_t i = 0;
        // for (const auto& currentNodeIndex : order)
        // {
        //     const INDEX noLabels = input.cardinality_[currentNodeIndex];
        //     mrf.add_unary_factor_with_relation(std::vector<REAL>(noLabels,0.0)); 
        //     originalNodeIndexTraversalPriority[currentNodeIndex] = i;
        //     i++;
        // }

        INDEX p = 0;
        for(INDEX i=0; i<input.number_of_variables_; ++i) {
            const INDEX noLabels = input.cardinality_[i];
            std::vector<REAL> pots(noLabels, 0.0);
            mrf.add_unary_factor_without_relation(pots.begin(), pots.end());
            origNodeTraversalPriority[order[i]] = i;
            // mrf.add_unary_factor(pots.begin(), pots.end());
        }

        for (std::size_t i = 0; i < order.size() - 1; i++ ) {
            mrf.add_unary_unary_relation(order[i], order[i + 1]);
        }

        REAL initial_lb = 0.0;
        for(INDEX i=0; i < input.number_of_cliques_; ++i) {
            if (input.clique_scopes_[i].size() == 1) {
                const INDEX var = input.clique_scopes_[i][0];
                //std::cout << "unary potential for variable " << var << ":\n";
                auto* f = mrf.get_unary_factor(var);
                assert(input.function_tables_[i].size() == input.cardinality_[var]);
                initial_lb += *std::min_element(input.function_tables_[i].begin(), input.function_tables_[i].end());
                for(INDEX x=0; x<input.function_tables_[i].size(); ++x) {
                    assert( (*f->GetFactor())[x] == 0.0);
                    (*f->GetFactor())[x] = input.function_tables_[i][x];
                }
            }
        }

        //std::cout << "initial lower bound unaries = " << initial_lb << "\n"; 

        // now the pairwise potentials. 
        for(INDEX i=0; i<input.number_of_cliques_; ++i) {
        if(input.clique_scopes_[i].size() == 2) {
            const INDEX var1 = input.clique_scopes_[i][0];
            const INDEX var2 = input.clique_scopes_[i][1];

            const INDEX dim1 = mrf.get_number_of_labels(var1);
            const INDEX dim2 = mrf.get_number_of_labels(var2);
            assert(var1<var2 && var2 < input.number_of_variables_);
            assert(input.function_tables_[i].size() == input.cardinality_[var1]*input.cardinality_[var2]);
            assert(input.function_tables_[i].size() == dim1*dim2);
            matrix<REAL> pairwise_cost(dim1,dim2);
            initial_lb += *std::min_element(input.function_tables_[i].begin(), input.function_tables_[i].end());
            //std::cout << "pairwise potential on (" << var1 << "," << var2 << "):\n";
            for(INDEX l1=0; l1<dim1; ++l1) {
                for(INDEX l2=0; l2<dim2; ++l2) {
                    pairwise_cost(l1,l2) = input.function_tables_[i][l1*dim2 + l2];
            //      std::cout << input.function_tables_[i][l1*dim2 + l2] << " ";
                }
            //   std::cout << "\n";
            }
            //std::cout << pairwise_cost;
            mrf.add_pairwise_factor_without_relation(var1, var2, pairwise_cost);
            if (origNodeTraversalPriority[var1] < origNodeTraversalPriority[var2])
                mrf.add_unary_pairwise_relation(var1, var2);
            else
                mrf.add_unary_pairwise_relation(var2, var1);
        }
        }

        /*
        initial_lb = 0.0;
        for(INDEX i=0; i<input.number_of_cliques_; ++i) {
        if(input.clique_scopes_[i].size() == 1) {
            const INDEX var = input.clique_scopes_[i][0];
            auto* f = mrf.get_unary_factor(var);
            initial_lb += f->LowerBound();
        }
        }
        std::cout << "initial lower bound unaries = " << initial_lb << "\n"; 

        // now the pairwise potentials. 
        for(INDEX i=0; i<input.number_of_cliques_; ++i) {
        if(input.clique_scopes_[i].size() == 2) {
            const INDEX var1 = input.clique_scopes_[i][0];
            const INDEX var2 = input.clique_scopes_[i][1];
            auto* f = mrf.get_pairwise_factor(var1,var2);
            initial_lb += f->LowerBound();
        }
        }

        std::cout << "initial lower bound = " << initial_lb << "\n"; 
        */
    }
}
}
#endif // LP_MP_HORIZON_TRACKING_UAI_CONSTRUCTOR_HXX
