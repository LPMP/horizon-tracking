#ifndef LP_MP_HORIZON_TRACKING_CONSTRUCTOR_HXX
#define LP_MP_HORIZON_TRACKING_CONSTRUCTOR_HXX

#include "vector.hxx"
#include "parse_rules.h"
#include "pegtl/parse.hh"
#include "tree_decomposition.hxx"
#include "arboricity.h"
// #include "mrf_problem_construction.hxx"
#include <string>

namespace LP_MP {

template<typename MRF_CONSTRUCTOR, typename MAX_FACTOR, typename PAIRWISE_MAX_FACTOR_MESSAGE>
class max_chain_constructor : public MRF_CONSTRUCTOR {
public:
using mrf_constructor = MRF_CONSTRUCTOR;
using max_factor_container = MAX_FACTOR;
using pairwise_max_factor_message_container = PAIRWISE_MAX_FACTOR_MESSAGE;

using mrf_constructor::mrf_constructor;

template<typename ITERATOR>
max_factor_container* add_max_chain(ITERATOR var_begin, ITERATOR var_end, 
                                    const tensor3_variable<REAL>& maxPairwisePotentials, 
                                    const tensor3_variable<REAL>& linearPairwisePotentials)
{
    const INDEX no_vars = std::distance(var_begin, var_end);
    std::vector<INDEX> no_labels;
    for(auto it = var_begin; it!=var_end; ++it) {
        const INDEX i = (*it);
        no_labels.push_back( this->GetNumberOfLabels(i) );
    }
    auto* f = new max_factor_container(maxPairwisePotentials, linearPairwisePotentials, no_labels);
    this->lp_->AddFactor(f);

    INDEX c=0;
    for(auto it = var_begin; it+1!=var_end; ++it, ++c) {
        const INDEX i = (*it);
        const INDEX j = (*it+1);
        auto* p = this->GetPairwiseFactor(i,j);
        this->lp_.template add_message<pairwise_max_factor_message_container>(p, f, c);
    }

    return f;
}
};

namespace UAIMaxPotInput {

   struct MaxPotInput {
      INDEX number_of_variables_;
      INDEX number_of_cliques_;
      std::vector<INDEX> cardinality_;
      std::vector<std::vector<INDEX>> clique_scopes_;
      std::vector<std::vector<REAL>> function_tables_;
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
      // LP_MP::UaiMrfInput::ParseProblem(filename, solver?);
      std::cout << "parsing " << filename << "\n";
      pegtl::file_parser problem(filename);
      std::vector<MaxPotInput> input;
      bool read_suc = problem.parse< grammar, action >(input);
      return read_suc;
   }

   bool ParseProblemGridAndDecomposeToChains(const std::string& filename)
   {
      // LP_MP::UaiMrfInput::ParseProblem(filename, solver?);
        std::cout << "parsing " << filename << "\n";
        pegtl::file_parser problem(filename);
        std::vector<MaxPotInput> input;
        bool read_suc = problem.parse< grammar, action >(input);
        assert(read_suc);
        assert(input.size() == 2); // One max potential and one linear potential field
        INDEX numNodes = input[0].number_of_variables_;
        INDEX xEdgeDistance = 1;
        INDEX yEdgeDistance = 0;
        INDEX horizLastNode = 0;
        INDEX horizLength = 0;
        INDEX vertLastNode = 0;
        INDEX vertLength = 0;
        for (const auto& currentCliqueScope : input[0].clique_scopes_)
        {
            if (currentCliqueScope.size() == 1)
                continue;

            assert(currentCliqueScope.size() <= 2);
            INDEX currentEdgeDistance = currentCliqueScope[1] - currentCliqueScope[0];
            if (yEdgeDistance == 0 && currentEdgeDistance != xEdgeDistance)
                yEdgeDistance = currentEdgeDistance;
                
            else if (currentEdgeDistance != xEdgeDistance)
                assert(yEdgeDistance == currentEdgeDistance);
            
            if (currentCliqueScope[0] == horizLastNode && currentEdgeDistance == xEdgeDistance)
            {
                horizLength++;
                horizLastNode = currentCliqueScope[1];
            }

            if (currentCliqueScope[0] == vertLastNode && currentEdgeDistance == yEdgeDistance)
            {
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
        for (INDEX currentChain = 0; currentChain < numChains; currentChain++)
        {
            if (currentChain < gridSizeY)
            {
                for (INDEX i = 0; i < gridSizeX; i++)
                    chains[currentChain].insert(i + currentChain * gridSizeX);
            }

            else
            {
                for (INDEX i = 0; i < gridSizeY; i++)
                    chains[currentChain].insert(i * gridSizeX + currentChain - gridSizeY);
            }
        }

        for (const auto& currentChain : chains)
        {
            std::vector<std::vector<INDEX>> functionTableSizes;

            for (auto currentNodeItr = currentChain.begin(); currentNodeItr != std::prev(currentChain.end()); ++currentNodeItr)
            {
                INDEX l1Size = input[0].cardinality_[*currentNodeItr];
                INDEX l2Size = input[0].cardinality_[*std::next(currentNodeItr, 1)];
                //Assuming equivalent max potential and linear potential graphs:
                assert(input[0].cardinality_[*currentNodeItr] == input[1].cardinality_[*currentNodeItr]);
                assert(input[0].cardinality_[*std::next(currentNodeItr, 1)] == input[1].cardinality_[*std::next(currentNodeItr, 1)]);
                functionTableSizes.push_back(std::vector<INDEX>{l1Size, l2Size});
            }

            tensor3_variable<REAL> linearPairwisePotentials(functionTableSizes.begin(), functionTableSizes.end());
            tensor3_variable<REAL> maxPairwisePotentials(functionTableSizes.begin(), functionTableSizes.end());
            
            for (INDEX currentPotentialsIndex = 0; currentPotentialsIndex < input.size(); currentPotentialsIndex++)
            { 
                INDEX cliqueIndexChain = 0;
                for (INDEX cliqueIndex = 0; cliqueIndex < input[currentPotentialsIndex].clique_scopes_.size(); cliqueIndex++)
                {
                    const auto& currentCliqueScope = input[currentPotentialsIndex].clique_scopes_[cliqueIndex];

                    // node potentials:
                    if (currentCliqueScope.size() == 1)
                    {
                        // Check if the current node is present in the current chain, then transfer the unaries to pairwise terms
                        if (currentChain.count(currentCliqueScope[0]) > 0)
                        {
                            INDEX chainDelta = *std::next(currentChain.begin(), 1) - *currentChain.begin();
                            INDEX edgeIndexToRight = (currentCliqueScope[0] - *currentChain.begin()) / chainDelta;
                            if (edgeIndexToRight < currentChain.size() - 1)
                            {
                                // Send to the edges on the right:
                                INDEX numStatesCurrentNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0]];
                                INDEX numStatesNextNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0] + chainDelta];
                                for (INDEX l1 = 0; l1 < numStatesCurrentNode; l1++)
                                {
                                    for (INDEX l2 = 0; l2 < numStatesNextNode; l2++)
                                    {
                                        if (currentPotentialsIndex == 0)
                                            linearPairwisePotentials(edgeIndexToRight, l1, l2) = input[currentPotentialsIndex].function_tables_[cliqueIndex][l1];
                                        else
                                            maxPairwisePotentials(edgeIndexToRight, l1, l2) = input[currentPotentialsIndex].function_tables_[cliqueIndex][l1];
                                    }
                                }
                            }
                            else
                            {
                                // Send to the left edges for the corner nodes:
                                INDEX edgeIndexToLeft = edgeIndexToRight - 1;
                                INDEX numStatesCurrentNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0]];
                                INDEX numStatesPrevNode = input[currentPotentialsIndex].cardinality_[currentCliqueScope[0] - chainDelta];
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
                        continue;   
                    }

                    // edge potentials:
                    if (currentChain.count(currentCliqueScope[0]) == 0 || 
                        currentChain.count(currentCliqueScope[1]) == 0) 
                        continue;   // Current clique is not present in the current chain.

                    INDEX dim23InputIndex = 0;
                    for (INDEX l1 = 0; l1 < input[currentPotentialsIndex].cardinality_[currentCliqueScope[0]]; l1++)
                    {
                        for (INDEX l2 = 0; l2 < input[currentPotentialsIndex].cardinality_[currentCliqueScope[1]]; l2++)
                        {
                            if (currentPotentialsIndex == 0)
                                linearPairwisePotentials(cliqueIndexChain, l1, l2) += input[currentPotentialsIndex].function_tables_[cliqueIndex][dim23InputIndex];
                            else
                                maxPairwisePotentials(cliqueIndexChain, l1, l2) += input[currentPotentialsIndex].function_tables_[cliqueIndex][dim23InputIndex];

                            dim23InputIndex++;
                        }
                    }
                    cliqueIndexChain++;
                }
            }

            max_chain_constructor newChain;
            newChain.add_max_chain(currentChain.begin(), currentChain.end(), maxPairwisePotentials, linearPairwisePotentials);
        }

        return read_suc;
   }
}
}


#endif // LP_MP_HORIZON_TRACKING_CONSTRUCTOR_HXX
