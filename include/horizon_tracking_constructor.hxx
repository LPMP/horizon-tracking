#ifndef LP_MP_HORIZON_TRACKING_CONSTRUCTOR_HXX
#define LP_MP_HORIZON_TRACKING_CONSTRUCTOR_HXX

#include "vector.hxx"
#include "parse_rules.h"
#include "pegtl/parse.hh"
#include "tree_decomposition.hxx"
#include "arboricity.h"
#include "mrf_problem_construction.hxx"
#include "tree_decomposition.hxx"
#include <string>

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
max_chain_factor_container* add_max_chain(ITERATOR var_begin, ITERATOR var_end,
                                    TENSORSIZE potSizeBegin, TENSORSIZE potSizeEnd, 
                                    const tensor3_variable<REAL>& maxPairwisePotentials,
                                    factor_tree<FMC>* t = nullptr)
{
    std::vector<INDEX> numLabels;
    for(auto it = var_begin; it!=var_end; ++it) {
        const INDEX i = (*it);
        numLabels.push_back( this->GetNumberOfLabels(i) );
    }
    tensor3_variable<REAL> linearPairwisePotentials(potSizeBegin, potSizeEnd);
    for(std::size_t n1=0; n1<maxPairwisePotentials.dim1(); ++n1) {
        for(std::size_t i=0; i<maxPairwisePotentials.dim2(n1); ++i) {
            for(std::size_t j=0; j<maxPairwisePotentials.dim3(n1); ++j) {
                linearPairwisePotentials(n1, i, j) = 0.0;
            }
        }
    }
    auto* f = this->lp_->template add_factor<max_chain_factor_container>(maxPairwisePotentials, linearPairwisePotentials, numLabels, false);

    INDEX c=0;
    for(auto it = var_begin; std::next(it, 1)!=var_end; ++it, ++c) {
        const INDEX i = (*it);
        const INDEX j = *std::next(it, 1);
        auto* p = this->GetPairwiseFactor(i,j);
        auto* m = this->lp_->template add_message<pairwise_max_factor_message_container>(p, f, c, c, c+1);

        if(t != nullptr) {
            t->add_message(m, Chirality::right); 
        }
    }

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

template<typename FMC, typename MAX_CHAIN_FACTOR, typename MAX_POTENTIAL_FACTOR, typename MAX_CHAIN_MAX_POTENTIAL_MESSAGE>
class max_disjoint_chain_constructor {
public:
using max_chain_factor_container = MAX_CHAIN_FACTOR;
using max_potential_factor_container = MAX_POTENTIAL_FACTOR;
using max_chain_max_potential_message_container = MAX_CHAIN_MAX_POTENTIAL_MESSAGE;

max_chain_factor_container* add_max_chain(const tensor3_variable<REAL>& linearPairwisePotentials, 
                                    const tensor3_variable<REAL>& maxPairwisePotentials,
                                    const std::vector<INDEX>& numLabels,
                                    factor_tree<FMC>* t = nullptr)
{
    auto* f = this->lp_->template add_factor<max_chain_factor_container>(maxPairwisePotentials, linearPairwisePotentials, numLabels, false);
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

        // first input the unaries, as pairwise potentials need them to be able to link to them
        // add unary factors with cost zero for each variables. There are models where unaries are not explicitly added.
        for(INDEX i=0; i<input.number_of_variables_; ++i) {
        const INDEX noLabels = input.cardinality_[i];
        mrf.AddUnaryFactor(i,std::vector<REAL>(noLabels,0.0));
        }

        REAL initial_lb = 0.0;
        for(INDEX i=0; i<input.number_of_cliques_; ++i) {
        if(input.clique_scopes_[i].size() == 1) {
            const INDEX var = input.clique_scopes_[i][0];
            //std::cout << "unary potential for variable " << var << ":\n";
            auto* f = mrf.GetUnaryFactor(var);
            assert(input.function_tables_[i].size() == input.cardinality_[var]);
            initial_lb += *std::min_element(input.function_tables_[i].begin(), input.function_tables_[i].end());
            for(INDEX x=0; x<input.function_tables_[i].size(); ++x) {
                //std::cout << input.function_tables_[i][x] << " ";
                assert( (*f->GetFactor())[x] == 0.0);
                (*f->GetFactor())[x] = input.function_tables_[i][x];
            }
            //std::cout << "\n";
        }
        }

        //std::cout << "initial lower bound unaries = " << initial_lb << "\n"; 

        // now the pairwise potentials. 
        for(INDEX i=0; i<input.number_of_cliques_; ++i) {
        if(input.clique_scopes_[i].size() == 2) {
            const INDEX var1 = input.clique_scopes_[i][0];
            const INDEX var2 = input.clique_scopes_[i][1];
            const INDEX dim1 = mrf.GetNumberOfLabels(var1);
            const INDEX dim2 = mrf.GetNumberOfLabels(var2);
            assert(var1<var2 && var2 < input.number_of_variables_);
            assert(input.function_tables_[i].size() == input.cardinality_[var1]*input.cardinality_[var2]);
            assert(input.function_tables_[i].size() == dim1*dim2);
            matrix<REAL> pairwise_cost(dim1,dim2);
            initial_lb += *std::min_element(input.function_tables_[i].begin(), input.function_tables_[i].end());
            //std::cout << "pairwise potential on (" << var1 << "," << var2 << "):\n";
            for(INDEX l1=0; l1<dim1; ++l1) {
                for(INDEX l2=0; l2<dim2; ++l2) {
                    pairwise_cost(l1,l2) = input.function_tables_[i][l2*dim1 + l1];
            //      std::cout << input.function_tables_[i][l2*dim1 + l1] << " ";
                }
            //   std::cout << "\n";
            }
            //std::cout << pairwise_cost;
            mrf.AddPairwiseFactor(var1,var2,pairwise_cost); // or do we have to transpose the values?
        }
        }

        /*
        initial_lb = 0.0;
        for(INDEX i=0; i<input.number_of_cliques_; ++i) {
        if(input.clique_scopes_[i].size() == 1) {
            const INDEX var = input.clique_scopes_[i][0];
            auto* f = mrf.GetUnaryFactor(var);
            initial_lb += f->LowerBound();
        }
        }
        std::cout << "initial lower bound unaries = " << initial_lb << "\n"; 

        // now the pairwise potentials. 
        for(INDEX i=0; i<input.number_of_cliques_; ++i) {
        if(input.clique_scopes_[i].size() == 2) {
            const INDEX var1 = input.clique_scopes_[i][0];
            const INDEX var2 = input.clique_scopes_[i][1];
            auto* f = mrf.GetPairwiseFactor(var1,var2);
            initial_lb += f->LowerBound();
        }
        }

        std::cout << "initial lower bound = " << initial_lb << "\n"; 
        */
    }

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

        auto* f = chain_constructor.add_max_chain(linearPairwisePotentials, maxPairwisePotentials, numLabels, &tree);
        max_chain_potentials.push_back(f);


        chain_constructor.add_max_potential(max_chain_potentials.begin(), max_chain_potentials.end(), &tree);
        tree.init();
        s.GetLP().add_tree(tree);

        return read_suc;
   }

   template<typename SOLVER>
   bool ParseProblemGridAndDecomposeToChains(const std::string& filename, SOLVER& s)
   {
        using FMC = typename SOLVER::FMC;
        auto& chain_constructor = s.template GetProblemConstructor<0>();

        std::cout << "parsing " << filename << "\n";
        pegtl::file_parser problem(filename);
        std::vector<MaxPotInput> input;
        bool read_suc = problem.parse< grammar, action >(input);
        assert(read_suc);
        assert(input.size() == 2); // One max potential and one linear potential field
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

        std::vector<typename std::remove_reference_t<decltype(chain_constructor)>::max_chain_factor_container*> max_chain_potentials;
        factor_tree<FMC> tree;
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

            tensor3_variable<REAL> maxPairwisePotentials(functionTableSizes.begin(), functionTableSizes.end());
            
            // Populate max potentials:
            for (INDEX currentPotentialsIndex = 1; currentPotentialsIndex < input.size(); currentPotentialsIndex++)
            { 
                INDEX cliqueIndexChain = 0;
                // process the pairwise potentials first:
                for (INDEX cliqueIndex = 0; cliqueIndex < input[currentPotentialsIndex].clique_scopes_.size(); cliqueIndex++)
                {
                    const auto& currentCliqueScope = input[currentPotentialsIndex].clique_scopes_[cliqueIndex];
                    if (currentCliqueScope.size() == 1)
                        continue; // Add unaries later.
                    
                    if (currentChain.count(currentCliqueScope[0]) == 0 || 
                        currentChain.count(currentCliqueScope[1]) == 0) 
                        continue;   // Current clique is not present in the current chain.

                    INDEX dim23InputIndex = 0;
                    for (INDEX l1 = 0; l1 < input[currentPotentialsIndex].cardinality_[currentCliqueScope[0]]; l1++)
                    {
                        for (INDEX l2 = 0; l2 < input[currentPotentialsIndex].cardinality_[currentCliqueScope[1]]; l2++)
                        {
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
                    
                    if (currentChain.count(currentCliqueScope[0]) == 0) 
                        continue;   // Current clique is not present in the current chain.

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
                                maxPairwisePotentials(edgeIndexToRight, l1, l2) += input[currentPotentialsIndex].function_tables_[cliqueIndex][l1];
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
                                maxPairwisePotentials(edgeIndexToLeft, l1, l2) += input[currentPotentialsIndex].function_tables_[cliqueIndex][l1];
                            }
                        }
                    }
                }
            }

            auto* f = chain_constructor.add_max_chain(currentChain.begin(), currentChain.end(), functionTableSizes.begin(), functionTableSizes.end(), maxPairwisePotentials, &tree);
            max_chain_potentials.push_back(f);

        }

        auto* f = chain_constructor.add_max_potential(max_chain_potentials.begin(), max_chain_potentials.end(), &tree);
        s.GetLP().add_tree(tree);
        s.GetLP().construct_decomposition();

        return read_suc;
   }

   
}
}


#endif // LP_MP_HORIZON_TRACKING_CONSTRUCTOR_HXX
