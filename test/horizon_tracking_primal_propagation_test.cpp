#include "horizon_tracking_test_helper.hxx"
#include "test.h"
#include "data_max_potential_grids_test.hxx"
#include "graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "horizon_tracking.h"
#include "LP_FWMAP.hxx"

using namespace LP_MP;

std::vector<std::string> solver_options_primal_prop = {
   {"max potential grid test"},
   {"--maxIter"}, {"10"},
   {"--timeout"}, {"60"}, // half minute
   {"--roundingReparametrization"}, {"damped_uniform"},
   {"--standardReparametrization"}, {"damped_uniform"},
   {"-v"}, {"2"}
};

int main(int argc, char** argv)
{
    using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING_CHAINS>, StandardVisitor>;
    solver_type solver(solver_options_primal_prop);
    compute_lower_bound_chains(solver, grid_uai_input_medium, grid_uai_input_medium_lb, false);

    solver.GetLP().ComputePass(1);
    auto& constructor = solver.template GetProblemConstructor<0>();
    std::vector<std::size_t> labeling;
    for(std::size_t i=0; i<constructor.get_number_of_variables(); ++i) {
	    auto* f = constructor.get_unary_factor(i);
	    std::size_t minimal_label = std::min_element(f->GetFactor()->begin(), f->GetFactor()->end()) - f->GetFactor()->begin();
	    labeling.push_back(minimal_label);
    }

	solver.GetLP().write_back_reparametrization();
    for (auto i = 0; i < solver.GetLP().GetNumberOfFactors(); i++) {
		solver.GetLP().GetFactor(i)->init_primal(); //TODO: DO OR DO NOT?
	}

    for(std::size_t i=0; i<constructor.get_number_of_variables(); ++i) {
	    auto* f = constructor.get_unary_factor(i);
	    f->GetFactor()->primal() = labeling[i];
	    f->propagate_primal_through_messages();
	    // check consistency of pairwise factors
	    for(auto* m : f->template get_messages<typename FMC_HORIZON_TRACKING_CHAINS::UnaryPairwiseMessageLeftContainer>()) {
		    assert(m->CheckPrimalConsistency());
	    }
	    for(auto* m : f->template get_messages<typename FMC_HORIZON_TRACKING_CHAINS::UnaryPairwiseMessageRightContainer>()) {
		    assert(m->CheckPrimalConsistency());
	    }
	    for(std::size_t p=0; p<constructor.get_number_of_pairwise_factors(); ++p) {
		    const auto idx = constructor.get_pairwise_variables(p);
		    if(idx[0] <= i && idx[1] <= i) {
			    auto* p_f = constructor.get_pairwise_factor(p);
			    assert(p_f->GetFactor()->primal()[0] == labeling[idx[0]]);
			    assert(p_f->GetFactor()->primal()[1] == labeling[idx[1]]);
			    auto msgs_to_chain = p_f->get_messages<typename FMC_HORIZON_TRACKING_CHAINS::pairwise_max_message_container>();
			    assert(msgs_to_chain.size() == 1);
			    assert(msgs_to_chain[0]->CheckPrimalConsistency()); 
		    }
	    }
    }

    assert(solver.GetLP().CheckPrimalConsistency());
}



