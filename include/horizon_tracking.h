#ifndef LP_MP_HORIZON_TRACKING_H
#define LP_MP_HORIZON_TRACKING_H

#include "factors_messages.hxx"
#include "LP_MP.h"
#include "simplex_factor.hxx"
#include "simplex_marginalization_message.hxx"
#include "horizon_tracking_factors.hxx"
#include "mrf_problem_construction.hxx"
#include "horizon_tracking_constructor.hxx"

#include "parse_rules.h"

namespace LP_MP {

struct FMC_HORIZON_TRACKING {
   constexpr static const char* name = "Horizon tracking problem";

   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_HORIZON_TRACKING, 0, true >;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_HORIZON_TRACKING, 1, false >;
   using max_chain_container = FactorContainer<max_potential_on_chain, FMC_HORIZON_TRACKING, 2, false>;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_HORIZON_TRACKING, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_HORIZON_TRACKING, 1 >;
   using pairwise_max_message_container = MessageContainer<pairwise_max_factor_message, 1, 2, message_passing_schedule::left, 1, variableMessageNumber, FMC_HORIZON_TRACKING, 2>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, max_chain_container >;
   using MessageList = meta::list< UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, pairwise_max_message_container >;

   using mrf_constructor = StandardMrfConstructor<FMC_HORIZON_TRACKING,0,1,0,1>;
   using constructor = max_chain_constructor<mrf_constructor, max_chain_container, pairwise_max_message_container>;
   using ProblemDecompositionList = meta::list<constructor>;
};

}

#endif // LP_MP_HORIZON_TRACKING_H
