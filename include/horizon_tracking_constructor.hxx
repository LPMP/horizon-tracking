#ifndef LP_MP_HORIZON_TRACKING_CONSTRUCTOR_HXX
#define LP_MP_HORIZON_TRACKING_CONSTRUCTOR_HXX

#include ""

namespace LP_MP {

template<typename MRF_CONSTRUCTOR, typename MAX_FACTOR, typename PAIRWISE_MAX_FACTOR_MESSAGE>
class horizon_tracking_constructor : public MRF_CONSTRUCTOR {
public:
using mrf_constructor = MRF_CONSTRUCTOR;
using max_factor_container = MAX_FACTOR;
using pairwise_max_factor_message_container = PAIRWISE_MAX_FACTOR_MESSAGE;

using mrf_constructor::mrf_constructor;

template<typename ITERATOR>
max_factor_container* add_max_factor_on_pairwise(ITERATOR pairwise_begin, ITERATOR pairwise_end, LP_tree* t)
{
    const INDEX no_max_vars = std::distance(pairwise_begin, pairwise_end);
    vector<INDEX> no_labels(no_vars);
    for(auto it = pairwise_begin; it!=pairwise_end; ++it) {
        const INDEX i = (*it)[0];
        const INDEX i = (*it)[1];
        assert(i<j);
        no_labels[i] = this->GetNumberOfLabels(i)*this->GetNumberOfLabels(j);
    }
    auto* f = new max_factor_container(no_labels.begin(), no_labels.end());
    this->lp_->AddFactor(f);

    INDEX c=0;
    for(auto it = pairwise_begin; it!=pairwise_end; ++it) {
        const INDEX i = (*it)[0];
        const INDEX i = (*it)[1];
        auto* p = this->GetPairwiseFactor(i,j);
        auto* m = new pairwise_max_Factor_message_container(c++, p, f);
        this->lp_->AddMessage(m);

        if(t != nullptr) {
            t->AddMessage(m);
        }
    } 

    return f;
}



};

}

#endif // LP_MP_HORIZON_TRACKING_CONSTRUCTOR_HXX
