#ifndef LP_MP_HORIZON_TRACKING_FACTORS_HXX
#define LP_MP_HORIZON_TRACKING_FACTORS_HXX

#include "vector.hxx"

namespace LP_MP {

    class max_factor {
        public:
            template<typename ITERATOR>
            max_factor(ITERATOR size_begin, ITERATOR size_end)
            :
                potential(size_begin, size_end),
                primal(std::distance(size_begin, size_end))
        {}
            
            REAL EvaluatePrimal() const
            {
                REAL cost = 0.0;
                for(INDEX i=0; i<potential.size(); ++i) {
                    cost += potential[i][ primal[i] ];
                }
                return cost;
            }

            REAL LowerBound() const
            {
                REAL lb = -std::numeric_limits<REAL>::infinity();
                for(INDEX i=0; i<potential.size(); ++i) {
                    const REAL best_local = std::min(potential[i].begin(), potential[i].end());
                    lb = std::max(best_local, lb);
                }

                return lb;
            }

            auto export_variables() { std::tie(potential); }

        private:
            two_dim_variable_array<REAL> potential;
            vector<INDEX> primal;

    };

    class pairwise_max_factor_message {
        public:
            pairwise_max_factor_message(const INDEX _entry) : entry(_entry) {}

            template<typename FACTOR, typename MSG>
            void RepamRight(FACTOR& r, const MSG& msgs)
            {
                for(INDEX i=0; i<r[entry].size(); ++i) {
                    r[entry][i] += msgs[i];
                }
            }

            template<typename FACTOR, typename MSG>
            void RepamLeft(FACTOR& l, const MSG& msgs)
            {
                INDEX c=0;
                for(i=0; i<l.dim1(); ++i) {
                    for(j=0; j<l.dim2(); ++j) {
                        l.cost(i,j) += msgs[c++];
                    }
                } 
            }

            template<typename LEFT_FACTOR, typename MSG>
            void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
            {
                vector<REAL> m(l.size());
                INDEX c=0;
                for(i=0; i<l.dim1(); ++i) {
                    for(j=0; j<l.dim2(); ++j) {
                        m[c++] = l(i,j);
                    }
                }
                msg -= omega*m; 
            }

        private:
            const INDEX entry;
    };

}

#endif // LP_MP_HORIZON_TRACKING_FACTORS_HXX
