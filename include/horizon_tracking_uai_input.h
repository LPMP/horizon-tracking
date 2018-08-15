#ifndef LP_MP_HORIZON_TRACKING_UAI_INPUT_H
#define LP_MP_HORIZON_TRACKING_UAI_INPUT_H

#include <vector>
#include "mrf_input.h"

namespace LP_MP {

    struct horizon_tracking_input {
        mrf_input mrf;
        std::vector<mrf_input> bottleneck_potentials;
    };

    namespace horizon_tracking_uai_input {

        horizon_tracking_input parse_file(const std::string& filename);
        horizon_tracking_input parse_string(const std::string& filename);

    } // namespace horizon_tracking_input

} // namespace LP_MP 

#endif // LP_MP_HORIZON_TRACKING_UAI_INPUT_H
