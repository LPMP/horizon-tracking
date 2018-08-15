#include "horizon_tracking_uai_input.h"
#include "mrf_uai_input_impl.hxx"
#include <cassert>

namespace LP_MP {

namespace horizon_tracking_uai_input {

    // import basic parsers
    using Parsing::opt_whitespace;
    using Parsing::mand_whitespace;
    using Parsing::opt_invisible;
    using Parsing::mand_invisible;
    using Parsing::positive_integer;
    using Parsing::real_number; 

    struct mrf_init : pegtl::string<'M','A','R','K','O','V'> {};
    struct mrf_init_line : pegtl::seq<opt_whitespace, mrf_init, opt_whitespace, pegtl::eol> {};

    struct bottleneck_potential_init : pegtl::string<'M','A','X','-','P','O','T','E','N','T','I','A','L','S'> {};
    struct bottleneck_potential_init_line : pegtl::seq<opt_whitespace, bottleneck_potential_init, opt_whitespace, pegtl::eol> {};

    struct mrf : pegtl::seq<
                 mrf_init_line,
                 pegtl::star<pegtl::not_at<bottleneck_potential_init_line>, pegtl::any>
                 > {};

    struct bottleneck_potential : pegtl::seq<
                                  bottleneck_potential_init_line,
                                  pegtl::star<pegtl::not_at<bottleneck_potential_init_line>, pegtl::any>
                                  > {};

    struct grammar : pegtl::seq<
                     mrf,
                     pegtl::star<bottleneck_potential>
                     > {};

    template< typename Rule >
        struct action
        : pegtl::nothing< Rule > {};

    template<> struct action< mrf > {
        template<typename INPUT>
            static void apply(const INPUT& in, horizon_tracking_input& input)
            {
                input.mrf = mrf_uai_input::parse_string<mrf_init>(in.string());
            }
    };

    template<> struct action< bottleneck_potential > {
        template<typename INPUT>
            static void apply(const INPUT& in, horizon_tracking_input& input)
            {
                input.bottleneck_potentials.push_back( mrf_uai_input::parse_string<bottleneck_potential_init>(in.string()) );
            }
    };

    horizon_tracking_input parse_file(const std::string& filename)
    {
        horizon_tracking_input input;
        pegtl::file_parser problem(filename);

        const bool read_success = problem.parse<grammar, action>(input);
        assert(read_success);
        return input; 
    }

    horizon_tracking_input parse_string(const std::string& input_string)
    {
        horizon_tracking_input input;
        const bool read_success = pegtl::parse<grammar, action>(input_string, "", input);
        assert(read_success);
        return input; 
    }

} // namespace horizon_tracking_input

} // namespace LP_MP 
