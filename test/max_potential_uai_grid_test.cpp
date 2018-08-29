#include "test_helper.hxx"
#include "data_max_potential_grids_test.hxx"

using namespace LP_MP;

int main()
{
    // Test on 2x2 Grid:
    {
        TestUAI(solver_options_small, grid_uai_input_small, 0.242203);
    }

    // Test on 3x3 Grid:
    {
        TestUAI(solver_options_small, grid_uai_input_3x3, 0.590735);
    }

    // Test on 5x5 Grid:
    {
        TestUAI(solver_options_medium, grid_uai_input_medium, 4.307381);
    }
}
