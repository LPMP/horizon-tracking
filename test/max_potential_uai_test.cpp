#include "test.h"
#include "horizon_tracking_constructor.hxx"

using namespace LP_MP;

int main()
{
    std::string path = "Chain-Small-5x1.uai";
    bool success = LP_MP::UAIMaxPotInput::ParseProblem(path);
    test(success == 1);
}
