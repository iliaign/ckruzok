#include "ecsrt_v_ist.hpp"

#include <cmath>

double ecsrt_v_ist(const double E, const double e)
{
    return std::atan2(std::sqrt(1.0 - e * e) * std::sin(E), std::cos(E) - e);
}
