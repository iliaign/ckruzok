#include "mean_motion.hpp"

#include <cmath>

double mean_motion(const double a, const double mu)
{
    return std::sqrt(mu / a) / a;
}
