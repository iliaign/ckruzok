
#include "ean_mean.hpp"
#include <cmath>

double ean_mean(const double E, const double e)
{
    return E - e * std::sin(E);
}
Ы