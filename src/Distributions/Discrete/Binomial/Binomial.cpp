#include "Binomial.h"
#include "basicfxn.h"
#include<algorithm>

namespace DiceForge {
Binomial::Binomial(uint_t n, real_t p)
    : n(n), p(p) {
    if (n < 0 || p < 0 || p > 1) {
        throw std::invalid_argument("Expected n > 0 and 0 <= p <= 1");
    }

    pmfs = std::vector<real_t>(n+1);
    cdfs = std::vector<real_t>(n+1);
    real_t cumulative = 0.0;

    for (int i = 0; i <= n; i++) {
        pmfs[i] = nCr(n, i) * std::pow(p, i) * std::pow((1-p), n-i);
        cumulative += pmfs[i];
        cdfs[i] = cumulative;
    }
    cdfs[n] = 1.0; // Incase of numerical errors.
}

int_t Binomial::next(real_t r) {
    auto it = std::lower_bound(cdfs.begin(), cdfs.end(), r);
    return static_cast<int_t>(std::distance(cdfs.begin(), it));
}

real_t Binomial::variance() const {
    return n * p * (1-p);
}

real_t Binomial::expectation() const {
    return n * p;
}

int_t Binomial::minValue() const {
    return 0;
}

int_t Binomial::maxValue() const {
    return n;
}

real_t Binomial::pmf(int_t k) const {
    if (k > n || k < 0)
    {
        return 0;
    }
    return pmfs[k];
}

real_t Binomial::cdf(int_t k) const {
    if(k<0) return 0.0;
    if(k>=n) return 1.0;
    return cdfs[k];
}
} // namespace DiceForge
