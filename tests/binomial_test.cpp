#include <gtest/gtest.h>
#include "diceforge_distributions.h" 
#include "diceforge_generators.h" 
#include <numeric>
#include <vector>

TEST(BinomialDistribution, MeanandVariance) {
    const int n = 1000, p = 0.5;
    const int samples = 1000000;
    double tolerance_frac = 0.00005;
    DiceForge::Binomial dist(n, p);
    DiceForge::MT64 rng(time(NULL));
    long long mean = 0.0, var = 0.0; 
    for (int i = 0; i < samples; i++) {
        int val = dist.next(rng.next_unit());
        mean += val;
        var += (val*val);
    }
    mean /= samples;
    var /= samples;
    var -= (mean*mean);
    EXPECT_NEAR(mean, n*p, n*p*tolerance_frac);
    EXPECT_NEAR(var, n*p*(1-p), n*p*(1-p)*tolerance_frac);
}

TEST(BinomialDistribution, ChiSquareTest) {
    const int n = 100;
    const double p = 0.5;
    const int samples = 100000;
    DiceForge::Binomial dist(n, p);
    DiceForge::MT64 rng(time(NULL));
    std::vector<int> observed(n + 1, 0);
    for (int i = 0; i < samples; i++) {
        int val = dist.next(rng.next_unit());
        if (val >= 0 && val <= n) {
            observed[val]++;
        }
    }
    double chi_square = 0.0;
    int valid_bins = 0;
    for (int k = 0; k <= n; ++k) {
        double expected = samples * dist.pmf(k);
        if (expected >= 5.0) {  
            double diff = observed[k] - expected;
            chi_square += (diff * diff) / expected;
            valid_bins++;
        }
    }
    int df = valid_bins - 1;
    double critical_value = 135.81; 
    EXPECT_LT(chi_square, critical_value);
}


