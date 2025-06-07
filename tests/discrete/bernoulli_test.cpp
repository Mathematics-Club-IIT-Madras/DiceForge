#include <gtest/gtest.h>
#include "diceforge_distributions.h" 
#include "diceforge_generators.h" 
#include <numeric>
#include <vector>

TEST(BernoulliDistribution, MeanAndVariance) {
    const double p = 0.3;
    const int samples = 1000000;
    double tolerance_frac = 0.005;  
    
    DiceForge::Bernoulli dist(p);
    DiceForge::MT64 rng(time(NULL));
    
    double mean = 0.0, var = 0.0;
    
    for (int i = 0; i < samples; i++) {
        int val = dist.next(rng.next_unit());
        mean += val;
        var += val;  
    }
    
    mean /= samples;
    var /= samples;
    var -= (mean * mean);
    
    EXPECT_NEAR(mean, p, p * tolerance_frac);
    EXPECT_NEAR(var, p * (1 - p), p * (1 - p) * tolerance_frac);
}

TEST(BernoulliDistribution, ChiSquareTest) {
    const double p = 0.3;
    const int samples = 100000;
    
    DiceForge::Bernoulli dist(p);
    DiceForge::MT64 rng(time(NULL));
    
    int one_freq = 0, zero_freq = 0;
    std::vector<int> observed(2, 0); 
    for (int i = 0; i < samples; i++) {
        int val = dist.next(rng.next_unit());
        if(val == 0) zero_freq++;
        else one_freq++;
    }
    double chi_square = 0.0;
    int valid_bins = 0;
    double expected = samples * (1-p);
    double diff = zero_freq - expected;
    chi_square += (diff*diff) / expected;
    expected = samples * p;
    diff = one_freq - expected;
    chi_square += (diff*diff)/expected;
    //int df = valid_bins - 1;  
    double critical_value = 10.828;
    EXPECT_LT(chi_square, critical_value);
}
