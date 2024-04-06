#ifndef DF_CAUCHY_H
#define DF_CAUCHY_H

#include "distribution.h"

namespace DiceForge {
    /// @brief DiceForge::Cauchy - A Continuous Probability Distribution (Cauchy) 
    class Cauchy : public Continuous {
        private:
            real_t x0, gamma, inv_gamma;
        public:
            /// @brief Initializes the Cauchy distribution about location x = x0 with scale gamma
            /// @param x0 centre of the distribution
            /// @param gamma scale factor of the distribution
            /// @note gamma > 0
            Cauchy(real_t x0 = 0, real_t gamma = 1);
            /// @brief Returns the next value of the random variable described by the distribution
            /// @param r A random real number uniformly distributed between 0 and 1
            real_t next(real_t r);
            /// @brief Returns the theoretical variance of the distribution
            /// @note The variation of a Cauchy distribution is undefined
            /// @returns NaN
            real_t variance() const override final;
            /// @brief Returns the theoretical expectation value of the distribution
            /// @note The expectation value of a Cauchy distribution is undefined
            /// @returns NaN
            real_t expectation() const override final;
            /// @brief Returns the minimum possible value of the random variable described by the distribution
            /// @returns negative infinity
            real_t minValue() const override final;
            /// @brief Returns the maximum possible value of the random variable described by the distribution
            /// @returns positive infinity
            real_t maxValue() const override final;
            /// @brief Probability density function of the Cauchy distribution
            real_t pdf(real_t x) const override final;
            /// @brief Cumulative distribution function of the Cauchy distribution
            real_t cdf(real_t x) const override final;
            /// @brief Returns x0 (centre of the distribution) 
            real_t get_x0() const;
            /// @brief Returns gamma (scale factor of the distribution) 
            real_t get_gamma() const;
    };

    /// @brief Fits the given sample points (x, y=pdf(x)) to a Cauchy distribution using interquartile estimations followed by
    /// non-linear least squares regression using modified Gauss-Newton
    /// @param x list of x coordinates
    /// @param y list of corresponding y coordinates where y = pdf(x)
    /// @param max_iter maximum iterations to attempt to fit the data (higher to try for better fits)
    /// @param epsilon minimum acceptable error tolerance while attempting to fit the data (smaller to try for better fits)
    /// @return A Cauchy distribution fit to the given sample points
    Cauchy fitToCauchy(std::vector<real_t> x, std::vector<real_t> y, int max_iter = 10000, real_t epsilon = 1e-6);
}


#endif
