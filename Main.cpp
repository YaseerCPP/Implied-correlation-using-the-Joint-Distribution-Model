#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// Define constants
const double PI = 3.14159265358979323846;

// Helper functions for the normal distribution
double normalPdf(double x) {
    return (1.0 / sqrt(2 * PI)) * exp(-0.5 * x * x);
}

// Pricing function for a simple European call option
double callOptionPrice(double S, double K, double T, double r, double sigma) {
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    return S * normalPdf(d1) - K * exp(-r * T) * normalPdf(d2);
}

// Objective function to minimize (squared difference between model and market prices)
double objectiveFunction(double correlation, const std::vector<double>& marketPrices,
                         const std::vector<double>& S1, const std::vector<double>& K1,
                         const std::vector<double>& S2, const std::vector<double>& K2,
                         double T, double r, double sigma1, double sigma2) {
    double sumSquaredError = 0.0;

    for (size_t i = 0; i < marketPrices.size(); ++i) {
        double modelPrice1 = callOptionPrice(S1[i], K1[i], T, r, sigma1);
        double modelPrice2 = callOptionPrice(S2[i], K2[i], T, r, sigma2);
        // Here, we simplify and assume the model price depends linearly on correlation
        double modelPrice = (modelPrice1 + modelPrice2) / 2; // This is a simplified assumption

        sumSquaredError += (modelPrice - marketPrices[i]) * (modelPrice - marketPrices[i]);
    }

    return sumSquaredError;
}

// Simple optimization using grid search
double optimizeCorrelation(const std::vector<double>& marketPrices,
                           const std::vector<double>& S1, const std::vector<double>& K1,
                           const std::vector<double>& S2, const std::vector<double>& K2,
                           double T, double r, double sigma1, double sigma2) {
    double bestCorrelation = 0.0;
    double minError = std::numeric_limits<double>::max();
    
    for (double correlation = -1.0; correlation <= 1.0; correlation += 0.01) {
        double error = objectiveFunction(correlation, marketPrices, S1, K1, S2, K2, T, r, sigma1, sigma2);
        if (error < minError) {
            minError = error;
            bestCorrelation = correlation;
        }
    }

    return bestCorrelation;
}

int main() {
    // Example market data
    std::vector<double> marketPrices = {10.0, 15.0, 12.0}; // Example market option prices
    std::vector<double> S1 = {100.0, 105.0, 110.0}; // Asset 1 prices
    std::vector<double> K1 = {95.0, 100.0, 105.0};  // Asset 1 strike prices
    std::vector<double> S2 = {200.0, 210.0, 220.0}; // Asset 2 prices
    std::vector<double> K2 = {190.0, 200.0, 210.0}; // Asset 2 strike prices

    double T = 1.0; // Time to maturity
    double r = 0.05; // Risk-free rate
    double sigma1 = 0.20; // Volatility of asset 1
    double sigma2 = 0.25; // Volatility of asset 2

    double impliedCorrelation = optimizeCorrelation(marketPrices, S1, K1, S2, K2, T, r, sigma1, sigma2);

    std::cout << "Implied Correlation: " << impliedCorrelation << std::endl;

    return 0;
}
