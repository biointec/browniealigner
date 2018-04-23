/***************************************************************************
 *   Copyright (C) 2014 - 2016 Jan Fostier (jan.fostier@intec.ugent.be)    *
 *   This file is part of Brownie                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <chrono>
#include <ctime>
#include <fstream>
#include <vector>
#include <map>
#define MAX_TIMERS 16

/**
 * Utility class with timers
 */
class Util
{
private:
        static int currentTimer;
        static std::chrono::time_point<std::chrono::system_clock> startTime[MAX_TIMERS];

public:
        /**
         * Create a string with a human readable version of a time period
         * @param time Time period (expressed in s)
         * @return String with a human readable version of a time period
         */
        static std::string humRead(double time);

        /**
         * Start a chronometer
         */
        static void startChrono();

        /**
         * Stop the chronometer
         * @return The time in
         */
        static double stopChrono();

        /**
         * Stop the chronometer and return a human readable string
         * @return A human readable string containg the elapsed time
         */
        static std::string stopChronoStr() {
                return humRead(stopChrono());
        }

        /**
         * Get a string containing the date and time
         * @return string containing date and time
         */
        static std::string getDateTime();

        /**
         * Compute the sensitivity
         * @param TN True negatives
         * @param FP False positives
         * The sensitivity
         */
        static double getSpecificity(double TN, double FP) {
                return ((TN+FP) == 0) ? 1.0 : TN / (TN + FP);
        }

        /**
         * Compute the specificity
         * @param TP True positives
         * @param FN False negatives
         * The specificity
         */
        static double getSensitivity(double TP, double FN) {
                return ((TP+FN) == 0) ? 1.0 : TP / (TP + FN);
        }

        static bool fileExists(const std::string& filename) {
                std::ifstream file(filename.c_str(), std::ios::in);
                bool OK = file.good();
                file.close();
                return OK;
        }

        /**
         * Compute the probability p(k) from a Poisson distribution with mean mu
         * @param k Number of observations
         * @param mu Average
         * @return The probability p(k)
         */
        static double poissonPDF(unsigned int k, double mu);

        /**
         * Compute the probability ratio p(k, mu1) / p(k, mu2)
         * @param k Number of observations
         * @param mu1 Expected number of observation (mean of distribution)
         * @param mu2 Expected number of observation (mean of distribution)
         * @return The probability p(k)
         */
        static double poissonPDFratio(unsigned int k, double mu1, double mu2);

        /**
         * Compute the probability p(k) from a negative bionomial(mu, var)
         * @param k Number of observations
         * @param mu Average
         * @param var Variance
         * @return The probability p(k)
         */
        static double negbinomialPDF(unsigned int k, double mu, double var);

        /**
         * Compute the log(probability p(k)) from a negative bionomial(mu, var)
         * @param k Number of observations
         * @param mu Average
         * @param var Variance
         * @return The probability p(k)
         */
        static double logNegbinomialPDF(unsigned int k, double mu, double var);

        /**
         * Compute the probability ratio p(k, mu1) / p(k, mu2)
         * @param k Number of observations
         * @param mu1 Mean of the distribution
         * @param var1 Variance of the first distribution
         * @param mu2 Mean of the second distribution
         * @param var2 Variance of the second distribution
         * @return The probability p(k)
         */
        static double negbinomialPDFratio(unsigned int k,
                                          double mu1, double var1,
                                          double mu2, double var2);

        /**
         * Compute the probability p(k) for the geometric distribution
         * @param k Number of observations
         * @param mu Mean of the distribution
         * @param mu2 Variance of the distribution
         * @return The probability p(k)
         */
        static double geometricPDF(unsigned int k, double mu);

        /**
         * Compute the log(probability p(k)) for the geometric distribution
         * @param k Number of observations
         * @param mu Mean of the distribution
         * @param mu2 Variance of the distribution
         * @return The probability p(k)
         */
        static double logGeometricPDF(unsigned int k, double mu);

        /**
         * Compute the probability ratio p(mu1) / p(k, mu2)
         * @param k Number of observations
         * @param mu1 Mean of the geometric distribution
         * @param mu2 Mean of the negative bionomial distribution
         * @param var2 Variance of the negative bionomial distribution
         * @return The probability p(k)
         */
        static double geometricnegbinomialPDFratio(unsigned int k, double mu1,
                                                   double mu2, double var2);

        /**
         * Compute the percentage of two size_t numbers
         * @param nom Nominator
         * @param den Denominator
         * @return The percentage
         */
        static double toPercentage(size_t nom, size_t den) {
                if (den == 0)
                        return 0;
                return 100.0 * double(nom) / double(den);
        }

        static void binomialMixtureEM(const std::map<unsigned int, double>& data,
                                      std::vector<double>& mu,
                                      std::vector<double>& sigma2,
                                      std::vector<double>& MC,
                                      int maxIterations = 20);
};

#endif
