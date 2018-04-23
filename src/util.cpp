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

#include "util.h"
#include "global.h"
#include <ctime>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace std::chrono;

time_point<system_clock> Util::startTime[MAX_TIMERS];
int Util::currentTimer = 0;

string Util::humRead(double time)
{
        uint64_t timeInt = uint64_t(time);

        uint64_t days = timeInt / 86400;
        timeInt -= days * 86400;
        uint64_t hours = timeInt / 3600;
        timeInt -= hours * 3600;
        uint64_t min = timeInt / 60;
        timeInt -= min * 60;
        uint64_t sec = timeInt;
        uint64_t ms = uint64_t((time - timeInt) * 1000);

        ostringstream iss;
        if (days > 0) {
                iss << days << "d " << hours << "h " << min << "min";
        } else if (hours > 0) {
                iss << hours << "h " << min << "min " << sec << "s";
        } else if (min > 0) {
                iss << min << "min " << sec << "s";
        } else if (sec > 0) {
                iss << sec << "s " << ms << "ms";
        } else {
                iss << ms << "ms";
        }

        return iss.str();
}

void Util::startChrono()
{
        // make sure we don't use too many timers
        assert(currentTimer < MAX_TIMERS);
        startTime[currentTimer] = system_clock::now();
        currentTimer++;
}

double Util::stopChrono()
{
        // make sure stopChrono isn't called too often
        currentTimer--;
        assert(currentTimer >= 0);

        chrono::duration<double> elapsed = system_clock::now() - startTime[currentTimer];
        return (elapsed.count());
}

string Util::getDateTime()
{
        time_t time = system_clock::to_time_t(system_clock::now());
        return string(ctime(&time));
}

double Util::poissonPDF(unsigned int k, double mu)
{
        return exp(k*log(mu)-mu-lgamma(k+1));
}

double Util::poissonPDFratio(unsigned int k, double mu1, double mu2)
{
        return exp(k*log(mu1) - mu1 - k*log(mu2) + mu2);
}

double Util::negbinomialPDF(unsigned int k, double mu, double sigma2)
{
        // make sure that sigma2 is bigger than mu
        if ((sigma2 - mu) < 1e-3)
                sigma2 = mu + 1e-3;

        double p = (sigma2 - mu)/sigma2;
        double r = mu*mu/(sigma2 - mu);
        return exp(lgamma(k + r) - lgamma(k+1) - lgamma(r) + r*log(1-p) + k*log(p));
}

double Util::logNegbinomialPDF(unsigned int k, double mu, double sigma2)
{
        // make sure that sigma2 is bigger than mu
        if ((sigma2 - mu) < 1e-3)
                sigma2 = mu + 1e-3;

        double p = (sigma2 - mu)/sigma2;
        double r = mu*mu/(sigma2 - mu);
        return lgamma(k + r) - lgamma(k+1) - lgamma(r) + r*log(1-p) + k*log(p);
}

double Util::negbinomialPDFratio(unsigned int k, double mu1, double sigma21,
                                 double mu2, double sigma22)
{
        double p1 = (sigma21 - mu1)/sigma21;
        double r1 = mu1*mu1/(sigma21 - mu1);
        double p2 = (sigma22 - mu2)/sigma22;
        double r2 = mu2*mu2/(sigma22 - mu2);
        return exp(lgamma(k + r1) - lgamma(r1) + r1*log(1-p1) + k*log(p1) -
                   lgamma(k + r2) + lgamma(r2) - r2*log(1-p2) - k*log(p2));
}

double Util::geometricPDF(unsigned int k, double mu)
{
        double p = 1.0 / mu;
        return p * pow(1 - p, k-1);
}

double Util::logGeometricPDF(unsigned int k, double mu)
{
        double p = 1.0 / mu;
        return log(p) + (k-1)*log(1.0-p);
}

double Util::geometricnegbinomialPDFratio(unsigned int k, double mu1,
                                          double mu2, double var2)
{
        double p1 = 1.0 / mu1;
        double p2 = (var2 - mu2)/var2;
        double r2 = mu2*mu2/(var2 - mu2);
        return exp(log(p1) + (k-1)*log(1.0-p1) - lgamma(k + r2) + lgamma(k+1) +
                   lgamma(r2) - r2*log(1.0-p2) - k*log(p2));
}

void Util::binomialMixtureEM(const map<unsigned int, double>& data,
                             vector<double>& mu, vector<double> &var,
                             vector<double>& MC, int maxIterations)
{
        // sanity check
        assert(var.size() == mu.size());

        // shortcuts
        int numComponents = mu.size();

        // initialize weights
        map<unsigned int, double> *weight = new map<unsigned int, double>[numComponents];

        for (int itCount = 0; itCount < maxIterations; itCount++) {

                // compute weights
                for (const auto& element : data) {
                        unsigned int x = element.first;

                        // compute the weights corresponding to the negative binomials
                        for (int j = 0; j < numComponents; j++) {
                                double nom = 1.0;
                                for (int k = 0; k < numComponents; k++) {
                                        if (k == j)
                                                continue;
                                        nom += MC[k]/MC[j] * negbinomialPDFratio(x, mu[k], var[k], mu[j], var[j]);
                                }
                                weight[j][x] = 1.0 / nom;
                        }
                }

                // compute mean
                for (int j = 0; j < numComponents; j++) {
                        double sum = 0.0, count = 0.0;
                        for (const auto& element : data) {
                                unsigned int x = element.first;
                                double y = element.second;

                                count += weight[j][x] * y;
                                sum += weight[j][x] * x * y;
                        }

                        mu[j] = sum / count;
                        MC[j] = count;
                }

                // compute variance
                for (int j = 0; j < numComponents; j++) {
                        double sum = 0.0, count = 0.0;
                        for (const auto& element : data) {
                                unsigned int x = element.first;
                                double y = element.second;

                                sum += weight[j][x] * y * (x - mu[j]) * (x - mu[j]);
                                count += weight[j][x] * y;
                        }

                        var[j] = sum / count;

                        // make sure the variance is at least the average
                        // (overdispersed Poisson model)
                        if ((var[j] - mu[j]) < FLOAT_SMALL)
                                var[j] = mu[j] + FLOAT_SMALL;
                }

                // display results
                /*for (int j = 0; j < numComponents; j++) {
                        cout << "Average of component " << j << ": " << mu[j] << endl;
                        cout << "Mixing coefficient of component " << j << ": " << MC[j] << endl;
                        cout << "Variance of component " << j << ": " << var[j] << endl;
                }*/
        }

        delete [] weight;
}
