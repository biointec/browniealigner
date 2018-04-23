/***************************************************************************
 *   Copyright (C) 2015 - 2016 Jan Fostier (jan.fostier@intec.ugent.be)    *
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

#include "global.h"
#include "tkmer.h"
#include "kmercounttable.h"
#include "util.h"

using namespace std;

// ============================================================================
// KMER SPECTRUM CLASS
// ============================================================================

ostream &operator<<(std::ostream &out, const KmerSpectrum& kms)
{
        out << "Kmer spectrum with " << kms.numErrComp << " error and "
            << kms.numTrueComp << " true components:" << endl;
        for (size_t j = 0; j < kms.numErrComp; j++)
                out << "\tError component " << j+1 << ": [" << kms.spectrumMu[j]
                    << ", " << kms.spectrumVar[j] << ", " << kms.spectrumMC[j]
                    << "]" << endl;
        for (size_t j = kms.numErrComp; j < kms.numTrueComp + kms.numErrComp; j++)
                out << "\tTrue component " << j+1 - kms.numErrComp << ": ["
                    << kms.spectrumMu[j] << ", " << kms.spectrumVar[j] << ", "
                    << kms.spectrumMC[j] << "]" << endl;

        out << "Estimated genome size based on this spectrum: "
            << kms.getEstimatedGenomeSize() << endl;
        out << "Coverage cutoff based on this spectrum: "
            << kms.covCutoff << endl;

        return out;
}

double KmerSpectrum::evalSpec(unsigned int x, unsigned int mult) const
{
        if (mult > 0) {                 // we're evaluating the true spectrum
                double mu, var, MC;
                if (mult <= numTrueComp) {
                        int i = numErrComp+mult-1;
                        mu = spectrumMu[i];
                        var = spectrumVar[i];
                        MC = spectrumMC[i];
                } else {
                        mu = mult * spectrumMu[numErrComp];
                        var = mult * spectrumVar[numErrComp];
                        MC = spectrumMC[numErrComp];
                }

                return MC * Util::negbinomialPDF(x, mu, var);
        }

        double P = 0.0;
        for (size_t i = 0; i < numErrComp; i++)
                P += spectrumMC[i] * Util::negbinomialPDF(x, spectrumMu[i], spectrumVar[i]);
        return P;
}

double KmerSpectrum::evalSpecLog(unsigned int x, unsigned int mult) const
{
        if (mult > 0) {                 // we're evaluating the true spectrum
                double mu, var, MC;
                if (mult <= numTrueComp) {
                        int i = numErrComp+mult-1;
                        mu = spectrumMu[i];
                        var = spectrumVar[i];
                        MC = spectrumMC[i];
                } else {
                        mu = mult * spectrumMu[numErrComp];
                        var = mult * spectrumVar[numErrComp];
                        MC = spectrumMC[numErrComp];
                }

                return log(MC) + Util::logNegbinomialPDF(x, mu, var);
        }

        double P = 0.0;
        for (size_t i = 0; i < numErrComp; i++)
                P += log(spectrumMC[i]) + Util::logNegbinomialPDF(x, spectrumMu[i], spectrumVar[i]);
        return P;
}

void KmerSpectrum::fitKmerSpectrum(double avgKmerCov)
{
        numErrComp = 2;         // TODO: export as a parameter
        numTrueComp = 2;        // TODO: export as a parameter

        spectrumMu.clear(); spectrumVar.clear(); spectrumMC.clear();

        for (size_t i = 0; i < numErrComp; i++) {
                spectrumMu.push_back(2.0 + i);
                spectrumVar.push_back(1.1 * spectrumMu.back());
                spectrumMC.push_back(1.0);
        }

        for (size_t i = 0; i < numTrueComp; i++) {
                spectrumMu.push_back((i+1) * avgKmerCov);
                spectrumVar.push_back(1.1 * spectrumMu.back());
                spectrumMC.push_back(1.0);
        }

        xMax = (numTrueComp + 1) * avgKmerCov;

        map<unsigned int, double> data;
        for (const auto& it : spectrum) {
                if (it.first > xMax)
                        break;
                data[it.first] = it.second;
        }

        Util::binomialMixtureEM(data, spectrumMu, spectrumVar, spectrumMC, 50);

        yMax = 1.2 * spectrumMC[numErrComp] * Util::negbinomialPDF(spectrumMu[numErrComp], spectrumMu[numErrComp], spectrumVar[numErrComp]);

        // find the minimum in between the error model and the unique peak
        covCutoff = ceil(spectrumMu[0]);
        for ( ; covCutoff < spectrumMu[numErrComp]; covCutoff++)
                if (getOddsRatio(covCutoff, 0, 1) < 0.5)
                        break;
}

size_t KmerSpectrum::getEstimatedGenomeSize() const
{
        // estimate the genome size
        size_t genomeSize = 0;
        for (size_t mult = 1; mult <= numTrueComp; mult++)
                genomeSize += mult * spectrumMC[mult - 1 + numErrComp];
        return genomeSize;
}

void KmerSpectrum::writeSpectrum(const string& filename) const
{
        ofstream ofs(filename.c_str());
        for (const auto& it : spectrum) {
                double fit = 0.0;
                for (size_t i = 0; i < numErrComp + numTrueComp; i++)
                        fit += spectrumMC[i] * Util::negbinomialPDF(it.first, spectrumMu[i], spectrumVar[i]);
                ofs << it.first << "\t" << it.second << "\t" << fit << endl;
        }
        ofs.close();
}

void KmerSpectrum::loadSpectrum(const string& filename)
{
        spectrum.clear();

        ifstream ifs(filename.c_str());
        while (true) {
                double x, y, dummy;
                ifs >> x >> y >> dummy;
                if (!ifs.good())
                        break;
                spectrum[x] = y;
        }
        ifs.close();
}

void KmerSpectrum::writeGNUPlotFile(const string& filename) const
{
        ofstream ofs(filename.c_str());
        ofs << "set output \"spectrum.ps\"" << endl;
        ofs << "set terminal postscript landscape" << endl;
        ofs << "set xrange [0:" << xMax << "]" << endl;
        ofs << "set yrange [0:" << yMax << "]" << endl;
        ofs << "plot \"spectrum.txt\" using 1:2 title \'k-mer spectrum\' with lines, ";
        ofs << "\"spectrum.txt\" using 1:3 title \'model fit\' with lines lt 3 lw 3 lc 3" << endl;
        ofs.close();
}

void KmerSpectrum::writeSpectrumFit(const string& filename) const
{
        ofstream ofs(filename.c_str());
        ofs << xMax << "\t" << yMax << endl;
        ofs << numErrComp << "\t" << numTrueComp << endl;
        for (size_t i = 0; i < numErrComp + numTrueComp; i++)
                ofs << spectrumMu[i] << "\t" << spectrumVar[i] << "\t" << spectrumMC[i] << endl;
        ofs << covCutoff << endl;
        ofs.close();
}

void KmerSpectrum::loadSpectrumFit(const string& filename)
{
        spectrumMu.clear(); spectrumVar.clear(); spectrumMC.clear();
        ifstream ifs(filename.c_str());
        ifs >> xMax >> yMax;
        ifs >> numErrComp >> numTrueComp;
        for (size_t i = 0; i < numErrComp + numTrueComp; i++) {
                double mu, var, MC;
                ifs >> mu >> var >> MC;
                spectrumMu.push_back(mu);
                spectrumVar.push_back(var);
                spectrumMC.push_back(MC);
        }
        ifs >> covCutoff;
        ifs.close();
}

// ============================================================================
// KMER COUNT TABLE CLASS
// ============================================================================

bool KmerCountTable::insert(const Kmer& kmer, KmerCount kmerCount)
{
        // chose a representative kmer
        Kmer reprKmer = kmer.getRepresentative();

        // reverse-complement the KmerCount if necessary
        if (kmer != reprKmer)
                kmerCount.revComp();

        // insert KmerCount in table
        auto it = kmerCountTable.insert(pair<Kmer, KmerCount>(reprKmer, kmerCount));

        return it.second;
}

NodeID KmerCountTable::incrKmerCount(const Kmer& kmer)
{
        // get the representative kmer
        Kmer reprKmer = kmer.getRepresentative();

        // insert value in table
        auto it = kmerCountTable.find(reprKmer);

        if (it != kmerCountTable.end()) {
                it->second.incrementCount();
                return (kmer == reprKmer) ? it->second.getNodeID() : -it->second.getNodeID();
        } else {
                numUniqueKmers++;
                return 0;
        }
}

void KmerCountTable::getKmerSpectrum(KmerSpectrum& spectrum)
{
        for (auto it : kmerCountTable) {
                unsigned x = it.second.getCount();
                spectrum[x]++;
        }

        spectrum[1] = numUniqueKmers;
}
