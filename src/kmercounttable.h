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

#ifndef KMER_SPECTRUM_H
#define KMER_SPECTRUM_H

#include <atomic>
#include <limits>
#include <vector>
#include <map>

#include <google/sparse_hash_map>

// ============================================================================
// KMER SPECTRUM CLASS
// ============================================================================

class KmerSpectrum {

private:
        std::map<unsigned int, size_t> spectrum;
        size_t xMax, yMax;

        // k-mer spectrum components
        size_t numErrComp, numTrueComp;
        std::vector<double> spectrumMu, spectrumVar, spectrumMC;

        double covCutoff;

public:
        /**
         * Default constructor
         */
        KmerSpectrum() : xMax(0), yMax(0), numErrComp(2),
                numTrueComp(2), covCutoff(5) {}

        /**
         * Operator [] overloading
         * @param key key
         * @return reference to the value
         */
        size_t& operator[](unsigned int key) {
                return spectrum[key];
        }

        /**
         * Fit a mixture model of negative binomials to the spectrum
         * @param avgKmerCov Estimate of the average true kmer coverage
         */
        void fitKmerSpectrum(double avgKmerCov);

        /**
         * Get the coverage cutoff
         * @return Coverage cutoff
         */
        double getCovCutoff() const {
                return covCutoff;
        }

        /**
         * Get the average k-mer coverage
         * @return The average k-mer coverage
         */
        double getAvgKmerCov() const {
                return spectrumMu[numErrComp];
        }

        /**
         * Evaluate the fit in point x given multiplicity mult
         * @param x Point in which to evaluate the spectrum
         * @param mult multiplicity
         * @return Spectrum evaluated in x given multiplicity mult
         */
        double evalSpec(unsigned int x, unsigned int mult) const;

        /**
         * Evaluate the logarithm of the fit in point x given multiplicity mult
         * @param x Point in which to evaluate the spectrum
         * @param mult multiplicity
         * @return Spectrum evaluated in x given multiplicity mult
         */
        double evalSpecLog(unsigned int x, unsigned int mult) const;

        /**
         * Get the odds ratio of observing kmerCov with mult1 over mult2
         * @param kmerCov Kmer coverage
         * @param mult1 Multiplicity 1 of the coverage (0 = error model)
         * @param mult2 Multiplicity 2 of the coverage (0 = error model)
         * @return The odds ratio
         */
        double getOddsRatio(unsigned int kmerCov, unsigned int mult1,
                            unsigned int mult2) const {
                return evalSpec(kmerCov, mult1) / evalSpec(kmerCov, mult2);
        }

        /**
         * Get the estimated genome size
         * @return The estimated genome size
         */
        size_t getEstimatedGenomeSize() const;

        /**
         * Write the spectrum to disk
         * @param filename File name
         */
        void writeSpectrum(const std::string& filename) const;

        /**
         * Load the spectrum from disk
         * @param filename File name
         */
        void loadSpectrum(const std::string& filename);

        /**
         * Write the GNUPlot auxiliary file
         * @param filename File name
         */
        void writeGNUPlotFile(const std::string& filename) const;

        /**
         * Write the spectrum fit file
         * @param filename File name
         */
        void writeSpectrumFit(const std::string& filename) const;

        /**
         * Load the spectrum from disk
         * @param filename File name
         */
        void loadSpectrumFit(const std::string& filename);

        /**
         * Operator<< overloading
         * @param out Output stream (input)
         * @param kms Kmer spectrum
         * @return Output stream (output)
         */
        friend std::ostream &operator<<(std::ostream &out, const KmerSpectrum& kms);
};

// ============================================================================
// KMER COUNT CLASS
// ============================================================================

class KmerCount {

private:
        std::atomic<uint16_t> kmerCount;     // times a kmer was observed
        NodeID nodeID;          // node identifier to which the kmer belongs

public:
        /**
         * Default constructor
         */
        KmerCount(NodeID nodeID_) : kmerCount(0), nodeID(nodeID_) {}

        /**
         * Copy constructor
         * @param rhs Right hand side
         */
        KmerCount(const KmerCount& rhs) {
                kmerCount = rhs.kmerCount.load();
                nodeID = rhs.nodeID;
        }

        /**
         * Assignment operator
         * @param rhs Right hand side
         */
        KmerCount& operator=(const KmerCount& rhs) {
                kmerCount = rhs.kmerCount.load();
                nodeID = rhs.nodeID;
                return *this;
        }

        /**
         * Increment the kmer count
         */
        void incrementCount() {
                uint16_t curVal, newVal;
                do {
                        curVal = kmerCount;
                        if (curVal == std::numeric_limits<uint16_t>::max())
                                return;
                        newVal = curVal + 1;
                } while (!std::atomic_compare_exchange_weak(&kmerCount, &curVal, newVal));
        }

        /**
         * The the kmer count
         * @return The kmer count
         */
        size_t getCount() const {
                return kmerCount;
        }

        /**
         * Reverse complement the KmerCount
         */
        void revComp() {
                nodeID = -nodeID;
        }

        /**
         * Get the node identifier
         * @return The node identifier
         */
        NodeID getNodeID() const {
                return nodeID;
        }
};

// ============================================================================
// KMER COUNT TABLE CLASS
// ============================================================================

class KmerCountTable {

private:
        std::atomic<size_t> numUniqueKmers;     // number of unique k-mers
        // actual container
        google::sparse_hash_map<Kmer, KmerCount, KmerHash> kmerCountTable;

public:
        /**
         * Default constructor
         */
        KmerCountTable() : numUniqueKmers(0) {}

        /**
         * Insert a < Kmer, KmerCount > element in the table
         * @param kmer Kmer key
         * @param kmerCount KmerCount value
         * @return True if the element was inserted, false otherwise
         */
        bool insert(const Kmer& kmer, KmerCount kmerCount);

        /**
         * Increment the count for a certain kmer
         * @param kmer Kmer key
         * @return The nodeID to which the kmer belongs, zero if kmer is not found
         */
        NodeID incrKmerCount(const Kmer& kmer);

        /**
         * Clear table
         */
        void clear() {
                numUniqueKmers = 0;
                kmerCountTable.clear();
        }

        /**
         * Get the size of the table
         * @return The size of the table
         */
        size_t size() const {
                return kmerCountTable.size();
        }

        /**
         * Pre-allocate memory to hold about hint elements
         * @param hint Number of elements
         */
        void resize(size_t hint) {
                kmerCountTable.resize(hint);
        }

        /**
         * Get the kmer spectrum from the KmerCountTable
         * @param spectrum Spectrum to generate (output)
         */
        void getKmerSpectrum(KmerSpectrum& spectrum);
};

#endif
