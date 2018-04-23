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

#ifndef KMERTABLE_H
#define KMERTABLE_H

#include "global.h"

#include <google/sparse_hash_set>
#include <mutex>
#include <condition_variable>

// ============================================================================
// TYPEDEFS
// ============================================================================

typedef google::sparse_hash_set<RKmer, RKmerHash> RKmerHashTable;

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class Settings;
class ReadFile;
class LibraryContainer;


// ============================================================================
// MIXING FUNCTION
// ============================================================================

class MixingLSB {

private:
        std::vector<KmerLSB> mixingF;   // mixing (permutation) table
        std::vector<KmerLSB> imixingF;  // inverse mixing table
public:
        /**
         * Default constructor
         */
        MixingLSB();

        /**
         * Mixing function for a kmer least significant bytes
         * @param input Unmixed kmer LSB input
         * @return Mixed kmer LSB
         */
        KmerLSB mix(KmerLSB input) const;

        /**
         * Inverse mixing function for a kmer least significant bytes
         * @param input Unmixed kmer LSB input
         * @return Mixed kmer LSB
         */
        KmerLSB invmix(KmerLSB input) const;
};

// ============================================================================
// KMER TABLE
// ============================================================================

class KmerTable {

private:
        const Settings& settings;               // reference to the settings object
        RKmerHashTable **tableThread;           // kmer hash table per thread
        RKmerHashTable **tables;                // kmer hash table
        MixingLSB mixFunction;                  // kmer lsb mixing function

        std::vector<std::vector<Kmer> > sharedKmerBuf;  // shared kmer buffer
        std::vector<std::mutex> sharedKmerBufMutex;     // shared kmer buffer mutex

        size_t numThreadReady;                  // number of worker threads ready
        std::mutex terminateMutex;              // kmer buffer mutex
        std::condition_variable terminateCV;    // read buffer full condition

        /**
         * Get the identifier of the thread that needs to process a kmer
         * @param kmer kmer to handle
         * @return [0 ... numThreads-1]
         */
        size_t getThreadIDForKmer(const Kmer& kmer) const;

        /**
         * Parse one read and generate the kmers
         * @param read Input read to process
         * @param kmerBuffer Output kmer buffers
         * @return True upon success, false otherwise
         */
        size_t parseRead(std::string &read,
                         std::vector<Kmer> *kmerBuffer);

        /**
         * Parse a buffer of reads and store kmers in temporary buffers per thread
         * @param readBuffer Input read buffer
         * @param kmerBuffer Output kmer buffers
         * @param myKmerBuf Vector to store local kmers
         */
        void parseReads(size_t thisThread,
                        std::vector<std::string>& readBuffer,
                        std::vector<Kmer>* kmerBuffer,
                        std::vector<Kmer>& myKmerBuf);

        /**
         * Actually store kmers in the tables
         * @param thisThread Identifier for this thread
         * @param kmerBuffer Kmers to store
         */
        void storeKmersInTable(size_t thisThread,
                               const std::vector<Kmer>& kmerBuffer);

        /**
         * Entry routine for worker thread
         * @param myID Unique threadID
         * @param input Pointer to the library container
         */
        void workerThread(size_t myID, LibraryContainer* inputs);

public:
        /**
         * Default constructor
         * @param settings Settings object
         */
        KmerTable(const Settings& settings) : settings(settings),
                tableThread(NULL), tables(NULL) {}

        /**
         * Destructor
         */
        ~KmerTable();

        /**
         * Read the input files specified in the command line
         * @param inputs Input libraries
         */
        void parseInputFiles(LibraryContainer &inputs);

        /**
         * Clear the table
         */
        void clear();

        /**
         * Find a kmer in the table
         * @param kmer Kmer to look for
         * @return pair< bool, bool >(found, flag)
         */
        std::pair<bool, bool> find(const Kmer &kmer) const;

        /**
         * Get the total number of kmers
         * @return The total number of kmers
         */
        size_t getNumKmers() const;

        /**
         * Get the total number of kmers with a coverage bigger than 1
         * @return The total number of kmers
         */
        size_t getNumKmersCovGTOne() const;

        /**
         * Write all kmers
         */
        void writeAllKmers(const std::string& filename);

        /**
         * Write the kmers with a coverage greater than one to disc
         */
        void writeKmersWithCovGTOne(const std::string& filename);

#ifdef DEBUG
        /**
         * Validate the first stage
         */
        void validateStage1();
#endif
};

#endif
