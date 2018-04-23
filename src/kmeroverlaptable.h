/***************************************************************************
 *   Copyright (C) 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)     *
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

#ifndef KMEROVERLAPTABLE_H
#define KMEROVERLAPTABLE_H

#include "kmeroverlap.h"

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class LibraryContainer;
class Settings;

// ============================================================================
// KMER OVERLAP TABLE
// ============================================================================

class KmerOverlapTable {

private:
        const Settings &settings;       // reference to the settings object
        google::sparse_hash_map<Kmer, KmerOverlap, KmerHash> table; // actual table

        /**
         * Get the unique kmer extending a given kmer to the left
         * @param kmer Kmer to be extended (input)
         * @param leftKmer Kmer that left-overlaps with kmer (output)
         * @return True if a unique left-overlapping kmer is found
         */
        bool getLeftUniqueKmer(const KmerOverlapRef &kmer,
                               KmerOverlapRef &leftKmer) const;

        /**
         * Get the unique kmer extending a given kmer to the right
         * @param kmer Kmer to be extended (input)
         * @param rightKmer Kmer that right-overlaps with kmer (output)
         * @return True if a unique right-overlapping kmer is found
         */
        bool getRightUniqueKmer(const KmerOverlapRef &kmer,
                                KmerOverlapRef &rightKmer) const;

        /**
         * Find a kmer in the table
         * @param kmer Kmer to look for
         * @return KmerRef containing iterator to the kmer and reversed flag
         */
        KmerOverlapRef find(const Kmer &kmer) const;

        /**
         * Insert a kmer in the table
         * @param kmer Kmer to insert
         * @return KmerRef containing iterator to the kmer and reversed flag
         */
        KmerOverlapRef insert(const Kmer &kmer);

        /**
         * Convert a deque of overlapping kmers to a string
         * @param kmerSeq A deque of overlapping kmers
         * @param output An stl string (output)
         */
        static void convertKmersToString(const std::deque<KmerOverlapRef> &kmerSeq,
                                         std::string &output);

        /**
         * Parse one read and generate the kmers
         * @param read Input read to process
         * @param kmerBuffer Output kmer buffers to be inserted
         * @return True upon success, false otherwise
         */
        void parseRead(std::string &read,
                       std::vector<std::pair<Kmer, KmerOverlap> >& kmerBuffer) ;

        /**
         * Parse a buffer of reads and store kmers in temporary buffers per thread
         * @param readBuffer Input read buffer
         * @param kmerBuffer Output kmer buffers
         */
        void parseReads(size_t thisThread,
                        std::vector<std::string>& readBuffer,
                        std::vector<std::pair<Kmer, KmerOverlap> >& kmerBuffer) ;

        /**
         * Entry routine for worker thread
         * @param myID Unique threadID
         */
        void workerThread(size_t myID, LibraryContainer* inputs);
public:
        /**
         * Default constructor
         * @param settings Settings object
         */
        KmerOverlapTable(const Settings& settings) : settings(settings) {}

        /**
         * Get the number of elements in the table
         * @return The number of elements
         */
        size_t size() const {
                return table.size();
        }

        /**
         * Clear the kmer table
         */
        void clear() {
                table.clear();
        }

        /**
         * Load the kmers from disc
         */
        void loadKmersFromDisc(const std::string& filename);

        /**
         * Find overlap between the kmers stored in the overlap table
         * @param inputs Input libraries
         */
        void parseInputFiles(LibraryContainer &inputs);

        /**
         * Extract nodes from graph
         * @param nodeFilename Filename for the nodes
         * @param arcFilename Filename for the arcs
         * @param metaDataFilename Filename for the metadata
         */
        void extractNodes(const std::string& nodeFilename,
                          const std::string& arcFilename,
                          const std::string& metaDataFilename);

#ifdef DEBUG
        /**
         * Validate the kmer table
         */
        void validateStage2();
#endif
};

#endif
