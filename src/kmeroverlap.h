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

#ifndef KMEROVERLAP_H
#define KMEROVERLAP_H

#include "tkmer.h"

#include <google/sparse_hash_map>
#include <deque>
#include <atomic>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class KmerOverlap;

// ============================================================================
// TYPEDEFS
// ============================================================================

// shortcut notation for a const iterator
typedef google::sparse_hash_map<Kmer, KmerOverlap, KmerHash>::const_iterator KmerOverlapIt;

// shortcut notation for a <Key, Data> pair
typedef std::pair<Kmer, KmerOverlap> KmerOverlapPair;

// ============================================================================
// KMER OVERLAP
// ============================================================================

class KmerOverlap {

private:
        const static unsigned char leftMask[4];
        const static unsigned char rightMask[4];
        const static unsigned char cLeftMask;
        const static unsigned char cRightMask;

public:
        std::atomic<unsigned char> bf;

        /**
         * Default constructor
         */
        KmerOverlap() : bf(0) {}

        /**
         * Default constructor
         */
        KmerOverlap(const KmerOverlap& rhs) : bf(rhs.bf.load()) {}

        /**
         * Create a kmer metadata from a given bitfield
         * @param metaData Given bitfield
         */
        KmerOverlap(uint8_t metaData) : bf(metaData) {}

        /**
         * Marks the left extension of a kmer with a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         */
        void markLeftOverlap(char nucleotide) {
                bf.fetch_or(leftMask[Nucleotide::charToNucleotide(nucleotide)]);
        }

        /**
         * Marks the right extension of a kmer with a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         */
        void markRightOverlap(char nucleotide) {
                bf.fetch_or(rightMask[Nucleotide::charToNucleotide(nucleotide)]);
        }

        /**
         * Unmarks the left extension of a kmer with a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         */
        void unmarkLeftOverlap(char nucleotide) {
                bf.fetch_and(~leftMask[Nucleotide::charToNucleotide(nucleotide)]);
        }

        /**
         * Unmarks the right extension of a kmer with a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         */
        void unmarkRightOverlap(char nucleotide) {
                bf.fetch_and(~rightMask[Nucleotide::charToNucleotide(nucleotide)]);
        }

        /**
         * Checks the left extension of a kmer with a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         * @return True of false
         */
        bool hasLeftOverlap(char nucleotide) const {
                return (bf & leftMask[Nucleotide::charToNucleotide(nucleotide)]) != 0;
        }

        /**
         * Checks the right extension of a kmer with a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         * @return True of false
         */
        bool hasRightOverlap(char nucleotide) const {
                return (bf & rightMask[Nucleotide::charToNucleotide(nucleotide)]) != 0;
        }

        /**
         * Check for a unique left extension of a kmer
         * @param nucleotide Unique left nucleotide (output)
         * @return True if a unique exists, false otherwise
         */
        bool hasUniqueLeftOverlap(char &nucleotide) const {
                unsigned int count = 0;
                for (unsigned int i = 0; i < 4; i++) {
                        if (bf & leftMask[i]) {
                                nucleotide = Nucleotide::nucleotideToChar(i);
                                count++;
                        }
                }

                return (count == 1);
        }

        /**
         * Check for a unique right extension of a kmer
         * @param nucleotide Unique right nucleotide (output)
         * @return True if a unique exists, false otherwise
         */
        bool hasUniqueRightOverlap(char &nucleotide) const {
                unsigned int count = 0;
                for (unsigned int i = 0; i < 4; i++) {
                        if (bf & rightMask[i]) {
                                nucleotide = Nucleotide::nucleotideToChar(i);
                                count++;
                        }
                }

                return (count == 1);
        }

        /**
         * Get the number of left-overlapping kmers
         * @return A number [0..4]
         */
        unsigned int getNumLeftOverlap() const {
                unsigned int count = 0;
                for (unsigned int i = 0; i < 4; i++)
                        if (bf & leftMask[i])
                                count++;
                return count;
        }

        /**
         * Get the number of right-overlapping kmers
         * @return A number [0..4]
         */
        unsigned int getNumRightOverlap() const {
                unsigned int count = 0;
                for (unsigned int i = 0; i < 4; i++)
                        if (bf & rightMask[i])
                                count++;
                return count;
        }

        /**
         * Get the left overlap encoded as bits
         * @return The left overlap
         */
        uint8_t getLeftOverlap() const {
                return bf >> 4;
        }

        /**
         * Get the right overlap encoded as bits
         * @return The right overlap
         */
        uint8_t getRightOverlap() const {
                return bf & cRightMask;
        }

        /**
         * Check if a kmer a left dead-end
         * @return True or false
         */
        bool isLeftDeadEnd() const {
                return !(bf & cLeftMask);
        }

        /**
         * Check if a kmer a left dead-end
         * @return True or false
         */
        bool isRightDeadEnd() const {
                return !(bf & cRightMask);
        }

        /**
         * Get reverse complement of the kmerMD
         * @return A kmerMD containing the reverse complement info
         */
        KmerOverlap getReverseComplement() const {
                KmerOverlap copy(bf);
                copy.bf = (((copy.bf & 0xaa) >> 1) | ((copy.bf & 0x55) << 1));
                copy.bf = (((copy.bf & 0xcc) >> 2) | ((copy.bf & 0x33) << 2));
                copy.bf = (((copy.bf & 0xf0) >> 4) | ((copy.bf & 0x0f) << 4));
                return copy;
        }

        /**
         * Operator ==
         * @param rhs Right hand side
         * @return True of false
         */
        bool operator==(const KmerOverlap& rhs) const {
                return bf == rhs.bf;
        }

        /**
         * Operator !=
         * @param rhs Righ hand size
         * @return True of false
         */
        bool operator!=(const KmerOverlap& rhs) const {
                return !operator==(rhs);
        }

        /**
         * Write a kmeroverlap to file
         * @param ofs Open output file stream
         */
        void write(std::ofstream& ofs) const {
                unsigned char v = bf;
                ofs.write((char*)&v, sizeof(v));
        }

        /**
         * Read a kmeroverlap from file
         * @param ifs Open input file stream
         */
        void read(std::ifstream& ifs) {
                unsigned char v;
                ifs.read((char*)&v, sizeof(v));
                bf = v;
        }
};

// ============================================================================
// KMER OVERLAP REFERENCE
// ============================================================================

// a kmer overlap reference is a pair of two values: a) the iterator that points
// to that kmer or its reverse complement in the table.  b) a boolean to
// indicate whether the iterator points the reverse complement kmer or not

class KmerOverlapRef : public std::pair<KmerOverlapIt, bool> {

public:
        /**
         * Default constructor
         */
        KmerOverlapRef() {}

        /**
         * Constructor
         * @param it Iterator to the table
         * @param reverse True if the iterator points to the reverse complement
         */
        KmerOverlapRef(KmerOverlapIt it, bool reverse) :
                std::pair<KmerOverlapIt, bool>(it, reverse) {}

        /**
         * Returns the nucleotide at the left side of the kmer
         * @return ASCII encoding of 'A', 'C', 'G' and 'T'
         */
        char peekNucleotideLeft() const {
                return (second) ?
                        Nucleotide::getComplement(first->first.peekNucleotideRight()) :
                        first->first.peekNucleotideLeft();
        }

        /**
         * Returns the nucleotide at the right side of the kmer
         * @return ASCII encoding of 'A', 'C', 'G' and 'T'
         */
        char peekNucleotideRight() const {
                return (second) ?
                        Nucleotide::getComplement(first->first.peekNucleotideLeft()) :
                        first->first.peekNucleotideRight();
        }

        /**
         * Mark the left extention of the kmer by a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         */
        void markLeftOverlap(char nucleotide) const {
                // shortcut notation
                KmerOverlap &metaData = const_cast<KmerOverlap&>(first->second);

                if (second)
                        metaData.markRightOverlap(Nucleotide::getComplement(nucleotide));
                else
                        metaData.markLeftOverlap(nucleotide);
        }

        /**
         * Mark the right extention of the kmer by a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         */
        void markRightOverlap(char nucleotide) const {
                // shortcut notation
                KmerOverlap &metaData = const_cast<KmerOverlap&>(first->second);

                if (second)
                        metaData.markLeftOverlap(Nucleotide::getComplement(nucleotide));
                else
                        metaData.markRightOverlap(nucleotide);
        }

        /**
         * Unmark the left extention of the kmer by a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         */
        void unmarkLeftOverlap(char nucleotide) const {
                // shortcut notation
                KmerOverlap &metaData = const_cast<KmerOverlap&>(first->second);

                if (second)
                        metaData.unmarkRightOverlap(Nucleotide::getComplement(nucleotide));
                else
                        metaData.unmarkLeftOverlap(nucleotide);
        }

        /**
         * Unmark the right extention of the kmer by a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         */
        void unmarkRightOverlap(char nucleotide) const {
                // shortcut notation
                KmerOverlap &metaData = const_cast<KmerOverlap&>(first->second);

                if (second)
                        metaData.unmarkLeftOverlap(Nucleotide::getComplement(nucleotide));
                else
                        metaData.unmarkRightOverlap(nucleotide);
        }

        /**
         * Check if a kmer can be left-extended by a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         * @return true of false
         */
        bool hasLeftOverlap(char nucleotide) const {
                if (second)
                        return first->second.hasRightOverlap(Nucleotide::getComplement(nucleotide));
                return first->second.hasLeftOverlap(nucleotide);
        }

        /**
         * Check if a kmer can be right-extended by a certain nucleotide
         * @param nucleotide Nucleotide under consideration
         * @return true of false
         */
        bool hasRightOverlap(char nucleotide) const {
                if (second)
                        return first->second.hasLeftOverlap(Nucleotide::getComplement(nucleotide));
                return first->second.hasRightOverlap(nucleotide);
        }

        /**
         * Get the number of left-overlapping kmers
         * @return A number [0..4]
         */
        unsigned int getNumLeftOverlap() const {
                if (second)
                        return first->second.getNumRightOverlap();
                else
                        return first->second.getNumLeftOverlap();
        }

        /**
         * Get the number of right-overlapping kmers
         * @return A number [0..4]
         */
        unsigned int getNumRightOverlap() const {
                if (second)
                        return first->second.getNumLeftOverlap();
                else
                        return first->second.getNumRightOverlap();
        }

        /**
         * Get the left overlap encoded as bits
         * @return The left overlap
         */
        uint8_t getLeftOverlap() const {
                if (second)
                        return first->second.getReverseComplement().getLeftOverlap();
                else
                        return first->second.getLeftOverlap();
        }

        /**
         * Get the right overlap encoded as bits
         * @return The right overlap
         */
        uint8_t getRightOverlap() const {
                if (second)
                        return first->second.getReverseComplement().getRightOverlap();
                else
                        return first->second.getRightOverlap();
        }

        /**
         * Check if a kmer is extendable to the left in a unique way
         * @param nucleotide Nucleotide that extends the kmer to the left
         * @return True of false
         */
        bool hasLeftUniqueOverlap(char &nucleotide) const {
                if (second) {
                        bool value = first->second.hasUniqueRightOverlap(nucleotide);
                        if (value)
                                nucleotide = Nucleotide::getComplement(nucleotide);
                        return value;
                }

                return first->second.hasUniqueLeftOverlap(nucleotide);
        }

        /**
         * Check if a kmer is extendable to the right in a unique way
         * @param nucleotide Nucleotide that extends the kmer to the right
         * @return True of false
         */
        bool hasRightUniqueOverlap(char &nucleotide) const {
                if (second) {
                        bool value = first->second.hasUniqueLeftOverlap(nucleotide);
                        if (value)
                                nucleotide = Nucleotide::getComplement(nucleotide);
                        return value;
                }

                return first->second.hasUniqueRightOverlap(nucleotide);
        }

        /**
         * Check if a kmer a left dead-end
         * @return True or false
         */
        bool isLeftDeadEnd() const {
                if (second)
                        return first->second.isRightDeadEnd();
                return first->second.isLeftDeadEnd();
        }

        /**
         * Check if a kmer a right dead-end
         * @return True or false
         */
        bool isRightDeadEnd() const {
                if (second)
                        return first->second.isLeftDeadEnd();
                return first->second.isRightDeadEnd();
        }

        /**
         * Get the kmer from the reference
         * @return The kmer
         */
        const Kmer getKmer() const {
                if (second)
                        return first->first.getReverseComplement();
                return first->first;
        }

        /**
         * Set a kmer as processed
         * @param Value True or false
         */
        void setProcessed(bool value) {
                const_cast<Kmer&>(first->first).setFlag1(value);
        }

        /**
         * Check if a kmer is processed
         * @return True or false
         */
        bool isProcessed() const {
                return first->first.getFlag1();
        }

        /**
         * Set a kmer as unique
         * @param Value True or false
         */
        void setUnique(bool value) {
                const_cast<Kmer&>(first->first).setFlag2(!value);
        }

        /**
         * Check if a kmer is unique
         * @return True or false
         */
        bool isUnique() const {
                return !first->first.getFlag2();
        }
};

#endif
