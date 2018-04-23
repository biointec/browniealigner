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

#ifndef TKMER_H
#define TKMER_H

#include "global.h"
#include "nucleotide.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <google/sparse_hash_set>

//============================================================================
// CLASS PROTOTYPES
// ============================================================================

class TString;

//============================================================================
// KMER TEMPLATE CLASS
// ============================================================================

template<size_t numBytes>
class TKmer {

private:
        static size_t k;        // the value of k (always odd)
        static size_t kMSB;     // the most significant occupied byte (k / 4)
        static size_t kMSLL;    // the most significant occupied uint64_t (k / 32)

        static uint8_t leftBit;
        static uint8_t rightBit;
        static uint64_t metaMask;

        uint8_t buf[numBytes];

        /**
         * Clear kmer
         */
        void clear() {
                memset(buf, 0, numBytes);
        }

public:
        /**
         * Default constructor
         */
        TKmer() {
                clear();
        }

        /**
         * Create a kmer from a c-string
         * @param str Null-terminated ASCII string of size >= k
         */
        TKmer(const char *str);

        /**
         * Create a kmer from a c++ string
         * @param str c++ string of size >= k
         * @param offset Offset marking the starting position
         */
        TKmer(const std::string& str, size_t offset = 0);

        /**
         * Create a kmer from a tight string
         * @param tString Tight string
         * @param offset Offset marking the starting position
         */
        TKmer(const TString &tString, size_t offset = 0);

        /**
         * Create a kmer from an input file stream
         * @param ifs Opened input file stream
         */
        TKmer(std::ifstream &ifs);

        /**
         * Create a kmer from a reduced kmer and lsb
         * @param rKmer Reduced kmer
         * @param lsb Least significant bytes (input)
         */
        TKmer(const TKmer<numBytes - 2>& kmer, KmerLSB lsb);

        /**
         * Create a reduced kmer from a kmer
         * @param rKmer Reduced kmer
         * @param lsb Least significant bytes (output)
         */
        TKmer(const TKmer<numBytes + 2>& kmer, KmerLSB &lsb);

        /**
         * Replace the kmer by its reverse word
         */
        void reverse();

        /**
         * Get reverse of the kmer
         * @return A kmer containing the reverse of this kmer
         */
        TKmer getReverse() const {
                TKmer copy = *this;
                copy.reverse();
                return copy;
        }

        /**
         * Replace the kmer by its complement word
         */
        void complement();

        /**
         * Get complement of the kmer
         * @return A kmer containing the complement of this kmer
         */
        TKmer getComplement() const {
                TKmer copy = *this;
                copy.complement();
                return copy;
        }

        /**
         * Calculate reverse complement of the kmer
         */
        void reverseComplement();

        /**
         * Get reverse complement of the kmer
         * @return A kmer containing the reverse complement of this kmer
         */
        TKmer getReverseComplement() const {
                TKmer copy = *this;
                copy.reverseComplement();
                return copy;
        }

        /**
         * Get the represensative kmer (smallest of kmer and reverse complement)
         * @return Representative kmer
         */
        TKmer getRepresentative() const {
                TKmer kmerRC = getReverseComplement();
                return (kmerRC < *this) ? kmerRC : *this;
        }

        /**
         * Set flag 1 value
         * @param value True of false
         */
        void setFlag1(bool value) {
                if (value)
                        buf[kMSB] |= leftBit;
                else
                        buf[kMSB] &= ~leftBit;
        }

        /**
         * Set flag 2 value
         * @param value True of false
         */
        void setFlag2(bool value) {
                if (value)
                        buf[kMSB] |= rightBit;
                else
                        buf[kMSB] &= ~rightBit;
        }

        /**
         * Get flag 1 value
         * @return True of false
         */
        bool getFlag1() const {
                return (buf[kMSB] & leftBit) != 0;
        }

        /**
         * Get flag 2 value
         * @return True of false
         */
        bool getFlag2() const {
                return (buf[kMSB] & rightBit) != 0;
        }

        /**
         * Push a nucleotide on the right side of the kmer
         * @param c ASCII encoding of 'A', 'C', 'G' and 'T'
         */
        void pushNucleotideRight(char c);

        /**
         * Push a nucleotide on the left side of the kmer
         * @param c ASCII encoding of 'A', 'C', 'G' and 'T'
         */
        void pushNucleotideLeft(char c);

        /**
         * Returns the nucleotide at the left side of the kmer
         * @return ASCII encoding of 'A', 'C', 'G' and 'T'
         */
        char peekNucleotideLeft() const {
                return Nucleotide::nucleotideToChar(buf[0]);
        }

        /**
         * Returns the nucleotide at the right side of the kmer
         * @return ASCII encoding of 'A', 'C', 'G' and 'T'
         */
        char peekNucleotideRight() const {
                int byteID = (k-1) / 4;
                int bitID = 2*((k-1) % 4);
                return Nucleotide::nucleotideToChar(buf[byteID] >> bitID);
        }

        /**
         * Operator '==' overloading
         * @param rhs Right hand side kmer
         * @return True if they're equal (i.e. string is equal)
         */
        bool operator==(const TKmer<numBytes> &rhs) const;

        /**
         * Operator '!=' overloading
         * @param rhs Right hand side kmer
         * @return True if they differ (i.e. string differs)
         */
        bool operator!=(const TKmer<numBytes> &rhs) const {
                return !(*this == rhs);
        }

        /**
         * Operator '<' overloading
         * @param rhs Right hand side kmer
         * @return True or false
         */
        bool operator<(const TKmer<numBytes> &rhs) const;

        /**
         * Operator '>' overloading
         * @param rhs Right hand side kmer
         * @return True or false
         */
        bool operator>(const TKmer<numBytes> &rhs) const;

        /**
         * Get a hash value for the kmer (code from Thomas Wang)
         * @return A hash value for the kmer
         */
        size_t getHash() const;

        /**
         * Convert kmer to a string
         * @return stl string
         */
        std::string str() const;

        /**
         * Specify kmer word size
         * @param wordSize Word size
         */
        static void setWordSize(size_t wordSize) {
                assert(wordSize % 2 == 1);
                assert(wordSize <= 4*numBytes-1);

                k = wordSize;
                kMSB = k / 4;
                kMSLL = k / 32;
                leftBit = 1 << (2*(k % 4) + 1);
                rightBit = 1 << 2*(k % 4);
                metaMask = uint64_t(3) << 2*(k % 32);
        }

        /**
         * Get the kmer word size
         */
        static size_t getK() {
                return k;
        }

        /**
         * Write a kmer to file
         * @param ofs Openen output file stream
         */
        void write(std::ofstream& ofs) const {
                ofs.write((char*)buf, kMSB+1);
        }

        /**
         * Write a kmer to file
         * @param ofs Openen output file stream
         */
        void writeBytes() const {
                for (int i = 0; i < numBytes; i++)
                        std::cout << (int)buf[i] << " ";
                std::cout << std::endl;
        }

        /**
         * Write a kmer to file
         * @param ofs Open output file stream
         */
        void writeNoFlags(std::ofstream& ofs) const {
                const size_t llSize = (numBytes + 7) / 8;
                uint64_t work[llSize];
                memcpy(work, buf, numBytes);
                work[kMSLL] &= ~metaMask;
                ofs.write((char*)work, kMSB+1);
        }

        /**
         * Operator<< overloading
         * @param out Output stream to add kmer to
         * @param kmer Right-hand side kmer
         * @return Output stream with the kmer added to it
         */
        friend std::ostream &operator<<(std::ostream &out, const TKmer<numBytes> &kmer) {
                for (size_t byteID = 0, offset = 0, i = 0; i < TKmer<numBytes>::k; i++) {
                        out << Nucleotide::nucleotideToChar(kmer.buf[byteID] >> offset);
                        offset += 2;
                        if (offset == 8) {
                                byteID++;
                                offset = 0;
                        }
                }

                return out;
        }

        friend class TKmer<numBytes - 2>;
        friend class TKmer<numBytes + 2>;

} ATTRIBUTE_PACKED;

// ============================================================================
// HASH FUNCTION
// ============================================================================

template<size_t numBytes>
struct TKmerHash {
                size_t operator()(const TKmer<numBytes> &kmer) const {
                        return kmer.getHash();
        }
};

// ============================================================================
// TKMER CLASS
// ============================================================================

template<size_t numBytes> size_t TKmer<numBytes>::k = 4*numBytes - 1;
template<size_t numBytes> size_t TKmer<numBytes>::kMSB = k / 4;
template<size_t numBytes> size_t TKmer<numBytes>::kMSLL = k / 32;

template<size_t numBytes> uint8_t TKmer<numBytes>::leftBit = 1 << (2*(k % 4) + 1);
template<size_t numBytes> uint8_t TKmer<numBytes>::rightBit = 1 << (2*(k % 4));
template<size_t numBytes> uint64_t TKmer<numBytes>::metaMask = uint64_t(3) << (2*(k % 32));

template<size_t numBytes>
TKmer<numBytes>::TKmer(const char* str)
{
        assert(strlen(str) >= k);

        const size_t llSize = (numBytes + 7) / 8;
        uint64_t work[llSize];

        for (size_t i = 0; i < kMSLL; i++, str += 32)
                work[i] = Nucleotide::pack32(str);

        work[kMSLL] = Nucleotide::pack32(str, k % 32);

        for (size_t i = kMSLL+1; i < llSize; i++)
                work[i] = 0;

        memcpy(buf, work, numBytes);
}

template<size_t numBytes>
TKmer<numBytes>::TKmer(const std::string& str, size_t offset)
{
        assert(str.size() >= (k + offset));

        const size_t llSize = (numBytes + 7) / 8;
        uint64_t work[llSize];

        const char *cstr = str.c_str() + offset;
        for (size_t i = 0; i < kMSLL; i++, cstr += 32)
                work[i] = Nucleotide::pack32(cstr);

        work[kMSLL] = Nucleotide::pack32(cstr, k % 32);

        for (size_t i = kMSLL+1; i < llSize; i++)
                work[i] = 0;

        memcpy(buf, work, numBytes);
}

template<size_t numBytes>
TKmer<numBytes>::TKmer(const TKmer<numBytes+KMERBYTEREDUCTION>& kmer, KmerLSB &lsb)
{
        memcpy(buf, kmer.buf + KMERBYTEREDUCTION, numBytes);
        // FIXME: Check for generic KMERBYTEREDUCTION
        lsb = kmer.buf[0];
        lsb += KmerLSB(kmer.buf[1]) << 8;
}

template<size_t numBytes>
TKmer<numBytes>::TKmer(const TKmer<numBytes-KMERBYTEREDUCTION>& kmer, KmerLSB lsb)
{
        memcpy(buf + KMERBYTEREDUCTION, kmer.buf, numBytes - KMERBYTEREDUCTION);
        buf[0] = lsb;
        buf[1] = lsb >> 8;
}

template<size_t numBytes>
std::string TKmer<numBytes>::str() const
{
        std::string str;
        str.resize(k);

        // fill in the string
        for (size_t p = 0, bitID = 0, byteID = 0; p < k; p++) {
                str[p] = Nucleotide::nucleotideToChar(buf[byteID] >> bitID);
                bitID += 2;
                if (bitID == 8) {
                        byteID++;
                        bitID = 0;
                }
        }

        return str;
}

template<size_t numBytes>
TKmer<numBytes>::TKmer(std::ifstream& ifs)
{
        ifs.read((char*)buf, kMSB+1);
        memset(buf+kMSB+1, 0, numBytes-kMSB-1);
}

template<size_t numBytes>
void TKmer<numBytes>::pushNucleotideRight(char c)
{
        const size_t llSize = (numBytes + 7) / 8;
        uint64_t work[llSize];
        work[llSize - 1] = 0;

        memcpy(work, buf, numBytes);

        // get & remove the metaData
        uint64_t metaData = work[kMSLL] & metaMask;
        work[kMSLL] &= ~metaMask;

        uint64_t leftBits = 0, rightBits = 0;
        for (ssize_t i = llSize - 1; i >= 0; i--) {
                rightBits = work[i] & 0x3;
                work[i] = (work[i] >> 2) | leftBits;
                leftBits = rightBits << 62;
        }

        // set the left nucleotide
        work[kMSLL] |= Nucleotide::charToNucleotide(c) << (2*(k % 32) - 2);

        // restore the metadata
        work[kMSLL] &= ~metaMask;
        work[kMSLL] |= metaData;

        memcpy(buf, work, numBytes);
}

template<size_t numBytes>
void TKmer<numBytes>::pushNucleotideLeft(char c)
{
        const size_t llSize = (numBytes + 7) / 8;
        uint64_t work[llSize];

        memcpy(work, buf, numBytes);

        // get & remove the metaData
        uint64_t metaData = work[kMSLL] & metaMask;
        work[kMSLL] &= ~metaMask;

        uint64_t leftBits = 0, rightBits = Nucleotide::charToNucleotide(c);
        for (size_t i = 0; i < llSize; i++) {
                leftBits = work[i] & (uint64_t(0x3) << 62);
                work[i] = (work[i] << 2) | rightBits;
                rightBits = leftBits >> 62;
        }

        // restore the metadata
        work[kMSLL] &= ~metaMask;
        work[kMSLL] |= metaData;

        memcpy(buf, work, numBytes);
}

template<size_t numBytes>
void TKmer<numBytes>::reverse()
{
        const size_t llSize = (numBytes + 7) / 8;
        uint64_t work[llSize];

        memcpy(work, buf, numBytes);

        // get the metaData
        uint64_t metaData = work[kMSLL] & metaMask;

        // invert all the words (64 bits galore)
        for (size_t i = 0; i <= kMSLL; i++) {
                uint64_t &w = work[i];

                w = (((w & 0xccccccccccccccccull) >> 2) |
                     ((w & 0x3333333333333333ull) << 2));
                w = (((w & 0xf0f0f0f0f0f0f0f0ull) >> 4) |
                     ((w & 0x0f0f0f0f0f0f0f0full) << 4));
                w = (((w & 0xff00ff00ff00ff00ull) >> 8) |
                     ((w & 0x00ff00ff00ff00ffull) << 8));
                w = (((w & 0xffff0000ffff0000ull) >> 16) |
                     ((w & 0x0000ffff0000ffffull) << 16));
                w = (((w & 0xffffffff00000000ull) >> 32) |
                     ((w & 0x00000000ffffffffull) << 32));
        }

        // swap all words
        for (size_t i = 0; i < (kMSLL + 1) / 2; i++) {
                uint64_t temp = work[i];
                work[i] = work[kMSLL-i];
                work[kMSLL-i] = temp;
        }

        // shift the words to the right
        uint64_t leftBits = 0, rightBits = 0;
        int numLeftBits = 2*(k % 32);
        int numRightBits = 64 - numLeftBits;
        uint64_t rightMask = (uint64_t(1) << numRightBits) - 1;
        for (ssize_t i = kMSLL; i >= 0; i--) {
                rightBits = work[i] & rightMask;
                work[i] = (work[i] >> numRightBits) | leftBits;
                leftBits = rightBits << numLeftBits;
        }

        // restore the metaData
        work[kMSLL] |= metaData;

        memcpy(buf, work, numBytes);
}

template<size_t numBytes>
void TKmer<numBytes>::complement()
{
        const size_t llSize = (numBytes + 7) / 8;
        uint64_t work[llSize];

        memcpy(work, buf, numBytes);

        // get the metaData
        uint64_t metaData = work[kMSLL] & metaMask;

        for (size_t i = 0; i <= kMSLL; i++)
                work[i] = ~work[i];

        work[kMSLL] &= (uint64_t(1) << 2*(k % 32)) - 1;
        work[kMSLL] |= metaData;

        memcpy(buf, work, numBytes);
}

template<size_t numBytes>
void TKmer<numBytes>::reverseComplement()
{
        const size_t llSize = (numBytes + 7) / 8;
        uint64_t work[llSize];

        memcpy(work, buf, numBytes);

        // get the metaData
        uint64_t metaData = work[kMSLL] & metaMask;

        // invert all the words (64 bits galore)
        for (size_t i = 0; i <= kMSLL; i++) {
                uint64_t &w = work[i];

                w = (((w & 0xccccccccccccccccull) >> 2) |
                     ((w & 0x3333333333333333ull) << 2));
                w = (((w & 0xf0f0f0f0f0f0f0f0ull) >> 4) |
                     ((w & 0x0f0f0f0f0f0f0f0full) << 4));
                w = (((w & 0xff00ff00ff00ff00ull) >> 8) |
                     ((w & 0x00ff00ff00ff00ffull) << 8));
                w = (((w & 0xffff0000ffff0000ull) >> 16) |
                     ((w & 0x0000ffff0000ffffull) << 16));
                w = (((w & 0xffffffff00000000ull) >> 32) |
                     ((w & 0x00000000ffffffffull) << 32));
                w = ~w;
        }

        // swap all words
        for (size_t i = 0; i < (kMSLL + 1) /2; i++) {
                uint64_t temp = work[i];
                work[i] = work[kMSLL - i];
                work[kMSLL - i] = temp;
        }

        // shift the words to the right
        uint64_t leftBits = 0, rightBits = 0;
        int numLeftBits = 2*(k % 32);
        int numRightBits = 64 - numLeftBits;
        uint64_t rightMask = (uint64_t(1) << numRightBits) - 1;
        for (ssize_t i = kMSLL; i >= 0; i--) {
                rightBits = work[i] & rightMask;
                work[i] = (work[i] >> numRightBits) | leftBits;
                leftBits = rightBits << numLeftBits;
        }

        // restore the metaData
        work[kMSLL] |= metaData;

        memcpy(buf, work, numBytes);
}

template<size_t numBytes>
bool TKmer<numBytes>::operator==(const TKmer<numBytes> &rhs) const
{
        const size_t llSize = (numBytes + 7) / 8;

        uint64_t w1[llSize], w2[llSize];
        w1[kMSLL] = w2[kMSLL] = 0;

        memcpy(w1, buf, numBytes);
        w1[kMSLL] &= ~metaMask;

        memcpy(w2, rhs.buf, numBytes);
        w2[kMSLL] &= ~metaMask;

        for (size_t i = 0; i <= kMSLL; i++)
                if (w1[i] != w2[i])
                        return false;

        return true;
}

template<size_t numBytes>
bool TKmer<numBytes>::operator<(const TKmer<numBytes> &rhs) const
{
        const size_t llSize = (numBytes + 7) / 8;

        uint64_t w1[llSize], w2[llSize];
        w1[kMSLL] = w2[kMSLL] = 0;

        memcpy(w1, buf, numBytes);
        w1[kMSLL] &= ~metaMask;

        memcpy(w2, rhs.buf, numBytes);
        w2[kMSLL] &= ~metaMask;

        for (ssize_t i = kMSLL; i >= 0; i--) {
                if (w1[i] != w2[i])
                        return (w1[i] < w2[i]);
        }

        return false;
}

template<size_t numBytes>
bool TKmer<numBytes>::operator>(const TKmer<numBytes> &rhs) const
{
        const size_t llSize = (numBytes + 7) / 8;

        uint64_t w1[llSize], w2[llSize];
        w1[kMSLL] = w2[kMSLL] = 0;

        memcpy(w1, buf, numBytes);
        w1[kMSLL] &= ~metaMask;

        memcpy(w2, rhs.buf, numBytes);
        w2[kMSLL] &= ~metaMask;

        for (ssize_t i = kMSLL; i >= 0; i--)
                if (w1[i] != w2[i])
                        return (w1[i] > w2[i]);

        return false;
}

template<size_t numBytes>
size_t TKmer<numBytes>::getHash() const
{
        const size_t llSize = (numBytes + 7) / 8;

        uint64_t work[llSize];
        work[kMSLL] = 0;

        memcpy(work, buf, numBytes);
        work[kMSLL] &= ~metaMask;

        size_t hash = 0;
        for (size_t i = 0; i <= kMSLL; i++) {
                uint64_t &w = work[i];

                w = ~w + (w << 21); // key = (key << 21) - key - 1;
                w = w ^ (w >> 24);
                w = (w + (w << 3)) + (w << 8); // key * 265
                w = w ^ (w >> 14);
                w = (w + (w << 2)) + (w << 4); // key * 21
                w = w ^ (w >> 28);
                w = w + (w << 31);
                hash = hash ^ size_t(w);
        }

        return hash;
}

//=============================================================================
// KMER ITERATOR
// ============================================================================

class KmerIt {

private:
        /**
         * Sets offset to the next valid kmer in the string.  Sets offset to
         * getEndPosition() if no valid kmer remains. Checks all k characters
         * for validity
         */
        void findFirstValidKmer() {
                for (size_t i = offset, v = 0; i < str.size(); i++) {
                        char c = str[i];
                        if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
                                v++;
                        else
                                v = 0;

                        if (v == Kmer::getK()) { // we've found a valid kmer
                                offset = i + 1 - Kmer::getK();
                                kmer = Kmer(str, offset);
                                return;
                        }
                }

                // no valid kmer was found
                offset = getEndPosition();
        }

        /**
         * Get the offset position reflecting the end-of-string
         * @return The end-of-string position
         */
        size_t getEndPosition() const {
                if (str.size() < Kmer::getK())
                        return 0;
                else
                        return str.size() + 1 - Kmer::getK();
        }

        /**
         * Increment the offset by a certain amount
         * @param amount Amount to increase the offset by
         */
        void incrementOffset(size_t amount = 1) {
                offset += amount;
                if (offset > getEndPosition())
                        offset = getEndPosition();
        }

        /**
         * Decrement the offset by a certain amount
         * @param amount Amount to decrease the offset by
         */
        void decrementOffset(size_t amount = 1) {
                if (int(offset)- amount < 0)
                        offset = 0;
                else
                        offset -= amount;

        }

        const std::string& str;         // const-ref to the string
        size_t offset;                  // current position in the string
        Kmer kmer;                      // current kmer

public:
        /**
         * Default constructor
         * @param tStr Tight string reference
         * @param offset Offset in the string
         */
        KmerIt(const std::string& str_) : str(str_), offset(0) {
                findFirstValidKmer();
        }

        /**
         * Prefix increment operator (move to the next kmer)
         * @return Reference to the object after incrementing
         */
        KmerIt& operator++() {
                // advance to the next character
                incrementOffset();

                // get out early
                if (!isValid())
                        return *this;

                // see if the next character is valid
                char next = str[offset + Kmer::getK() - 1];
                if (next == 'A' || next == 'C' || next == 'G' || next == 'T') {
                        kmer.pushNucleotideRight(next);
                } else {
                        // advance by k characters
                        incrementOffset(Kmer::getK());

                        // find the first valid kmer
                        findFirstValidKmer();

                        // if no valid kmer was found skip entire read
                        if (isValid())
                                kmer = Kmer(str, offset);
                }

                return *this;
        }


             /**
         * Prefix decrement operator (move to the previous kmer)
         * @return Reference to the object after decrementing
         */
        KmerIt& operator--() {
                // advance to the previous character
                decrementOffset();

                // get out early
                if (!isValid())
                        return *this;

                // see if the previous character is valid
                char previous = str[offset ];
                if (previous == 'A' || previous == 'C' || previous == 'G' || previous == 'T') {
                        kmer.pushNucleotideLeft(previous);
                } else {
                        // go back by k characters
                        decrementOffset(Kmer::getK());

                        // find the first valid kmer
                        findFirstValidKmer();

                        // if no valid kmer was found skip entire read
                        if (isValid())
                                kmer = Kmer(str, offset);
                }

                return *this;
        }


        /**
         * Postfix increment operator (move to the next kmer)
         * @return Copy of the object before incrementing
         */
        KmerIt operator++(int) {
                KmerIt copy(*this);
                operator++();
                return copy;
        }
        /**
         * Postfix decrement operator (move to the previous kmer)
         * @return Copy of the object before decrementing
         */
        KmerIt operator--(int) {
                KmerIt copy(*this);
                operator--();
                return copy;
        }
        /**
         * Get the current kmer
         * @return The current kmer
         */
        Kmer getKmer() const {
                return kmer;
        }

        /**
         * Is the current kmer valid
         * @return true or false
         */
        bool isValid() const {
                return offset != getEndPosition();
        }

        /**
         * Get the position offset of the current kmer in the string
         * @return The offset
         */
        size_t getOffset() const {
                return offset;
        }

        /**
         * Does the kmer have left overlap
         * @return true or false
         */
        bool hasLeftOverlap() const {
                if (!isValid() || (offset == 0))
                        return false;

                char c = str[offset-1];
                return (c == 'A' || c == 'C' || c == 'G' || c == 'T');
        }

        /**
         * Get the character left to the current k-mer
         * @return The character left to the current k-mer
         */
        char getLeftOverlap() const {
                assert(hasLeftOverlap());
                return str[offset-1];
        }

        /**
         * Does the kmer have right overlap
         * @return true or false
         */
        bool hasRightOverlap() const {
                if (!isValid() || ((offset+1) == getEndPosition()))
                        return false;

                char c = str[offset + Kmer::getK()];
                return (c == 'A' || c == 'C' || c == 'G' || c == 'T');
        }

        /**
         * Get the character right to the current k-mer
         * @return The character right to the current k-mer
         */
        char getRightOverlap() const {
                assert(hasRightOverlap());

                return str[offset + Kmer::getK()];
        }
};

#endif
