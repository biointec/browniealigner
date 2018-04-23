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

#ifndef NUCLEOTIDE_H
#define NUCLEOTIDE_H

#include "global.h"
#include <string>
#include <algorithm>

// ============================================================================
// NUCLEOTIDE CLASS
// ============================================================================

class Nucleotide {

private:
        static const NucleotideID charToNucleotideLookup[4];
        static const char charMask;
        static const char nucleotideToCharLookup[4];
        static const NucleotideID nucleotideMask;

public:
        /**
         * Convert a character into a two bit encoding
         * @param c ASCII encoding of 'A', 'C', 'G' and 'T'
         * @return 0 for 'A', 1 for 'C', 2 for 'G' and 3 for 'T'
         */
        static NucleotideID charToNucleotide(char c) {
                return charToNucleotideLookup[(c >> 1) & charMask];
        }

        /**
         * Convert a two bit encoded nucleotide into a character
         * @param n Two bit encoded nucleotide
         * @returns 'A' for 0, 'C' for 1, 'G' for 2 and T for 3
         */
        static char nucleotideToChar(NucleotideID n) {
                return nucleotideToCharLookup[n & nucleotideMask];
        }

        /**
         * Pack four nucleotides of a string into a quad (byte)
         * @param str Input string (at least of size 4)
         * @returns packed 8 bits encoding (byte)
         */
        static uint8_t packQuad(const char *str) {
                uint8_t quad = 0;
                quad |= uint8_t(Nucleotide::charToNucleotide(*str++));
                quad |= uint8_t(Nucleotide::charToNucleotide(*str++)) << 2;
                quad |= uint8_t(Nucleotide::charToNucleotide(*str++)) << 4;
                quad |= uint8_t(Nucleotide::charToNucleotide(*str++)) << 6;
                return quad;
        }

        /**
         * Convert n (max 4) nucleotides of a string into a quad (byte)
         * @param str Input string (at least of size n)
         * @param n Number of characters
         * @returns packed 8 bits encoding (byte)
         */
        static uint8_t packQuad(const char *str, size_t n) {
                uint8_t quad = 0;
                for (size_t i = 0; i < n; i++)
                        quad |= uint8_t(Nucleotide::charToNucleotide(*str++)) << 2*i;
                return quad;
        }

        /**
         * Pack n (max 32) nucleotides of a string into a 64 bit uint
         * @param str Input string (at least of size 32)
         * @returns packed 32 bits encoding (byte)
         */
        static uint64_t pack32(const char *str, size_t n = 32) {
                uint64_t quad = 0;
                for (size_t i = 0; i < n; i++)
                        quad |= uint64_t(Nucleotide::charToNucleotide(*str++)) << 2*i;
                return quad;
        }

        /**
         * Convert packed four nucleotides into a string
         * @param quad Packed four nucleotides
         * @param str Output string (at least of size 4) (output)
         * @returns packed 8 bits encoding
         */
        static void unpackQuad(uint8_t quad, char *str) {
                *str++ = nucleotideToChar(quad);
                *str++ = nucleotideToChar(quad >> 2);
                *str++ = nucleotideToChar(quad >> 4);
                *str++ = nucleotideToChar(quad >> 6);
        }

        /**
         * Convert packed four nucleotides into a string
         * @param quad Packed four nucleotides
         * @param n Number of characters
         * @param str Output string (at least of size n) (output)
         * @returns packed 8 bits encoding
         */
        static void unpackQuad(uint8_t quad, size_t n, char *str) {
                for (size_t i = 0; i < n; i++)
                        *str++ = nucleotideToChar(quad >> 2*i);
        }

        /**
         * Convert n (max 4) characters of a string into an unsigned 8 bit
         * @param str Input string (at least of size n)
         * @param n Number of characters
         * @returns packed 8 bits encoding
         */
        static uint8_t packString(const char *str, size_t n) {
                uint8_t quad = 0;
                for (size_t i = 0; i < n; i++)
                        quad |= Nucleotide::charToNucleotide(*str++) << 2*i;
                return quad;
        }

        /**
        * Convert a nucleotide into its reverse complement
         * @param n ASCII encoding of 'A', 'C', 'G' and 'T'
         * @return ASCII encoding of the reverse complement
         */
        static char getComplement(char c) {
                uint8_t n = (uint8_t)charToNucleotide(c);
                return nucleotideToChar(3-n);
        }

        /**
         * Convert a two bit encoded nucleotide into its reverse complement
         * @param n Two bit encoded nucleotide
         * @return Two bit encoded complement
         */
        static uint8_t getComplement(uint8_t n) {
                return 3 - (n & nucleotideMask);
        }

        /**
         * Reverse an stl string
         * @param str The string to be reversed (input/output)
         */
        static void reverse(std::string& str) {
                std::reverse(str.begin(), str.end());
        }

        /**
         * Complement each nucleotide in an stl string
         * @param str The string to be complemented (input/output)
         */
        static void complement(std::string& str) {
                for (size_t i = 0; i < str.size(); i++)
                        str[i] = getComplement(str[i]);
        }

        /**
         * Reverse complement an stl string
         * @param str The string to be reverse complemented (input/output)
         */
        static void revCompl(std::string& str) {
                reverse(str);
                complement(str);
        }

        /**
         * Get reverse of an stl string
         * @param str The string to be reversed
         * @return The reverse of a string
         */
        static std::string getReverse(const std::string& str) {
                std::string R = str;
                reverse(R);
                return R;
        }

        /**
         * Get complement of an stl string
         * @param str The string to be complemented
         * @return The complement of a string
         */
        static std::string getComplement(const std::string& str) {
                std::string C = str;
                complement(C);
                return C;
        }

        /**
         * Get reverse complement of an stl string
         * @param str The string to be reverse complemented
         * @return The reverse complement of a string
         */
        static std::string getRevCompl(const std::string& str) {
                std::string RC = str;
                reverse(RC);
                complement(RC);
                return RC;
        }

        /**
         * Fill in the NW alignment matrix of two sequences
         * @param seq1 Input sequence 1 (stl string)
         * @param seq2 Input sequence 2 (stl string)
         * @param M Pre-allocated alignment matrix (input / output)
         *          Note: 0th col and 0th row must have already been filled in!
         */
        static void getNWAlignMatrix(const std::string& seq1,
                                     const std::string& seq2, int **M);

        /**
         * Calculate the NW alignment given the alignment matrix
         * @param seq1 Input sequence 1 (stl string)
         * @param seq2 Input sequence 2 (stl string)
         * @param M Input alignment matrix M
         * @param align1 Aligned sequence 1 ('-' for gaps) (output)
         * @param align2 Aligned sequence 2 ('-' for gaps) (output)
         */
        static void getNWAlignment(const std::string& seq1,
                                   const std::string& seq2,
                                   int **M,
                                   std::string &align1,
                                   std::string &align2);

        /**
         * Calculate the maximum NW alignment score for a given sequence
         * @param seq Input sequence (stl string)
         * @return The maximum alignment score
         */
        static int maxNWScore(const std::string& seq);
};

#endif
