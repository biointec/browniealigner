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

#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include "tstring.h"
#include <iostream>

using namespace std;

// ============================================================================
// ALIGNMENT RESULT
// ============================================================================

class AlnRes {

//private:
public:
        /**
         * s1len and s2len can be shorter than the corresponding sequence length
         * if the alignment contains trailing gaps in either of the sequences.
         * In that case s1len/s2len contain the sequence length that were
         * aligned to characters or internal gaps in the other sequence.
         */

        int s1len;      // number of characters in s1 that were actually aligned
        int s2len;      // number of characters in s2 that were actually aligned
        int score;      // alignment score
public:
        AlnRes() : s1len(0), s2len(0), score(0) {}

        AlnRes(int s1len, int s2len, int score) :
                s1len(s1len), s2len(s2len), score(score) {}
};

// ============================================================================
// MATRIX CLASS
// ============================================================================

class Matrix {

private:
        size_t numRow;          // number of actively used rows
        size_t numCol;          // number of actively used columns
        vector<int> data;       // actual data

public:
        /**
         * Default constructor (initializes to 100x100 matrix)
         */
        Matrix() : numRow(0), numCol(0) {}

        /**
         * Default constructor
         * @param numRow Number of rows
         * @param numCol Number of columns
         */
        Matrix(size_t numRow, size_t numCol) : numRow(numRow), numCol(numCol) {
                data.resize(numRow * numCol);
        }

        /**
         * Resize the current matrix, allocate memory only if necessary
         * @param numRow Number of rows
         * @param numCol Number of columns
         */
        void resize(size_t numRow, size_t numCol) {
                this->numRow = numRow;
                this->numCol = numCol;
                if (data.size() < (numRow * numCol))
                        data.resize(numRow * numCol);
        }

        /**
         * Operator () overloading to access element (i, j) of the matrix
         * @param i Row index
         * @param j Column index
         * @return Copy of the element at position (i, j)
         */
        int operator() (int i, int j) const {
                return data[i * numCol + j];    // row-major
        }

        /**
         * Operator () overloading to access element (i, j) of the matrix
         * @param i Row index
         * @param j Column index
         * @return Reference to element at position (i, j)
         */
        int& operator() (int i, int j) {
                return data[i * numCol + j];
        }
};

// ============================================================================
// ALIGNMENT CLASS
// ============================================================================

class NWAligner {

private:
        size_t currSize;// size of matrix (as a linear array)
        int numRow;

        int maxIndel;   // maximum number of indels
        int M;          // match score
        int I;          // mismatch penalty
        int G;          // gap score
        int *matrix;    // alignment matrix

        Matrix D;       // dense score matrix

        /**
         * Allocate memory to align sequences with lengths l1 and l2
         * @param l1 sequence length 1
         * @param l2 sequence length 2
         */
        void reserveBanded(size_t l1, size_t l2);

        /**
         * Return the alignment score when aligning two characters
         * @param a character a in ACTG+N alphabet, N = wildcard
         * @param b character b in ACTG+N alphabet, N = wildcard
         * @return alignment score (diagonal)
         */
        int S(char a, char b) const;

public:
        /**
         * Default constructor
         * @param maxIndel Maximum number of insertion or deletions
         * @param M Match score
         * @param I Mismatch penalty
         * @param G Gap score
         */
        NWAligner(int maxIndel = 3, int M = 1, int I = -1, int G = -3);

        /**
         * Destructor
         */
        ~NWAligner() {
                delete [] matrix;
        }

        /**
         * Perform the alignment between two sequences
         * @param s1 First string to align
         * @param s2 Second string to align
         * @return The alignment score (higher is better)
         */
        AlnRes align(const string &s1, const string &s2);

        /**
         * Get the gap score (negative number)
         * @return The gap score
         */
        int getGapScore() const {
                return G;
        }

        /**
         * Print matrix to stdout
         */
        void printAlignment(const AlnRes& alnRes,
                            const string &s1, const string &s2) const;
         /**
          * aligns two argument according to the alignment result
          */
        void applyAlignmentToArguments(const AlnRes& alnRes,
                            string &s1, string &s2);

        void applyAlignmentToArgumentsBanded(string& s1, string& s2);

        /**
         * Remove trailing gaps from the alignment
         * @param alnRes Alignment result input
         * @return Alignment result with gaps removed (if any)
         */
        AlnRes trimTrailingGaps(const AlnRes& alnRes) const;

        /**
         * Print matrix to stdout
         */
        void printMatrix(const AlnRes& alnRes) const;

        /**
         *
         *
         */
        //AlnRes trimTrailingGapsS2(const string& s1, const string& s2) const;

        int operator() (int i, int j) const {
                int k = max(i, j);
                int l = maxIndel - i + j;
                return matrix[k * (2 * maxIndel + 1) + l];
        }

        int& operator() (int i, int j) {
                int k = max(i, j);
                int l = maxIndel - i + j;
                return matrix[k * (2 * maxIndel + 1) + l];
        }

        /**
         * Perform the banded alignment between two sequences
         * @param s1 First string
         * @param s2 Second string
         * @return The alignment score (higher is better)
         */
        int alignBanded(const string &s1, const string &s2);

        /**
         * Global alignment but don't penalize trailing gaps in either s1 OR s2
         * @param s1 First string
         * @param s2 Second string
         * @return The alignment result
         */
        AlnRes alnGlobFreeEndGap(const string &s1, const string &s2);

        /**
         * Global alignment but don't penalize trailing gaps in s1
         * @param s1 First string
         * @param s2 Second string
         * @return The alignment result
         */
        AlnRes alnGlobFreeEndGapS2(const string &s1, const string &s2);

        /**
         * Get the maximal attainable score
         * @param l length of the sequence
         * @return maximal attainable score
         */
        int getMaxScore(size_t l) const {
                return M * l;
        }

        /**
         * Print matrix to stdout
         */
        void printAlignmentBanded(const string &s1, const string &s2) const;
};

#endif
