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

#include "nucleotide.h"

using namespace std;

// ============================================================================
// NUCLEOTIDE CLASS
// ============================================================================

const NucleotideID Nucleotide::charToNucleotideLookup[4] = {0, 1, 3, 2};
const char Nucleotide::charMask = 3;
const char Nucleotide::nucleotideToCharLookup[4] = {65, 67, 71, 84};
const NucleotideID Nucleotide::nucleotideMask = 3;

int max(int a, int b, int c)
{
        int t = (a > b) ? a : b;
        return (t > c) ? t : c;
}

void Nucleotide::getNWAlignMatrix(const string& seq1, const string& seq2,
                                  int **M)
{
        for (size_t i = 1; i <= seq1.size(); i++) {
                NucleotideID c1 = Nucleotide::charToNucleotide(seq1[i-1]);
                for (size_t j = 1; j <= seq2.size(); j++) {
                        NucleotideID c2 = Nucleotide::charToNucleotide(seq2[j-1]);
                        int s1 = M[i-1][j-1] + scoreMatrix[c1][c2];
                        int s2 = M[i-1][j] + indelScore;
                        int s3 = M[i][j-1] + indelScore;
                        M[i][j] = max(s1, s2, s3);
                }
        }
}

void Nucleotide::getNWAlignment(const string& seq1, const string& seq2,
                                int **M, string &align1, string &align2)
{
        align1.clear();
        align2.clear();
        align1.reserve(seq1.size() + seq2.size());
        align2.reserve(seq1.size() + seq2.size());

        // backtracking from NW matrix
        size_t i = seq1.size(), j = seq2.size();
        while (i > 0 && j > 0) {
                NucleotideID c1 = Nucleotide::charToNucleotide(seq1[i-1]);
                NucleotideID c2 = Nucleotide::charToNucleotide(seq2[j-1]);

                if (M[i][j] == M[i-1][j-1] + scoreMatrix[c1][c2]) {
                        align1.push_back(seq1[--i]);
                        align2.push_back(seq2[--j]);
                } else if (M[i][j] == M[i-1][j] + indelScore) {
                        align1.push_back(seq1[--i]);
                        align2.push_back('-');
                } else if (M[i][j] == M[i][j-1] + indelScore) {
                        align1.push_back('-');
                        align2.push_back(seq2[--j]);
                } else {
                        // should never get here !
                        assert(false);
                }
        }

        while (i-- > 0)
                align1.push_back('-');
        while (j-- > 0)
                align2.push_back('-');

        reverse(align1);
        reverse(align2);
}

int Nucleotide::maxNWScore(const std::string& seq)
{
        int score = 0;
        for (size_t i = 0; i < seq.size(); i++) {
                NucleotideID c = Nucleotide::charToNucleotide(seq[i]);
                score += scoreMatrix[c][c];
        }

        return score;
}
