/***************************************************************************
 *   Copyright (C) 2014 - 2017 Jan Fostier (jan.fostier@ugent.be)          *
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

#include "alignment.h"
#include <string>
#include <stdlib.h>
#include <stdexcept>
#include <sstream>

using namespace std;

// ============================================================================
// ALIGNMENT CLASS
// ============================================================================

void NWAligner::reserveBanded(size_t l1, size_t l2)
{
        // check the dimensions of s1 and s2
        int thisMaxDim = max(l1, l2);
        int thisMinDim = min(l1, l2);
        if ((thisMaxDim - thisMinDim) > maxIndel) {
                ostringstream oss;
                oss << "Cannot align sequences with length " << l1 << " and "
                    << l2 << " as the maximum indel is " << maxIndel << "\n";
                throw runtime_error(oss.str());
        }

        size_t reqSize = (2*maxIndel+1) * (thisMaxDim+1);

        // reallocate memory if necessary
        if (reqSize > currSize) {
                currSize = reqSize;
                delete [] matrix;
                matrix = new int[currSize];
        }
}

int NWAligner::S(char a, char b) const
{
        if (a == b)
                return M;
        if ((a == 'N') || (b == 'N'))
                return M;
        return I;
}

AlnRes NWAligner::alnGlobFreeEndGap(const string& s1, const string& s2)
{
        NWAligner& F = *this;   // shorthand notation

        // perform a normal, end-to-end alignment
        alignBanded(s1, s2);

        // remove trailing gaps
        int i = s1.size();
        int j = s2.size();

        // remove gaps at the end of s2, if any
        while ((i > 0) && (j < i + maxIndel) && (F(i, j) == F(i-1, j) + G))
                i--;
        if (i < s1.size())
                return AlnRes(i, j, F(i, j));

        // remove gaps at the end of s1, if any
        while ((j > 0) && (j > i - maxIndel) && (F(i, j) == F(i, j-1) + G))
                j--;
        if (j < s2.size())
                return AlnRes(i, j, F(i, j));

        // there were no trailing gaps, return the global aln score
        return AlnRes(i, j, F(i, j));
}

AlnRes NWAligner::alnGlobFreeEndGapS2(const string& s1, const string& s2)
{
        NWAligner& F = *this;   // shorthand notation

        // quick exit in case of too many gaps  FIXME
        if (abs((int)s1.length() - (int)s2.length()) > maxIndel)
                return AlnRes(0, 0, G*(s1.length() + s2.length()));

        // perform a normal, end-to-end alignment
        alignBanded(s1, s2);

        // remove trailing gaps
        int i = s1.size();
        int j = s2.size();

        // remove gaps at the end of s2, if any
        while ((i > 0) && (j < i + maxIndel) && (F(i, j) == F(i-1, j) + G))
                i--;
        if (i < s1.size())
                return AlnRes(i, j, F(i, j));

        // there were no trailing gaps, return the global aln score
        return AlnRes(i, j, F(i, j));
}

AlnRes NWAligner::align(const string& s1, const string& s2)
{
        // reserve memory (+1 for first/row column with gap penalties)
        D.resize(s1.length() + 1, s2.length() + 1);

        // initialize the borders of the matrix
        for (int i = 0; i <= s1.length(); i++)
                D(i, 0) = i * G;

        // initialize the borders of the matrix
        for (int j = 0; j <= s2.length(); j++)
                D(0, j) = j * G;

        // initialize the rest of the bulk of the matrix
        for (int i = 1; i <= (int)s1.length(); i++) {
                for (int j = 1; j <= (int)s2.length(); j++) {
                        int diag = D(i-1, j-1) + S(s1[i-1], s2[j-1]);
                        int up   = D(i-1, j) + G;
                        int left = D(i, j-1) + G;

                        D(i, j) = max(diag, max(up, left));
                }
        }

        return AlnRes( s1.length(), s2.length(), D(s1.length(), s2.length()) );
}

AlnRes NWAligner::trimTrailingGaps(const AlnRes& alnRes) const
{
        // remove trailing gaps
        int i = alnRes.s1len;
        int j = alnRes.s2len;

        // remove gaps at the end of s2, if any
        while ( (i > 0) && (D(i, j) == D(i-1, j) + G) )
                i--;
        if (i < alnRes.s1len)
                return AlnRes(i, j, D(i, j));

        // remove gaps at the end of s1, if any
        while ( (j > 0) && (D(i, j) == D(i, j-1) + G) )
                j--;
        if (j < alnRes.s2len)
                return AlnRes(i, j, D(i, j));

        // there were no trailing gaps, return the global aln score
        return alnRes;
}

int NWAligner::alignBanded(const string& s1, const string& s2)
{
        NWAligner& F = *this;   // shorthand notation

        // allocate memory if necessary
        reserveBanded(s1.length(), s2.length());

        // initialize the borders of the matrix
        for (int i = 0; i <= maxIndel; i++) {
                F(i, 0) = i * G;
                F(0, i) = i * G;
        }

        // initialize the rest of the bulk of the matrix
        for (int i = 1; i <= (int)s1.length(); i++)
        {
                int startj = max(1, i - maxIndel);
                int endj = min((int)s2.length(), i + maxIndel);

                for (int j = startj; j <= endj; j++)
                {
                        int diag = F(i-1, j-1) + S(s1[i-1], s2[j-1]);
                        int up   = (j < i + maxIndel) ? F(i-1, j) + G : diag-1;
                        int left = (j > i - maxIndel) ? F(i, j-1) + G : diag-1;

                        F(i, j) = max(diag, max(up, left));
                }
        }

        return F(s1.length(), s2.length());
}

NWAligner::NWAligner(int maxIndel, int M, int I, int G) :
        currSize(100), maxIndel(maxIndel), M(M), I(I), G(G)
{
        matrix = new int[currSize];
}

void NWAligner::printMatrix(const AlnRes& alnRes) const
{
        for (int i = 0; i <= alnRes.s1len; i++) {
                for (int j = 0; j <= alnRes.s2len; j++)
                        cout << D(i, j) << "\t";
                cout << "\n";
        }

        /*for (int l = 0; l < 2*maxIndel+1; l++) {
                for (int k = 0; k < maxDim+1; k++)
                        cout << matrix[k * (2 * maxIndel + 1) + l] << "\t";
                cout << "\n";
        }*/
}

void NWAligner::printAlignment(const AlnRes& alnRes,
                               const string& s1, const string& s2) const
{
        string al1, al2;

        int i = alnRes.s1len;
        int j = alnRes.s2len;
        while (i > 0 || j > 0) {
                if ((i > 0) && (D(i, j) == D(i-1, j) + G)) {
                        al1.push_back(s1[i-1]);
                        al2.push_back('-');
                        i--;
                } else if ((j > 0) && (D(i, j) == D(i, j-1) + G)) {
                        al1.push_back('-');
                        al2.push_back(s2[j-1]);
                        j--;
                } else {
                        al1.push_back(s1[i-1]);
                        al2.push_back(s2[j-1]);
                        i--;
                        j--;
                }
        }

        reverse(al1.begin(), al1.end());
        reverse(al2.begin(), al2.end());

        cout << al1 << "\n";
        for (int i = 0; i < max((int)al1.size(), (int)al2.size()); i++)
                if (al1[i] == al2[i] || al1[i] == 'N' || al2[i] == 'N')
                        cout << "|";
                else
                        cout << "*";
        cout << "\n" << al2 << "\n";
}

void NWAligner::printAlignmentBanded(const string& s1, const string& s2) const
{
        const NWAligner& F = *this;   // shorthand notation

        string al1, al2;

        int i = s1.size();
        int j = s2.size();
        while (i > 0 || j > 0) {
                if ((i > 0) && (j > 0) && (F(i, j) == F(i-1, j-1) + S(s1[i-1], s2[j-1]))) {
                        al1.push_back(s1[i-1]);
                        al2.push_back(s2[j-1]);
                        i--;
                        j--;
                } else if ((i > 0) && (j < i + maxIndel) && (F(i, j) == F(i-1, j) + G)) {
                        al1.push_back(s1[i-1]);
                        al2.push_back('-');
                        i--;
                } else {
                        al1.push_back('-');
                        al2.push_back(s2[j-1]);
                        j--;
                }
        }

        reverse(al1.begin(), al1.end());
        reverse(al2.begin(), al2.end());

        cout << al1 << "\n";
        for (int i = 0; i < max((int)al1.size(), (int)al2.size()); i++)
                if (al1[i] == al2[i] || al1[i] == 'N' || al2[i] == 'N')
                        cout << "|";
                else
                        cout << "*";
        cout << "\n" << al2 << "\n";
}
