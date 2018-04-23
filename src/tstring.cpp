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

#include "tstring.h"
#include "tkmer.h"
#include <cstring>
#include <sstream>

const uint8_t TString::charToNucleotideLookup[4] = {0, 1, 3, 2};
const char TString::charMask = 3;
const char TString::nucleotideToCharLookup[4] = {65, 67, 71, 84};
const uint8_t TString::nucleotideMask = 3;

using namespace std;

TString::TString(string str) : length(0), buf(NULL)
{
        setSequence(str);
}

TString::TString(ifstream& ifs)
{
        ifs.read((char*)&length, sizeof(length));
        if (!ifs.good())
                length = 0;

        // store the actual string
        size_t numBytes = (length + 3) / 4;
        buf = new uint8_t[numBytes];

        ifs.read((char*)buf, numBytes);
}

void TString::read(ifstream& ifs)
{
        size_t oldNumBytes = (length + 3) / 4;

        ifs.read((char*)&length, sizeof(length));
        if (!ifs.good()) {
                length = 0;
                return;
        }

        // store the actual string
        size_t numBytes = (length + 3) / 4;
        if (numBytes > oldNumBytes) {
                delete [] buf;
                buf = new uint8_t[numBytes];
        }

        ifs.read((char*)buf, numBytes);
}

void TString::setSequence(const std::string& str)
{
        size_t strSize = str.size();

        size_t numBytes = (strSize + 3) / 4;
        if (((length + 3) / 4) != numBytes) {
                delete [] buf;
                buf = new uint8_t[numBytes];
        }

        memset(buf, 0, numBytes);

        // store the length of the string
        length = strSize;

        for (size_t i = 0, byteID = 0, bitOffset = 0; i < length; i++) {
                buf[byteID] |= Nucleotide::charToNucleotide(str[i]) << bitOffset;

                bitOffset += 2;
                if (bitOffset == 8)  {
                        bitOffset = 0;
                        byteID++;
                }
        }
}

string TString::getSequence() const
{
        ostringstream oss;
        oss << *this;

        return oss.str();
}

string TString::substr(size_t offset, size_t len) const
{
        if (offset >= length)
                return string();

        len = min(len, length - offset);
        string result;
        result.reserve(len);

        size_t byteID = offset / 4, byteOff = 2 * (offset % 4);

        for (size_t i = offset; i < offset + len; i++) {
                result.push_back(Nucleotide::nucleotideToChar(buf[byteID] >> byteOff));

                byteOff += 2;
                if (byteOff == 8)  {
                        byteOff = 0;
                        byteID++;
                }
        }

        return result;
}

void TString::complement()
{
        const size_t numBytes = (length + 3) / 4;
        const size_t llSize = (numBytes + 7) / 8;
        uint64_t *work = new uint64_t[llSize];

        memcpy(work, buf, numBytes);

        for (size_t i = 0; i < llSize; i++)
                work[i] = ~work[i];

        work[llSize-1] &= (uint64_t(1) << 2*(length % 32)) - 1;
        memcpy(buf, work, numBytes);
        delete [] work;
}

void TString::reverse()
{
        const size_t numBytes = (length + 3) / 4;
        const size_t llSize = (numBytes + 7) / 8;
        uint64_t *work = new uint64_t[llSize];

        memcpy(work, buf, numBytes);

        // invert all the words (64 bits galore)
        for (size_t i = 0; i < llSize; i++) {
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
        for (size_t i = 0; i < llSize / 2; i++) {
                uint64_t temp = work[i];
                work[i] = work[llSize-1-i];
                work[llSize-1-i] = temp;
        }

        // shift the words to the right
        uint64_t leftBits = 0, rightBits = 0;
        int numLeftBits = 2*(length % 32);
        int numRightBits = 64 - numLeftBits;
        uint64_t rightMask = (uint64_t(1) << numRightBits) - 1;
        for (ssize_t i = llSize-1; i >= 0; i--) {
                rightBits = work[i] & rightMask;
                work[i] = (work[i] >> numRightBits) | leftBits;
                leftBits = rightBits << numLeftBits;
        }

        memcpy(buf, work, numBytes);
        delete [] work;
}

void TString::reverseComplement()
{
        const size_t numBytes = (length + 3) / 4;
        const size_t llSize = (numBytes + 7) / 8;
        uint64_t *work = new uint64_t[llSize];

        memcpy(work, buf, numBytes);

        // invert all the words (64 bits galore)
        for (size_t i = 0; i < llSize; i++) {
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
        for (size_t i = 0; i < llSize / 2; i++) {
                uint64_t temp = work[i];
                work[i] = work[llSize-1-i];
                work[llSize-1-i] = temp;
        }

        // shift the words to the right
        uint64_t leftBits = 0, rightBits = 0;
        int numLeftBits = 2*(length % 32);
        int numRightBits = 64 - numLeftBits;
        uint64_t rightMask = (uint64_t(1) << numRightBits) - 1;
        for (ssize_t i = llSize-1; i >= 0; i--) {
                rightBits = work[i] & rightMask;
                work[i] = (work[i] >> numRightBits) | leftBits;
                leftBits = rightBits << numLeftBits;
        }

        memcpy(buf, work, numBytes);
        delete [] work;
}

bool TString::gotoNextChar(size_t& byteID, size_t& byteOff) const
{
        if (byteOff == 0) {
                if (byteID == 0)
                        return false;
                byteOff = 8;
                byteID--;
        }

        byteOff -= 2;
        return true;
}

void TString::initOffsets(size_t& byteID, size_t& byteOff) const
{
        byteID = (length + 3) / 4 - 1;
        byteOff = ((length-1) % 4) * 2;
}

void TString::append(const TString& tString)
{
        size_t lBytes = (length + 3) / 4;
        size_t rBytes = (tString.length + 3) / 4;

        size_t dstByteID = length / 4;
        size_t dstBitOff = (length % 4) << 1;

        length += tString.length;
        size_t tBytes = (length + 3) / 4;

        // create new space to store the string
        uint8_t *dstBuf = new uint8_t[tBytes];

        // copy the first string
        memcpy(dstBuf, buf, lBytes);
        memset(dstBuf + lBytes, 0, tBytes - lBytes);

        // copy the second string
        size_t cDstBitOff = 8 - dstBitOff;

        for (size_t srcByteID = 0; srcByteID < rBytes - 1; srcByteID++) {
                dstBuf[dstByteID++] |= tString.buf[srcByteID] << dstBitOff;
                dstBuf[dstByteID] |= tString.buf[srcByteID] >> cDstBitOff;
        }

        dstBuf[dstByteID++] |= tString.buf[rBytes-1] << dstBitOff;
        if (dstByteID < tBytes)
                dstBuf[dstByteID] |= tString.buf[rBytes-1] >> cDstBitOff;

        delete [] buf;
        buf = dstBuf;
}

char TString::operator[](int index) const
{
        size_t byteID = index / 4;
        size_t bitID = 2*(index % 4);

        return Nucleotide::nucleotideToChar(buf[byteID] >> bitID);
}

std::ostream &operator<<(std::ostream &out, const TString &tString)
{
        size_t byteID = 0, byteOff = 0;

        for (size_t i = 0; i < tString.length; i++) {
                out << Nucleotide::nucleotideToChar(tString.buf[byteID] >> byteOff);

                byteOff += 2;
                if (byteOff == 8)  {
                        byteOff = 0;
                        byteID++;
                }
        }

        return out;
}
