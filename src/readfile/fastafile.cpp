/***************************************************************************
 *   Copyright (C) 2010-215 Jan Fostier (jan.fostier@intec.ugent.be)       *
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

#include "fastafile.h"
#include <cstdlib>

using namespace std;

// ============================================================================
// FASTA FILE
// ============================================================================

bool FastAFile::getNextRead(string &read)
{
        // empty output strings
        read.clear();

        // read lines until '>' encountered (skip ; lines)
        string dummy;
        while (rfHandler->good() && rfHandler->peekCharacter() != '>')
                rfHandler->getLine();

        // read description lines
        if (rfHandler->good())
                rfHandler->getLine();

        // read the actual read
        while (rfHandler->good() && rfHandler->peekCharacter() != '>') {
                const char *result = rfHandler->getLine();
                read.append(result);
                if (!read.empty() && read[read.size() - 1] == '\n')
                        read.erase(read.size() - 1);
        }

        return !read.empty();
}

bool FastAFile::getNextRecord(ReadRecord& record)
{
        record.clear();

        // read lines until '>' encountered (skip ; lines)
        string dummy;
        while (rfHandler->good() && rfHandler->peekCharacter() != '>') {
                const char *result = rfHandler->getLine();
                record.preRead.append(result);
        }

        // read description lines
        if (rfHandler->good()) {
                const char *result = rfHandler->getLine();
                record.preRead.append(result);
        }

        // read the actual read
        while (rfHandler->good() && rfHandler->peekCharacter() != '>') {
                const char *result = rfHandler->getLine();
                record.read.append(result);
                if (!record.read.empty() && record.read.back() == '\n')
                        record.read.pop_back();
        }

        record.postRead = '\n';

        return !record.read.empty();
}

void FastAFile::writeRecord(const ReadRecord& record)
{
        rfHandler->writeLine(record.preRead);

        size_t offset = 0;
        while (offset < record.read.size()) {
                int remaining = record.read.size() - offset;
                rfHandler->writeLine(record.read.substr(offset, min(70, remaining)));
                offset += 70;
                if (offset < record.read.size())
                        rfHandler->writeChar('\n');
        }

        rfHandler->writeLine(record.postRead);
}


