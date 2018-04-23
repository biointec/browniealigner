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

#include "fastqfile.h"
#include <cstdlib>

using namespace std;

// ============================================================================
// FASTQ FILE
// ============================================================================

bool FastQFile::getNextRead(string &read)
{
        // empty output strings
        read.clear();

        // read sequence identifier
        char c = rfHandler->getCharacter();

        // end of file might be reached
        if (!good())
                return false;
        if (c != '@')
                throw ios::failure("File doesn't appear to be in fastQ format");
        rfHandler->getLine();

        // read the actual read
        read = rfHandler->getLine();
        if (!read.empty() && read[read.size() - 1] == '\n')
                read.erase(read.size() - 1);

        // read the + line
        rfHandler->getLine();
        // read the quality scores
        rfHandler->getLine();

        return !read.empty();
}

bool FastQFile::getNextRecord(ReadRecord& record)
{
        // empty output strings
        record.clear();

        // read sequence identifier
        char c = rfHandler->peekCharacter();

        // end of file might be reached
        if (!good())
                return false;
        if (c != '@')
                throw ios::failure("File doesn't appear to be in fastQ format");
                
        record.preRead = rfHandler->getLine();

        // read the actual read
        record.read = rfHandler->getLine();
        if (!record.read.empty() && record.read.back() == '\n') {
                record.read.pop_back();
                record.postRead = '\n';
        }

        record.qualityOff = record.postRead.size();
        // read the + line
        record.postRead.append(rfHandler->getLine());
        // read the quality scores
        record.postRead.append(rfHandler->getLine());

        return !record.read.empty();
}
