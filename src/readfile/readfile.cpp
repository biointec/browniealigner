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

#include <cstring>
#include <cassert>
#include "readfile.h"

using namespace std;

#define BUFFER_SIZE 4096

// ============================================================================
// REGULAR READFILE HANLDER
// ============================================================================

bool RegularReadFileHandler::open(const string& filename, ReadFileMode mode)
{
        this->mode = mode;
        if (mode == READ)
                fh = fopen(filename.c_str(), "r");
        else
                fh = fopen(filename.c_str(), "w");

        if (fh == NULL)
                throw ios_base::failure("Cannot open file: " + filename);
        return true;
}

void RegularReadFileHandler::reset()
{
        if (mode == READ)
                fseek(fh, 0, SEEK_SET);
}

void RegularReadFileHandler::close()
{
        fclose(fh);
}

// ============================================================================
// GZIPPED READFILE HANDLER
// ============================================================================

#ifdef HAVE_ZLIB

bool GZipReadFileHandler::open(const std::string& filename, ReadFileMode mode)
{
        this->mode = mode;
        const char* m = (mode == READ) ? "r" : "w";
        ifs = gzopen(filename.c_str(), m);
        if (ifs == Z_NULL)
                throw ios_base::failure("Cannot open file: " + filename);
        return true;
}

void GZipReadFileHandler::reset() {
        gzrewind(ifs);
}

void GZipReadFileHandler::close()
{
        gzclose(ifs);
        ifs = Z_NULL;
}

#endif

// ============================================================================
// READFILE
// ============================================================================

ReadFile::ReadFile(bool gzipped) : rfHandler(NULL)
{
#ifndef HAVE_ZLIB
        if (gzipped)
                throw ios_base::failure("Velvet was compiled without "
                                        "zlib support");
        rfHandler = new RegularReadFileHandler();
#else
        if (gzipped)
                rfHandler = new GZipReadFileHandler();
        else
                rfHandler = new RegularReadFileHandler();
#endif

}

void ReadFile::writeRecord(const ReadRecord& record)
{
        rfHandler->writeLine(record.preRead);
        rfHandler->writeLine(record.read);
        rfHandler->writeLine(record.postRead);
}
