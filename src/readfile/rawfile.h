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

#ifndef RAWFILE_H
#define RAWFILE_H

#include "readfile.h"

class RawFile : public ReadFile
{

public:
        /**
         * Default constructor
         * @param gzipped True if the file is gzipped
         */
        RawFile(bool gzipped) : ReadFile(gzipped) {};

        /**
         * Get the next read from a file
         * @param read String containing the read (output)
         */
        bool getNextRead(std::string &read);

        /**
         * Get the next record from a file
         * @param record Record to store output
         */
        bool getNextRecord(ReadRecord& record);
};

#endif
