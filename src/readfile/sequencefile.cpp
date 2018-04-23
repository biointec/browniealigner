/***************************************************************************
 *   Copyright (C) 2010 Jan Fostier (jan.fostier@intec.ugent.be)           *
 *   Original Velvet code by Daniel Zerbino (zerbino@ebi.ac.uk)            *
 *                                                                         *
 *   This file is part of Velvet 2.0                                       *
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

#include "sequencefile.h"

using namespace std;

// ============================================================================
// SEQUENCE FILE
// ============================================================================

bool SequenceFile::getNextRead(string &read, string &description)
{
        description = rfHandler->getLine();
        if (description.empty())
                return false;
        if (description[0] != '>')
                throw ios_base::failure("File doesn't appear to a "
                                        "Velvet sequence file");
        description.erase(description.begin());

        if (!rfHandler->good())
                return false;

        read = rfHandler->getLine();

        return !read.empty();
}

void SequenceFile::writeRead(const string& read, const string& description)
{
        rfHandler->writeChar('>');
        rfHandler->writeLine(description);
        rfHandler->writeLine(read);
}
