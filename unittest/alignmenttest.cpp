/***************************************************************************
 *   Copyright (C) 2015 Jan Fostier (jan.fostier@intec.ugent.be)           *
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

#include <gtest/gtest.h>
#include "../src/alignment.h"

using namespace std;

TEST(Alignment, AlignmentTest)
{
        NWAligner align(3, 1, -1, -3);

        string s1 = "ACGTACGTAC";
        string s2 = "ACGTACGTAC";

        AlnRes res = align.align(s1, s2);

        //align.printAlignment(s1, s2);

        ASSERT_EQ(res.score, 10);

        s2 = "ACGTACGCAC";  // one mismatch
        res = align.align(s1, s2);
        //align.printAlignment(s1, s2);

        ASSERT_EQ(res.score, 8);

        s2 = "ACGTACCGTAC";  // one insertion
        res = align.align(s1, s2);
        //align.printAlignment(s1, s2);

        ASSERT_EQ(res.score, 7);

        s2 = "ACGTAGTAC";  // one deletion
        res = align.align(s1, s2);
        //align.printAlignment(s1, s2);

        ASSERT_EQ(res.score, 6);
}
