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

#include "dsnode.h"

Arc* DSNode::arcs = NULL;

bool DSNode::deleteLeftArc(NodeID targetID)
{
        int i = leftID;
        for ( ; i < leftID + arcInfo.p.numLeft; i++)
                if (arcs[i].getNodeID() == targetID)
                        break;

        // arc is not found, get out
        if (i == leftID + arcInfo.p.numLeft)
                return false;

        // shift remaing arcs one position down
        for (i++; i < leftID + arcInfo.p.numLeft; i++)
                arcs[i - 1] = arcs[i];

        arcInfo.p.numLeft--;
        arcs[leftID + arcInfo.p.numLeft].deleteArc();

        return true;
}

bool DSNode::deleteRightArc(NodeID targetID)
{
        ArcID i = rightID;
        for ( ; i < rightID + arcInfo.p.numRight; i++)
                if (arcs[i].getNodeID() == targetID)
                        break;

        // arc is not found, get out
        if (i == rightID + arcInfo.p.numRight)
                return false;

        // shift remaing arcs one position down
        for (i++; i < rightID + arcInfo.p.numRight; i++)
                arcs[i - 1] = arcs[i];

        arcInfo.p.numRight--;
        arcs[rightID + arcInfo.p.numRight].deleteArc();

        return true;
}
