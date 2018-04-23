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

#ifndef KMERNPP_H
#define KMERNPP_H

#include "global.h"

// ============================================================================
// NODE POSITION PAIR CLASS
// ============================================================================

class NodePosPair : public std::pair<NodeID, NodePosition>
{
public:
        /**
         * Default constructor
         */
        NodePosPair() :
                std::pair<NodeID, NodePosition>(0, 0) {}

        /**
         * Default constructor
         * @param id Node identifier
         * @param pos Position identifier
         */
        NodePosPair(NodeID id, NodePosition pos) :
                std::pair<NodeID, NodePosition>(id, pos) {}

        /**
         * Get the node identifier
         * @return Node identifier
         */
        NodeID getNodeID() const {
                return first;
        }

        /**
         * Check whether the object points to valid position
         * @return True of false
         */
        bool isValid() const {
                return first != 0;
        }

        /**
         * Get the offset position
         * @return Offset position
         */
        NodePosition getPosition() const {
                return second;
        }
};

#endif
