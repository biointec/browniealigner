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

#ifndef NODEENDSTABLE_H
#define NODEENDSTABLE_H

#include "global.h"
#include "tkmer.h"
#include <google/sparse_hash_map>

// ============================================================================
// NODE END METADATA
// ============================================================================

class NodeEndMD {

private:
        NodeID nodeID;  // identifier for the node

public:
        /**
         * Default constructor
         * @param nodeID Identifier for the node
         * @param side Node end specifier
         */
        NodeEndMD(NodeID nodeID) : nodeID(nodeID) {};

        /**
         * Get the node identifier
         * @return The node identifier
         */
        NodeID getNodeID() const {
                return nodeID;
        }

} ATTRIBUTE_PACKED;

// ============================================================================
// TYPEDEFS
// ============================================================================

// shortcut notation for a <Key, Data> pair
typedef std::pair<Kmer, NodeEndMD> Value;

// shortcut notation for a const iterator
typedef google::sparse_hash_map<Kmer, NodeEndMD, KmerHash>::const_iterator NETableIt;

// ============================================================================
// NODE END REFERENCE
// ============================================================================

// a node end reference is a pair of two values: a) the iterator that points to
// that node end or its reverse complement in the table.  b) a boolean to
// indicate whether the iterator points the reverse complement or not
class NodeEndRef : public std::pair<NETableIt, bool> {

public:
        /**
         * Default constructor
         */
        NodeEndRef() {}

        /**
         * Constructor
         * @param it Iterator to the table
         * @param reverse True if the iterator points to the reverse complement
         */
        NodeEndRef(NETableIt it, bool reverse) :
                std::pair<NETableIt, bool>(it, reverse) {}

        /**
         * Get the node identifier
         * @return The node identifier
         */
        NodeID getNodeID() const {
                NodeID id = first->second.getNodeID();
                return (second) ? -id : id;
        }

        /**
         * Get the kmer from the reference
         * @return The kmer
         */
        const Kmer getKmer() const {
                if (second)
                        return first->first.getReverseComplement();
                return first->first;
        }
};

// ============================================================================
// NODE END TABLE
// ============================================================================

class NodeEndTable {

private:
        google::sparse_hash_map<Kmer, NodeEndMD, KmerHash> table;    // actual table
        bool doubleStranded;    // double stranded reads or not

public:
        /**
         * Default constructor
         * @param doubleStranded Double stranded or not
         * @param size Initial size of the map
         */
        NodeEndTable(bool doubleStranded_, size_t size = 0) :
                doubleStranded(doubleStranded_) { table.resize(size); }

        /**
         * Insert a Kmer in the table
         * @param kmer Kmer to be inserted
         * @param nodeID Identifier for the node
         * @return True if the kmer is inserted, false otherwise
         */
        bool insert(const Kmer &kmer, NodeID nodeID);

        /**
         * Find a kmer in the table
         * @param kmer Kmer to look for
         * @return NodeEndRef containing iterator to the kmer and reversed flag
         */
        NodeEndRef find(const Kmer &kmer) const;

        /**
         * Get the past-end iterator
         * @return The past-end iterator
         */
        NETableIt end() const {
                return table.end();
        }
};

#endif
