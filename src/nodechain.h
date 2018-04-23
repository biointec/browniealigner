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

#ifndef NODECHAIN_H
#define NODECHAIN_H

#include "global.h"

#include <vector>
#include <set>
#include <deque>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class NodeChain;

// ============================================================================
// << DECLARATIONS
// ============================================================================

std::ostream &operator<<(std::ostream &out, const std::vector<NodeID> &path);

// ============================================================================
// SORT DECLARATIONS
// ============================================================================

bool sortChainByOcc(const NodeChain& left, const NodeChain& right);
bool sortChainByNodeID(const NodeChain& left, const NodeChain& right);

// ============================================================================
// NODE CHAIN CLASS
// ============================================================================

class NodeChain : public std::vector<NodeID>
{
private:
        size_t count;   // number of times the nodechain was observed

public:
        /**
         * Default constructor
         */
        NodeChain() : count(0) {}

        /**
         * Constructor from an input vector
         * @param input Input vector
         */
        NodeChain(const std::vector<NodeID>& input) : count(1) {
                for (const auto& it : input)
                        push_back(it);
        }

        /**
         * Increment the count by one
         */
        void incrCount() {
                count++;
        }

        /**
         * Get the count
         * @return count
         */
        size_t getCount() const {
                return count;
        }

        /**
         * Set the count
         * @param target Target count
         */
        void setCount(size_t target) {
                count = target;
        }

        /**
         * Operator < overloading
         * @param rhs Right hand side
         */
        bool operator<(const NodeChain& rhs) const;

        /**
         * Operator < overloading
         * @param rhs Right hand side
         */
        bool operator==(const NodeChain& rhs) const;

        /**
         * Get the reverse complement of the node chain
         * @return The reverse complementary chain
         */
        NodeChain getReverseComplement() const;

        /**
         * Get the representative node chain
         * @return The representative node chain
         */
        NodeChain getRepresentative() const {
                NodeChain RC = getReverseComplement();
                return (RC < *this) ? RC : *this;
        }

        /**
         * Operator<< overloading
         * @param out Output file stream (input)
         * @param ncc NodeChainContainer to display
         * @return Output file stream
         */
        friend std::ostream &operator<<(std::ostream &out, const NodeChain& nc);
};

// ============================================================================
// NODE CHAIN POSITION CLASS
// ============================================================================

class NodeChainPos : public std::pair<size_t, size_t>
{
private:
        bool reverse;

public:
        /**
         * Default constructor
         * @param first First element (nodeChain identifier)
         * @param second Second element (position in nodeChain)
         */
        NodeChainPos(size_t first, size_t second, bool reverse_ = false) :
                pair<size_t, size_t>(first, second), reverse(reverse_) {}


        /**
         * Get the identifier of the nodeChain
         * @return The identifier of the nodeChain
         */
        size_t getChainID() const {
                return first;
        }

        /**
         * Get the position of the nodeChain
         * @return The position of the nodeChain
         */
        size_t getChainPos() const {
                return second;
        }

        /**
         * Does the NodeChainPos point to a reverse complement
         * @return true or false
         */
        bool getReverse() const {
                return reverse;
        }

        /**
         * Prefix++ overloading
         * @return Reference to this object
         */
        NodeChainPos& operator++() {
                if (reverse)
                        second--;
                else
                        second++;
                return *this;
        }

        /**
         * Postfix++ overloading
         * @return Current object copy before applying operator
         */
        NodeChainPos operator++(int) {
                NodeChainPos result(*this);
                ++(*this);
                return result;
        }
};

// ============================================================================
// NODE CHAIN CONTAINER CLASS
// ============================================================================

class NodeChainContainer : public std::vector<NodeChain>
{
private:

        std::multimap<NodeID, NodeChainPos> index;
        std::map<std::pair<int64_t, int64_t>, size_t> PER;

        /**
         * Build the index
         */
        void buildIndex();

        /**
         * Get nodeChain section starting at node with identifier nodeID
         * @param nodeID Identifier of the target node
         * @param id Identifier of the nodeChain
         * @param pos Position within the nodeChain
         * @return A nodeChain section starting at a given node
         */
        NodeChain getNodeChainSection(NodeID nodeID, size_t id, size_t pos) const;

        /**
         * Replace the oldPattern by the newPattern in the container, and this
         * in both strands. The seedPattern is necessarily the first fraction
         * of the oldPattern, and will all instances of seedPattern will be replaced
         * @param seedPattern Pattern used to find (partial) occurrences of oldPattern
         * @param oldPattern Pattern to be replaced (or a fraction thereof)
         * @param newPattern New pattern
         * @return Then number of nodechains there were updated
         */
        size_t updateNodeChains(const NodeChain& seedPattern,
                                const NodeChain& oldPattern,
                                const NodeChain& newPattern);

public:
        /**
         * Default constructor
         */
        NodeChainContainer() {}

        /**
         * Default constructor
         * @param input A vector of node chains
         */
        NodeChainContainer(const std::vector<NodeChain>& input);

        /**
         * Load the container from disk
         * @param filename File names of the input files
         */
        void addContainers(const std::vector<std::string>& filenames);

        /**
         * Get all nodeChains starting at node with identifier nodeID
         * @param leftID Identifier of the left node
         * @param rightID Identifier of the right node
         * @return A vector of all nodeChains starting at this node
         */
        std::vector<NodeChain> getNodeChainTree(NodeID leftID, NodeID rightID) const;

        /**
         * Get all reductions starting at node with identifier nodeID
         * @param leftID Identifier of the left node
         * @param rightID Identifier of the right node
         * @return A vector of all reductions starting at this node
         */
        std::vector<NodeChain> getNonBranchingPath(NodeID leftID, NodeID rightID) const;

        void findOcc(const NodeChain& target, std::vector<NodeChainPos>& occ) const;

        /**
         * Update all node chains to reflect the changes imposed by a reduction
         * @param reduction Reduction under consideration
         * @return The number of node chains that were altered
         */
        size_t processReduction(const NodeChain& reduction);

        /**
         * Update all node chains to reflect the changes imposed by a concatenation
         * @param concatenation Concatenation under consideration
         * @return The number of node chains that were altered
         */
        size_t processConcatentation(const NodeChain& concatenation);

        void insertIndex(NodeID nodeID, NodeChainPos pos);
        void deleteIndex(NodeID nodeID, NodeChainPos oldPos);
        void updateIndex(NodeID nodeID, NodeChainPos oldPos, NodeChainPos newPos);
        size_t haveConsensus(std::vector<std::pair<NodeID, size_t> >& consensus) const;

        void printIndex() const;

        //std::vector<NodeID>::iterator begin() { return nodeChains.begin(); }

        /**
         * Check whether a reduction is consistent with all node chains
         * @param reduction Reduction to validate
         * @return True or false
         */
        bool isNonBranchingPath(const NodeChain& reduction) const;

        /**
         * Smooth nodeChains passing through nodeID such that they all agree
         */
        void smoothPath(NodeID nodeID);

        /**
         * Get the nodeID given a nodeChainPos
         * @param pos The nodeChainPos
         * @return nodeID at corresponding location
         */
        NodeID getNodeID(const NodeChainPos& pos) const {
                if (pos.getReverse())
                        return -at(pos.getChainID())[pos.getChainPos()];
                else
                        return at(pos.getChainID())[pos.getChainPos()];
        }

        /**
         * Get the count given a nodeChainPos
         * @param pos The nodeChainPos
         * @return count at corresponding location
         */
        size_t getCount(const NodeChainPos& pos) const {
                return at(pos.getChainID()).getCount();
        }

        /**
         * Get the index of a specific nodechain
         * @param nc Nodechain
         */
        size_t getChainIdx(const NodeChain& nc) const;

        /**
         * Check whether a node chain position is valid
         * @param pos The node chain position
         * @return true or false
         */
        bool validNCP(const NodeChainPos& pos) const {
                if (size() <= pos.getChainID())
                        return false;
                if (at(pos.getChainID()).size() <= pos.getChainPos())
                        return false;
                return true;
        }

        void removeSubChain(const NodeChainPos& pos);

        void printPER() const;

        /**
         * Operator<< overloading
         * @param out Output file stream (input)
         * @param ncc NodeChainContainer to display
         * @return Output file stream
         */
        friend std::ostream &operator<<(std::ostream &out, const NodeChainContainer& ncc);
};

#endif
