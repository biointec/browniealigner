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

#ifndef SSNODE_H
#define SSNODE_H

#include "global.h"
#include "dsnode.h"
#include "nucleotide.h"
#include "tkmer.h"

// ============================================================================
// NODE CLASS
// ============================================================================

class SSNode {

private:
        static DSNode *nodes;   // pointer to the double stranded nodes

        NodeID nodeID;          // identiffier of the node
        DSNode *dsNode;         // reference to the double stranded node

public:
        /**
         * Default constructor
         */
        SSNode() : nodeID(0), dsNode(NULL) {}

        /**
         * Constructor
         * @param dsNode Double stranded node
         * @param ID Unique identifier of the node
         */
        SSNode(DSNode* dsNode, NodeID nodeID) : nodeID(nodeID), dsNode(dsNode) {
                assert(nodeID != 0);
        }

        /**
         * Constructor
         * @param it Iterator pointing to this node
         */
        SSNode(const ArcIt& it) : nodeID(it->getNodeID()), dsNode(nodes + abs(nodeID)) {
                assert(nodeID != 0);
        }

        /**
         * Set the flag
         * @param flag True of false
         */
        void setFlag(bool flag) {
                dsNode->setFlag(flag);
        }

        /**
         * Get the flag
         * @return True of false
         */
        bool getFlag() const {
                return dsNode->getFlag();
        }

        /**
         * Get the avarge kmer coverage
         * @return The average kmer coverage
         */
        double getAvgKmerCov(){
                return (double)dsNode->getKmerCov()/(double)dsNode->getMarginalLength();
        }

        /**
         * Get the expected multiplicity for this node
         * @param avgKmerCov Global average kmer coverage
         */
        int getExpMult(double avgKmerCov) {
                return round(getAvgKmerCov() / avgKmerCov);
        }

        /**
         * Set the read start coverage
         * @param target The target read start coverage
         */
        void setReadStartCov(Coverage target) {
                dsNode->setReadStartCov(target);
        }

        /**
         * Get the read start coverage
         * @return The read start coverage
         */
        size_t getReadStartCov() const {
                return dsNode->getReadStartCov();
        }

        /**
         * Atomically increment the read start coverage
         */
        void incReadStartCov() {
                dsNode->incReadStartCov();
        }

        /**
         * Set the kmer coverage
         * @param target The kmer coverage
         */
        void setKmerCov(Coverage target) {
                dsNode->setKmerCov(target);
        }

        /**
         * Get the kmer coverage
         * @return The kmer coverage
         */
        size_t getKmerCov() const {
                return dsNode->getKmerCov();
        }

        /**
         * Atomically increment the kmer coverage
         */
        void incKmerCov() {
                dsNode->incKmerCov();
        }

        /**
         * Invalidate this node
         */
        void invalidate() {
                dsNode->invalidate();
        }

        /**
         * Check if a node is invalidated
         * @return True or false
         */
        bool isValid() const {
                if (nodeID == 0)
                        return false;
                return dsNode->isValid();
        }

        /**
         * Get the identifier of this node
         */
        NodeID getNodeID() const {
                return nodeID;
        }

        /**
         * Operator '==' overloading
         * @param rhs Right hand side SSNode
         * @return True if they're equal
         */
        bool operator==(const SSNode &rhs) const {
                if (dsNode != rhs.dsNode)
                        return false;
                return (nodeID == rhs.nodeID);
        }

        /**
         * Operator '!=' overloading
         * @param rhs Right hand side kmer
         * @return True if they're different
         */
        bool operator!=(const SSNode &rhs) const {
                return !(*this == rhs);
        }

        /**
         * Get the length of the node (# nucleotides in DNA string)
         * @return The length of the node
         */
        size_t getLength() const {
                return dsNode->getLength();
        }

        /**
         * The marginal length == length - k + 1
         * @return The marginal length of the node
         */
        size_t getMarginalLength() const {
                return dsNode->getMarginalLength();
        }

        /**
         * Get the number of left arcs
         * @return The number of left arcs
         */
        uint8_t getNumLeftArcs() const {
                return (nodeID > 0) ?
                        dsNode->getNumLeftArcs() : dsNode->getNumRightArcs();
        }

        /**
         * Get the number of right arcs
         * @return The number of right arcs
         */
        uint8_t getNumRightArcs() const {
                return (nodeID > 0) ?
                        dsNode->getNumRightArcs() : dsNode->getNumLeftArcs();
        }

        /**
         * Get the number of left arcs
         * @param The number of left arcs
         */
        void setNumLeftArcs(uint8_t numArcs) const {
                if (nodeID > 0)
                        dsNode->setNumLeftArcs(numArcs);
                else
                        dsNode->setNumRightArcs(numArcs);
        }

        /**
         * Get the number of right arcs
         * @param The number of right arcs
         */
        void setNumRightArcs(uint8_t numArcs) const {
                if (nodeID > 0)
                        dsNode->setNumRightArcs(numArcs);
                else
                        dsNode->setNumLeftArcs(numArcs);
        }

        void swapRightArcsSign() {
                if (nodeID > 0)
                        dsNode->swapRightArcsSign();
                else
                        dsNode->swapLeftArcsSign();
        }

        void copyRightArcs(SSNode &source) {
                int numRightArcs = source.getNumRightArcs();
                source.setNumRightArcs(0);
                setNumRightArcs(numRightArcs);
                setFirstRightArcID(source.getFirstRightArcID());
                // if they have opposite sign
                if ((nodeID < 0) != (source.getNodeID() < 0))
                        swapRightArcsSign();
        }

        /**
         * Get an iterator pointing to the first left arc
         * @return an iterator pointing to the first left arc
         */
        ArcIt leftBegin() const {
                return (nodeID > 0) ?
                        dsNode->leftBegin(false) : dsNode->rightBegin(true);
        }

        /**
         * Get an iterator pointing past the last left arc
         * @retur nan iterator pointing to the last left arc
         */
        ArcIt leftEnd() const {
                return (nodeID > 0) ?
                        dsNode->leftEnd(false) : dsNode->rightEnd(true);
        }

        /**
         * Get an iterator pointing to the first right arc
         * @return an iterator pointing to the first left arc
         */
        ArcIt rightBegin() const {
                return (nodeID > 0) ?
                        dsNode->rightBegin(false) : dsNode->leftBegin(true);
        }

        /**
         * Get an iterator pointing past the last right arc
         * @retur nan iterator pointing to the last right arc
         */
        ArcIt rightEnd() const {
                return (nodeID > 0) ?
                        dsNode->rightEnd(false) : dsNode->leftEnd(true);
        }

        /**
         * Delete the left arcs
         */
        void deleteAllLeftArcs() {
                if (nodeID > 0)
                        dsNode->deleteLeftArcs();
                else
                        dsNode->deleteRightArcs();
        }

        /**
         * Delete the right arcs
         */
        void deleteAllRightArcs() {
                if (nodeID > 0)
                        dsNode->deleteRightArcs();
                else
                        dsNode->deleteLeftArcs();
        }

        /**
         * Get the identifier for the left right arc
         * @return The identifier for the left right arc
         */
        ArcID getFirstLeftArcID() {
                if (nodeID > 0)
                        return dsNode->getFirstRightArcID();
                return dsNode->getFirstLeftArcID();
        }

        /**
         * Get the identifier for the first right arc
         * @return The identifier for the first right arc
         */
        ArcID getFirstRightArcID() {
                if (nodeID > 0)
                        return dsNode->getFirstRightArcID();
                return dsNode->getFirstLeftArcID();
        }

        /**
         * Set the identifier for the left right arc
         * @param The identifier for the left right arc
         */
        void setFirstLeftArcID(ArcID target) {
                if (nodeID > 0)
                        return dsNode->setFirstRightArcID(target);
                return dsNode->setFirstLeftArcID(target);
        }

        /**
         * Set the identifier for the first right arc
         * @param The identifier for the first right arc
         */
        void setFirstRightArcID(ArcID target) {
                if (nodeID > 0)
                        return dsNode->setFirstRightArcID(target);
                return dsNode->setFirstLeftArcID(target);
        }

        /**
         * Delete a specific left arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteLeftArc(NodeID targetID) {
                if (nodeID > 0)
                        return dsNode->deleteLeftArc(targetID);
                return dsNode->deleteRightArc(-targetID);
        }

        /**
         * Delete a specific right arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteRightArc(NodeID targetID) {
                if (nodeID > 0)
                        return dsNode->deleteRightArc(targetID);
                return dsNode->deleteLeftArc(-targetID);
        }

        /**
         * Get a specific left arc
         * @param targetID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* getLeftArc(NodeID targetID) const {
                if (nodeID > 0)
                        return dsNode->getLeftArc(targetID);
                return dsNode->getRightArc(-targetID);
        }

        /**
         * Get a specific right arc
         * @param targetID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* getRightArc(NodeID targetID) const {
                if (nodeID > 0)
                        return dsNode->getRightArc(targetID);
                return dsNode->getLeftArc(-targetID);
        }

        /**
         * Get the nodeID a specific left arc
         * @param nucleotide Nucleotide under consideration
         * @return Node identifier
         */
       NodeID getLeftArcNodeID(char nucleotide) const {
                for (ArcIt it = leftBegin(); it != leftEnd(); it++) {
                        if (SSNode(it).peekNucleotideMarginalRight() == nucleotide)
                                return it->getNodeID();
                }

                return 0;
        }

        /**
         * Get a specific right arc
         * @param nucleotide Nucleotide under consideration
         * @return Pointer to the specific arc, NULL if not found
         */
        NodeID getRightArcNodeID(char nucleotide) const {
                for (ArcIt it = rightBegin(); it != rightEnd(); it++)
                        if (SSNode(it).peekNucleotideMarginalLeft() == nucleotide)
                                return it->getNodeID();

                return 0;
        }

        /**
         * Replace a left arc with a different one
         * @param origID Identifier for the original target node
         * @param newID Identifier for the new target node
         */
        void replaceLeftArc(NodeID origID, NodeID newID) {
                if (nodeID > 0) {
                        if (dsNode->getLeftArc(origID) == NULL)
                                std::cout << "Paniek ! " << std::endl;
                } else
                        if (dsNode->getRightArc(-origID) == NULL)
                                std::cout << "Paniek 2 ! " << std::endl;
                if (nodeID > 0) {
                        dsNode->getLeftArc(origID)->setNodeID(newID);
                } else {
                        dsNode->getRightArc(-origID)->setNodeID(-newID);
                }
        }

        /**
         * Replace a right arc with a different one
         * @param origID Identifier for the original target node
         * @param newID Identifier for the new target node
         */
        void replaceRightArc(NodeID origID, NodeID newID) {
                if (nodeID > 0)
                        dsNode->getRightArc(origID)->setNodeID(newID);
                else
                        dsNode->getLeftArc(-origID)->setNodeID(-newID);
        }

        /**
         * Get the sequence of this node
         * @return stl string containing the sequence
         */
        std::string getSequence() const {
                std::string seq = dsNode->getSequence();
                if (nodeID < 0)
                        Nucleotide::revCompl(seq);

                return seq;
        }

        /**
         * Get a subsequence of this node
         * @param offset Start offset
         * @param len Length of node
         * @return stl string containing the sequence
         */
        std::string substr(size_t offset, size_t len) const {
                if (nodeID > 0)
                        return dsNode->substr(offset, len);
                else
                        return Nucleotide::getRevCompl(dsNode->substr(getLength() - len - offset, len));
        }

        /**
         * Get a nucleotide at a specified position ('-' for out-of bounds)
         * @param pos Position in the sequence
         * @return Nucleotide at specified position
         */
        char getNucleotide(NodePosition pos) const {
                // check for out-of-bounds
                if (pos >= getLength())
                        return '-';
                if (nodeID < 0)
                        return Nucleotide::getComplement(dsNode->getNucleotide(getLength() - pos - 1));
                else
                        return dsNode->getNucleotide(pos);
        }

        /**
         * Set the sequence of this node
         * @param str String containing only 'A', 'C', 'G' and 'T'
         */
        void setSequence(const std::string& str) {
                if (nodeID > 0)
                        dsNode->setSequence(str);
                else
                        dsNode->setSequence(Nucleotide::getRevCompl(str));
        }

        /**
         * Get the left kmer of this node
         * @return The left kmer of this node
         */
        Kmer getLeftKmer() const {
                const std::string& seq = getSequence();
                return Kmer(seq);
        }

        /**
         * Get the right kmer of this node
         * @return The right kmer of this node
         */
        Kmer getRightKmer() const {
                const std::string& seq = getSequence();
                return Kmer(seq, seq.size() - Kmer::getK());
        }

        /**
         * Get the leftmost nucleotide of this node
         * @return The leftmost nucleotide
         */
        char peekNucleotideLeft() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideLeft();
                return Nucleotide::getComplement(dsNode->peekNucleotideRight());
        }

        /**
         * Get the rightmost nucleotide of this node
         * @return The rightmost nucleotide
         */
        char peekNucleotideRight() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideRight();
                return Nucleotide::getComplement(dsNode->peekNucleotideLeft());
        }

        /**
         * Get the nucleotide at position k - 1
         * @return The nucleotide at position k - 1
         */
        char peekNucleotideMarginalLeft() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideMarginalLeft();
                return Nucleotide::getComplement(dsNode->peekNucleotideMarginalRight());
        }

        /**
         * Get the nucleotide at position size - k
         * @return The nucleotide at position size - k
         */
        char peekNucleotideMarginalRight() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideMarginalRight();
                return Nucleotide::getComplement(dsNode->peekNucleotideMarginalLeft());
        }

        void inheritRightArcs(SSNode& target) {
                // make sure the node has currently no right arcs
                assert(getNumRightArcs() == 0);

                // update the arc information for the connected nodes
                for (ArcIt it = target.rightBegin(); it != target.rightEnd(); it++) {
                        SSNode rightNode(it);
                        rightNode.replaceLeftArc(target.getNodeID(), nodeID);
                }

                // copy the arcs
                copyRightArcs(target);
        }

        /**
         * Set the static node pointer
         * @param nodes Pointer to the double stranded nodes
         */
        static void setNodePointer(DSNode *nodes_) {
                nodes = nodes_;
        }
};

#endif
