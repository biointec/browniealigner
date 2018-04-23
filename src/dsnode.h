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

#ifndef DSNODE_H
#define DSNODE_H

#include "global.h"
#include "arc.h"
#include "tkmer.h"
#include "tstring.h"

#include <map>
#include <set>
#include <atomic>

// ============================================================================
// ARC ITERATOR CLASS
// ============================================================================

class ArcIt {

private:
        Arc* arcPtr;    // pointer to the positive arc
        bool reversed;  // true if this is an arc from a RC node
        Arc arc;        // arc copy, with nodeID sign swapped if necessary

        void setArc() {
                arc = *arcPtr;
                if (reversed)
                        arc.setNodeID(-arc.getNodeID());
        }

public:
        /**
         * Default constructor
         * @param id Identifier of the arc
         */
        ArcIt(Arc *arcPtr, bool reversed) : arcPtr(arcPtr), reversed(reversed) {
                setArc();
        }

        /**
         * Overloading of == operator
         * @return true of false
         */
        bool operator==(const ArcIt& it) {
                // make sure we're not comparing apples to pears
                assert(reversed == it.reversed);
                return arcPtr == it.arcPtr;
        }

        /**
         * Overloading the != operator
         * @return true of false
         */
        bool operator!=(const ArcIt &it) {
                // make sure we're not comparing apples to pears
                assert(reversed == it.reversed);
                return arcPtr != it.arcPtr;
        }

        /**
         * Overloading of postfix ++ operator
         * @return Copy of the iterator before the ++ operation
         */
        ArcIt operator++(int notused) {
                ArcIt copy = *this;
                arcPtr++;
                setArc();
                return copy;
        }

        /**
         * Overloading of prefix ++ operator
         * @return Reference to the iterator after ++ operator
         */
        ArcIt& operator++() {
                arcPtr++;
                setArc();
                return *this;
        }

        /**
         * Deference operator
         * @return a reference to the arc
         */
        const Arc& operator*() const {
                return arc;
        }

        /**
         * Deference operator
         * @return a reference to the arc
         */
        const Arc* operator->() const {
                return &arc;
        }
};

// ============================================================================
// DOUBLE STRANDED NODE CLASS
// ============================================================================

class DSNode {

private:
        static Arc* arcs;

        typedef union {
                struct Packed {
                        uint8_t numLeft:3;              // [0...4]
                        uint8_t numRight:3;             // [0...4]
                        uint8_t invalid:1;
                        uint8_t flag:1;
                } p;
                uint8_t up;
        } Bitfield;

        TString sequence;       // DNA sequence
        ArcID leftID;           // ID of the first left arc or merged node
        ArcID rightID;          // ID of the first right arc or merged node
        Bitfield arcInfo;       // number of arcs at each node

        std::atomic<Coverage> readStartCov;
        std::atomic<Coverage> kmerCov;
        std::atomic<bool> myFlag;

public:
        /**
         * Set the static arc pointer
         * @param arcPtr The static arc pointer
         */
        static void setArcsPointer(Arc *arcPtr) {
                arcs = arcPtr;
        }

        /**
         * Default constructor
         */
        DSNode() : leftID(0), rightID(0), readStartCov(0), kmerCov(0), myFlag(false) {
                arcInfo.up = 0;
        }

        /**
         * Set the read start coverage
         * @param target The target read start coverage
         */
	void setReadStartCov(Coverage target) {
                readStartCov = target;
        }

        /**
         * Get the read start coverage
         * @return The read start coverage
         */
        Coverage getReadStartCov() const {
                return readStartCov;
        }

        /**
         * Atomically increment the read start coverage
         */
        void incReadStartCov() {
                readStartCov++;
        }

        /**
         * Set the kmer coverage
         * @param target The kmer coverage
         */
        void setKmerCov(Coverage target) {
                kmerCov = target;
        }

        /**
         * Get the kmer coverage
         * @return The kmer coverage
         */
        Coverage getKmerCov() const {
                return kmerCov;
        }

        /**
         * Atomically increment the kmer coverage
         */
        void incKmerCov() {
                kmerCov++;
        }

        /**
         * Set the flag
         * @param flag True of false
         */
        void setFlag(bool flag) {
                //arcInfo.p.flag = (flag) ? 1 : 0;
                myFlag = flag;
        }

        /**
         * Get the flag
         * @return True of false
         */
        bool getFlag() const {
                //return arcInfo.p.flag;
                return myFlag;
        }

        void swapRightArcsSign() {
                for (int i = 0; i < arcInfo.p.numRight; i++)
                        arcs[rightID + i].setNodeID(-arcs[rightID + i].getNodeID());
        }

        void swapLeftArcsSign() {
                for (int i = 0; i < arcInfo.p.numLeft; i++)
                        arcs[leftID + i].setNodeID(-arcs[leftID + i].getNodeID());
        }

        /**
         * Get the identifier for the first left arc
         * @return The identifier for the first left arc
         */
        ArcID getFirstLeftArcID() const {
                return leftID;
        }

        /**
         * Get the identifier for the first right arc
         * @return The identifier for the first right arc
         */
        ArcID getFirstRightArcID() const {
                return rightID;
        }

        /**
         * Set the identifier for the first left arc
         * @param The identifier for the first left arc
         */
        void getFirstLeftArcID(ArcID target) {
                leftID = target;
        }

        /**
         * Get the identifier for the first right arc
         * @param The identifier for the first right arc
         */
        void setFirstRightArcID(ArcID target) {
                rightID = target;
        }

        /**
         * Invalidate a node (= mark as deleted)
         */
        void invalidate() {
                arcInfo.p.invalid = 1;
        }

        /**
         * Check if a node is valid
         * @return True or false
         */
        bool isValid() const {
                return (arcInfo.p.invalid == 0);
        }

        /**
         * Get the length of the node (# nucleotides in DNA string)
         * @return The length of the node
         */
        size_t getLength() const {
                return sequence.getLength();
        }

        /**
         * The marginal length == length - k + 1
         * @return The marginal length of the node
         */
        NodeLength getMarginalLength() const {
                return getLength() - Kmer::getK() + 1;
        }

        /**
         * Set the number of left arcs
         * @param numLeft The number of left arcs
         */
        void setNumLeftArcs(uint8_t numLeft) {
                arcInfo.p.numLeft = numLeft;
        }

        /**
         * Set the number of right arcs
         * @param numright The number of right arcs
         */
        void setNumRightArcs(uint8_t numRight) {
                arcInfo.p.numRight = numRight;
        }

        /**
         * Get the number of left arcs
         * @return The number of left arcs
         */
        uint8_t getNumLeftArcs() const {
                return arcInfo.p.numLeft;
        }

        /**
         * Get the number of right arcs
         * @return The number of right arcs
         */
        uint8_t getNumRightArcs() const {
                return arcInfo.p.numRight;
        }

        /**
         * Delete all left arcs
         */
        void deleteLeftArcs() {
                for (int i = 0; i < arcInfo.p.numLeft; i++)
                        arcs[leftID + i].deleteArc();
                arcInfo.p.numLeft = 0;
        }

        /**
         * Delete all right arcs
         */
        void deleteRightArcs() {
                for (int i = 0; i < arcInfo.p.numRight; i++)
                        arcs[rightID + i].deleteArc();
                arcInfo.p.numRight = 0;
        }

        /**
         * Get a specific left arc
         * @param nodeID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* getLeftArc(NodeID nodeID) {
                for (int i = 0; i < arcInfo.p.numLeft; i++)
                        if (arcs[leftID + i].getNodeID() == nodeID)
                                return arcs + leftID + i;
                return NULL;
        }

        /**
         * Get a specific right arc
         * @param nodeID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* getRightArc(NodeID nodeID) {
                for (int i = 0; i < arcInfo.p.numRight; i++)
                        if (arcs[rightID + i].getNodeID() == nodeID)
                                return arcs + rightID + i;
                return NULL;
        }

        /**
         * Delete a specific left arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteLeftArc(NodeID targetID);

        /**
         * Delete a specific right arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteRightArc(NodeID targetID);

        /**
         * Get an iterator pointing to the first left arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the first left arc
         */
        ArcIt leftBegin(bool reversed = false) const {
                return ArcIt(arcs + leftID, reversed);
        }

        /**
         * Get an iterator pointing past the last left arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the last left arc
         */
        ArcIt leftEnd(bool reversed = false) const {
                return ArcIt(arcs + leftID + arcInfo.p.numLeft, reversed);
        }

        /**
         * Get an iterator pointing to the first right arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the first left arc
         */
        ArcIt rightBegin(bool reversed = false) const {
                return ArcIt(arcs + rightID, reversed);
        }

        /**
         * Get an iterator pointing past the last right arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the last right arc
         */
        ArcIt rightEnd(bool reversed = false) const {
                return ArcIt(arcs + rightID + arcInfo.p.numRight, reversed);
        }

        /**
         * Set the ID of the first outgoing left arc
         * @param leftLeftID Identifier to the first outgoing left arc
         */
        void setFirstLeftArcID(uint32_t firstLeftID) {
                leftID = firstLeftID;
        }

        /**
         * Set the ID of the first outgoing right arc
         * @param firstRightID Identifier to the first outgoing right arc
         */
        void setFirstRightArcID(uint32_t firstRightID) {
                rightID = firstRightID;
        }

        /**
         * Set the sequence of this node
         * @param str String containing only 'A', 'C', 'G' and 'T'
         */
        void setSequence(const std::string& str) {
                sequence.setSequence(str);
        }

        /**
         * Get the sequence of this node
         * @return The sequence of this node
         */
        std::string getSequence() const {
                return sequence.getSequence();
        }

        /**
         * Get a subsequence of this node
         * @param offset Start offset
         * @param len Length of node
         * @return stl string containing the sequence
         */
        std::string substr(size_t offset, size_t len) const {
                return sequence.substr(offset, len);
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
                return sequence[pos];
        }

        /**
         * Get the tight sequence of this node
         * @return The tight sequence
         */
        const TString& getTSequence() const {
                return sequence;
        }

        /**
         * Get the leftmost nucleotide of this node
         * @return The leftmost nucleotide
         */
        char peekNucleotideLeft() const {
                return sequence.peekNucleotideLeft();
        }

        /**
         * Get the rightmost nucleotide of this node
         * @return The rightmost nucleotide
         */
        char peekNucleotideRight() const {
                return sequence.peekNucleotideRight();
        }

        /**
         * Get the nucleotide at position k - 1
         * @return The nucleotide at position k - 1
         */
        char peekNucleotideMarginalLeft() const {
                return sequence.peekNucleotideMarginalLeft();
        }

        /**
         * Get the nucleotide at position size - k
         * @return The nucleotide at position size - k
         */
        char peekNucleotideMarginalRight() const {
                return sequence.peekNucleotideMarginalRight();
        }

        /**
         * Get the leftmost kmer of this node
         * @return The leftmost kmer
         */
        Kmer getLeftKmer() const {
                return Kmer(sequence, 0);
        }

        /**
         * Get the rightmost kmer of this node
         * @return The rightmost kmer
         */
        Kmer getRightKmer() const {
                return Kmer(sequence, sequence.getLength() - Kmer::getK());
        }

        /**
         * Write a node to file
         * @param ofs Open output file stream
         */
        void write(std::ofstream& ofs) const {
                double kmerCov = getKmerCov();
                Coverage readStCov = getReadStartCov();

                ofs.write((char*)&kmerCov, sizeof(kmerCov));
                ofs.write((char*)&readStCov, sizeof(readStCov));
                ofs.write((char*)&leftID, sizeof(leftID));
                ofs.write((char*)&rightID, sizeof(rightID));
                ofs.write((char*)&arcInfo, sizeof(arcInfo));

                sequence.write(ofs);
        }

        /**
         * Load a node from file
         * @param ifs Open input file stream
         */
        void read(std::ifstream& ifs) {
                double kmerCov;
                Coverage readStCov;

                ifs.read((char*)&kmerCov, sizeof(kmerCov));
                ifs.read((char*)&readStCov, sizeof(readStCov));
                ifs.read((char*)&leftID, sizeof(leftID));
                ifs.read((char*)&rightID, sizeof(rightID));
                ifs.read((char*)&arcInfo, sizeof(arcInfo));

                setKmerCov(kmerCov);
                setReadStartCov(readStCov);

                sequence.read(ifs);
        }
};

#endif
