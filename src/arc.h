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

#ifndef ARC_H
#define ARC_H

#include "global.h"
#include <atomic>

// ============================================================================
// ARC CLASS
// ============================================================================

class Arc {

private:
        NodeID nodeID;                  // ID of node to which the arc points
        std::atomic<Coverage> cov;      // arc coverage
        std::atomic<bool> flag;         // arc flag

#ifdef DEBUG
        bool trueArc;                   // does the arc exist in the genome?
#endif

public:
        /**
         * Default constructor
         */
        Arc() : nodeID(0), cov(0), flag(false) {
#ifdef DEBUG
                trueArc = false;
#endif
        };

        /**
         * Copy constructor
         */
        Arc(const Arc& rhs) : nodeID(rhs.nodeID), cov(rhs.cov.load()), flag(rhs.flag.load()) {
#ifdef DEBUG
                trueArc = rhs.trueArc;
#endif
        }

        /**
         * Assignment operator
         */
        Arc& operator=(const Arc& rhs) {
                nodeID = rhs.nodeID;
                cov = rhs.cov.load();
                flag = rhs.flag.load();
#ifdef DEBUG
                trueArc = rhs.trueArc;
#endif
                return *this;
        }

        /**
         * Set the target nodeID the arc is pointing to
         * @param targetNodeID The ID of the node the arc is pointing to
         */
        void setNodeID(NodeID targetNodeID) {
                nodeID = targetNodeID;
        }

        /**
         * Get the target nodeID and the side of the target node
         * @return targetNodeID The ID of the node the arc is pointing to
         */
        NodeID getNodeID() const {
                return nodeID;
        }

        /**
         * Set the coverage of the arc
         * @param coverage Coverage of the arc
         */
        void setCoverage(Coverage coverage) {
                cov = coverage;
        }

        /**
         * Get the arc coverage
         * @return The arc coverage
         */
        Coverage getCoverage() const {
                return cov;
        }

        /**
         * Atomically increment of the coverage
         */
        void incReadCov() {
                cov++;
        }

        /**
         * Delete arc (mark as invalid)
         */
        void deleteArc() {
                nodeID = 0;
        }

        /**
         * Check if arc is valid (== not deleted)
         * @return True of false
         */
        bool isValid() const {
                return nodeID != 0;
        }

        /**
         * Write a node to file
         * @param ofs Open output file stream
         */
        void write(std::ofstream& ofs) const {
                ofs.write((char*)&nodeID, sizeof(nodeID));
                ofs.write((char*)&cov, sizeof(cov));
        }

        /**
         * Load a node from file
         * @param ifs Open input file stream
         */
        void read(std::ifstream& ifs) {
                ifs.read((char*)&nodeID, sizeof(nodeID));
                ifs.read((char*)&cov, sizeof(cov));
        }

#ifdef DEBUG
        /**
         * Is the arc a true arc?
         * @return True of false
         */
        bool getTrueArc() const {
                return trueArc;
        }

        /**
         * Set the true arc flag
         * @param flag True of false
         */
        void setTrueArc(bool flag) {
                trueArc = flag;
        }
#endif

        /**
         * Get the arc flag
         * @return True of false
         */
        bool getFlag() const {
                return flag;
        }

        /**
         * Set the arc flag
         * @param flag True of false
         */
        void setFlag(bool flag_) {
                flag = flag_;
        }
};

#endif
