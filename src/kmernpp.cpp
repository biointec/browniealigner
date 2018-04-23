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

#include "kmernpp.h"
#include "graph.h"
#include "settings.h"
#include "tstring.h"

using namespace std;

// ============================================================================
// KMER - NODE POSITON PAIR TABLE
// ============================================================================

void DBGraph::buildKmerNPPTable()
{
        // count the number of k-mers in the graph
        size_t numKmers = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                const DSNode& node = nodes[id];
                if (!node.isValid())
                        continue;
                numKmers += node.getMarginalLength();
        }

        kmerNPPTable.clear();
        kmerNPPTable.resize(numKmers);

        // populate the table with kmers
        for (NodeID id = 1; id <= numNodes; id++) {
                const DSNode &node = nodes[id];
                if (!node.isValid())
                        continue;
                const TString& tStr = node.getTSequence();
                Kmer kmer(tStr);
                NodePosition pos = 0;
                insertNPP(kmer, NodePosPair(id, pos++));

                for (size_t i = Kmer::getK(); i < tStr.getLength(); i++) {
                        kmer.pushNucleotideRight(tStr[i]);
                        insertNPP(kmer, NodePosPair(id, pos++));
                }
        }
}

void DBGraph::buildKmerCountTable(KmerCountTable& table)
{
        // count the number of k-mers in the graph
        size_t numKmers = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                const DSNode& node = nodes[id];
                if (!node.isValid())
                        continue;
                numKmers += node.getMarginalLength();
        }

        table.clear();
        table.resize(numKmers);

        // populate the table with kmers
        for (NodeID id = 1; id <= numNodes; id++) {
                const DSNode &node = nodes[id];
                if (!node.isValid())
                        continue;

                // process the first kmer
                const TString& tStr = node.getTSequence();
                Kmer kmer(tStr);
                table.insert(kmer, KmerCount(id));

                // process the rest of the kmers
                for (size_t i = Kmer::getK(); i < tStr.getLength(); i++) {
                        kmer.pushNucleotideRight(tStr[i]);
                        table.insert(kmer, KmerCount(id));
                }
        }
}

void DBGraph::destroyKmerNPPTable()
{
        kmerNPPTable.clear();
}

bool DBGraph::insertNPP(const Kmer& kmer, NodePosPair npp)
{
        // chose a representative kmer
        Kmer reprKmer = settings.isDoubleStranded() ?
                kmer.getRepresentative() : kmer;

        bool reverse = (kmer != reprKmer);

        // translate the nodeID and position to that representative
        if (reverse)
                revCompNPP(npp);

        // insert value in table
        auto it = kmerNPPTable.insert(pair<Kmer, NodePosPair>(reprKmer, npp));

        return it.second;
}

NodePosPair DBGraph::findNPP(Kmer const &kmer) const
{
        // choose a representative kmer
        Kmer reprKmer = settings.isDoubleStranded() ?
                kmer.getRepresentative() : kmer;

        bool reverse = (kmer != reprKmer);

        // find the kmer in the table
        auto it = kmerNPPTable.find(reprKmer);

        // if it is not found, get out
        if (it == kmerNPPTable.end())
                return NodePosPair(0, 0);

        // if is found, check whether it needs to be reverse-complemented
        NodePosPair npp = it->second;
        if (reverse)
                revCompNPP(npp);

        return npp;
}

void DBGraph::revCompNPP(NodePosPair& npp) const
{
        NodeID nodeID = npp.getNodeID();
        NodePosition pos = npp.getPosition();
        const SSNode node = getSSNode(nodeID);
        npp = NodePosPair(-nodeID, node.getMarginalLength() - 1 - pos);
}

bool DBGraph::consecutiveNPP(NodePosPair& left, NodePosPair& right) const
{
        // return false if one of the npps is invalid
        if (!left.isValid() || !right.isValid())
                return false;

        // left and right belong to the same node?
        if (left.getNodeID() == right.getNodeID())
                if ((left.getPosition()+1) == right.getPosition())
                        return true;

        // left and right belong to the connected nodes?
        SSNode leftNode = getSSNode(left.getNodeID());
        if (leftNode.getRightArc(right.getNodeID()) == NULL)
                return false;

        // if so, make sure the left kmer is the last kmer in the left node
        if (left.getPosition() != leftNode.getMarginalLength() - 1)
                return false;

        // if so, make sure the right kmer is the first kmer in the right node
        return right.getPosition() == 0;
}
