/***************************************************************************
 *   Copyright (C) 2017 Jan Fostier (jan.fostier@ugent.be)                 *
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

#include <stack>
#include <queue>
#include <climits>

#include "graph.h"
#include "graphaln.h"

using namespace std;

// ============================================================================
// GRAPH ALIGNMENT CLASS
// ============================================================================

std::ostream &operator<<(std::ostream &out, const GraphAlignment& ga)
{
        for (auto it : ga.path)
                out << it.nodeID << " ";
        out << "(score: " << ga.getScore() << ")";
        return out;
}

std::string GraphAlignment::getSequence(const DBGraph& dbg) const
{
        string result;

        for (const auto& it : path) {
                string nodeSeq = dbg.getSSNode(it.nodeID).getSequence();
                result.append(nodeSeq.substr(it.nodeBegin, it.nodeEnd-it.nodeBegin));
        }

        return result;
}

double GraphAlignment::getAvgKmerCov(const DBGraph& dbg) const
{
        double totKmer = 0, totLen = 0;
        for (auto it : path) {
                double len = it.nodeEnd - it.nodeBegin;
                totKmer += dbg.getSSNode(it.nodeID).getKmerCov() * len;
                totLen += len;
        }

        return totKmer / totLen;
}

// ============================================================================
// GRAPH ALIGNER CLASS
// ============================================================================

void GraphAligner::visualizeNodeAln(const NodeAlignment& nodeAln,
                                    const GraphAlignment& ga, const string& P) const
{
        SSNode node = dbg.getSSNode(nodeAln.nodeID);

        /*for (int i = 0; i < nodeAln.pattBegin; i++)
                cout << " ";*/
        size_t nodeLen = nodeAln.nodeEnd - nodeAln.nodeBegin;
        size_t readLen = nodeAln.pattEnd - nodeAln.pattBegin;
        cout << node.substr(nodeAln.nodeBegin, nodeLen) << " (id: "
             << nodeAln.nodeID << ", node: " << nodeAln.score << "/"
             << aligner.getMaxScore(readLen) << ", path: "
             << ga.getScore() << "/" << getMaxAttainableScore(P, ga)
             << ", depth: "  << ga.getDepth()
             << ", pathlen: " << ga.getAlignedLength() << ") " << endl;
}

int GraphAligner::getMaxAttainableScore(const std::string& P,
                                        const GraphAlignment& ga,
                                        size_t shortestPath) const
{
        // don't overflow the gapPenalty computation
        if (shortestPath == numeric_limits<size_t>::max())
                return numeric_limits<int>::min();

        // if the shortest path to the dst is longer than the unaligned pattern
        size_t unalnPattLen = P.length() - ga.getAlignedLength();
        int gapPenalty = (shortestPath > unalnPattLen) ?
                aligner.getGapScore() * (shortestPath - unalnPattLen) : 0;

        // upper bound to the score that can be attained
        return ga.getScore() + aligner.getMaxScore(unalnPattLen) + gapPenalty;
}

NodeAlignment GraphAligner::alnToNodeOpenEnd(const string& patt,
                                            size_t pattBegin,
                                            const NodePosPair& nodeBegin)
{
        SSNode node = dbg.getSSNode(nodeBegin.getNodeID());

        size_t nodeLen = node.getLength() - nodeBegin.getPosition();
        size_t pattLen = patt.length() - pattBegin;

        // first guess: alignment length equals shortest of both strings
        size_t nodeAlLen = min(nodeLen, pattLen);
        size_t pattAlLen = nodeAlLen;

        string nodeAlStr, pattAlStr;
        AlnRes alnRes;

        // while alignment ends with gaps...
        do {
                nodeAlStr = node.substr(nodeBegin.getPosition(), nodeAlLen);
                pattAlStr = patt.substr(pattBegin, pattAlLen);
                alnRes = aligner.align(nodeAlStr, pattAlStr);
                alnRes = aligner.trimTrailingGaps(alnRes);
                //aligner.printAlignment(alnRes, nodeAlStr, pattAlStr);

                // extend node/pattern length if alignment end with gaps
                nodeAlLen += min(nodeLen - nodeAlLen, pattAlLen - alnRes.s2len);
                pattAlLen += min(pattLen - pattAlLen, nodeAlLen - alnRes.s1len);

        } while ((nodeAlLen > nodeAlStr.size()) || (pattAlLen > pattAlStr.size()));

        float relScore = (float)alnRes.score / (float)aligner.getMaxScore(alnRes.s2len);

        return NodeAlignment(node.getNodeID(), nodeBegin.getPosition(),
                             nodeBegin.getPosition() + alnRes.s1len, pattBegin,
                             pattBegin + alnRes.s2len, alnRes.score, relScore);
}

NodeAlignment GraphAligner::nodeAlnTo(const string& patt, size_t pattBegin,
                                      const NodePosPair& nodeBegin,
                                      const NodePosPair& nodeEnd)
{
        // input assertions
        assert(nodeBegin.getNodeID() == nodeEnd.getNodeID());
        assert(nodeEnd.getPosition() >= nodeBegin.getPosition());

        SSNode node = dbg.getSSNode(nodeBegin.getNodeID());

        size_t pattLen = patt.length() - pattBegin;

        // the node alignment length is fixed
        size_t nodeAlLen = nodeEnd.getPosition() - nodeBegin.getPosition();
        // first guess: pattern alignment length equals node length
        size_t pattAlLen = min(pattLen, nodeAlLen);

        string nodeAlStr = node.substr(nodeBegin.getPosition(), nodeAlLen);
        string pattAlStr;
        AlnRes alnRes;

        // while alignment ends with gaps...
        do {
                pattAlStr = patt.substr(pattBegin, pattAlLen);
                alnRes = aligner.align(nodeAlStr, pattAlStr);
                alnRes = aligner.trimTrailingGaps(alnRes);
                //aligner.printAlignment(alnRes, nodeAlStr, pattAlStr);

                // extend pattern length if alignment end with gaps in pattern
                pattAlLen += min(pattLen - pattAlLen, nodeAlLen - alnRes.s1len);

        } while (pattAlLen > pattAlStr.size());

        // reintroduce trailing gaps in alignment in case patt can not be
        // extended - we want end-to-end alignment of node string)
        alnRes.score += (nodeAlLen - alnRes.s1len) * aligner.getGapScore();
        alnRes.s1len = nodeAlLen;

        float relScore = (float)alnRes.score / (float)aligner.getMaxScore(alnRes.s2len);

        return NodeAlignment(node.getNodeID(), nodeBegin.getPosition(),
                             nodeEnd.getPosition(), pattBegin,
                             pattBegin + alnRes.s2len, alnRes.score, relScore);
}

NodeAlignment GraphAligner::nodeAln(const string& patt,
                                    size_t pattBegin, size_t pattEnd,
                                    const NodePosPair& nodeBegin,
                                    const NodePosPair& nodeEnd)
{
        // input assertions
        assert(nodeBegin.getNodeID() == nodeEnd.getNodeID());
        assert(pattEnd >= pattBegin);
        assert(nodeEnd.getPosition() >= nodeBegin.getPosition());

        SSNode node = dbg.getSSNode(nodeBegin.getNodeID());

        // alignment lengths are fixed
        size_t nodeAlLen = nodeEnd.getPosition() - nodeBegin.getPosition();
        size_t pattAlLen = pattEnd - pattBegin;

        string nodeAlStr = node.substr(nodeBegin.getPosition(), nodeAlLen);
        string pattAlStr = patt.substr(pattBegin, pattAlLen);

        AlnRes alnRes = aligner.align(nodeAlStr, pattAlStr);

        //aligner.printAlignment(alnRes, nodeAlStr, pattAlStr);

        float relScore = (float)alnRes.score / (float)aligner.getMaxScore(pattAlLen);

        return NodeAlignment(node.getNodeID(), nodeBegin.getPosition(),
                             nodeEnd.getPosition(), pattBegin, pattEnd,
                             alnRes.score, relScore);
}


SearchRes GraphAligner::DFSAln_guided(const string& P, const NodePosPair& srcNPP,
                                      GraphAlignment& bestGA, bool markovFilter, bool branchAndBoundFileter,  vector<NodeID>seedNodeChain)
{
        //
        seedNodeChain.pop_back();
        GraphAlignment currGA;
        size_t numVisits = 0;
        stack<QueueEl> pq;
        pq.push( QueueEl(alnToNodeOpenEnd(P, 0, srcNPP), 1) );
        vector<GraphAlignment> bestResults;
        while (!pq.empty() && (numVisits++ < maxVisits))
        {
                // pick the first element from the priority queue
                NodeAlignment currNA = pq.top().first;
                size_t depth = pq.top().second;
                pq.pop();
                currGA.updatePath(currNA, depth);
                //visualizeNodeAln(currNA, currGA, P);

                if (P.length() == currGA.getAlignedLength()) {  // fully aligned
                        if (currGA.getScore() > bestGA.getScore()){
                                bestGA = currGA;

                                bestResults.clear();
                                bestResults.push_back(bestGA);
                        }
                        else{
                                if  (currGA.getScore()== bestGA.getScore()){
                                        bestResults.push_back(currGA);
                                }
                        }
                        continue;
                }
                // explore all right neighbors of the current node
                vector<NodeAlignment> v;
                SSNode node = dbg.getSSNode(currNA.nodeID);

                // explore this path only if it has potential to improve bestGA
                if (branchAndBoundFileter && getMaxAttainableScore(P, currGA) < bestGA.getScore())
                        continue;       // branch-and-bound

                vector<NodeID> nodeChain;
                set<NodeID> eligibleNodes;
                if (markovFilter){
                        currGA.getNodeChain(nodeChain );
                        nodeChain.insert(nodeChain.begin(),seedNodeChain.begin(),seedNodeChain.end());
                        //by default all the right nodes are eligibleNodes, except we are sure some of them are not
                        eligibleNodes = mch.getPotentialNextNodes( nodeChain);
                     
                }
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {

                        NodePosPair nextBegin = getNodeBegin(it->getNodeID());
                        if (markovFilter && eligibleNodes.find( it->getNodeID()) == eligibleNodes.end())
                                continue;
                        v.push_back(alnToNodeOpenEnd(P, currNA.pattEnd, nextBegin));
                }
                // greedily prioritize towards node that appears to be best
                sort(v.begin(), v.end());
                for (const auto& nextNA : v)
                        pq.push( QueueEl(nextNA, depth + 1) );
        }
        if (bestResults.size() >1 )
                        return MULTIPLE;
        else
                return (pq.empty() ) ? OPTIMAL : EXHAUSTED;
                
        

}


//jan version
SearchRes GraphAligner::DFSAln(const string& P, const NodePosPair& srcNPP,
                                 GraphAlignment& bestGA)
{
        GraphAlignment currGA;
        size_t numVisits = 0;

        stack<QueueEl> pq;
        pq.push( QueueEl(alnToNodeOpenEnd(P, 0, srcNPP), 1) );

        while (!pq.empty() && (numVisits++ < maxVisits))
        {
                // pick the first element from the priority queue
                NodeAlignment currNA = pq.top().first;
                size_t depth = pq.top().second;
                pq.pop();
                currGA.updatePath(currNA, depth);
                //visualizeNodeAln(currNA, currGA, P);

                if (P.length() == currGA.getAlignedLength()) {  // fully aligned
                        if (currGA.getScore() > bestGA.getScore())
                                bestGA = currGA;
                        continue;
                }

                // explore all right neighbors of the current node
                vector<NodeAlignment> v;
                SSNode node = dbg.getSSNode(currNA.nodeID);
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        // explore this path only if it has potential to improve bestGA
                        if (getMaxAttainableScore(P, currGA) <= bestGA.getScore())
                                continue;       // branch-and-bound

                        NodePosPair nextBegin = getNodeBegin(it->getNodeID());
                        v.push_back(alnToNodeOpenEnd(P, currNA.pattEnd, nextBegin));
                }

                // greedily prioritize towards node that appears to be best
                sort(v.begin(), v.end());
                for (const auto& nextNA : v)
                        pq.push( QueueEl(nextNA, depth + 1) );
        }

        return pq.empty() ? OPTIMAL : EXHAUSTED;
}



SearchRes GraphAligner::DFSAlnTo(const string& P,
                                   const NodePosPair& srcBegin,
                                   const NodePosPair& dstEnd,
                                   const NodeProp<size_t>& lenToDst,
                                   GraphAlignment& bestGA)
{
        GraphAlignment currGA;  // graph alignment currently under investigation
        size_t numVisits = 0;   // number of nodes visited
        stack<QueueEl> pq;      // priority queue (stack = DFS)

        // add the alignment to the end of the source node
        NodePosPair srcEnd = getNodeEnd(srcBegin.getNodeID());
        pq.push( QueueEl(nodeAlnTo(P, 0, srcBegin, srcEnd), 1) );

        // at the destination node, also check the alignment to dstEnd
        if (srcBegin.getNodeID() == dstEnd.getNodeID())
                pq.push( QueueEl(nodeAln(P, 0, P.length(), srcBegin, dstEnd), 1) );

        while (!pq.empty() && (numVisits++ < maxVisits))
        {
                // pop and apply the top element from the priority queue
                NodeAlignment currNA = pq.top().first;
                size_t currDepth = pq.top().second;
                pq.pop();
                currGA.updatePath(currNA, currDepth);
                visualizeNodeAln(currNA, currGA, P);

                if (currNA.getNodeEnd() == dstEnd) {    // dstEnd reached
                        if (currGA.getScore() > bestGA.getScore())
                                bestGA = currGA;
                        continue;
                }

                // explore all right neighbors of the current node
                vector<NodeAlignment> v;
                SSNode node = dbg.getSSNode(currNA.nodeID);
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        NodeID nextID = it->getNodeID();

                        // at dst node, also check the alignment to dstEnd
                        if (nextID == dstEnd.getNodeID())
                                v.push_back(nodeAln(P, currNA.pattEnd, P.length(),
                                                    getNodeBegin(nextID), dstEnd));

                        // explore node only if it has potential to improve bestGA
                        if (getMaxAttainableScore(P, currGA, lenToDst[nextID]) <= bestGA.getScore())
                                continue;       // branch-and-bound



                        // push the alignment to the entire node
                        v.push_back(nodeAlnTo(P, currNA.pattEnd, getNodeBegin(nextID),
                                              getNodeEnd(nextID)));
                }

                // greedily prioritize towards node that appears to be best
                sort(v.begin(), v.end());
                for (const auto& nextNA : v)
                        pq.push( QueueEl(nextNA, currDepth + 1) );
        }

        return pq.empty() ? OPTIMAL : EXHAUSTED;
}

void GraphAligner::findShortestPathToNPP(const NodePosPair& dstEnd,
                                         size_t maxLen,
                                         NodeProp<size_t>& lenToDst,
                                         vector<NodeID>& modNodes) const
{
        // we perform a Dijkstra walk to the left
        priority_queue<PathDFS, vector<PathDFS>, PathDFSComp> heap;
        heap.push( PathDFS(dstEnd.getNodeID(), 0, dstEnd.getPosition()) );

        while( !heap.empty() ) {
                // get the top node
                PathDFS currTop = heap.top();
                heap.pop();

                const NodeID& currID = currTop.nodeID;

                // if node was visited before, get out: this is possible because
                // we don't update distances in the PQ but rather insert doubles
                // Note: the first visit to a node always represents the
                // shortest path from the source node to that node
                if (currTop.length >= lenToDst[currID])
                        continue;

                lenToDst[currID] = currTop.length;
                modNodes.push_back(currID);

                // if the path is too long
                if (currTop.length > maxLen)
                        continue;

                // add the left neighbors of the current node to the PQ
                SSNode curr = dbg.getSSNode(currID);
                for (ArcIt it = curr.leftBegin(); it != curr.leftEnd(); it++) {
                        SSNode left = dbg.getSSNode(it->getNodeID());
                        size_t leftLength = currTop.length + left.getMarginalLength();

                        // if left was not visited before
                        if (leftLength < lenToDst[it->getNodeID()])
                                heap.push( PathDFS(it->getNodeID(), currID, leftLength) );
                }
        }

        // shortest path from the dst node to dstNPP is the loop that first
        // contains the entire dst node (by definition)
        SSNode dst = dbg.getSSNode(dstEnd.getNodeID());
        size_t dstLen = numeric_limits<size_t>::max();
        for (ArcIt it = dst.rightBegin(); it != dst.rightEnd(); it++)
                dstLen = min(lenToDst[it->getNodeID()], dstLen);
        if (dstLen < numeric_limits<size_t>::max())
                dstLen += dst.getMarginalLength();
        lenToDst[dst.getNodeID()] = dstLen;
}

void GraphAligner::findParallelPath(NodeID nodeID, GraphAlignment& bestGraphAln)
{
        bestGraphAln.clear();
        GraphAlignment currGraphAln;

        SSNode node = dbg.getSSNode(nodeID);
        assert(node.isValid());

        string str = node.getSequence();

        // consider all left neighbors ...
        for (ArcIt lIt = node.leftBegin(); lIt != node.leftEnd(); lIt++) {
                SSNode left = dbg.getSSNode(lIt->getNodeID());

                // ... and their right neigbors as starting points
                for (ArcIt rIt = left.rightBegin(); rIt != left.rightEnd(); rIt++) {
                        // skip the original node
                        if (rIt->getNodeID() == nodeID)
                                continue;

                        DFSAln(str, NodePosPair(rIt->getNodeID(), 0), currGraphAln);

                        // if the current path improves the score, update current best
                        if (currGraphAln.getScore() < bestGraphAln.getScore())
                                bestGraphAln = currGraphAln;
                }
        }
}
void GraphAlignment::getNodeChain(vector<NodeID> &visitedNodes) {

        for (NodeAlignment nodeAl :path)
                visitedNodes.push_back( nodeAl.nodeID);

}

